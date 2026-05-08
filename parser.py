from collections import defaultdict
import csv
import re


A_ONLY_SAM = "output/A_only.sam"
PAIRED_SAM = "output/paired.sam"
OUT_CSV = "output/regression_table.csv"

MIN_FRAG = 0
MAX_FRAG = 500


def parse_optional_tags(fields):
    tags = {}

    for field in fields[11:]:
        parts = field.split(":", 2)
        if len(parts) != 3:
            continue

        tag, tag_type, value = parts

        if tag_type == "i":
            value = int(value)
        elif tag_type == "f":
            value = float(value)

        tags[tag] = value

    return tags


def is_unmapped(flag):
    return bool(flag & 4)


def is_reverse_strand(flag):
    return bool(flag & 16)


def is_first_in_pair(flag):
    return bool(flag & 64)


def is_second_in_pair(flag):
    return bool(flag & 128)


def is_secondary(flag):
    return bool(flag & 256)


def strand_from_flag(flag):
    return "-" if is_reverse_strand(flag) else "+"


def reference_length_from_cigar(cigar):
    if cigar == "*":
        return None

    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(length)
        if op in {"M", "D", "N", "=", "X"}:
            total += length

    return total


def read_sam_records(sam_path):
    records = []

    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")

            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            pnext = int(fields[7]) if fields[7] != "0" else 0
            tlen = int(fields[8])
            seq = fields[9]

            if is_unmapped(flag):
                continue

            tags = parse_optional_tags(fields)

            records.append({
                "qname": qname,
                "flag": flag,
                "rname": rname,
                "pos": pos,
                "mapq": mapq,
                "cigar": cigar,
                "pnext": pnext,
                "tlen": tlen,
                "seq": seq,
                "strand": strand_from_flag(flag),
                "ref_len": reference_length_from_cigar(cigar),
                "AS": tags.get("AS"),
                "XS": tags.get("XS"),
                "NM": tags.get("NM"),
            })

    return records


def parse_a_only_features(a_only_sam):

    records = read_sam_records(a_only_sam)
    records_by_qname = defaultdict(list)

    for rec in records:
        records_by_qname[rec["qname"]].append(rec)

    a_features = {}

    for qname, recs in records_by_qname.items():
        recs_sorted = sorted(
            recs,
            key=lambda r: r["AS"] if r["AS"] is not None else -10**9,
            reverse=True
        )

        candidate_count = len(recs_sorted)

        for rank, rec in enumerate(recs_sorted, start=1):
            key = (qname, rec["pos"], rec["strand"])

            a_as = rec["AS"]
            a_xs = rec["XS"]

            if a_as is not None and a_xs is not None:
                a_score_gap = a_as - a_xs
            else:
                a_score_gap = None

            a_mapq_clean = None if rec["mapq"] == 255 else rec["mapq"]

            a_features[key] = {
                "qname": qname,
                "A_pos": rec["pos"],
                "A_strand": rec["strand"],
                "A_AS": a_as,
                "A_XS": a_xs,
                "A_score_gap": a_score_gap,
                "A_MAPQ": rec["mapq"],
                "A_MAPQ_clean": a_mapq_clean,
                "A_NM": rec["NM"],
                "A_candidate_count": candidate_count,
                "A_rank": rank,
                "A_ref_len": rec["ref_len"],
            }

    return a_features


def parse_paired_candidates(paired_sam):
    records = read_sam_records(paired_sam)
    records_by_qname = defaultdict(list)

    for rec in records:
        records_by_qname[rec["qname"]].append(rec)

    candidate_pairs = []

    for qname, recs in records_by_qname.items():
        first_mates = [r for r in recs if is_first_in_pair(r["flag"])]
        second_mates = [r for r in recs if is_second_in_pair(r["flag"])]

        used_second_indices = set()

        for a in first_mates:
            match_index = None

            for j, b in enumerate(second_mates):
                if j in used_second_indices:
                    continue

                same_candidate_pair = (
                    a["pos"] == b["pnext"] and
                    a["pnext"] == b["pos"] and
                    a["tlen"] == -b["tlen"]
                )

                if same_candidate_pair:
                    match_index = j
                    break

            if match_index is None:
                continue

            b = second_mates[match_index]
            used_second_indices.add(match_index)

            candidate_pairs.append({
                "qname": qname,
                "A_pos": a["pos"],
                "A_strand": a["strand"],
                "B_pos": b["pos"],
                "B_strand": b["strand"],
                "A_paired_AS": a["AS"],
                "B_paired_AS": b["AS"],
                "TLEN": a["tlen"],
                "primary": (not is_secondary(a["flag"])) and (not is_secondary(b["flag"])),
                "A_ref_len_paired": a["ref_len"],
                "B_ref_len_paired": b["ref_len"],
            })

    return candidate_pairs


def expected_b_center_from_a(pair):
    
    expected_fragment_length = (MIN_FRAG + MAX_FRAG) / 2

    a_pos = pair["A_pos"]
    b_pos = pair["B_pos"]

    a_len = pair["A_ref_len_paired"] or 0
    b_len = pair["B_ref_len_paired"] or 0

    if b_pos >= a_pos:
        return a_pos + expected_fragment_length - b_len
    else:
        return a_pos - expected_fragment_length + a_len


def build_regression_rows():
    a_features = parse_a_only_features(A_ONLY_SAM)
    paired_candidates = parse_paired_candidates(PAIRED_SAM)

    rows = []

    for pair in paired_candidates:
        key = (pair["qname"], pair["A_pos"], pair["A_strand"])

        if key not in a_features:
            continue

        features = a_features[key]

        expected_center = expected_b_center_from_a(pair)

        target_z = abs(pair["B_pos"] - expected_center)

        bowtie_fixed_width = MAX_FRAG - MIN_FRAG

        row = {
            "qname": pair["qname"],
            "A_pos": pair["A_pos"],
            "A_strand": pair["A_strand"],
            "B_pos": pair["B_pos"],
            "B_strand": pair["B_strand"],

            "A_AS": features["A_AS"],
            "A_XS": features["A_XS"],
            "A_score_gap": features["A_score_gap"],
            "A_MAPQ": features["A_MAPQ"],
            "A_MAPQ_clean": features["A_MAPQ_clean"],
            "A_NM": features["A_NM"],
            "A_candidate_count": features["A_candidate_count"],
            "A_rank": features["A_rank"],

            "expected_B_center": expected_center,
            "target_z": target_z,
            "target_width": 2 * target_z,

            "bowtie_fixed_width": bowtie_fixed_width,

            "TLEN_for_evaluation_only": pair["TLEN"],
            "primary": pair["primary"],
        }

        rows.append(row)

    return rows


def write_csv(rows, out_csv):
    if not rows:
        print("No rows were created.")
        return

    fieldnames = list(rows[0].keys())

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows):
    print(f"Number of regression rows: {len(rows)}")

    if not rows:
        return

    avg_target_width = sum(r["target_width"] for r in rows) / len(rows)
    bowtie_width = rows[0]["bowtie_fixed_width"]

    print(f"Bowtie fixed window width: {bowtie_width}")
    print(f"Average minimum adaptive width needed: {avg_target_width:.2f}")
    print(f"Width reduction relative to Bowtie fixed width: {bowtie_width - avg_target_width:.2f}")

    print("\nFirst 10 rows:")
    for row in rows[:10]:
        print(row)


if __name__ == "__main__":
    rows = build_regression_rows()
    write_csv(rows, OUT_CSV)
    print_summary(rows)
    print(f"\nWrote regression table to: {OUT_CSV}")
