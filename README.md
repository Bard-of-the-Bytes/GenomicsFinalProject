# GenomicsFinalProject
Bowtie 2 adaptive window sizing project for Computational Genomics.

Input Files:
- data/reads1.fq. Mate A reads in FASTQ format.
- data/reads2.fq. Mate B reads in FASTQ format.
- data/ref.fa. Reference genome in FASTA format. This is the Bowtie 2 lambda phage reference.

Output Files:
- output/A_only.sam. SAM output from aligning mate A alone, providing information such as alignment score, edit distance, strand, and position.
- output/paired.sam. SAM output from paired-end alignment using both mates, providing Bowtie 2’s reported candidate pairings.
- output/regression_table.csv. Intermediate CSV produced by parser.py and read by regression.py.
- index/ref_index.*.bt2. Bowtie 2 index files produced from ref.fa.

Scripts:
- parser.py. Python parser which combines the A-only SAM file and paired-end SAM file into a regression table.
- regression.py. Python linear regression script which reads the regression table, fits the linear regression model, and prints the final adaptive-window results.

How to Run the Project from Terminal.
- Ensure Bowtie 2, Python, and all relevant dependencies are installed.
- cd ~/GenomicsFinalProject-main (enter the project directory).
- bowtie2-build data/ref.fa index/ref_index (build the Bowtie 2 index from ref.fa, producing index/ref_index.*.bt2 files).
- bowtie2 -x index/ref_index -U data/reads1.fq -a -S output/A_only.sam (align mate A, producing output/A_only.sam).
- bowtie2 -x index/ref_index -1 data/reads1.fq -2 data/reads2.fq -a --no-mixed --no-discordant -I 0 -X 500 -S output/paired.sam (align both mates as paired ends, producing output/paired.sam).
- python3 parser.py (build the regression table, producing output/regression_table.csv).
- python3 regression.py (fit the final regression model, producing terminal output).
