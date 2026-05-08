import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_absolute_error


CSV_PATH = "output/regression_table.csv"


def main():
    df = pd.read_csv(CSV_PATH)

    print(f"Loaded {len(df)} rows")

    bowtie_fixed_width = df["bowtie_fixed_width"].iloc[0]
    print(f"Bowtie fixed window width: {bowtie_fixed_width}")

    df["A_strand_minus"] = (df["A_strand"] == "-").astype(int)

    feature_cols = ["A_AS", "A_NM", "A_strand_minus"]
    target_col = "target_z"

    X = df[feature_cols]
    y = df[target_col]

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.25,
        random_state=42
    )

    model = LinearRegression()
    model.fit(X_train, y_train)

    pred_train = model.predict(X_train)
    pred_test = model.predict(X_test)

    pred_train = np.maximum(0, pred_train)
    pred_test = np.maximum(0, pred_test)

    print("\nLinear regression model:")
    print(f"Intercept: {model.intercept_:.3f}")
    for name, coef in zip(feature_cols, model.coef_):
        print(f"{name}: {coef:.3f}")

    print("\nPrediction quality:")
    print(f"Train R^2: {r2_score(y_train, pred_train):.3f}")
    print(f"Test R^2:  {r2_score(y_test, pred_test):.3f}")
    print(f"Test MAE:  {mean_absolute_error(y_test, pred_test):.3f}")

    oracle_avg_width = np.mean(2 * y)
    print("\nOracle lower bound:")
    print(f"Average minimum width needed: {oracle_avg_width:.2f}")
    print(
        f"Oracle width reduction vs Bowtie: "
        f"{100 * (bowtie_fixed_width - oracle_avg_width) / bowtie_fixed_width:.2f}%"
    )


    train_underprediction = y_train - pred_train

    print("\nAdaptive rule evaluation with safety margins:")
    print("Quantile\tMargin\tRetention\tAvg adaptive width\tWidth reduction")

    for q in [0.50, 0.75, 0.80, 0.90, 0.95, 0.975, 0.99]:
        margin = np.quantile(train_underprediction, q)

        margin = max(0, margin)

        z_hat = pred_test + margin

        retained = z_hat >= y_test
        retention = retained.mean()

        avg_adaptive_width = np.mean(2 * z_hat)

        width_reduction = (
            bowtie_fixed_width - avg_adaptive_width
        ) / bowtie_fixed_width

        print(
            f"{q:.3f}\t\t"
            f"{margin:.2f}\t"
            f"{retention:.3f}\t\t"
            f"{avg_adaptive_width:.2f}\t\t\t"
            f"{100 * width_reduction:.2f}%"
        )


if __name__ == "__main__":
    main()
