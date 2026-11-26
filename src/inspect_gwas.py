import pandas as pd


def inspect_gwas(path):

    print("Inspecting GWAS summary statistics")

    df = pd.read_csv(
        path,
        sep=r"\s+",
        compression="gzip"
    )

    print("\nDATASET DIMENSIONS:")
    print(df.shape)

    print("\nMISSING VALUES:")
    print(df.isna().sum())

    print("\nDUPLICATE SNP IDENTIFIERS:")
    print(df["SNP"].duplicated().sum())

    print("\nSAMPLE SIZE:")
    print(df["N"].describe())

    print("\nCHI-SQUARE STATISTICS:")
    print(df["CHISQ"].describe())

    print("\nZ-SCORE STATISTICS:")
    print(df["Z"].describe())

    print("\nCHI-SQUARE VS Z-SCORE")

    df["Z_SQUARED"] = df["Z"] ** 2

    df["CHI2_DIFFERENCE"] = (
        df["CHISQ"] - df["Z_SQUARED"]
    )

    print(
        "Mean absolute difference:",
        df["CHI2_DIFFERENCE"].abs().mean()
    )

    print(
        "Maximum absolute difference:",
        df["CHI2_DIFFERENCE"].abs().max()
    )

    print(
        "Correlation between CHISQ and Z squared:",
        df["CHISQ"].corr(df["Z_SQUARED"])
    )

    print(
        "Variants with difference > 0.01:",
        (df["CHI2_DIFFERENCE"].abs() > 0.01).sum()
    )

    print("\nALLELE COUNTS:")
    print(df["A1"].value_counts())
    print(df["A2"].value_counts())

    print("\nINVALID NUMERIC VALUES:")
    print("Non-finite CHISQ:", (~df["CHISQ"].map(
        lambda x: pd.notna(x) and abs(x) != float("inf")
    )).sum())

    return df


if __name__ == "__main__":

    path = "data/GWAS/UKB_460K.blood_WHITE_COUNT.sumstats.gz"

    inspect_gwas(path)