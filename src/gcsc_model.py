# Basic prototype

import numpy as np
import pandas as pd
import statsmodels.api as sm

# X: individuals x genes
# twas_chi2: Series indexed by gene
# gene_set: list of genes in selected set
# N: GWAS sample size

def simplified_gcsc(X, twas_chi2, gene_set, N):

    # Calculate gene co-regulation matrix
    R = X.corr()

    # Calculate squared correlations
    R2 = R ** 2

    # All-gene co-regulation scores
    c0 = R2.sum(axis=1)

    # Gene-set co-regulation scores
    cS = R2[gene_set].sum(axis=1)

    # Prepare regression data
    genes = twas_chi2.index.intersection(R.index)

    y = twas_chi2.loc[genes] - 1
    predictors = pd.DataFrame({
        "baseline": N * c0.loc[genes],
        "gene_set": N * cS.loc[genes]
    })

    # Fit simplified regression
    model = sm.OLS(
        y,
        predictors
    ).fit()

    # Estimate gene-set heritability component
    tau_S = model.params["gene_set"]
    M_set = len(gene_set)

    h2_set = M_set * tau_S

    return {
        "tau_0": model.params["baseline"],
        "tau_S": tau_S,
        "h2_set": h2_set,
        "p_value": model.pvalues["gene_set"],
        "model": model
    }