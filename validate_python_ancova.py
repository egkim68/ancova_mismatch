"""Recomputes the eight analyses with statsmodels on the data sets exported by export_for_python.R and
compares them with the p values computed in R. Needs numpy, pandas, scipy, statsmodels."""
import sys
import numpy as np, pandas as pd, scipy, statsmodels
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

ALPHA = 0.05
data = pd.read_csv("check/datasets.csv"); rp = pd.read_csv("check/r_pvalues.csv").set_index("case")
print(f"Python {sys.version.split()[0]} | numpy {np.__version__} | scipy {scipy.__version__} | statsmodels {statsmodels.__version__}")
rows = []
for case, d in data.groupby("case"):
    d = d.copy(); d["A"] = d["A"].astype("category"); d["xc"] = d["x"] - d["x"].mean()
    p = {}
    p["Oneway"] = anova_lm(smf.ols("y ~ C(A)", d).fit()).loc["C(A)", "PR(>F)"]
    p["R_default"] = anova_lm(smf.ols("y ~ C(A) * x", d).fit()).loc["C(A)", "PR(>F)"]
    p["TypeII"] = anova_lm(smf.ols("y ~ C(A) * x", d).fit(), typ=2).loc["C(A)", "PR(>F)"]
    p["ANCOVA"] = anova_lm(smf.ols("y ~ x + C(A)", d).fit(), typ=2).loc["C(A)", "PR(>F)"]
    p["TypeIII_raw"] = anova_lm(smf.ols("y ~ C(A, Sum) * x", d).fit(), typ=3).loc["C(A, Sum)", "PR(>F)"]
    fc = smf.ols("y ~ C(A, Sum) * xc", d).fit()
    p["TypeIII_ctr"] = anova_lm(fc, typ=3).loc["C(A, Sum)", "PR(>F)"]
    p["TypeIII_ctr_HC3"] = anova_lm(fc, typ=3, robust="hc3").loc["C(A, Sum)", "PR(>F)"]
    p_int = anova_lm(smf.ols("y ~ C(A) * x", d).fit(), typ=2).loc["C(A):x", "PR(>F)"]
    p["Drop_int"] = p["ANCOVA"] if p_int > ALPHA else p["TypeIII_raw"]
    rows.append(pd.Series({k: abs(v - rp.loc[case, k]) for k, v in p.items()}, name=case))
diff = pd.DataFrame(rows)
print(f"{len(diff)} data sets; largest absolute difference in p over all analyses: {diff.values.max():.2e}")
print(diff.max().apply(lambda v: f"{v:.2e}").to_string())
