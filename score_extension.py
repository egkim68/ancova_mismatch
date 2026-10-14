"""Scores H1-H7 (weighted, then held-out) and E1/E2 of Addendum A on a results CSV.
Usage: python score_extension.py output/ancova_ext/ancova_results.csv
Weights: layouts with 2, 3 and 4 levels carry the shares .40, .35, .25 of the main run; 5 levels carry zero."""
import sys
import numpy as np, pandas as pd

M = ['Oneway', 'R_default', 'TypeII', 'ANCOVA', 'TypeIII_raw', 'TypeIII_ctr', 'TypeIII_ctr_HC3', 'Drop_int']
R = ['rej_' + m for m in M]
SHARE = {2: .40, 3: .35, 4: .25}

def wmean(sub, col_or_fn):
    """Weighted mean over layouts: mean within each layout, combined with SHARE over the layouts present."""
    vals, ws = [], []
    for a, w in SHARE.items():
        g = sub[sub.a == a]
        if len(g): vals.append(col_or_fn(g)); ws.append(w)
    ws = np.array(ws) / sum(ws)
    return np.average(np.array(vals), axis=0, weights=ws)

def harmful(g): return (g[R] > 0.10).mean().values * 100

def score(null, label):
    ok = {}
    inter, ci, far, uv = null.slope_range > 0, null.cov_imbalance > 0, null.d0 >= 2, null.sd_ratio > 1
    h = dict(zip(M, wmean(null, harmful)))
    h2 = dict(zip(M, wmean(null[~inter & ci], harmful)))
    h3 = dict(zip(M, wmean(null[inter & far], harmful)))
    he = dict(zip(M, wmean(null[~uv], harmful))); hu = dict(zip(M, wmean(null[uv], harmful)))
    dis = float(wmean(null, lambda g: np.array([g.flip_any.mean() * 100]))[0])
    others = [v for k, v in h.items() if k != 'TypeIII_ctr_HC3']
    ok['H1'] = h['R_default'] >= 50 and h['R_default'] - max(v for k, v in h.items() if k != 'R_default') >= 10
    ok['H2'] = h2['Oneway'] >= 50 and h2['R_default'] >= 50 and h2['TypeII'] <= 20 and h2['ANCOVA'] <= 20
    ok['H3'] = h3['TypeIII_raw'] >= 60 and h3['TypeIII_raw'] - h3['TypeIII_ctr'] >= 30
    ok['H4'] = he['TypeIII_ctr'] <= 5 and hu['TypeIII_ctr'] >= 15
    ok['H5'] = h['TypeIII_ctr_HC3'] <= 10 and h['TypeIII_ctr_HC3'] < min(others)
    ok['H6'] = dis >= 25
    ok['H7'] = h['Drop_int'] >= 30
    print(f"\n[{label}] null designs (layouts 2-4): {int((null.a <= 4).sum())}")
    print('  harmful %:', {k: round(v, 1) for k, v in h.items()})
    print('  no-inter + cov imbalance:', {k: round(v, 1) for k, v in h2.items() if k in ('Oneway', 'R_default', 'TypeII', 'ANCOVA')})
    print('  inter + far:', {k: round(v, 1) for k, v in h3.items() if k in ('TypeIII_raw', 'TypeIII_ctr')})
    print('  equal var ctr %.1f / unequal var ctr %.1f; HC3 %.1f / %.1f; disagreement %.1f' %
          (he['TypeIII_ctr'], hu['TypeIII_ctr'], he['TypeIII_ctr_HC3'], hu['TypeIII_ctr_HC3'], dis))
    return ok

d = pd.read_csv(sys.argv[1]); null = d[d.a_effect == 0].copy()
print(f"{len(d)} designs, {len(null)} null designs; layouts: {null.a.value_counts().sort_index().to_dict()}")
res = {'weighted': score(null, 'weighted'), 'held-out': score(null[null.design_id % 4 == 0], 'held-out')}
print('\nE1 (H1-H7):'); print(pd.DataFrame(res))
e2 = {}
for a in sorted(null.a.unique()):
    g = null[null.a == a]; h = pd.Series(harmful(g), index=M)
    e2[a] = bool(h.idxmax() == 'R_default' and h['TypeIII_ctr_HC3'] < h['TypeIII_ctr'])
    print(f"layout {a}: harmful% R_default {h['R_default']:.1f} (max: {h.idxmax()}), ctr {h['TypeIII_ctr']:.1f}, HC3 {h['TypeIII_ctr_HC3']:.1f}, n={len(g)}")
print('E2 by layout:', e2)
if 'cov_shape' in null:
    print('\nDescriptive, harmful % by covariate shape (HC3 / ctr / raw / R_default):')
    for s, g in null.groupby('cov_shape'): print(' ', s, len(g), dict(zip(['HC3', 'ctr', 'raw', 'R_default'], harmful(g)[[6, 5, 4, 1]].round(1))))
    print('Descriptive, harmful % by smallest group (HC3 / ctr):')
    for s, g in null.groupby(pd.cut(null.min_group, [0, 3, 5, 9, 1e9], labels=['3', '4-5', '6-9', '10+'])): print(' ', s, len(g), harmful(g)[[6, 5]].round(1))
