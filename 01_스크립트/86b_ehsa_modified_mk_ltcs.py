# 86b_ehsa_modified_mk_ltcs.py - Emerging hot spot trend test variants: classic Mann-Kendall, Hamed-Rao autocorrelation-corrected Mann-Kendall, and Mann-Kendall with BH-FDR
#  across districts. Output: results_std3/S27_ehsa_mk_variants.csv, S27_ehsa_mk_summary.csv
import os, sys, math, numpy as np, pandas as pd, pymannkendall as pmk
from scipy.stats import kendalltau
from statsmodels.stats.multitest import multipletests
R = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "04_분석", "LTCS", "results_std3")
eb = pd.read_csv(f"{R}/S12c_yearly_EB_SIR.csv", encoding="utf-8-sig"); s6 = pd.read_csv(f"{R}/S6_spatial_RR.csv")
area = dict(zip(s6.fbd_region, s6.area)); L = open(f"{R}/nb_LTCS.graph").read().split("\n"); nb = {}
for line in L[1:]:
    t = line.split()
    if len(t) >= 2: nb[int(t[0])] = [int(x) for x in t[2:]]
EBm = eb.pivot(index="year", columns="fbd_region", values="ebSIR").dropna(axis=1)
regs = list(EBm.columns); pos = {area[r]: k for k, r in enumerate(regs)}
nbo = [[pos[a] for a in nb[area[r]] if a in pos] for r in regs]
YRS = list(EBm.index); T = len(YRS); thr = math.ceil(0.9 * T); X = EBm.values
def gistar(x):
    nn = len(x); xb = x.mean(); s = math.sqrt((x ** 2).sum() / nn - xb ** 2); z = np.zeros(nn)
    for i in range(nn):
        nbr = [i] + nbo[i]; k = len(nbr); w = 1.0 / k
        z[i] = (x[nbr].sum() * w - xb) / (s * math.sqrt((nn * (k * w * w) - 1.0) / (nn - 1)))
    return z
Z = np.vstack([gistar(X[t]) for t in range(T)])
yrs = np.array(YRS, float)
def last_run(b):
    c = 0
    for v in b[::-1]:
        if v: c += 1
        else: break
    return c
def esri(zv, tau, p):
    hot = zv > 1.96; cold = zv < -1.96; nh, nc = hot.sum(), cold.sum(); up = (p < 0.05) and tau > 0; dn = (p < 0.05) and tau < 0
    if hot[-1]:
        if nh >= thr: return "Intensifying hot" if up else ("Diminishing hot" if dn else "Persistent hot")
        if nh == 1: return "New hot"
        if last_run(hot) >= 2 and last_run(hot) == nh: return "Consecutive hot"
        return "Oscillating hot" if cold.any() else "Sporadic hot"
    if cold[-1]:
        if nc >= thr: return "Intensifying cold" if dn else ("Diminishing cold" if up else "Persistent cold")
        if nc == 1: return "New cold"
        if last_run(cold) >= 2 and last_run(cold) == nc: return "Consecutive cold"
        return "Oscillating cold" if hot.any() else "Sporadic cold"
    if nh >= thr: return "Historical hot"
    if nc >= thr: return "Historical cold"
    return "No pattern"
m = Z.shape[1]
orig = [kendalltau(yrs, Z[:, k]) for k in range(m)]
hr = [pmk.hamed_rao_modification_test(Z[:, k]) for k in range(m)]
p_orig = np.array([o[1] for o in orig]); p_fdr = multipletests(p_orig, method="fdr_bh")[1]
rows = []
for k in range(m):
    tau = orig[k][0]
    rows.append(dict(fbd_region=regs[k], EHSA_original=esri(Z[:, k], tau, p_orig[k]), EHSA_hamed_rao=esri(Z[:, k], hr[k].Tau, hr[k].p),
                     EHSA_mk_bh_fdr=esri(Z[:, k], tau, p_fdr[k]), mk_tau=round(tau, 3), mk_p=round(p_orig[k], 4), hr_p=round(hr[k].p, 4), mk_q_bh=round(p_fdr[k], 4)))
V = pd.DataFrame(rows)
s11 = pd.read_csv(f"{R}/S11_ehsa.csv", encoding="utf-8-sig")
chk = V.merge(s11[["fbd_region", "EHSA"]], on="fbd_region"); agree = (chk.EHSA == chk.EHSA_original).mean()
print("reproduce S11 agreement:", round(agree, 4))
V.to_csv(f"{R}/S27_ehsa_mk_variants.csv", index=False, encoding="utf-8-sig")
cats = sorted(set(V.EHSA_original) | set(V.EHSA_hamed_rao) | set(V.EHSA_mk_bh_fdr))
S = pd.DataFrame({"category": cats, "original_MK": [(V.EHSA_original == c).sum() for c in cats], "hamed_rao_MK": [(V.EHSA_hamed_rao == c).sum() for c in cats], "MK_BH_FDR": [(V.EHSA_mk_bh_fdr == c).sum() for c in cats]})
S.loc[len(S)] = ["changed_vs_original", 0, int((V.EHSA_hamed_rao != V.EHSA_original).sum()), int((V.EHSA_mk_bh_fdr != V.EHSA_original).sum())]
S.to_csv(f"{R}/S27_ehsa_mk_summary.csv", index=False, encoding="utf-8-sig"); print(S.to_string(index=False))
