#!/usr/bin/env python3
# 90_ltcs_summary_numbers.py - Robustness verdicts and every number reported in the manuscript (unrounded values plus one half-up rounding helper), and a SHA-256
#  signature of inputs, code and result tables. No model fitting. Output: numbers_ltcs.json, results_std3/S90_robustness_verdict.csv, result_signature_LTCS.json
import os, glob, json, math, hashlib
import numpy as np, pandas as pd
from decimal import Decimal, ROUND_HALF_UP
HERE = os.path.dirname(os.path.abspath(__file__)); SK = os.path.dirname(os.path.dirname(HERE)); R = os.path.join(HERE, "results_std3")
def rnd(x, nd=2):
    if x is None or (isinstance(x, float) and math.isnan(x)): return None
    q = Decimal(1).scaleb(-nd); return float(Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP))
rd = lambda f: pd.read_csv(f) if os.path.exists(f) else None
M3 = rd(f"{R}/S3_main_IRR.csv"); FIT = rd(f"{R}/S3_main_fit.csv"); X0 = rd(f"{R}/S60_extra_X0_main_fullprec.csv")
HYP = M3[M3.role == "hypothesis"]["cov"].tolist(); FORCED = M3[M3.role == "forced"]["cov"].tolist()
FP = X0.set_index("cov") if X0 is not None else None
chk = {}
if FP is not None:
    for v in M3["cov"]: chk[v] = abs(round(math.exp(FP.beta[v]), 3) - M3.set_index("cov").IRR[v])
S7 = rd(f"{R}/S7_multiverse_summary.csv"); S7i = S7.set_index("cov") if S7 is not None else None
ALT = []
for f in sorted(glob.glob(f"{R}/S41_part_*.csv")):
    d = pd.read_csv(f)
    if "analysis" in d and not d.analysis.astype(str).str.startswith(("G2_noise", "H1_nco", "bw_twoway")).all(): ALT.append(d.assign(src="88g")[["analysis", "cov", "IRR", "lo", "hi", "cred", "src"]])
for f in sorted(glob.glob(f"{R}/S60_extra_X*.csv")):
    d = pd.read_csv(f); an = d.analysis.iloc[0]
    if an.startswith(("X0_", "X3_", "X7_")): continue
    ALT.append(d.assign(src="88X")[["analysis", "cov", "IRR", "lo", "hi", "cred", "src"]])
for f in sorted(glob.glob(f"{R}/S43_part_*.csv")):
    d = pd.read_csv(f); ALT.append(d.assign(src="88k")[["analysis", "cov", "IRR", "lo", "hi", "cred", "src"]])
SENSV = ["s577", "s2", "sc", "s8", "y1524", "x2021", "y1324", "sexE", "supply"]
for v in SENSV:
    f = os.path.join(HERE, f"results_{v}", "S3_main_IRR.csv")
    if os.path.exists(f): ALT.append(pd.read_csv(f).assign(analysis=f"sens_{v}", src="sens")[["analysis", "cov", "IRR", "lo", "hi", "cred", "src"]])
A = pd.concat(ALT, ignore_index=True) if ALT else pd.DataFrame(columns=["analysis", "cov", "IRR", "lo", "hi", "cred", "src"])
A = A[~A.analysis.astype(str).str.startswith(("G2_noise", "H1_nco", "bw_twoway"))]
# MI Rubin(88X X7)
MI = [pd.read_csv(f) for f in sorted(glob.glob(f"{R}/S60_extra_X7_mi_r*.csv"))]; RUB = {}
if MI:
    Q = pd.concat(MI)
    for v, b in Q.groupby("cov"):
        m = len(b); qb = b.beta.mean(); W = (b.sd ** 2).mean(); Bv = b.beta.var(ddof=1); T = W + (1 + 1 / m) * Bv
        RUB[v] = dict(m=m, beta=qb, se=math.sqrt(T), IRR=math.exp(qb), lo=math.exp(qb - 1.96 * math.sqrt(T)), hi=math.exp(qb + 1.96 * math.sqrt(T)), cred=bool(abs(qb) > 1.96 * math.sqrt(T)))
    pd.DataFrame(RUB).T.reset_index(names="cov").to_csv(f"{R}/S60_X7_mi_rubin.csv", index=False)
def ev(r): r = 1 / r if r < 1 else r; return r + math.sqrt(r * (r - 1))
rows = []
from statsmodels.stats.multitest import multipletests
for v in HYP:
    m = M3.set_index("cov").loc[v]
    if FP is not None: b, sd = FP.beta[v], FP.sd[v]; irr, lo, hi = math.exp(b), math.exp(b - 1.959964 * sd), math.exp(b + 1.959964 * sd); p2 = 2 * (1 - 0.5 * (1 + math.erf(abs(b / sd) / math.sqrt(2))))
    else: irr, lo, hi, b, sd = m.IRR, m.lo, m.hi, math.log(m.IRR), float("nan"); p2 = 2 * min(m.P_gt1, 1 - m.P_gt1)
    lim = lo if irr > 1 else hi
    a = A[A["cov"] == v]; same = (a.cred.astype(bool) & (np.sign(np.log(a.IRR)) == np.sign(math.log(irr))))
    mv = S7i.loc[v] if S7i is not None and v in S7i.index else None
    rows.append(dict(cov=v, label=m.get("label_en", v), beta=b, sd=sd, IRR=irr, lo=lo, hi=hi, P_gt1=m.P_gt1, cred_main=bool(m.cred), p_two_sided_approx=p2,
                     E_point=ev(irr), E_limit=(ev(lim) if (lim > 1) == (irr > 1) else 1.0),
                     mv_cred_frac=(mv.cred_frac if mv is not None else np.nan), mv_sign_agree=(mv.sign_agree if mv is not None else np.nan), mv_robust=(bool(mv.robust) if mv is not None else False),
                     n_alt=int(a.analysis.nunique()), n_alt_cred_same=int(a[same].analysis.nunique()), alt_frac=(a[same].analysis.nunique() / a.analysis.nunique() if len(a) else np.nan),
                     kh_kept=int(sum(1 for k in ("KH_typeII", "KH_typeIII", "KH_typeIV") if len(a[(a.analysis == k) & same]))),
                     mi_IRR=RUB.get(v, {}).get("IRR"), mi_lo=RUB.get(v, {}).get("lo"), mi_hi=RUB.get(v, {}).get("hi"), mi_cred=RUB.get(v, {}).get("cred")))
T = pd.DataFrame(rows); T["q_BH"] = multipletests(T.p_two_sided_approx, method="fdr_bh")[1]
T["verdict"] = np.where(~T.cred_main, "무연관", np.where(T.mv_robust & (T.alt_frac >= 0.80), "강건", "부분"))
T.to_csv(os.path.join(R, "S90_robustness_verdict.csv"), index=False)
print(T[["cov", "IRR", "lo", "hi", "cred_main", "mv_cred_frac", "mv_sign_agree", "alt_frac", "n_alt", "q_BH", "E_point", "E_limit", "verdict"]].round(3).to_string())
# numbers json
NUM = {"_meta": {"tag": "LTCS", "issue": 72, "version": "std3", "created": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"), "round_function": "rnd(x, nd) = Decimal(repr(float(x))).quantize(10**-nd, ROUND_HALF_UP); every displayed manuscript value uses this one function", "raw_values": "반올림 전 원값", "main_source": "S60_extra_X0_main_fullprec.csv(81 S3와 같은 식 재적합, 3자리 일치 확인) + S3_main_*.csv"}}
INFO = json.load(open(os.path.join(HERE, "panel_INFO_LTCS_std3.json")))
NUM["data"] = dict(n_rows=INFO["n_rows"], n_districts=INFO["n_sgg"], n_years=INFO["n_years"], year_min=INFO["years"][0], year_max=INFO["years"][1], total_cases=INFO["total_cases"], total_pop65=INFO["total_offset"],
                   crude_rate_per_1000_65plus=INFO["total_cases"] / INFO["total_offset"] * 1e3, n_corrected_cells_village=44, n_filled_cells_hospitals=4, vif_dropped=INFO.get("vif_dropped"), national_rate_per_100k_by_year=INFO["national_rate_by_year"])
fit = FIT.iloc[0].to_dict() if FIT is not None else {}
NUM["main_fit"] = {k: (float(v) if isinstance(v, (int, float, np.floating, np.integer)) else v) for k, v in fit.items()}
hp = rd(f"{R}/S3_main_hyper.csv"); NUM["main_hyper"] = {r.iloc[0]: dict(mean=r.iloc[1], lo=r.iloc[2], hi=r.iloc[3]) for _, r in hp.iterrows()} if hp is not None else {}
NUM["likelihood_forced_only"] = rd(f"{R}/S1_likelihood.csv").set_index("family")[["WAIC", "DIC"]].to_dict("index")
for nm, f in [("likelihood6", "S1c_likelihood5.csv"), ("refit_twice", "S1b_refit_twice.csv"), ("residual_moran", "S8_residual_moran.csv"), ("negcontrol", "S3_negcontrol.csv"), ("moran_year", "S6_moran_year.csv"),
               ("window_IRR", "S4_window_IRR.csv"), ("interaction", "S4_interaction.csv"), ("covid", "S5_covid.csv"), ("typeI", "S9_typeI_IRR.csv"), ("ehsa_summary", "S11_ehsa_summary.csv"), ("persistence", "S12_persistence.csv"),
               ("variance_partition", "S42_variance_partition.csv"), ("ppc", "S42_ppc.csv"), ("kh_struct", "S43_sensitivity_fit.csv"), ("plus05_moran", "S60_P1_moran_constant.csv"), ("extra_fits", "S60_extra_fit_A.csv"), ("kh", "S43_sensitivity_fit.csv"), ("evalues_88g", "S41_evalues.csv"), ("retention_88g", "S41_retention_by_cov.csv"), ("nco_88g", "S41_nco_detail.csv"), ("ehsa_mk", "S27_ehsa_mk_summary.csv"), ("window_grid", "S4c_window_grid_table.csv"), ("likelihood6_coef", "S1c_likelihood5_coef.csv"), ("repro_meta", "S0_repro_meta.csv"), ("mi_rubin", "S60_X7_mi_rubin.csv")]:
    d = rd(f"{R}/{f}"); NUM[nm] = d.to_dict("records") if d is not None else None
sc = rd(os.path.join(HERE, "HM_setcompare_LTCS_std3.csv")); NUM["set_compare"] = sc.to_dict("records") if sc is not None else None
NUM["best_set"] = open(os.path.join(HERE, "HM_best_LTCS_std3.txt")).read().splitlines()
NUM["hypothesis"] = T.to_dict("records")
NUM["forced"] = M3[M3.role == "forced"].to_dict("records")
rr = rd(f"{R}/S6_spatial_RR.csv")
if rr is not None:
    o = rr[rr.observed]; NUM["spatial"] = dict(n_hotspot=int((o["class"] == "hotspot").sum()), n_coldspot=int((o["class"] == "coldspot").sum()), RR_min=float(o.RR.min()), RR_max=float(o.RR.max()), RR_p90_p10=float(o.RR.quantile(.9) / o.RR.quantile(.1)),
                                  top10=o.sort_values("RR", ascending=False).head(10)[["fbd_region", "RR", "lo", "hi"]].to_dict("records"), bottom10=o.sort_values("RR").head(10)[["fbd_region", "RR", "lo", "hi"]].to_dict("records"))
NUM["sensitivity_versions"] = {}
for v in SENSV:
    f = os.path.join(HERE, f"results_{v}", "S3_main_IRR.csv")
    if os.path.exists(f): d = pd.read_csv(f); NUM["sensitivity_versions"][v] = d[d.role == "hypothesis"][["cov", "IRR", "lo", "hi", "cred"]].to_dict("records")
NUM["multiverse"] = S7.to_dict("records") if S7 is not None else None
NUM["alt_analyses_n"] = int(A.analysis.nunique()); NUM["alt_analyses_list"] = sorted(A.analysis.astype(str).unique().tolist())
NUM["main_vs_fullprec_maxdiff_IRR3"] = max(chk.values()) if chk else None
from statistics import NormalDist
ZA, ZB = NormalDist().inv_cdf(0.975), NormalDist().inv_cdf(0.80)
PW = {"method": "minimum detectable effect per 1 SD at 80% power, two-sided alpha 0.05: IRR_MDE = exp((z_0.975 + z_0.80) x posterior SD of beta), posterior SD from the main model (full-precision refit X0, NB BYM2+RW1+iid, 229 districts x 14 years); normal approximation to the posterior",
      "z_alpha": ZA, "z_beta": ZB, "n_cells": int(INFO["n_rows"]), "variables": {}}
src = FP if FP is not None else None
for v in M3["cov"]:
    sdv = float(src.sd[v]) if src is not None else (math.log(M3.set_index("cov").hi[v]) - math.log(M3.set_index("cov").lo[v])) / (2 * 1.959964)
    PW["variables"][v] = dict(role=M3.set_index("cov").role[v], posterior_sd=sdv, crI_halfwidth_log=1.959964 * sdv, MDE_IRR_up=math.exp((ZA + ZB) * sdv), MDE_IRR_down=math.exp(-(ZA + ZB) * sdv))
hs = [PW["variables"][v]["MDE_IRR_up"] for v in HYP]; PW["hypothesis_MDE_IRR_up_range"] = [min(hs), max(hs)]; PW["hypothesis_MDE_IRR_up_median"] = float(np.median(hs))
NUM["power"] = PW
PN = pd.read_csv(os.path.join(HERE, "panel_LTCS_std3.csv"), encoding="utf-8-sig")
COV = [c for c in M3["cov"]] + [c for c in NUM.get("negcontrol_vars", [])]
T1 = {}
for nm, (a_, b_) in {"2011-2015": (2011, 2015), "2016-2020": (2016, 2020), "2021-2024": (2021, 2024), "2011-2024": (2011, 2024)}.items():
    d = PN[PN.year.between(a_, b_)]; g = d.groupby("fbd_region").agg(c=("cases", "sum"), e=("E", "sum")); sir = g.c / g.e
    T1[nm] = dict(years=[a_, b_], n_cells=int(len(d)), cases_grade12=float(d.cases.sum()), pop65_person_years=float(d.offset.sum()), crude_rate_per_1000_65plus=float(d.cases.sum() / d.offset.sum() * 1e3),
                  cases_per_year_mean=float(d.groupby("year").cases.sum().mean()), pop65_per_year_mean=float(d.groupby("year").offset.sum().mean()),
                  SIR_district=dict(median=float(sir.median()), q1=float(sir.quantile(.25)), q3=float(sir.quantile(.75)), p10=float(sir.quantile(.1)), p90=float(sir.quantile(.9)), p90_p10=float(sir.quantile(.9) / sir.quantile(.1)), min=float(sir.min()), max=float(sir.max())),
                  covariates={v: dict(mean=float(d[v + "_raw"].mean()), sd=float(d[v + "_raw"].std(ddof=1)), label=(M3.set_index("cov").label_en[v] if "label_en" in M3 else v)) for v in M3["cov"] if v + "_raw" in d})
NUM["table1"] = dict(note="cell-level means/SDs of raw covariate values (district-years in window); SIR = district sum(cases)/sum(E) within window, E = sex-age indirectly standardised expected count (yearly sum E = sum cases)", windows=T1)
def clean(o):
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, list): return [clean(v) for v in o]
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating, float)): return None if (isinstance(o, float) and math.isnan(o)) or (isinstance(o, np.floating) and np.isnan(o)) else float(o)
    if isinstance(o, np.bool_): return bool(o)
    return o
json.dump(clean(NUM), open(os.path.join(HERE, "numbers_ltcs.json"), "w"), ensure_ascii=False, indent=1)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
DB = os.environ["ST_REF_DATASET"]
IN = [DB, os.path.join(SK, "02_패널", "EXT_canonical_pool_LTCS.csv"), os.path.join(SK, "02_패널", "LTCS_Estd3.csv"), os.path.join(HERE, "panel_LTCS_std3.csv"), os.path.join(HERE, "sets_LTCS_std3.json"), os.path.join(HERE, "config_LTCS_std3.json"), os.path.join(HERE, "hypothesis_table_LTCS.csv"), os.path.join(HERE, "graph_593.graph"), os.path.join(SK, "00_shp", "final.shp")] + sorted(glob.glob(os.path.join(HERE, "sets_LTCS_*.json")))
CODE = [os.path.join(SK, "01_스크립트", f) for f in ["80_build_disease_panel_v.py", "_st_config.py", "81_fit_inla_ltcs.R", "81R_ltcs_refit_lik.R", "85_window_maps.R", "86_yearly_maps_ehsa.py", "86b_ehsa_modified_mk_ltcs.py", "87_window_grid.py", "88g_sens_prep.py", "88g_sens_suite.R", "88h_ppc_varpart.R", "88k_struct_kh.R", "88X_ltcs_extra.R", "89g_sens_combine.py", "99L_build_Estd_ltcs.py", "rslot.sh"]] + [os.path.join(HERE, f) for f in ["build_pool_from_sqlite_ltcs.py", "build_pool_ltcs_additions.py", "make_pool_sets_ltcs.py", "run_ltcs_chain.sh", "90_ltcs_summary_numbers.py", "gate/code/11_ltcs_cell_gate.py", "gate/code/06_village_corrections.py", "gate/code/10_ltcs_cell_gate.py"]]
OUTF = sorted(glob.glob(os.path.join(HERE, "results_std3", "*.csv")) + [f for v in SENSV for f in glob.glob(os.path.join(HERE, f"results_{v}", "*.csv"))]) + [os.path.join(HERE, "numbers_ltcs.json")]
SIG = {"created": NUM["_meta"]["created"], "inputs": {os.path.relpath(p, SK) if p.startswith(SK) else p: sha(p) for p in IN if os.path.exists(p)},
       "code": {os.path.relpath(p, SK): sha(p) for p in CODE if os.path.exists(p)}, "results": {os.path.relpath(p, SK): sha(p) for p in OUTF if os.path.exists(p)}}
json.dump(SIG, open(os.path.join(HERE, "result_signature_LTCS.json"), "w"), ensure_ascii=False, indent=1)
print("numbers_ltcs.json sha256", sha(os.path.join(HERE, "numbers_ltcs.json")), "| 서명 결과 파일", len(SIG["results"]), "| 주모형 3자리 대조 최대차", NUM["main_vs_fullprec_maxdiff_IRR3"])
