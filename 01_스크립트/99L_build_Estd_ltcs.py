#!/usr/bin/env python3
# 99L_build_Estd_ltcs.py - Outcome counts and expected counts. cases = residents aged 65+ certified at long-term care grade 1 or 2 (KOSIS DT_35006_N009), summed from
#  general-district rows to the 229 districts (code-based linkage; double counting of city and sub-district rows prevented; district sum = national total asserted).
#  Expected counts by indirect standardisation in 2 sexes x 5 age bands (65-69 ... 85+) with national rates (DT_35006_N036; 2013 rates used for 2011-2012)
#  and resident population by sex and age (ST_FBDZIP); sum of E = sum of cases by year; no added constant. Alternative E: sex-only, age-only, crude.
#  Missing 85+ age cells are set to 0 only after the district total is asserted to equal the sum of the age bands. ST_EFIX2015=1 also writes E with fixed 2015 national rates (ST_OUT_EFIX).
#  Output: 02_패널/LTCS_Estd3.csv; checks in 04_분석/LTCS/gate/E/
import os, sys, json, glob, sqlite3, importlib.util, unicodedata
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); PAN = os.path.join(ROOT, "02_패널"); os.makedirs(PAN, exist_ok=True); SEL = os.environ.get("ST_OUT_SEL", os.path.join(ROOT, "04_분석", "LTCS", "gate", "E")); os.makedirs(SEL, exist_ok=True)
RAW = os.environ["ST_KOSIS_RAW"]
DB = os.environ["ST_REF_DATASET"]
def _load(nm, f):
    s = importlib.util.spec_from_file_location(nm, os.path.join(HERE, f)); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
FBDZIP = os.environ.get("ST_FBDZIP"); cw = _load("cw", "40_region_crosswalk.py"); m50 = _load("m50", "50_build_master_panel.py")
NFC = lambda s: unicodedata.normalize("NFC", str(s))
assert os.environ.get("ST_FBDZIP"), "ST_FBDZIP(창원 2023 교정 인구 폴더)를 지정"
YEARS = list(range(2011, 2025)); SEXES = ["남자", "여자"]
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
base = pd.read_sql("select year, key, fbd_region from disease_panel where tag='PERIO'", con); base["fbd_region"] = base.fbd_region.map(NFC); base["key"] = base.key.map(NFC)
PKEYS = set(base.key); assert len(base) == 3206 and len(PKEYS) == 229
K2R = dict(zip(base.key, base.fbd_region))
BANDS18 = [f"{lo}-{lo+4}" for lo in range(0, 85, 5)] + ["85+"]; OLD = ["65-69", "70-74", "75-79", "80-84", "85+"]
SIDO_HDR = {k for k in cw.SIDO_STD if len(k) >= 3 and (k.endswith("시") or k.endswith("도"))}
FIVE = [f"{lo}-{lo+4}세" for lo in range(0, 85, 5)]; OLD85 = ["85-89세", "90-94세", "95-99세", "100세 이상"]
EXTRA = {("충청북도", "청원군"): "충청북도|청주시"}; rows = []; NMISS85 = []
for y in YEARS:
    d = pd.read_csv(os.path.join(FBDZIP, f"{y}.csv"), encoding="cp949", dtype=str)
    d["행정구역별"] = d["행정구역별"].map(NFC).str.strip(); d["연령별"] = d["연령별"].map(NFC).str.strip()
    for s in SEXES: d[s] = pd.to_numeric(d[s].str.replace(",", ""), errors="coerce")
    sido = None; city = None; blocks = (d["행정구역별"] != d["행정구역별"].shift()).cumsum()
    for _, g in d.groupby(blocks, sort=False):
        nm = g["행정구역별"].iat[0]
        if nm == "전국": continue
        if nm in SIDO_HDR:
            sido = nm; city = None
            if nm != "세종특별자치시": continue
            g = g.drop_duplicates("연령별")
        sg = nm.split("-")[-1]
        if sg.endswith("시"): city = sg
        elif sg.endswith("구") and city is not None and cw.sido_std(sido) not in ("서울시", "부산시", "대구시", "인천시", "광주시", "대전시", "울산시"): continue
        else: city = None
        key = EXTRA.get((cw.sido_std(sido), sg), cw.make_key(sido, sg))
        if key not in PKEYS and key.endswith("군") and key[:-1] + "시" in PKEYS: key = key[:-1] + "시"
        for s in SEXES:
            a = g.set_index("연령별")[s]
            if pd.isna(a.get("계", np.nan)): continue
            f = a.reindex(FIVE); o = a.reindex(OLD85); assert f.notna().all(), (y, nm, s)
            NMISS85.append(int(o.isna().sum())); assert abs(a["계"] - f.sum() - o.fillna(0).sum()) < 0.5, ("missing 85+ cells are not true zeros", y, nm, s, a["계"], f.sum(), o.sum())
            o = o.where(o.notna(), 0.0)
            rows.append({"year": y, "sex": s, "key": NFC(key), **{b: f[b + "세"] for b in BANDS18[:-1]}, "85+": o.sum()})
print(f"[99L] missing 85+ age cells: {sum(NMISS85)} in {len(NMISS85)} blocks; each block total equals the sum of its age bands (true zeros)")
POP = pd.DataFrame(rows).groupby(["year", "key", "sex"], as_index=False).sum(numeric_only=True)
assert set(POP.key) >= PKEYS, sorted(PKEYS - set(POP.key))
cwn = POP[(POP.key == NFC("경상남도|창원시")) & (POP.year == 2023)][BANDS18].sum().sum()
print(f"[99S] 창원 2023 인구 {cwn:,.0f} (정답 1,009,038)"); assert abs(cwn - 1009038) < 1
NATPOP = POP[POP.key.isin(PKEYS)].groupby(["year", "sex"])[BANDS18].sum()

CITY_SIDO = {"수원시": "경기도", "성남시": "경기도", "안양시": "경기도", "안산시": "경기도", "고양시": "경기도", "용인시": "경기도", "부천시": "경기도", "화성시": "경기도",
             "청주시": "충청북도", "천안시": "충청남도", "전주시": "전라북도", "포항시": "경상북도", "창원시": "경상남도"}
def read_tbl(tbl, years):
    out = []
    for y in years:
        d = json.load(open(os.path.join(RAW, tbl, f"{y}.json"))); D = pd.DataFrame(d); D["year"] = y
        D["v"] = pd.to_numeric(D.DT.astype(str).str.replace(",", ""), errors="coerce"); out.append(D)
    return pd.concat(out, ignore_index=True)
TRACE = []
def to_sgg(D, tbl, sexcol, catcol, cats, extra=None):
    D = D[D[catcol].isin(cats)].copy()
    if extra is not None: D = D[extra(D)]
    D["C1_NM"] = D.C1_NM.map(NFC).str.strip(); D["C1"] = D.C1.astype(str)
    if tbl.startswith("DT_35006"):
        nat = D.C1 == "15350AC400Z1"; hdr = D.C1.str.len().eq(12) & D.C1.str.startswith("15350AC400") & ~nat
        par = set(zip(D.year, D.C1.str[:12].where(D.C1.str.len().gt(12))))
        leaf = hdr & ~pd.Series([(y, c) in par for y, c in zip(D.year, D.C1)], index=D.index); hdr = hdr & ~leaf
        H = D[hdr].drop_duplicates(["year", "C1"]); hmap = {(y, c): n for y, c, n in zip(H.year, H.C1, H.C1_NM)}
        D["sido_raw"] = [hmap.get((y, c[:12])) if len(c) > 12 else (n if lf else None) for y, c, n, lf in zip(D.year, D.C1, D.C1_NM, leaf)]
    else:
        nat = D.C1 == "00"; hdr = D.C1.str.len().eq(2) & D.C1.str.isdigit() & ~nat & D.C1_NM.map(lambda n: n in cw.SIDO_STD or n.endswith("도") or n.endswith("특별시") or n.endswith("광역시") or n.endswith("자치시"))
        par = set(zip(D.year, D.C1.str[:2].where(D.C1.str.len().eq(5))))
        leaf = hdr & ~pd.Series([(y, c) in par for y, c in zip(D.year, D.C1)], index=D.index); hdr = hdr & ~leaf
        H = D[hdr].drop_duplicates(["year", "C1"]); hmap = {(y, c): n for y, c, n in zip(H.year, H.C1, H.C1_NM)}
        D["sido_raw"] = [hmap.get((y, c[:2])) if len(c) == 5 else (n if lf else None) for y, c, n, lf in zip(D.year, D.C1, D.C1_NM, leaf)]
    kid = ~nat & ~hdr
    orphan = kid & D.sido_raw.isna()
    D.loc[orphan, "sido_raw"] = D.loc[orphan, "C1_NM"].map(lambda n: next((v for c, v in CITY_SIDO.items() if n.replace(" ", "").startswith(c)), None))
    assert D.loc[kid, "sido_raw"].notna().all(), D.loc[kid & D.sido_raw.isna(), ["year", "C1", "C1_NM"]].drop_duplicates().values.tolist()
    K = D[kid].copy(); pr = K[["sido_raw", "C1_NM"]].drop_duplicates()
    mp = {(a, b): m50.map_region(a, b, PKEYS) for a, b in zip(pr.sido_raw, pr.C1_NM)}
    K["key"] = [mp[(a, b)]["key"] for a, b in zip(K.sido_raw, K.C1_NM)]; K["role"] = [mp[(a, b)]["role"] for a, b in zip(K.sido_raw, K.C1_NM)]
    K["flag"] = [mp[(a, b)]["flag"] for a, b in zip(K.sido_raw, K.C1_NM)]
    un = K[K.key.isna()][["year", "sido_raw", "C1_NM"]].drop_duplicates()
    assert un.empty, un.values.tolist()
    K["key"] = K.key.map(NFC)
    hasprim = set(map(tuple, K[K.role == "primary"][["key", "year"]].drop_duplicates().values))
    K["used"] = ~((K.role == "sub").values & np.array([(k, y) in hasprim for k, y in zip(K.key, K.year)]))
    reg = K.drop_duplicates(["year", "C1", "C1_NM"])
    for r in reg.itertuples():
        TRACE.append(dict(표=tbl, year=r.year, 원코드=r.C1, 원시도=r.sido_raw, 원지역=r.C1_NM, key=r.key, fbd_region=K2R.get(r.key), 역할=r.role, 처리=r.flag or ("일대일" if r.role == "primary" else r.role), 사용=bool(r.used)))
    S = K[K.used].groupby(["year", sexcol, "key"]).v.sum(min_count=1)
    N = D[nat].groupby(["year", sexcol]).v.sum(min_count=1)
    HS = D[hdr | leaf].groupby(["year", sexcol]).v.sum(min_count=1)
    return S, N, HS, int((~nat & ~hdr).groupby(D.year).sum().max()), int(K[K.used].groupby("year").key.nunique().min())

SPEC = {
    "LTCS": dict(tbl="DT_35006_N009", sexcol="C2_NM", catcol="C4_NM", cats=["1등급", "2등급"], sexmap={"남자": "남자", "여자": "여자", "합계": "계"}),
    "LTCC": dict(tbl="DT_35006_N009", sexcol="C2_NM", catcol="C4_NM", cats=["1등급", "2등급", "3등급", "4등급", "5등급", "인지지원등급"], sexmap={"남자": "남자", "여자": "여자", "합계": "계"}),
    "LTCA": dict(tbl="DT_35006_N005", sexcol="C2_NM", catcol="C4_NM", cats=["65세이상"], sexmap={"남자": "남자", "여자": "여자", "합계": "계"}),
}
CNT = {}; NATCHK = []; RAWN = {}
for T, sp in SPEC.items():
    D = read_tbl(sp["tbl"], YEARS)
    S, N, HS, nraw, nkey = to_sgg(D, sp["tbl"], sp["sexcol"], sp["catcol"], sp["cats"])
    S = S.rename(index=sp["sexmap"], level=1); N = N.rename(index=sp["sexmap"], level=1); HS = HS.rename(index=sp["sexmap"], level=1)
    tot = S.groupby(level=[0, 1]).sum()
    for (y, s), v in tot.items():
        NATCHK.append(dict(TAG=T, 표=sp["tbl"], year=y, sex=s, 시군구합=v, 전국=N.get((y, s), np.nan), 비=v / N.get((y, s), np.nan), 시도헤더합=HS.get((y, s), np.nan), 원표지역행=nraw, 합산키=nkey))
    CNT[T] = S
D6 = read_tbl("DT_11761_N006", range(2011, 2019)); D9 = read_tbl("DT_11761_N009", range(2019, 2025))
S6, N6, H6, r6, k6 = to_sgg(D6, "DT_11761_N006", "C3_NM", "C2_NM", ["1급", "2급", "3급"])
S9, N9, H9, r9, k9 = to_sgg(D9, "DT_11761_N009", "C3_NM", "C2_NM", ["심한장애"])
for S, N, HS, tb, nr, nk in ((S6, N6, H6, "DT_11761_N006", r6, k6), (S9, N9, H9, "DT_11761_N009", r9, k9)):
    for (y, s), v in S.groupby(level=[0, 1]).sum().items():
        NATCHK.append(dict(TAG="DISV", 표=tb, year=y, sex=s, 시군구합=v, 전국=N.get((y, s), np.nan), 비=v / N.get((y, s), np.nan), 시도헤더합=HS.get((y, s), np.nan), 원표지역행=nr, 합산키=nk))
CNT["DISV"] = pd.concat([S6, S9])
NC = pd.DataFrame(NATCHK)
NC["헤더_전국비"] = NC.시도헤더합 / NC.전국
NC["판정"] = np.where(NC.비.sub(1).abs() <= 0.001, "일치", np.where(NC.헤더_전국비.sub(1).abs() <= 0.001, "원표 내부 불일치(시군구 행≠시도 헤더)", "불일치"))
NC.to_csv(os.path.join(SEL, "E_전국합대조_LTCS.csv"), index=False, encoding="utf-8-sig")
print(NC.groupby(["TAG", "표"]).agg(비_min=("비", "min"), 비_max=("비", "max"), 지역행=("원표지역행", "max"), 키=("합산키", "min")).to_string())
print(NC[NC.판정 != "일치"].to_string()); assert (NC.판정 != "불일치").all()
TR = pd.DataFrame(TRACE); TR.to_csv(os.path.join(SEL, "E_지역처리표_LTCS.csv"), index=False, encoding="utf-8-sig")
print("[99S] 처리표", len(TR), "행; 처리 종류", TR.groupby(["표", "처리"]).size().to_dict())

def ltc_nat(tbl, yrs, catcol, cats, agecol, sexcol, sexes={"남자": "남자", "여자": "여자"}):
    D = read_tbl(tbl, yrs); D = D[D[catcol].isin(cats) & D[sexcol].isin(list(sexes))]
    AG = {"65~69세": "65-69", "70~74세": "70-74", "75~79세": "75-79", "80~84세": "80-84", "85세이상": "85+"}
    D["band"] = D[agecol].map(AG); D = D[D.band.notna()]
    return D.assign(sex=D[sexcol].map(sexes)).groupby(["year", "sex", "band"]).v.sum().unstack()[OLD]
NAT = {}
n36y = sorted(int(os.path.basename(f)[:4]) for f in glob.glob(os.path.join(RAW, "DT_35006_N036", "*.json")))
n37y = sorted(int(os.path.basename(f)[:4]) for f in glob.glob(os.path.join(RAW, "DT_35006_N037", "*.json")))
NAT["LTCS"] = ltc_nat("DT_35006_N036", [y for y in n36y if y <= 2024], "C5_NM", ["1등급", "2등급"], "C2_NM", "C3_NM")
NAT["LTCC"] = ltc_nat("DT_35006_N036", [y for y in n36y if y <= 2024], "C5_NM", ["1등급", "2등급", "3등급", "4등급", "5등급", "인지지원등급"], "C2_NM", "C3_NM")
D37 = read_tbl("DT_35006_N037", [y for y in n37y if y <= 2024]); print("[99S] N037 열", {c: sorted(D37[c].unique())[:8] for c in D37.columns if c.endswith("_NM") and c.startswith("C")})
c_sex = next(c for c in ["C1_NM", "C2_NM", "C3_NM", "C4_NM"] if c in D37 and set(D37[c]) >= {"남자", "여자"}); c_age = next(c for c in ["C1_NM", "C2_NM", "C3_NM", "C4_NM"] if c in D37 and "85세이상" in set(D37[c]))
NAT["LTCA"] = ltc_nat("DT_35006_N037", [y for y in n37y if y <= 2024], c_sex, ["남자", "여자"], c_age, c_sex)
CMP = []
for T in ["LTCS", "LTCC", "LTCA"]:
    A = NAT[T].sum(1); B = CNT[T].groupby(level=[0, 1]).sum()
    for (y, s), v in A.items(): CMP.append(dict(TAG=T, year=y, sex=s, 전국연령표=v, 시군구표=B.get((y, s), np.nan), 비=v / B.get((y, s), np.nan)))
CMP = pd.DataFrame(CMP); print("[99S] 전국 연령표 / 시군구 표 합 비:", CMP.groupby("TAG").비.agg(["min", "max"]).round(4).to_dict("index"))
def dis_band(a):
    a = str(a).replace(" ", "")
    if a in ("85-89세", "90-94세", "95-99세", "100세이상"): return "85+"
    if a.endswith("세") and "-" in a:
        lo, hi = a[:-1].split("-")
        if int(hi) - int(lo) == 4: return f"{lo}-{hi}"
    return None
d2 = read_tbl("DT_11761_N002", range(2011, 2019)); d2 = d2[d2.C2_NM.isin(["1급", "2급", "3급"]) & d2.C3_NM.isin(SEXES)].assign(sex=lambda x: x.C3_NM)
d7 = read_tbl("DT_11761_N007", range(2019, 2025)); c7s = next(c for c in ["C2_NM", "C3_NM"] if set(d7[c]) >= {"남자", "여자"}); c7g = next(c for c in ["C2_NM", "C3_NM"] if "심한장애" in set(d7[c]))
d7 = d7[(d7[c7g] == "심한장애") & d7[c7s].isin(SEXES)].assign(sex=lambda x: x[c7s])
DN = pd.concat([d2, d7]); DN["band"] = DN.C1_NM.map(dis_band); DN = DN[DN.band.notna()]
NAT["DISV"] = DN.groupby(["year", "sex", "band"]).v.sum().unstack()[BANDS18]
chk = NAT["DISV"].sum(1) / CNT["DISV"].groupby(level=[0, 1]).sum().reindex(NAT["DISV"].index)
print("[99S] DISV 전국 연령표 / 시군구 합:", round(chk.min(), 4), round(chk.max(), 4))
RATES = {}; RLOG = []
for T in ["LTCS", "LTCC", "LTCA", "DISV"]:
    bands = BANDS18 if T == "DISV" else OLD; R = NAT[T] / NATPOP[bands].reindex(NAT[T].index)
    full = []
    for y in YEARS:
        src = y if y in R.index.get_level_values(0) else min(R.index.get_level_values(0), key=lambda z: abs(z - y))
        x = R.xs(src, level="year").copy(); x["year"] = y; x["율_연도"] = src; full.append(x.reset_index())
        RLOG.append(dict(TAG=T, year=y, 율_연도=src, **{f"{s}_{b}": round(float(R.loc[(src, s), b]), 6) for s in SEXES for b in (bands if T != "DISV" else ["60-64", "75-79", "85+"])}))
    RATES[T] = pd.concat(full).set_index(["year", "sex"])
pd.DataFrame(RLOG).to_csv(os.path.join(SEL, "E_전국율_LTCS.csv"), index=False, encoding="utf-8-sig")
T = "LTCS"; bands = OLD
tot = CNT[T].xs("계", level=1) if "계" in CNT[T].index.get_level_values(1) else CNT[T].groupby(level=[0, 2]).sum()
P = base.copy(); P["cases"] = tot.reindex(pd.MultiIndex.from_arrays([P.year, P.key])).values; assert P.cases.notna().all()
NATC = NAT[T]
E = np.zeros(len(P)); Esex = np.zeros(len(P)); Eage = np.zeros(len(P)); off = np.zeros(len(P)); p85 = np.zeros(len(P)); pop80 = np.zeros(len(P))
agerate = {}
for y in YEARS:
    src = y if y in NATC.index.get_level_values(0) else min(NATC.index.get_level_values(0), key=lambda z: abs(z - y))
    agerate[y] = NATC.xs(src, level="year").sum(0) / NATPOP.xs(src, level="year")[bands].sum(0)
sexcnt = CNT[T][CNT[T].index.get_level_values(1).isin(SEXES)].groupby(level=[0, 1]).sum()
X0 = {}
for s in SEXES:
    X = P[["year", "key"]].merge(POP[POP.sex == s], on=["year", "key"], how="left", validate="one_to_one"); assert X[bands].notna().all().all(); X0[s] = X
    rr = RATES[T].xs(s, level="sex").loc[X.year.values, bands].values
    E += (X[bands].values * rr).sum(1); off += X[bands].values.sum(1); p85 += X["85+"].values; pop80 += X[["80-84", "85+"]].values.sum(1)
    natpop65 = NATPOP.xs(s, level="sex")[bands].sum(1)
    sr = np.array([sexcnt.loc[(y, s)] / natpop65.loc[y] for y in X.year.values]); Esex += X[bands].values.sum(1) * sr
    Eage += (X[bands].values * np.vstack([agerate[y].values for y in X.year.values])).sum(1)
totpop = sum(X0[s][BANDS18].values.sum(1) for s in SEXES)
rs = lambda v: v * P.year.map(P.groupby("year").cases.sum()).values / pd.Series(v).groupby(P.year.values).transform("sum").values
P["offset"] = off; P["E_raw"] = E; P["E_std"] = rs(E); P["E_sex"] = rs(Esex); P["E_std_age"] = rs(Eage)
P["E_crude"] = P.offset * P.year.map(P.groupby("year").cases.sum() / P.groupby("year").offset.sum())
P["pop85_share_of65"] = p85 / off * 100; P["pop85_share"] = p85 / totpop * 100; P["pop80_share_fbd5"] = pop80 / totpop * 100; P["pop65_share"] = off / totpop * 100; P["pop_total_fbd5"] = totpop
for c in ["E_std", "E_sex", "E_std_age", "E_crude"]:
    assert (P.groupby("year")[c].sum() - P.groupby("year").cases.sum()).abs().max() < 1e-6 * P.cases.sum() and (P[c] > 0).all(), c
print(f"[99L] LTCS 행 {len(P)}, 합 {int(P.cases.sum()):,}, 0칸 {(P.cases == 0).sum()}, E_sex/E_std {np.round(np.percentile(P.E_sex / P.E_std, [0, 50, 100]), 3)}, E_age/E_std {np.round(np.percentile(P.E_std_age / P.E_std, [0, 50, 100]), 3)}")
P[["year", "fbd_region", "offset", "cases", "E_std", "E_sex", "E_std_age", "E_crude", "pop65_share", "pop85_share", "pop85_share_of65", "pop80_share_fbd5", "pop_total_fbd5"]].to_csv(os.environ.get("ST_OUT_PANEL", os.path.join(PAN, "LTCS_Estd3.csv")), index=False, encoding="utf-8-sig")
if os.environ.get("ST_EFIX2015") == "1":
    Efix = np.zeros(len(P))
    for s in SEXES:
        X = X0[s]; r15 = RATES[T].xs(s, level="sex").loc[2015, bands].values.astype(float)
        Efix += (X[bands].values * r15).sum(1)
    P["E_fix2015"] = Efix; assert (P.E_fix2015 > 0).all()
    print(f"[99L] E_fix2015 합/cases 합 연도별: {(P.groupby('year').E_fix2015.sum() / P.groupby('year').cases.sum()).round(4).to_dict()}")
    P[["year", "fbd_region", "cases", "E_std", "E_fix2015"]].to_csv(os.environ["ST_OUT_EFIX"], index=False, encoding="utf-8-sig")
