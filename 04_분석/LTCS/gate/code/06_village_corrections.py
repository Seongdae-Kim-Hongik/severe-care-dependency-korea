# 06_village_corrections.py - Village sports facility counts: corrects row shifts in the source workbook (ministry register) found against neighbouring years, and sets
#  remaining collapsed or exploded cells to missing (interpolated within district later). Input: reference dataset (ST_REF_DATASET).
#  Output: 04_분석/LTCS/gate/tables/T07*.csv
import sqlite3, os, json, unicodedata
import pandas as pd, numpy as np
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = os.environ["ST_REF_DATASET"]
NFC = lambda s: unicodedata.normalize("NFC", s) if isinstance(s, str) else s
c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
R = pd.read_sql('select year,sido_std,fbd_region,sigungu_cd,population,"vil__계" as n,"vil__시설조성면적_m2" as area from region_year where year between 2011 and 2024', c)
R["fbd_region"] = R.fbd_region.map(NFC)
N = R.pivot(index="sigungu_cd", columns="year", values="n"); A = R.pivot(index="sigungu_cd", columns="year", values="area")
name = R.drop_duplicates("sigungu_cd").set_index("sigungu_cd").fbd_region
N0, A0 = N.copy(), A.copy(); LOG = []
def put(cd, y, n_new, a_new, step, evidence):
    LOG.append(dict(step=step, year=y, sigungu_cd=cd, fbd_region=name[cd], n_sqlite=N0.at[cd, y], area_sqlite=A0.at[cd, y], n_corrected=n_new, area_corrected=a_new, evidence=evidence))
    N.at[cd, y] = n_new; A.at[cd, y] = a_new
assert N0.at["37310", 2023] == 461
put("37310", 2023, np.nan, np.nan, "S1_군위2023_이중집계_결측", "2023년말 마을체육시설 xlsx: 대구 시트 군위군 행 = 달성군 행과 전 열 동일(48·238,090㎡·체력단련 317), 경북 시트 군위군 행 = 의성군 값(413·115,649㎡, 2022·2024 의성과 같은 면적). 빌더가 두 행을 합해 461")
chain = ["37320", "37330", "37340", "37350", "37360", "37370", "37380", "37390", "37400", "37410", "37420", "37430"]
cells = {cd: (N0.at[cd, 2023], A0.at[cd, 2023]) for cd in chain}
put("37320", 2023, 413.0, 115649.0, "S2_경북2023_한행밀림_재배치", "경북 시트 '군위군' 라벨 행의 값(413·115,649㎡) = 의성 2022·2024 면적과 같음")
for prev, cd in zip(chain[:-1], chain[1:]):
    put(cd, 2023, cells[prev][0], cells[prev][1], "S2_경북2023_한행밀림_재배치", f"경북 시트 '{name[prev]}' 라벨 행의 값이 {name[cd]}의 2022·2024 지문(계·면적)과 같음")
lab21 = "수원 고양 용인 성남 부천 안산 화성 남양주 안양 평택 의정부 파주 시흥 김포 광명 광주 군포 이천 오산 하남 양주 구리 안성 포천 의왕 여주 양평 동두천 과천 가평 연천".split()
lab22 = "수원 용인 고양 성남 화성 부천 남양주 안산 평택 안양 시흥 김포 파주 의정부 광주 광명 하남 군포 오산 양주 이천 구리 안성 의왕 포천 양평 여주 동두천 과천 가평 연천".split()
gcd = {name[cd].replace("경기도", "")[:-1]: cd for cd in name.index if name[cd].startswith("경기도")}
c21 = {lab: (N0.at[gcd[lab], 2021], A0.at[gcd[lab], 2021]) for lab in lab21}
for i, (la, lb) in enumerate(zip(lab21, lab22)):
    if la == lb: continue
    put(gcd[lb], 2021, c21[la][0], c21[la][1], "S3_경기2021_라벨순서_재배치", f"2021년말 xlsx 경기 시트 {i+1}번째 행(라벨 {la}) 값 = 2022년말 같은 위치 라벨 {lb}의 지문")
c20 = {k: (N0.at[gcd[k], 2020], A0.at[gcd[k], 2020]) for k in ["이천", "하남", "양주"]}
for tgt, src in [("하남", "이천"), ("양주", "하남"), ("이천", "양주")]:
    put(gcd[tgt], 2020, c20[src][0], c20[src][1], "S4_경기2020_세행회전_재배치", f"2020년말 xlsx 경기 시트 {src} 행 값 = {tgt}의 2019·2021(재배치) 지문")
assert (N0.at[gcd["가평"], 2019], N0.at[gcd["과천"], 2019]) == (25.0, 36.0)
put(gcd["가평"], 2019, 36.0, A0.at[gcd["가평"], 2019], "S5_경기2019_계열교차_되돌림", "가평 2018·2020=36, 과천 2018·2020=25; 2019만 서로 뒤바뀜")
put(gcd["과천"], 2019, 25.0, A0.at[gcd["과천"], 2019], "S5_경기2019_계열교차_되돌림", "같은 교차의 짝")
pop = R.pivot(index="sigungu_cd", columns="year", values="population")
yrs = list(range(2011, 2025))
for cd in N.index:
    s = N.loc[cd]; a = A.loc[cd]
    for i, y in enumerate(yrs):
        if pd.isna(s[y]) or y == 2015: continue
        obs = [t for t in yrs if not pd.isna(s[t]) and t != y]
        prv = [t for t in obs if t < y]; nxt = [t for t in obs if t > y]
        if not prv or not nxt: continue
        p, q = prv[-1], nxt[0]
        far = lambda u, v: abs(u - v) > 0.05 * max(abs(v), 1)
        sido = R.loc[R.sigungu_cd == cd, "sido_std"].iat[0]
        other = R[(R.sido_std == sido) & (R.sigungu_cd != cd) & (R.year.isin([y - 1, y + 1]))]
        foreign = bool(((N.loc[other.sigungu_cd.values, :].values == None).any()) if False else False)
        foreign = any((N.at[o, t] == s[y]) and (A.at[o, t] == a[y]) for o, t in zip(other.sigungu_cd, other.year) if t in N.columns)
        cnt_far = abs(s[y] - s[p]) >= 0.25 * max(s[p], 1) and abs(s[y] - s[q]) >= 0.25 * max(s[q], 1)
        rowmis = (not pd.isna(a[p])) and far(a[y], a[p]) and far(a[y], a[q]) and (a[p] == a[q]) and (cnt_far or foreign)
        rp, rq = s[y] / pop.at[cd, y], s[p] / pop.at[cd, p], 
        vp, vq, vy = s[p] / pop.at[cd, p], s[q] / pop.at[cd, q], s[y] / pop.at[cd, y]
        spike1 = vp > 0 and vq > 0 and ((vy >= 3 * vp and vy >= 3 * vq) or (vy <= vp / 3 and vy <= vq / 3))
        spike2 = False
        if len(nxt) >= 2:
            q2 = nxt[1]; vq2 = s[q2] / pop.at[cd, q2]
            spike2 = vp > 0 and vq2 > 0 and all((v <= vp / 3 and v <= vq2 / 3) or (v >= 3 * vp and v >= 3 * vq2) for v in (vy, vq))
        if rowmis or spike1 or spike2:
            why = ";".join(k for k, f in [("행내용불일치(면적 앞뒤와 5% 넘게 다름·앞뒤끼리 같음 + 계 25% 이상 차이 또는 다른 지역 지문과 일치)", rowmis), ("단년 3배 일탈", spike1), ("2년 3배 일탈 후 복귀", spike2)] if f)
            put(cd, y, np.nan, np.nan, "S6_잔여_결측처리(규칙)", f"{why}; 앞 {p}={s[p]:g}, 뒤 {q}={s[q]:g}")
            if spike2: put(cd, q, np.nan, np.nan, "S6_잔여_결측처리(규칙)", f"2년 3배 일탈 후 복귀의 둘째 해; 앞 {p}={s[p]:g}, 복귀 {nxt[1]}={s[nxt[1]]:g}")
L = pd.DataFrame(LOG).drop_duplicates(["year", "sigungu_cd", "step"])
L.to_csv(os.path.join(B, "tables/T07_마을체육시설_정정기록.csv"), index=False, encoding="utf-8-sig")
out = N.stack(future_stack=True).rename("vil_n_corrected").reset_index().rename(columns={"level_1": "year"})
out = out.merge(N0.stack(future_stack=True).rename("vil_n_sqlite").reset_index().rename(columns={"level_1": "year"}), on=["sigungu_cd", "year"])
out.to_csv(os.path.join(B, "tables/T07b_village_sports_facilities_corrected.csv"), index=False, encoding="utf-8-sig")
pd.set_option("display.width", 260); print(L.drop(columns=["evidence"]).to_string()); print(L.groupby("step").size())
