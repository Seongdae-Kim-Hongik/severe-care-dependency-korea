#!/bin/bash
# run_ltcs_chain.sh - Full model chain, part 1: VIF gate, 81 (principal model and multiverse blocks), 81R, sensitivity versions (80 + 81), 85, 86, 86b, 87, 88g prep,
#  88k, 88X, multiverse merge. Requires ST_GRAPH. R jobs run through rslot.sh with one INLA thread; return codes go to run_ltcs_rc.log.
SK="${ST_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"; SC="$SK/01_스크립트"; L="$SK/04_분석/LTCS"; T=LTCS; RS="$SC/rslot.sh"; GR="${ST_GRAPH:?set ST_GRAPH to the principal adjacency graph file}"
TL="$L/run_ltcs_times.log"; RCLOG="$L/run_ltcs_rc.log"; LG="$L/logs"; mkdir -p "$LG"
rcl(){ echo "$(date '+%F %T') rc=$2 $1" >> "$RCLOG"; }; tm(){ echo "- $1 $(date '+%F %T')" >> "$TL"; }
wl(){ while :; do l=$(sysctl -n vm.loadavg | awk '{print $2}'); awk -v l="$l" 'BEGIN{exit !(l<=12)}' && break; echo "[wait] load1 $l"; sleep 60; done; }
PIDS=(); NAMES=(); bgadd(){ PIDS+=("$1"); NAMES+=("$2"); }
waitall(){ local i rc; for i in "${!PIDS[@]}"; do wait "${PIDS[$i]}"; rc=$?; rcl "${NAMES[$i]}" "$rc"; done; PIDS=(); NAMES=(); }
chmod +x "$RS" "$0"; bash -n "$0" || { rcl precheck_bash 1; exit 1; }
for f in 81_fit_inla_ltcs.R 81R_ltcs_refit_lik.R 85_window_maps.R 88g_sens_suite.R 88h_ppc_varpart.R 88k_struct_kh.R 88k_struct_kh_anem.R 88X_ltcs_extra.R; do Rscript -e "invisible(parse('$SC/$f'))" || { rcl "precheck_parse_$f" 1; exit 1; }; done
rcl precheck 0; tm "start (prechecks passed)"
export TAG=$T ST_VER=std3 ST_THREADS=1 ST_GRAPH="$GR" SK_RSLOTS=3 NCO=LTCA,LTCC
cd "$SC"
python3 - "$L/panel_INFO_LTCS_std3.json" <<'PYEOF' || { rcl vif_gate 1; echo "VIF GATE FAILED"; exit 1; }
import json,sys; J=json.load(open(sys.argv[1])); o=J.get("vif_over_threshold",[]); print("vif_over_threshold",o,"vif_dropped",J.get("vif_dropped")); sys.exit(1 if o else 0)
PYEOF
rcl vif_gate 0
wl; ST_MV_DEFER=1 ST_MV_PARTS=4 $RS Rscript 81_fit_inla_ltcs.R $T > "$LG/81_main.log" 2>&1 & P81=$!
while [ ! -f "$L/results_std3/S7_multiverse_specs.csv" ]; do sleep 5; kill -0 $P81 2>/dev/null || break; done
for R in 1-50 51-100 101-150 151-200; do wl; ST_MV_RANGE=$R $RS Rscript 81_fit_inla_ltcs.R $T > "$LG/81_mv_$R.log" 2>&1 & bgadd $! "81_mv_$R"; sleep 3; done
MVP=("${PIDS[@]}"); MVN=("${NAMES[@]}"); PIDS=(); NAMES=()
wl; $RS Rscript 81R_ltcs_refit_lik.R $T > "$LG/81R.log" 2>&1 & bgadd $! 81R_refit_lik
for V in s577 s2 sc s8 y1524 x2021 y1324 sexE supply; do (cd "$SK" && ST_VER=$V python3 01_스크립트/80_build_disease_panel_v.py $T) > "$LG/80_$V.log" 2>&1; rcl "80_$V" $?; done
for V in s577 s2 sc s8 y1524 x2021 y1324 sexE supply; do wl; ( if [ $V = s577 ]; then unset ST_GRAPH; fi; ST_VER=$V ST_MV_N=2 ST_SENS_LITE=1 $RS Rscript 81_fit_inla_ltcs.R $T ) > "$LG/81_$V.log" 2>&1 & bgadd $! "81_$V"; sleep 2; done
wait $P81; rc=$?; rcl 81_main_S1-S6_S9 $rc; tm "81 principal model (S1-S6, S9) done rc=$rc"; [ $rc -eq 0 ] || { echo "81 FAILED"; exit 1; }
(cd "$SK" && python3 01_스크립트/87_window_grid.py $T std3) > "$LG/87.log" 2>&1; rcl 87_window_grid $?
wl; $RS Rscript 85_window_maps.R $T std3 > "$LG/85.log" 2>&1; rcl 85_window_maps $?
(cd "$SK" && python3 01_스크립트/86_yearly_maps_ehsa.py $T std3) > "$LG/86.log" 2>&1; rcl 86_yearly_maps_ehsa $?
(cd "$SK" && python3 01_스크립트/86b_ehsa_modified_mk_ltcs.py) > "$LG/86b.log" 2>&1; rcl 86b_ehsa_mk $?
tm "85 86 86b 87 done"
(cd "$SK" && python3 01_스크립트/88g_sens_prep.py $T std3 LTCA,LTCC) > "$LG/88gprep.log" 2>&1; rc=$?; rcl 88g_sens_prep $rc
for G in 1 2 3 4 5 6; do wl; GROUP=$G $RS Rscript 88g_sens_suite.R > "$LG/88g_$G.log" 2>&1 & bgadd $! "88g_$G"; sleep 2; done
wl; $RS Rscript 88h_ppc_varpart.R > "$LG/88h.log" 2>&1 & bgadd $! 88h
for G in nb struct kh review; do wl; $RS Rscript 88k_struct_kh.R $G > "$LG/88k_$G.log" 2>&1 & bgadd $! "88k_$G"; sleep 2; done
for X in A B C; do wl; EXG=$X $RS Rscript 88X_ltcs_extra.R > "$LG/88X_$X.log" 2>&1 & bgadd $! "88X_$X"; sleep 2; done
for i in "${!MVP[@]}"; do wait "${MVP[$i]}"; rcl "${MVN[$i]}" $?; done; tm "multiverse blocks done"
ST_MV_PARTS=4 $RS Rscript 81_fit_inla_ltcs.R $T > "$LG/81_merge.log" 2>&1; rc=$?; rcl 81_S7_merge $rc; tm "81 S7 merge done rc=$rc"
waitall
if grep -q 'KH_typeII FAILED\|KH_typeIV FAILED' "$LG/88k_kh.log"; then KH_INIT=0 $RS Rscript 88k_struct_kh_anem.R kh > "$LG/88k_kh_init0.log" 2>&1; rcl 88k_kh_init0 $?; fi
tm "sensitivity suite and versions done"
(cd "$SK" && TAG=$T ST_VER=std3 python3 01_스크립트/89g_sens_combine.py $T std3) > "$LG/89g.log" 2>&1; rcl 89g_sens_combine $?
(cd "$L" && python3 90_ltcs_summary_numbers.py) > "$LG/90.log" 2>&1; rcl 90_summary_numbers $?
tm "part 1 done"; echo ALL_DONE
