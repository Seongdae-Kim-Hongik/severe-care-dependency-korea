#!/bin/bash
# run_ltcs_chain_part2.sh - Full model chain, part 2 (after the multiverse merge has written S7_multiverse_summary.csv): 88g groups 1-6, 88h, 88k review, 89g and 90.
SK="${ST_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"; SC="$SK/01_스크립트"; L="$SK/04_분석/LTCS"; T=LTCS; RS="$SC/rslot.sh"; GR="${ST_GRAPH:?set ST_GRAPH to the principal adjacency graph file}"
TL="$L/run_ltcs_times.log"; RCLOG="$L/run_ltcs_rc.log"; LG="$L/logs"
rcl(){ echo "$(date '+%F %T') rc=$2 $1" >> "$RCLOG"; }; tm(){ echo "- $1 $(date '+%F %T')" >> "$TL"; }
wl(){ while :; do l=$(sysctl -n vm.loadavg | awk '{print $2}'); awk -v l="$l" 'BEGIN{exit !(l<14)}' && break; sleep 60; done; }
PIDS=(); NAMES=(); bgadd(){ PIDS+=("$1"); NAMES+=("$2"); }
waitall(){ local i rc; for i in "${!PIDS[@]}"; do wait "${PIDS[$i]}"; rc=$?; rcl "${NAMES[$i]}" "$rc"; done; PIDS=(); NAMES=(); }
bash -n "$0" || exit 1; test -s "$L/results_std3/S7_multiverse_summary.csv" || { rcl tail_precheck 1; exit 1; }
export TAG=$T ST_VER=std3 ST_THREADS=1 ST_GRAPH="$GR" SK_RSLOTS=4 NCO=LTCA,LTCC; cd "$SC"; tm "part 2 start (88g, 88h, 88k review)"
for G in 1 2 3 4 5 6; do wl; GROUP=$G $RS Rscript 88g_sens_suite.R > "$LG/88g_${G}_p2.log" 2>&1 & bgadd $! "88g_${G}_p2"; sleep 2; done
wl; $RS Rscript 88h_ppc_varpart.R > "$LG/88h_p2.log" 2>&1 & bgadd $! 88h_p2
wl; $RS Rscript 88k_struct_kh.R review > "$LG/88k_review_p2.log" 2>&1 & bgadd $! 88k_review_p2
waitall; tm "part 2 88g 88h 88k done"
(cd "$SK" && TAG=$T ST_VER=std3 python3 01_스크립트/89g_sens_combine.py $T std3) > "$LG/89g_p2.log" 2>&1; rcl 89g_sens_combine_p2 $?
(cd "$L" && python3 90_ltcs_summary_numbers.py) > "$LG/90_p2.log" 2>&1; rcl 90_summary_numbers_p2 $?
tm "part 2 done"; echo ALL_DONE
