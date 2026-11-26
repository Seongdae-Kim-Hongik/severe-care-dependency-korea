#!/bin/bash
# rslot.sh - Runs a command when one of SK_RSLOTS (default 3) CPU slots is free, at nice 10. Usage: rslot.sh Rscript 81_fit_inla_ltcs.R LTCS
N=${SK_RSLOTS:-3}; D=${TMPDIR:-/tmp}/ltcs_rslots; mkdir -p "$D"
while :; do
  for i in $(seq 1 "$N"); do
    if mkdir "$D/s$i" 2>/dev/null; then
      echo $$ > "$D/s$i/pid"; trap 'rm -rf "$D/s'$i'"' EXIT INT TERM
      nice -n 10 "$@"; exit $?
    fi
    p=$(cat "$D/s$i/pid" 2>/dev/null); [ -n "$p" ] && ! kill -0 "$p" 2>/dev/null && rm -rf "$D/s$i"
  done
  sleep 5
done
