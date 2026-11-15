#!/usr/bin/env bash
# Re-run every script and store its stdout in results/logs/<script>.txt
# Usage (from repository root):  bash scripts/run_all.sh
set -u
cd "$(dirname "$0")"
mkdir -p ../results/logs
export PYTHONIOENCODING=utf-8 MPLBACKEND=Agg
run() { name=$1; shift; echo ">> $name $*"; s=$(date +%s);
        python "$name.py" "$@" > "../results/logs/$name.txt" 2>&1; rc=$?;
        echo "   exit=$rc  $(( $(date +%s) - s )) s" | tee -a ../results/logs/_timings.txt; }
: > ../results/logs/_timings.txt
run blg_kp
run blg_dot
run blg_floquet_full
run blg_ac_fast 0 0.005 0.010 0.020 0.030 0.050
run blg_sweetspot
run blg_doubledot
run blg_dext_scan
run blg_d60
run blg_compare
run blg_d60_verify
run blg_converge
run blg_valley
