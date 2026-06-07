#!/usr/bin/env bash
# src/run_netmhcpan.sh — Step 2: run NetMHCpan-4.1 (EL + BA) on the surrogate peptide list.
#
# WHY the binary directly, NOT the IEDB predict_binding.py wrapper:
#   The IEDB Python wrapper HANGS (>1 h, never even spawns the netMHCpan binary). The binary in
#   -p (peptide-list) mode predicts all peptides in ~1 s with a single model load. It is the
#   identical NetMHCpan-4.1 EL+BA model. See PROGRESS_LOG (2026-06-07 Step 2 run) for the full
#   debug trail (ruled out drvfs, peptide count, CPU contention, stale lock).
#
# Output columns of interest: Score_EL, Rnk_EL (presentation head) and Score_BA, Rnk_BA,
#   Aff(nM) (binding-affinity head; Aff(nM) = predicted IC50 for the predicted-vs-measured test).
#
# Run from WSL (the tool is WSL-native):
#   MHCI_HOME=/home/max/mhc_i bash src/run_netmhcpan.sh
set -euo pipefail

MHCI_HOME="${MHCI_HOME:-/home/max/mhc_i}"
ALLELE="${ALLELE:-HLA-A02:01}"
NETMHCPAN="$MHCI_HOME/method/netmhcpan-4.1-executable/netmhcpan_4_1_executable/netMHCpan"

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN="$REPO/results/step2/input_peptides.txt"
OUT="$REPO/results/step2/pred_netmhcpan_el_ba_raw.txt"
ERR="$REPO/results/step2/pred.stderr.log"

command -v tcsh >/dev/null || { echo "ERROR: tcsh not on PATH (sudo apt install tcsh)"; exit 1; }
[ -f "$NETMHCPAN" ] || { echo "ERROR: netMHCpan not found at $NETMHCPAN (set MHCI_HOME)"; exit 1; }
[ -f "$IN" ]        || { echo "ERROR: peptide list missing: $IN (run build_surrogate_input.py)"; exit 1; }

echo "NetMHCpan-4.1 (EL+BA) | allele=$ALLELE | $(wc -l < "$IN") peptides"
tcsh "$NETMHCPAN" -p -BA -a "$ALLELE" -f "$IN" > "$OUT" 2> "$ERR"
n=$(grep -c 'HLA-A\*' "$OUT" || true)
echo "wrote $OUT  ($n prediction rows)"
