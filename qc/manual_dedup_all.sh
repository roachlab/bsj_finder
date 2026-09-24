#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# Manual re-run of the mmseqs2_dedup rule across ALL samples.
# Replicates workflow/rules/mmseqs2_dedup.smk exactly:
#   - same params  (min-seq-id 0.95, -c 0.70, --cov-mode 0)
#   - same edge cases (0 seqs -> empty outputs; 1 seq -> skip clustering)
#   - SSE4.1 static mmseqs (tools/mmseqs/bin/mmseqs) for this AVX2-less CPU
# Outputs go to qc/manual_dedup/ (does NOT touch results/), so we can diff
# against results/dedup_results/ afterwards.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

MMSEQS="$(pwd)/tools/mmseqs/bin/mmseqs"
[[ -x "$MMSEQS" ]] || { echo "missing $MMSEQS" >&2; exit 1; }
MIN_SEQ_ID=0.95
MIN_COV=0.70
COV_MODE=0
OUT=qc/manual_dedup
LOG=qc/manual_dedup_run.log
mkdir -p "$OUT"
: > "$LOG"

echo "[$(date)] starting manual dedup: $(basename "$MMSEQS")" | tee -a "$LOG"
for in_fasta in results/bsj_results/*_bsj_sequences.fasta.gz; do
  s="$(basename "$in_fasta" _bsj_sequences.fasta.gz)"
  o_fasta="$OUT/$s.bsj_dedup.fasta"
  o_tsv="$OUT/$s.bsj_cluster.tsv"
  workdir="$OUT/${s}_mmseqs2_tmp"

  SEQ_COUNT=$(gunzip -c "$in_fasta" | grep -c "^>" || echo 0)
  SEQ_COUNT=$(echo "$SEQ_COUNT" | tr -d '[:space:]')
  echo "[$(date)] $s input=$SEQ_COUNT" >> "$LOG"

  if [[ "$SEQ_COUNT" -eq 0 ]]; then
    gunzip -c "$in_fasta" > "$o_fasta"
    printf "query\ttarget\n" > "$o_tsv"
    echo "[$(date)] $s -> empty (0 input)" >> "$LOG"
    continue
  fi

  if [[ "$SEQ_COUNT" -eq 1 ]]; then
    gunzip -c "$in_fasta" > "$o_fasta"
    gunzip -c "$in_fasta" | grep "^>" | sed 's/^>//' | awk '{print $1"\t"$1}' > "$o_tsv"
    echo "[$(date)] $s -> single seq, skipped clustering" >> "$LOG"
    continue
  fi

  mkdir -p "$workdir/tmp"
  SEQDB="$workdir/seqDB"; CLUDB="$workdir/cluDB"; REPDB="$workdir/repDB"

  "$MMSEQS" createdb "$in_fasta" "$SEQDB" >>"$LOG" 2>&1
  "$MMSEQS" linclust "$SEQDB" "$CLUDB" "$workdir/tmp" \
      --min-seq-id "$MIN_SEQ_ID" -c "$MIN_COV" --cov-mode "$COV_MODE" >>"$LOG" 2>&1
  "$MMSEQS" result2repseq "$SEQDB" "$CLUDB" "$REPDB" >>"$LOG" 2>&1
  "$MMSEQS" convert2fasta "$REPDB" "$o_fasta" >>"$LOG" 2>&1
  "$MMSEQS" createtsv "$SEQDB" "$SEQDB" "$CLUDB" "$o_tsv" >>"$LOG" 2>&1

  OUT_COUNT=$(grep -c ">" "$o_fasta" || echo 0)
  rm -rf "$workdir"
  echo "[$(date)] $s unique=$OUT_COUNT" >> "$LOG"
done
echo "[$(date)] done" | tee -a "$LOG"
