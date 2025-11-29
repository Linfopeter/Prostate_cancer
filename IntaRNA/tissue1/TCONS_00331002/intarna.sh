#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# =======================
# IntaRNA pipeline (heuristic -> exact) for TCONS vs Homo_sapiens_ENST*
# =======================
# Behavior
# - Run inside a TCONS_* folder (one TCONS_*.fasta + many Homo_sapiens_ENST*.fa)
#   OR in a parent directory; then it will recurse into each TCONS_* folder.
# - Produces per-target:
#     heuristic_{TARGET}.csv, heuristic_{TARGET}.txt
#     exact_{TARGET}.csv,     exact_{TARGET}.txt
# - Produces per-folder summary: summary_exact.csv (recreated each run)
#
# Tunables (env vars):
: "${THREADS:=}"        # auto-detect below if empty
: "${PADDING:=25}"      # bases to expand around heuristic hit for exact mode
: "${OUTNUM:=3}"        # number of interactions to report (heuristic & exact)
: "${FORCE:=0}"         # set FORCE=1 to overwrite per-target outputs even if present

# --- detect threads portably ---
if [[ -z "${THREADS}" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    THREADS="$(nproc)"
  elif command -v sysctl >/dev/null 2>&1; then
    THREADS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
  else
    THREADS=4
  fi
fi

# --- need(): require a tool in PATH ---
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing '$1' in PATH" >&2; exit 127; }; }
need IntaRNA
need awk
need sed
need grep

# --- FASTA length (sum of non-header lines) ---
seq_len() {
  # $1 = fasta file
  awk 'BEGIN{len=0} /^>/ {next} {gsub(/[ \t\r\n]/,""); len+=length($0)} END{print len}' "$1"
}

# --- clamp v to [lo,hi] ---
clamp() {
  # v lo hi
  local v="$1" lo="$2" hi="$3"
  if (( v < lo )); then echo "$lo"; elif (( v > hi )); then echo "$hi"; else echo "$v"; fi
}

# --- parse first (best) CSV row into variables ---
# expects header including: id1,start1,end1,id2,start2,end2,E  (among others)
parse_best_row() {
  # $1 = csv
  awk -F';' '
    NR==1{
      for(i=1;i<=NF;i++){h[$i]=i}
      next
    }
    NR==2{
      # expose as shell assignments
      printf "BEST_ID1=%s\nBEST_S1=%s\nBEST_E1=%s\nBEST_ID2=%s\nBEST_S2=%s\nBEST_E2=%s\nBEST_E=%s\n",
             $h["id1"], $h["start1"], $h["end1"],
             $h["id2"], $h["start2"], $h["end2"], $h["E"]
      exit
    }' "$1"
}

# --- Run IntaRNA (one mode) producing CSV OR TXT ---
run_intarna_csv() {
  # $1 mode(H/M), $2 query.fa, $3 target.fa, $4 out_csv
  IntaRNA \
    -q "$2" -t "$3" \
    --threads "$THREADS" \
    -m "$1" \
    --acc C --accW 150 --accL 100 \
    --outMode C \
    --outCsvCols 'id2,id1,start2,end2,start1,end1,subseqDP,hybridDP,E,Etotal,ED2,ED1,Pu2,Pu1,seedStart2,seedEnd2,seedStart1,seedEnd1,seedE' \
    --outCsvSort E \
    --outNumber "$OUTNUM" \
    --out "$4"
}

run_intarna_txt() {
  # $1 mode(H/M), $2 query.fa, $3 target.fa, $4 out_txt, [optional qRegion] [optional tRegion]
  local mode="$1" q="$2" t="$3" out="$4" qR="${5-}" tR="${6-}"
  local args=(
    -q "$q" -t "$t"
    --threads "$THREADS"
    -m "$mode"
    --acc C --accW 150 --accL 100
    --outMode D
    --out "$out"
  )
  [[ -n "$qR" ]] && args+=( --qRegion "$qR" )
  [[ -n "$tR" ]] && args+=( --tRegion "$tR" )
  IntaRNA "${args[@]}"
}

# --- per-folder worker: run pipeline inside one TCONS_* directory ---
process_one_folder() {
  local d="$1"
  ( cd "$d"
    echo "[INFO] === Processing folder: $PWD ==="
    # find query
    local q
    q="$(ls -1 TCONS_*.fasta 2>/dev/null | head -n1 || true)"
    [[ -n "$q" ]] || { echo "[WARN] No TCONS_*.fasta found in $PWD, skipping." >&2; return 0; }
    local qlen; qlen="$(seq_len "$q")"

    # init / reset summary
    local summary="summary_exact.csv"
    : > "$summary"
    echo "target_id,query_id,E_exact,start1,end1,start2,end2,qRegion,tRegion,heuristic_csv,exact_csv" >> "$summary"

    # enumerate targets
    shopt -s nullglob
    local t
    local found_any=0
    for t in Homo_sapiens_ENST*.fa; do
      found_any=1
      local base="${t##*/}"
      base="${base%.fa}"

      local h_csv="heuristic_${base}.csv"
      local h_txt="heuristic_${base}.txt"
      local e_csv="exact_${base}.csv"
      local e_txt="exact_${base}.txt"

      if (( FORCE == 0 )) && [[ -s "$e_csv" && -s "$e_txt" ]]; then
        echo "[SKIP] exact exists for $base"
        # still add to summary (pull best row if present)
        if [[ -s "$e_csv" ]]; then
          # shellcheck disable=SC1090
          eval "$(parse_best_row "$e_csv" || true)"
          if [[ -n "${BEST_E:-}" ]]; then
            # unknown qRegion/tRegion here; leave blank
            echo "${BEST_ID1},${BEST_ID2},${BEST_E},${BEST_S1},${BEST_E1},${BEST_S2},${BEST_E2},,,"\
                 "${h_csv},${e_csv}" >> "$summary"
          fi
        fi
        continue
      fi

      echo "[INFO] Target: $t"

      # 1) HEURISTIC CSV
      run_intarna_csv "H" "$q" "$t" "$h_csv" || {
        echo "[WARN] Heuristic CSV failed for $t; skipping." >&2; continue;
      }

      # 1b) HEURISTIC TXT (pretty interaction print)
      run_intarna_txt "H" "$q" "$t" "$h_txt" || true

      # 2) derive regions from best heuristic hit (+/- PADDING)
      # shellcheck disable=SC1090
      eval "$(parse_best_row "$h_csv" || true)"
      if [[ -z "${BEST_E:-}" ]]; then
        echo "[WARN] No interactions in heuristic for $t"
        continue
      fi

      local tlen; tlen="$(seq_len "$t")"
      # Heuristic returns 1-based positions inclusive
      local q_start="$BEST_S2" q_end="$BEST_E2" t_start="$BEST_S1" t_end="$BEST_E1"

      # expand by PADDING and clamp to sequence lengths
      local qR_start qR_end tR_start tR_end
      qR_start="$(clamp $(( q_start - PADDING )) 1 "$qlen")"
      qR_end="$(clamp   $(( q_end   + PADDING )) 1 "$qlen")"
      tR_start="$(clamp $(( t_start - PADDING )) 1 "$tlen")"
      tR_end="$(clamp   $(( t_end   + PADDING )) 1 "$tlen")"

      local qRegion="${qR_start}-${qR_end}"
      local tRegion="${tR_start}-${tR_end}"

      # 3) EXACT CSV (restricted to regions)
      IntaRNA \
        -q "$q" -t "$t" \
        --threads "$THREADS" \
        -m M \
        --acc C --accW 150 --accL 100 \
        --qRegion "$qRegion" \
        --tRegion "$tRegion" \
        --outMode C \
        --outCsvCols 'id2,id1,start2,end2,start1,end1,subseqDP,hybridDP,E,Etotal,ED2,ED1,Pu2,Pu1,seedStart2,seedEnd2,seedStart1,seedEnd1,seedE' \
        --outCsvSort E \
        --outNumber "$OUTNUM" \
        --out "$e_csv"

      # 3b) EXACT TXT for the pretty interaction (same regions)
      run_intarna_txt "M" "$q" "$t" "$e_txt" "$qRegion" "$tRegion" || true

      # 4) append best exact hit to summary
      # shellcheck disable=SC1090
      eval "$(parse_best_row "$e_csv" || true)"
      if [[ -n "${BEST_E:-}" ]]; then
        echo "${BEST_ID1},${BEST_ID2},${BEST_E},${BEST_S1},${BEST_E1},${BEST_S2},${BEST_E2},${qRegion},${tRegion},${h_csv},${e_csv}" \
          >> "$summary"
      fi
    done
    shopt -u nullglob

    if (( found_any == 0 )); then
      echo "[WARN] No targets (Homo_sapiens_ENST*.fa) found in $PWD"
    else
      echo "[INFO] Wrote $(wc -l < "$summary") lines to $summary"
    fi
  )
}

# --- Decide scope: current folder or recurse over TCONS_* subfolders ---
if ls -1 TCONS_*.fasta >/dev/null 2>&1; then
  # looks like we are inside a TCONS_* directory
  process_one_folder "."
else
  # look for subdirs named TCONS_* that contain a TCONS_*.fasta
  shopt -s nullglob
  any=0
  for d in TCONS_*; do
    if [[ -d "$d" ]] && compgen -G "$d/TCONS_*.fasta" > /dev/null; then
      any=1
      process_one_folder "$d"
    fi
  done
  shopt -u nullglob
  if (( any == 0 )); then
    echo "[ERROR] No TCONS_* folders with TCONS_*.fasta found here: $PWD" >&2
    exit 1
  fi
fi
