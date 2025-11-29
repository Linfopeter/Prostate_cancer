#!/usr/bin/env bash
# =============================================================================
# IntaRNA — Heuristic → Exact por cada target ENST/Homo dentro de cada TCONS_*
# - CSV crudo de IntaRNA (SIN edición): todas las columnas definidas por IntaRNA
#   vía --outMode C --outCsvCols '' y separador ';'
# - TXT detallado (outMode D) para heurístico y exacto
# - Iteración robusta por carpetas TCONS_*
# =============================================================================
set -euo pipefail
export LC_ALL=C

# ---------- Parámetros ----------
THREADS="${THREADS:-20}"
PAD="${PAD:-50}"              # padding ±nt para delimitar ventana exacta
SEP=','                       # separador CSV (default IntaRNA)
OUTCOLS=''                    # '' => TODAS las columnas definidas por IntaRNA

HEUR_ARGS=( --mode H --outMode C --outSep "${SEP}" --outCsvCols "${OUTCOLS}" --outNumber 1 )
EXACT_ARGS=( --mode M --outMode C --outSep "${SEP}" --outCsvCols "${OUTCOLS}" --outNoGUend --outNoLP )

# ---------- Helpers ----------
fasta_len() { awk 'BEGIN{L=0} /^>/{next} {gsub(/[ \t\r]/,""); L+=length($0)} END{print L}' "$1"; }

calc_regions() {
  local s1="$1" e1="$2" s2="$3" e2="$4" tlen="$5" qlen="$6" pad="$7"
  (( s1>e1 )) && { local x=$s1; s1=$e1; e1=$x; }
  (( s2>e2 )) && { local y=$s2; s2=$e2; e2=$y; }
  local tL=$(( s1-pad )); (( tL<1 )) && tL=1
  local tR=$(( e1+pad )); (( tR>tlen )) && tR=$tlen
  local qL=$(( s2-pad )); (( qL<1 )) && qL=1
  local qR=$(( e2+pad )); (( qR>qlen )) && qR=$qlen
  printf -- "--tRegion=%d-%d --qRegion=%d-%d" "$tL" "$tR" "$qL" "$qR"
}

target_basename() {
  local base; base="$(basename "$1")"
  if   [[ "$base" =~ (ENST[0-9]+_[0-9]+) ]]; then echo "${BASH_REMATCH[1]}"
  elif [[ "$base" =~ (ENST[0-9]+)       ]]; then echo "${BASH_REMATCH[1]}"
  else echo "${base%.*}"; fi
}

list_targets() {
  shopt -s nullglob
  local arr=() t
  arr+=(Homo_sapiens_ENST*sequence.fa)
  arr+=(ENST*.fa)
  arr+=(ENST*.fasta)
  for t in "${arr[@]}"; do [[ -e "$t" ]] && printf '%s\n' "$t"; done
  shopt -u nullglob
}

write_detailed_txt() {
  local mode="$1"; local tff="$2"; local qf="$3"; local out_txt="$4"; local err_txt="$5"; shift 5
  local extra_args=( "$@" )
  IntaRNA -t "$tff" -q "$qf" --threads "$THREADS" --mode "$mode" --outMode D \
          "${extra_args[@]}" > "$out_txt" 2> "$err_txt" || true
}

run_one_target() {
  local qf="$1" tff="$2" threads="$3" pad="$4"
  local tname; tname="$(target_basename "$tff")"

  echo "[INFO] $(pwd) | Target=$tname -> Heuristic (CSV + TXT)"

  # Heurístico: CSV crudo (todas columnas) + TXT detallado
  if ! IntaRNA -t "$tff" -q "$qf" --threads "$threads" "${HEUR_ARGS[@]}" \
        > "heuristic_${tname}.csv" 2> "heuristic_${tname}.err"; then
    echo "[WARN] Heurístico falló (exit). Revisa heuristic_${tname}.err"
  fi
  write_detailed_txt "H" "$tff" "$qf" "heuristic_${tname}.txt" "heuristic_${tname}.txt.err"

  # Si no hubo hits, no intentamos exacto
  if [[ ! -s "heuristic_${tname}.csv" || $(wc -l < "heuristic_${tname}.csv") -lt 2 ]]; then
    echo "[WARN] $(pwd) | $tname: sin hits heurísticos. Saltando exacto."
    [[ -s "heuristic_${tname}.err" ]] && tail -n 5 "heuristic_${tname}.err" || true
    return 0
  fi

  # Tomar la PRIMERA fila de resultados heurísticos para centrar la ventana exacta
  # Notas:
  #  - NO reordenamos ni tocamos el CSV; leemos columnas por nombre.
  #  - start1/end1 = target; start2/end2 = query, según docs.
  local tlen qlen; tlen="$(fasta_len "$tff")"; qlen="$(fasta_len "$qf")"
  local s1 e1 s2 e2
  read -r s1 e1 s2 e2 < <(
    awk -v FS="${SEP}" '
      NR==1 {
        for(i=1;i<=NF;i++) h[$i]=i
        next
      }
      NR==2 {
        print $(h["start1"]), $(h["end1"]), $(h["start2"]), $(h["end2"])
      }' "heuristic_${tname}.csv"
  )
  if [[ -z "${s1:-}" || -z "${e1:-}" || -z "${s2:-}" || -z "${e2:-}" ]]; then
    echo "[WARN] $(pwd) | $tname: no pude parsear coordenadas del heurístico."
    return 0
  fi

  # Exacto centrado en ventana alrededor del hit heurístico
  local region_str region_args=()
  region_str="$(calc_regions "$s1" "$e1" "$s2" "$e2" "$tlen" "$qlen" "$pad")"
  read -r -a region_args <<< "$region_str"

  echo "[INFO] $(pwd) | Target=$tname -> Exact ${region_args[*]} (CSV + TXT)"

  if ! IntaRNA -t "$tff" -q "$qf" --threads "$threads" "${EXACT_ARGS[@]}" "${region_args[@]}" \
        > "exact_${tname}.csv" 2> "exact_${tname}.err"; then
    echo "[WARN] Exacto falló (exit). Revisa exact_${tname}.err"
  fi
  write_detailed_txt "M" "$tff" "$qf" "exact_${tname}.txt" "exact_${tname}.txt.err" "${region_args[@]}"
}

summarize_folder() {
  # Resumen no altera CSVs originales; sólo lee para contar filas y extraer mejor E
  local out="summary_exact.csv"
  echo "tcons_folder,target_csv,rows,best_E,best_id1,best_id2,best_start1,best_end1,best_start2,best_end2" > "$out"
  shopt -s nullglob
  local f rows
  for f in exact_*.csv; do
    rows=$(( $(wc -l < "$f") - 1 ))
    if (( rows <= 0 )); then
      printf "%s,%s,0,NA,NA,NA,NA,NA,NA,NA\n" "$(basename "$(pwd)")" "$(basename "$f")" >> "$out"
      continue
    fi
    # Encontrar la fila con E mínima leyendo la columna por nombre (sin ordenar archivo)
    awk -v FS="${SEP}" -v OFS="," -v fbase="$(basename "$f")" -v tdir="$(basename "$(pwd)")" '
      NR==1 { for(i=1;i<=NF;i++) h[$i]=i; bestE=""; best="" ; next }
      NR>1 {
        e = $(h["E"])+0
        if (best=="" || e < bestE) { bestE=e; best=$0 }
      }
      END {
        if (best=="") { print tdir, fbase, 0, "NA","NA","NA","NA","NA","NA","NA"; exit }
        split(best,a,FS)
        print tdir, fbase, NR-1, bestE, a[h["id1"]], a[h["id2"]], a[h["start1"]], a[h["end1"]], a[h["start2"]], a[h["end2"]]
      }' "$f" >> "$out"
  done
  shopt -u nullglob
}

run_folder() {
  local dir="$1"
  cd "$dir"

  local LOG="$PWD/intarna.log"; touch "$LOG"

  # Query
  local qf
  qf="$(ls TCONS_*.fasta 2>/dev/null | head -n1 || true)"
  if [[ -z "$qf" ]]; then
    echo "[WARN] $dir: sin query TCONS_*.fasta — omito." | tee -a "$LOG"
    cd - >/dev/null; return 0
  fi

  # Targets
  mapfile -t targets < <(list_targets)
  if ((${#targets[@]}==0)); then
    echo "[WARN] $dir: sin targets Homo/ENST." | tee -a "$LOG"
    cd - >/dev/null; return 0
  fi

  echo "[INFO] ===== Carpeta: $dir | Query: $qf | Targets: ${#targets[@]} =====" | tee -a "$LOG"

  local tff
  for tff in "${targets[@]}"; do
    { echo "[LOG] $(date -Iseconds) | target=$tff"; run_one_target "$qf" "$tff" "$THREADS" "$PAD"; } | tee -a "$LOG"
  done

  summarize_folder
  echo "[INFO] ===== Resumen escrito en $(pwd)/summary_exact.csv =====" | tee -a "$LOG"

  cd - >/dev/null
}

# ---------- Main ----------
root="${1:-.}"
mapfile -t FOLDERS < <(find "$root" -maxdepth 1 -type d -name 'TCONS_*' | sort)
for d in "${FOLDERS[@]}"; do
  run_folder "$d"
done

echo "[DONE] Refinamiento exacto por ventanas completado."
