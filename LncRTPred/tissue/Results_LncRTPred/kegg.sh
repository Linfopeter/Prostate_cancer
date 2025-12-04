#!/usr/bin/env bash
# ==============================================================================
# kegg.sh — KEGG (humano) por carpeta TCONS y archivo annotated*.csv
# Ejecutar DESDE: Results_LncRTPred
# Salida por archivo: TCONS_xxxx/annotated_*.kegg.tsv con:
#   gene  kegg_gene_id  pathway_id  pathway_name
# Requisitos: bash, awk, sort, tr, curl
# ==============================================================================

set -euo pipefail
export LC_ALL=C

API_SLEEP="${API_SLEEP:-0.25}"   # pausa entre llamadas a KEGG
MAX_RETRY=4                      # reintentos por llamada
BASE_DIR="$(pwd)"
MAP_PATH="$BASE_DIR/_pathway_map.tsv"

# ---------------- Helpers ----------------
curl_retry() {
  # Uso: curl_retry URL
  local url="$1"
  local attempt=1
  local delay=0.5
  while :; do
    # -f: fail on HTTP error / -sS: silent+show errors / --max-time: evita cuelgues
    if out="$(curl -fsS --max-time 10 --retry 0 "$url")"; then
      printf '%s' "$out"
      return 0
    fi
    if (( attempt >= MAX_RETRY )); then
      return 1
    fi
    sleep "$delay"
    delay=$(awk -v d="$delay" 'BEGIN{printf "%.2f", d*1.7+0.2}')
    attempt=$((attempt+1))
  done
}

ensure_pathway_map() {
  if [[ ! -s "$MAP_PATH" ]]; then
    echo "[i] Descargando catálogo KEGG humano → $MAP_PATH"
    curl_retry "https://rest.kegg.jp/list/pathway/hsa" \
      | awk -F'\t' '{ sub(/^path:/,"",$1); print $1"\t"$2 }' > "$MAP_PATH" || {
        echo "[!] No pude descargar el catálogo de pathways. ¿Red OK?" >&2
        exit 1
      }
  fi
}

# Obtiene el primer KEGG gene ID humano para un símbolo (hsa:NNNN)
get_hsa_gene_id() {
  # Buscamos por símbolo sin 'hsa:' para evitar 400; luego filtramos hsa:
  # Ejemplo de línea: hsa:2332\tFMR1; fragile X mental retardation 1
  local symbol="$1"
  curl_retry "https://rest.kegg.jp/find/genes/${symbol}" \
    | awk -F'\t' '$1 ~ /^hsa:/ {print $1; exit}'
}

# Dado un hsa:NNNN, lista pathways (path:hsa04612). Luego quitamos 'path:'.
list_pathways_for_hsa() {
  local hsa_id="$1"
  curl_retry "https://rest.kegg.jp/link/pathway/${hsa_id}" \
    | awk '{ sub(/^path:/,"",$2); print $2 }'
}

process_csv() {
  local dir="$1"
  local csv="$2"
  local base="${csv%.csv}"
  local out="${base}.kegg.tsv"
  local tmp_genes="$dir/.genes.tmp"
  local tmp_raw="$dir/.raw.tmp"

  # Extrae símbolos únicos de la columna 4 (gene_name)
  awk -F',' 'NR>1 && $4!="" {print $4}' "$csv" | tr -d '\r' | sort -u > "$tmp_genes"
  if [[ ! -s "$tmp_genes" ]]; then
    echo "[!] $csv: col4 vacía; omito."
    rm -f "$tmp_genes" "$tmp_raw"
    return 0
  fi

  : > "$tmp_raw"
  while IFS= read -r sym; do
    [[ -z "$sym" ]] && continue

    # Mapear símbolo → hsa:NNNN (robusto a 400)
    hsa_id="$(get_hsa_gene_id "$sym" || true)"
    if [[ -z "${hsa_id:-}" ]]; then
      # No hallado en KEGG
      continue
    fi

    # Obtener pathways de ese gen
    if paths="$(list_pathways_for_hsa "$hsa_id" || true)"; then
      if [[ -n "$paths" ]]; then
        # Escribir filas: gene \t hsa:NNNN \t hsa04612
        while IFS= read -r pid; do
          printf "%s\t%s\t%s\n" "$sym" "$hsa_id" "$pid" >> "$tmp_raw"
        done <<< "$paths"
      fi
    fi

    sleep "$API_SLEEP"
  done < "$tmp_genes"

  # Escribir salida con nombre de pathway
  {
    echo -e "gene\tkegg_gene_id\tpathway_id\tpathway_name"
    if [[ -s "$tmp_raw" ]]; then
      awk -F'\t' 'FNR==NR{map[$1]=$2; next} {name=(($3 in map)?map[$3]:"NA"); print $0"\t"name}' \
        "$MAP_PATH" "$tmp_raw"
    fi
  } > "$out"

  rm -f "$tmp_genes" "$tmp_raw"
  echo "[✓] $dir → $(basename "$out")"
}

process_dir() {
  local dir="$1"
  shopt -s nullglob
  local csvs=("$dir"/annotated*.csv)
  if (( ${#csvs[@]} == 0 )); then
    echo "[!] $dir: sin annotated*.csv; omito."
    return 0
  fi
  for csv in "${csvs[@]}"; do
    process_csv "$dir" "$csv"
  done
}

# ---------------- Main ----------------
# 1) asegurar mapa
ensure_pathway_map

# 2) recorrer carpetas TCONS_* secuencial (evita rate-limit 400)
mapfile -t DIRS < <(find "$BASE_DIR" -maxdepth 1 -type d -name 'TCONS_*' | sort)
if (( ${#DIRS[@]} == 0 )); then
  echo "ERROR: No hay carpetas TCONS_* aquí." >&2
  exit 1
fi

for d in "${DIRS[@]}"; do
  process_dir "$d"
done

echo "Listo."
