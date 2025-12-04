#!/usr/bin/env bash
# =============================================================================
# LncRTPred — Pipeline idempotente y portable (protein_coding, single-line, 4 shards)
# Requisitos: bash, awk, grep, sed, gzip/gunzip, wget, coreutils (sha*, sort, etc.)
# Estructura esperada:
#   ./Human/Input_Data/lncrnas_pbmcs.fasta
#   ./Human/predict_lncrna_mrna_human.py
#   ./Human/Models/*.sav
# Salidas:
#   ./Human/Input_Data/query_lncrna.fa
#   ./data/target_mrna.fa
#   ./Human/split_mrna_parts/mrna_parte_{1..4}.fa
#   ./Human/Output_Data/Predicted_Val.part{1..4}.csv
#   ./Human/Output_Data/Resultado_Final.csv
# =============================================================================
set -euo pipefail


# -------------------- Activate conda env --------------------
ENV_NAME="lncrtpred_env"

if command -v conda >/dev/null 2>&1; then
  # Initialize conda shell function and activate env
  eval "$(conda shell.bash hook)"
  if conda env list | grep -qE "^\s*${ENV_NAME}\s"; then
    echo "✅ Activating conda env: ${ENV_NAME}"
    conda activate "$ENV_NAME"
  else
    echo "⚙️  Creating conda env '${ENV_NAME}' (LncRTPred dependencies)…"
    conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda python=3.8 pandas=1.3 numpy=1.20 scikit-learn=0.24 lightgbm
    conda activate "$ENV_NAME"
  fi
else
  echo "❌ conda not found in PATH. Please load conda before running this script."
  exit 1
fi


# -------------------- Config --------------------
ENSEMBL_URL="https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$ROOT/data"
HUMAN_DIR="$ROOT/Human"
IN_DIR="$HUMAN_DIR/Input_Data"
OUT_DIR="$HUMAN_DIR/Output_Data"
SPLIT_DIR="$HUMAN_DIR/split_mrna_parts"
RAW_GZ="$DATA_DIR/Homo_sapiens.GRCh38.cdna.all.fa.gz"
TARGET_FA="$DATA_DIR/target_mrna.fa"            # protein_coding + single-line
QUERY_IN_FASTA="$IN_DIR/TCONS_00590288.fasta"
QUERY_FA="$IN_DIR/query_lncrna.fa"              # saneado + single-line
PY_SCRIPT="$HUMAN_DIR/predict_lncrna_mrna_human.py"
RESULT_FINAL="$OUT_DIR/Resultado_Final.csv"
SHARDS=4

# -------------------- Helpers --------------------
log(){ echo -e "[`date +'%F %T'`] $*"; }

require_file(){
  local f="$1"
  if [[ ! -s "$f" ]]; then
    echo "[ERR] Falta archivo requerido: $f" >&2
    exit 1
  fi
}

# single-line FASTA + filtro alfabeto A/C/G/T/U (case-insensitive) + UPPER
# Si encuentra caracteres fuera del alfabeto, aborta listando headers problemáticos.
fasta_singleline_validate(){
  # $1=infile, $2=outfile
  local in="$1" out="$2"
  awk '
    BEGIN{RS=">"; ORS=""; bad=0}
    NR>1{
      split($0, L, "\n")
      h=L[1]
      seq=""
      for(i=2;i<=length(L);i++){
        gsub(/[ \t\r]/,"",L[i])
        seq=seq L[i]
      }
      # Validar alfabeto
      if (seq !~ /^[ACGTUacgtu]*$/){
        if (!reported[h]++){
          print h "\n" >> "/dev/stderr"
        }
        bad=1
      }
      # Emitir en mayúsculas, single-line
      print ">" h "\n" toupper(seq) "\n" > "'$out'"
    }
    END{
      if (bad==1){
        exit 42
      }
    }
  ' "$in" 2> "$out.badalphabet.headers" || {
    if [[ $? -eq 42 ]]; then
      echo "[ERR] Se encontraron registros con caracteres fuera de A/C/G/T/U. Revise headers en: $out.badalphabet.headers" >&2
      rm -f "$out"
      exit 1
    else
      echo "[ERR] Fallo procesando $in" >&2
      rm -f "$out" "$out.badalphabet.headers"
      exit 1
    fi
  }
  rm -f "$out.badalphabet.headers"
}

# Conteo seguro de registros FASTA
fasta_count(){
  # $1=infile
  grep -c '^>' "$1"
}

# Crear shards por records (balanceados, ceiling)
split_fasta_into_shards(){
  # $1=infile, $2=outdir, $3=shards
  local in="$1" outdir="$2" shards="$3"
  mkdir -p "$outdir"
  local total recs_per part
  total=$(fasta_count "$in")
  if [[ "$total" -eq 0 ]]; then
    echo "[ERR] FASTA vacío: $in" >&2
    exit 1
  fi
  # ceiling
  recs_per=$(( (total + shards - 1) / shards ))
  log "Dividiendo $in → $shards shards (~$recs_per registros/archivo)"
  awk -v outdir="$outdir" -v perfile="$recs_per" '
    BEGIN{RS=">"; ORS=""; c=0; f=1; n=0}
    NR>1{
      rec=">" $0
      n++
      if (n==1 || ((n-1)%perfile)==0){
        if (c){ close(file) }
        file=sprintf("%s/mrna_parte_%d.fa", outdir, f++)
        c=1
      }
      print rec > file
    }' "$in"
  # QA: asegurar 1..shards existen (últimos pueden faltar si total<shards)
  local i any=0
  for ((i=1;i<=shards;i++)); do
    if [[ -s "$outdir/mrna_parte_$i.fa" ]]; then any=1; fi
  done
  if [[ "$any" -ne 1 ]]; then
    echo "[ERR] No se generaron shards en $outdir" >&2
    exit 1
  fi
}

# Concatena resultados Predicted_Val.* a Resultado_Final.csv con encabezado único y orden por shard
concat_results(){
  local out_csv="$1"; shift
  local parts=("$@")
  echo "lncrna,mrna,Decision_Tree(%),K-Nearest-Neighbours(%),Random_Forest(%),LightGBM(%),Cumulative_Model_Score(%)" > "$out_csv"
  local p
  for p in "${parts[@]}"; do
    if [[ ! -s "$p" ]]; then
      echo "[ERR] Falta resultado: $p" >&2
      exit 1
    fi
    tail -n +2 "$p" >> "$out_csv"
  done
  log "Resultado combinado: $out_csv (filas: $(($(wc -l < "$out_csv")-1)))"
}

# -------------------- Preparación dirs --------------------
mkdir -p "$DATA_DIR" "$OUT_DIR" "$SPLIT_DIR"

# -------------------- Verificaciones base --------------------
require_file "$PY_SCRIPT"
require_file "$QUERY_IN_FASTA"

# -------------------- 1) Descargar Ensembl cDNA all (reanudable) --------------------
if [[ ! -s "$RAW_GZ" ]]; then
  log "Descargando Ensembl cDNA all (GRCh38 r115)…"
  wget -c -O "$RAW_GZ" "$ENSEMBL_URL"
else
  log "Encontrado: $RAW_GZ (omito descarga)"
fi

# Integridad gzip
log "Verificando integridad gzip…"
gzip -t "$RAW_GZ" || { echo "[ERR] Archivo corrupto: $RAW_GZ" >&2; rm -f "$RAW_GZ"; exit 1; }

# -------------------- 2) Filtrar protein_coding + single-line → target_mrna.fa --------------------
if [[ ! -s "$TARGET_FA" ]]; then
  log "Extrayendo transcript_biotype:protein_coding y normalizando a single-line/UPPER…"
  tmp_raw="$TARGET_FA.tmp.raw"
  tmp_ok="$TARGET_FA.tmp.ok"
  EXCL_TSV="$DATA_DIR/target_mrna.excl.tsv"

  # 2a) Extrae SOLO transcript_biotype:protein_coding (más estricto que gene_biotype)
  #     y deja single-line (sin validar alfabeto todavía).
  awk '
    BEGIN{RS=">"; ORS=""}
    NR>1{
      split($0, L, "\n")
      h=L[1]; seq=""
      for(i=2;i<=length(L);i++){ gsub(/[ \t\r]/,"",L[i]); seq=seq L[i] }
      if (h ~ /transcript_biotype:protein_coding/){
        print ">" h "\n" toupper(seq) "\n"
      }
    }' <(gunzip -c "$RAW_GZ") > "$tmp_raw"

  # 2b) Valida alfabeto A/C/G/T/U. En vez de abortar, EXCLUYE los que fallen
  #     y deja registro auditable en $EXCL_TSV
  echo -e "transcript_id\treason" > "$EXCL_TSV"
  awk -v excl="$EXCL_TSV" '
    BEGIN{RS=">"; ORS=""}
    NR>1{
      split($0,L,"\n"); h=L[1]; seq=L[2]
      if (seq ~ /^[ACGTU]*$/){
        print ">" h "\n" seq "\n"
      } else {
        # reporta el id (hasta el primer espacio de la cabecera) y motivo
        split(h, F, " ")
        print F[1] "\tNON_ACGTU" >> excl
      }
    }' "$tmp_raw" > "$tmp_ok"

  mv -f "$tmp_ok" "$TARGET_FA"
  rm -f "$tmp_raw"
  log "target_mrna.fa listo: $(grep -c '^>' "$TARGET_FA") registros (excluidos: $(($(wc -l < "$EXCL_TSV")-1)))"
else
  log "Encontrado: $TARGET_FA (omito filtrado)"
fi


# -------------------- 3) Preparar query lncrna.fa (single-line validado) --------------------
if [[ ! -s "$QUERY_FA" ]]; then
  log "Normalizando query lncRNA → $QUERY_FA (single-line, A/C/G/T/U, UPPER)…"
  fasta_singleline_validate "$QUERY_IN_FASTA" "$QUERY_FA"
  log "query_lncrna.fa listo: $(fasta_count "$QUERY_FA") registros"
else
  log "Encontrado: $QUERY_FA (omito normalización)"
fi

# -------------------- 4) Partir target en 4 shards por records --------------------
# Re-generar shards sólo si faltan o están vacíos
need_split=0
for i in $(seq 1 $SHARDS); do
  [[ -s "$SPLIT_DIR/mrna_parte_$i.fa" ]] || { need_split=1; break; }
done
if [[ "$need_split" -eq 1 ]]; then
  rm -f "$SPLIT_DIR"/mrna_parte_*.fa || true
  split_fasta_into_shards "$TARGET_FA" "$SPLIT_DIR" "$SHARDS"
else
  log "Shards existentes en $SPLIT_DIR (omito split)"
fi

# QA rápida
log "QA conteos:"
log "  - query : $(fasta_count "$QUERY_FA")"
log "  - target: $(fasta_count "$TARGET_FA")"
for i in $(seq 1 $SHARDS); do
  if [[ -s "$SPLIT_DIR/mrna_parte_$i.fa" ]]; then
    log "  - mrna_parte_$i.fa: $(fasta_count "$SPLIT_DIR/mrna_parte_$i.fa")"
  fi
done

# -------------------- 5) Ejecutar LncRTPred por shard (secuencial) --------------------
# Nota: si deseas paralelizar, puedes lanzar cada bloque en subshell "&" y luego wait.
parts_csv=()
for i in $(seq 1 $SHARDS); do
  shard="$SPLIT_DIR/mrna_parte_$i.fa"
  [[ -s "$shard" ]] || { log "[SKIP] No existe $shard (posible total<shards)."; continue; }
  part_csv="$OUT_DIR/Predicted_Val.part$i.csv"
  parts_csv+=("$part_csv")

  if [[ -s "$part_csv" ]]; then
    log "[OK] Resultado shard $i ya existe: $part_csv (omito corrida)"
    continue
  fi

  log "=== Shard $i/$SHARDS ==="
  # idempotente: copia target al lugar que espera LncRTPred
  cp -f "$shard" "$IN_DIR/target_mrna.fa"

  # Ejecutar predictor
  ( cd "$HUMAN_DIR" && python "$PY_SCRIPT" )

  # Guardar resultado individual y limpiar target temporal
  if [[ ! -s "$OUT_DIR/Predicted_Val.csv" ]]; then
    echo "[ERR] No se generó Output_Data/Predicted_Val.csv para shard $i" >&2
    exit 1
  fi
  mv -f "$OUT_DIR/Predicted_Val.csv" "$part_csv"
  rm -f "$IN_DIR/target_mrna.fa"
  log "[OK] Guardado: $part_csv (filas: $(wc -l < "$part_csv"))"
done

# -------------------- 6) Concatenar resultados --------------------
if [[ "${#parts_csv[@]}" -eq 0 ]]; then
  echo "[ERR] No hay partes procesadas para concatenar." >&2
  exit 1
fi
concat_results "$RESULT_FINAL" "${parts_csv[@]}"

log "Listo ✅"
