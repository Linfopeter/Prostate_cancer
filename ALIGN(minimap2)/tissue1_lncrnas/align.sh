#!/usr/bin/env bash
# ==========================================================
# Map lncRNAs (FASTA) -> GRCh38 (minimap2 splice-aware)
# - Entorno conda autocontenido (minimap2/samtools/bedtools)
# - Auto-fix de samtools (pin de versiones y librerías)
# - Descarga idempotente del genoma (r115, primary assembly)
# - Salidas: SAM (intermedio), BAM+BAI, PAF, BED12, CSV
# - Rutas relativas; seguro ante espacios en nombres
# ==========================================================
set -euo pipefail

# -------------------------------
# Config
# -------------------------------
ENV_NAME="aln_env"
QUERY="./tissue1_filtered_xloc_tcons.fasta"
GENOME_URL="https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GENOME_GZ="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GENOME_FA="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
THREADS="$(nproc 2>/dev/null || echo 4)"

# Salidas
SAM="lncrna_vs_GRCh38.aln.sam"
BAM="lncrna_vs_GRCh38.sorted.bam"
PAF="lncrna_vs_GRCh38.paf"
BED12="lncrna_vs_GRCh38.bed12"
CSV="lncrna_vs_GRCh38.alignments.csv"

# -------------------------------
# Helpers
# -------------------------------
ensure_conda() {
  command -v conda >/dev/null 2>&1 || { echo "❌ conda no encontrado."; exit 1; }
  eval "$(conda shell.bash hook)"
}

ensure_env() {
  if conda env list | grep -qE "^\s*${ENV_NAME}\s"; then
    echo "✅ Entorno '${ENV_NAME}' existe."
  else
    echo "⚙️  Creando entorno '${ENV_NAME}' (versiones fijadas)…"
    conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
      minimap2=2.28 \
      samtools=1.19.2 \
      htslib=1.19.1 \
      bedtools=2.31.1 \
      ncurses=6.4 \
      libdeflate>=1.18 \
      zlib>=1.2.13 \
      bzip2 \
      xz
  fi
  conda activate "$ENV_NAME"
}

heal_samtools_if_needed() {
  if samtools --version >/dev/null 2>&1; then
    echo "✅ samtools OK."
    return
  fi
  echo "🩹 Reparando samtools y librerías dentro del entorno…"
  conda install -y -c conda-forge -c bioconda \
    samtools=1.19.2 htslib=1.19.1 ncurses=6.4 libdeflate>=1.18 zlib>=1.2.13 bzip2 xz
  # Segundo intento
  samtools --version >/dev/null 2>&1 || { echo "❌ samtools sigue fallando."; exit 1; }
  echo "✅ samtools reparado."
}

# -------------------------------
# 1) Entorno
# -------------------------------
ensure_conda
ensure_env
heal_samtools_if_needed
command -v minimap2 >/dev/null || { echo "❌ minimap2 no disponible."; exit 1; }
command -v bedtools  >/dev/null || { echo "❌ bedtools no disponible."; exit 1; }

# -------------------------------
# 2) Entrada
# -------------------------------
[[ -s "$QUERY" ]] || { echo "❌ FASTA de consulta vacío o ausente: $QUERY"; exit 1; }
echo "✅ Query: $QUERY"

# -------------------------------
# 3) Genoma (idempotente)
# -------------------------------
if [[ -s "$GENOME_FA" ]]; then
  echo "✅ Genoma listo: $GENOME_FA"
else
  if [[ ! -s "$GENOME_GZ" ]]; then
    echo "🌐 Descargando GRCh38 r115 (primary assembly)…"
    if command -v wget >/dev/null 2>&1; then
      wget -c "$GENOME_URL" -O "$GENOME_GZ"
    else
      curl -L -o "$GENOME_GZ" "$GENOME_URL"
    fi
  else
    echo "✅ Encontrado comprimido: $GENOME_GZ"
  fi
  echo "🗜️  Descomprimiendo (atómico)…"
  tmp="${GENOME_FA}.tmp"
  if command -v pigz >/dev/null 2>&1; then pigz -dc "$GENOME_GZ" > "$tmp"; else gunzip -c "$GENOME_GZ" > "$tmp"; fi
  mv -f "$tmp" "$GENOME_FA"
fi

# -------------------------------
# 4) Alineación (SAM -> BAM)  [sin secundarios]
# -------------------------------
if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
  echo "✅ BAM ya existe. Saltando mapeo."
else
  echo "🚀 minimap2 (splice:hq)…"
  minimap2 -t "$THREADS" -ax splice:hq -uf -k14 --secondary=no "$GENOME_FA" "$QUERY" > "$SAM"

  echo "🧹 SAM -> BAM ordenado + index…"
  samtools view -@ "$THREADS" -bS "$SAM" | samtools sort -@ "$THREADS" -o "$BAM" -
  samtools index "$BAM"
  rm -f "$SAM"
fi

# -------------------------------
# 4.1) Colapsar a la MEJOR alineación por qname
#       Criterios: AS desc, refspan desc, NM asc; excluir suplementarias
# -------------------------------
BEST_BAM="lncrna_vs_GRCh38.best.bam"
if [[ -s "$BEST_BAM" && -s "${BEST_BAM}.bai" ]]; then
  echo "✅ BEST_BAM ya existe. Saltando selección."
else
  echo "🧠 Seleccionando mejor alineación por qname (AS↓, refspan↓, NM↑)…"
  hdr="$(mktemp)"; body="$(mktemp)"; outsam="$(mktemp)"

  # Cabeceras
  samtools view -H "$BAM" > "$hdr"

  # Candidatos: excluir suplementarias (-F 0x800); secundarios ya estaban desactivados en minimap2
  # Elegir por qname usando AS, NM y refspan (derivado de CIGAR)
  samtools view -F 0x800 "$BAM" \
  | awk '
    BEGIN{OFS="\t"}
    function ref_span(cigar,   s,op,len,span,tmp) {
      span=0; tmp=cigar;
      while (match(tmp,/([0-9]+)([MIDNSHP=X])/)) {
        len=substr(tmp,RSTART,RLENGTH-1)+0;
        op =substr(tmp,RSTART+length(len),1);
        if (op=="M" || op=="D" || op=="N" || op=="=" || op=="X") span+=len;
        tmp=substr(tmp,RSTART+RLENGTH);
      }
      return span;
    }
    {
      qname=$1; cigar=$6;
      as=nm="";
      for(i=12;i<=NF;i++){
        if($i ~ /^AS:i:/){split($i,a,":"); as=a[3]+0}
        else if($i ~ /^NM:i:/){split($i,a,":"); nm=a[3]+0}
      }
      if(as=="") as=-999999;  # por si faltara (raro)
      if(nm=="") nm=999999;
      span=ref_span(cigar);

      # Guardar mejor por qname
      if(!(qname in seen) ||
         (as  > best_as[qname]) ||
         (as == best_as[qname] && span > best_span[qname]) ||
         (as == best_as[qname] && span == best_span[qname] && nm < best_nm[qname])) {
           best_as[qname]=as; best_span[qname]=span; best_nm[qname]=nm; rec[qname]=$0; seen[qname]=1;
      }
    }
    END{
      for(q in rec) print rec[q];
    }' > "$body"

  # Ensamblar SAM y convertir a BAM ordenado
  cat "$hdr" "$body" > "$outsam"
  samtools view -@ "$THREADS" -bS "$outsam" \
    | samtools sort -@ "$THREADS" -o "$BEST_BAM" -
  samtools index "$BEST_BAM"
  rm -f "$hdr" "$body" "$outsam"
fi

# -------------------------------
# 5) PAF (opcional, del mejor BAM)
#     Si quieres el PAF colapsado, re-emitimos desde el FASTA directamente
#     o lo omitimos. Aquí lo regenero sobre el conjunto completo (como antes).
# -------------------------------
if [[ -s "$PAF" ]]; then
  echo "✅ PAF ya existe."
else
  echo "🗂️  Escribiendo PAF (sin colapsar, informativo)…"
  minimap2 -t "$THREADS" -x splice:hq -uf -k14 --secondary=no "$GENOME_FA" "$QUERY" > "$PAF"
fi

# -------------------------------
# 6) BED12 y CSV (a partir del BAM colapsado)
# -------------------------------
BED12="lncrna_vs_GRCh38.best.bed12"
CSV="lncrna_vs_GRCh38.best.alignments.csv"

echo "🧭 BEST_BAM -> BED12…"
bedtools bamtobed -bed12 -i "$BEST_BAM" > "$BED12"

echo "📄 BED12 -> CSV…"
awk 'BEGIN{
  FS=OFS="\t";
  print "chrom,start,end,qname,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts"
}{ print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12 }' "$BED12" | sed 's/\t/,/g' > "$CSV"

echo "✅ Listo."
echo "BEST_BAM: $BEST_BAM"
echo "BED12:    $BED12"
echo "CSV:      $CSV"