#!/usr/bin/env bash
set -euo pipefail

# ================================
# lncRNA functional annotation (Ensembl r115 / Gencode v46)
# - Activates conda env "aln_env" and ensures bedtools present
# - Downloads GTF to annotation/ (resume-safe) and uses it
# - Deduplicates to ONE label per lncRNA with precedence:
#   antisense > sense_overlapping > sense_intronic > bidirectional > intergenic
# ================================

# ---- Config (override via env or CLI) ----
ENV_NAME="${ENV_NAME:-aln_env}"
LNC_BED12="${1:-lncrna_vs_GRCh38.best.bed12}"

OUT_DIR="annotation"
mkdir -p "$OUT_DIR"

GTF_URL="${GTF_URL:-https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz}"
GTF_GZ="${OUT_DIR}/$(basename "$GTF_URL")"
GTF="${OUT_DIR}/Homo_sapiens.GRCh38.115.gtf"

PROM_WIN="${PROM_WIN:-1000}"
OUT="${OUT_DIR}/lncrna.annotated.tsv"

# -------------------------------
# 0) Conda: activate env and ensure tools
# -------------------------------
ensure_conda() {
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
  else
    echo "❌ conda not found in PATH. Please load conda and re-run." >&2
    exit 1
  fi
}
ensure_env() {
  if conda env list | grep -qE "^\s*${ENV_NAME}\s"; then
    echo "✅ Using conda env: ${ENV_NAME}"
  else
    echo "⚙️  Creating env '${ENV_NAME}' with bedtools + basics…"
    conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda bedtools wget pigz coreutils curl gawk
  fi
  conda activate "$ENV_NAME"
}
ensure_tools() {
  local need=()
  for t in bedtools awk sort join wget curl pigz; do
    command -v "$t" >/dev/null 2>&1 || need+=("$t")
  done
  if (( ${#need[@]} )); then
    echo "📦 Installing missing tools in '${ENV_NAME}': ${need[*]}"
    local pkgs=()
    for t in "${need[@]}"; do
      case "$t" in
        awk) pkgs+=("gawk") ;;
        sort|join) pkgs+=("coreutils") ;;
        wget) pkgs+=("wget") ;;
        curl) pkgs+=("curl") ;;
        pigz) pkgs+=("pigz") ;;
        bedtools) pkgs+=("bedtools") ;;
      esac
    done
    conda install -y -n "$ENV_NAME" -c conda-forge -c bioconda "${pkgs[@]}"
  fi
}
ensure_conda; ensure_env; ensure_tools

# -------------------------------
# 1) Inputs
# -------------------------------
[[ -s "$LNC_BED12" ]] || { echo "❌ Missing or empty: $LNC_BED12"; exit 1; }

# -------------------------------
# 2) Get GTF into annotation/ (resume-safe)
# -------------------------------
if [[ ! -s "$GTF" ]]; then
  echo "🌐 Fetching GTF into ${OUT_DIR}/ …"
  if [[ ! -s "$GTF_GZ" ]]; then
    if command -v wget >/dev/null 2>&1; then
      wget -c "$GTF_URL" -O "$GTF_GZ"
    else
      curl -L -C - -o "$GTF_GZ" "$GTF_URL"
    fi
  else
    echo "✅ Found compressed GTF: $GTF_GZ"
  fi
  echo "🗜️  Decompressing GTF…"
  tmp="${GTF}.tmp"
  if command -v pigz >/dev/null 2>&1; then pigz -dc "$GTF_GZ" > "$tmp"; else gunzip -c "$GTF_GZ" > "$tmp"; fi
  mv -f "$tmp" "$GTF"
else
  echo "✅ Using existing GTF: $GTF"
fi

# -------------------------------
# 3) Helpers (GTF attribute parsing)
# -------------------------------
attr_awk='
function geta(s,k,  r,m) {
  r = k"[[:space:]]+\"([^\"]+)\"";
  if (match(s, r, m)) return m[1];
  return "";
}
'

# -------------------------------
# 4) Build BEDs: protein-coding genes/exons → introns → promoters
# -------------------------------
GENE_BED="${OUT_DIR}/genes.pc.bed6"
EXON_BED="${OUT_DIR}/exons.pc.bed6"
INTRON_BED="${OUT_DIR}/introns.pc.bed6"
PROM_BED="${OUT_DIR}/promoters.pc.bed6"

if [[ ! -s "$GENE_BED" || ! -s "$EXON_BED" ]]; then
  echo "🧩 Parsing GTF → BED6 (protein_coding genes/exons)…"
  LC_ALL=C awk -F'\t' -v OFS='\t' "$attr_awk"'
    $0 ~ /^#/ {next}
    {
      type=$3; chr=$1; start=$4-1; end=$5; strand=$7; attrs=$9;
      gid=geta(attrs,"gene_id"); gname=geta(attrs,"gene_name");
      gtype=geta(attrs,"gene_biotype"); if (gtype=="") gtype=geta(attrs,"gene_type");
      if (gtype!="protein_coding") next;

      if (type=="gene") {
        print chr, start, end, gid "|" gname, 0, strand > "'"$GENE_BED"'"
      }
      else if (type=="exon") {
        tid=geta(attrs,"transcript_id");
        print chr, start, end, gid "|" gname "|" tid, 0, strand > "'"$EXON_BED"'"
      }
    }' "$GTF"
  for b in "$GENE_BED" "$EXON_BED"; do LC_ALL=C sort -k1,1 -k2,2n -o "$b" "$b"; done
fi

if [[ ! -s "$INTRON_BED" ]]; then
  echo "🧩 Computing introns = genes - exons (strand-aware)…"
  bedtools subtract -s -a "$GENE_BED" -b "$EXON_BED" \
    | awk '($3>$2)' \
    | LC_ALL=C sort -k1,1 -k2,2n -o "$INTRON_BED"
fi

if [[ ! -s "$PROM_BED" ]]; then
  echo "🧩 Building ±${PROM_WIN} bp promoter windows around TSS…"
  awk -v OFS='\t' -v W="$PROM_WIN" '
    {
      chr=$1; start=$2; end=$3; name=$4; score=$5; strand=$6;
      if (strand=="+") { tss=start; p1=tss-W; p2=tss+W }
      else             { tss=end;   p1=tss-W; p2=tss+W }
      if (p1<0) p1=0;
      print chr, p1, p2, name, score, strand
    }' "$GENE_BED" \
    | LC_ALL=C sort -k1,1 -k2,2n -o "$PROM_BED"
fi

# -------------------------------
# 5) Collapse lncRNA BED12 → BED6
# -------------------------------
LNC_BED6="${OUT_DIR}/lnc.bed6"
if [[ ! -s "$LNC_BED6" ]]; then
  echo "🧭 Collapsing lncRNA BED12 → BED6 (span)…"
  awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6}' "$LNC_BED12" \
    | LC_ALL=C sort -k1,1 -k2,2n -o "$LNC_BED6"
fi

# -------------------------------
# 6) Overlaps (sense/antisense/promoter)
# -------------------------------
ANTI="${OUT_DIR}/lnc_vs_gene.opposite.tsv"
SENSE_EXON="${OUT_DIR}/lnc_vs_exon.same.tsv"
SENSE_INTRON="${OUT_DIR}/lnc_vs_intron.same.tsv"
BIDI_RAW="${OUT_DIR}/lnc_vs_prom.opposite.tsv"

echo "🔎 Intersecting overlaps…"
bedtools intersect -s  -a "$LNC_BED6" -b "$EXON_BED"   -wa -wb > "$SENSE_EXON"
bedtools intersect -s  -a "$LNC_BED6" -b "$INTRON_BED" -wa -wb > "$SENSE_INTRON"
bedtools intersect -S  -a "$LNC_BED6" -b "$GENE_BED"   -wa -wb > "$ANTI"
bedtools intersect -S  -a "$LNC_BED6" -b "$PROM_BED"   -wa -wb > "$BIDI_RAW"

# -------------------------------
# 7) Build labeled tables (BEST per lncRNA within each class)
# -------------------------------
best_by_overlap() {
  # usage: best_by_overlap CLASS TARGET < input > output
  local CLASS="$1" TARGET="$2"
  awk -v OFS="\t" -v CLASS="$CLASS" -v TARGET="$TARGET" '
    function max(a,b){return a>b?a:b}
    function min(a,b){return a<b?a:b}
    {
      lchr=$1; ls=$2; le=$3; lname=$4; lscore=$5; lstrand=$6;
      rchr=$7; rs=$8; re=$9; rname=$10; rscore=$11; rstrand=$12;

      # overlap length
      ov = min(le,re) - max(ls,rs); if (ov<0) ov=0;

      # parse rname: gid|gname|tid? (tid may be absent)
      n=split(rname, parts, /\|/);
      gid=parts[1]; gname=(n>=2?parts[2]:""); tid=(n>=3?parts[3]:"");

      gspan = (re-rs);
      key = lname;

      # keep best by: ov desc, gspan desc, gid asc
      if( !(key in seen) ||
          (ov    >  best_ov[key]) ||
          (ov   ==  best_ov[key] && gspan >  best_span[key]) ||
          (ov   ==  best_ov[key] && gspan == best_span[key] && gid < best_gid[key]) ) {
        best_ov[key]=ov; best_span[key]=gspan; best_gid[key]=gid;
        rec[key]=lname"\t"CLASS"\t"lchr"\t"ls"\t"le"\t"lstrand"\t"gid"\t"gname"\t"tid"\t"TARGET;
        seen[key]=1;
      }
    }
    END{ for(k in rec) print rec[k]; }'
}

tmp_dir="$(mktemp -d)"; trap 'rm -rf "$tmp_dir"' EXIT

echo "🏷️  Building labeled tables (best per lncRNA/class)…"
best_by_overlap "antisense"         "gene"     < "$ANTI"       | sort -u > "$tmp_dir/antisense.tsv"
best_by_overlap "sense_overlapping" "exon"     < "$SENSE_EXON" | sort -u > "$tmp_dir/sense_overlapping.tsv"

# Intronic: exclude anything that already overlapped exon (same strand)
cut -f1 "$tmp_dir/sense_overlapping.tsv" | sort -u > "$tmp_dir/_has_exon.lst" || true
awk 'NR==FNR{e[$1]=1;next} !($4 in e)' "$tmp_dir/_has_exon.lst" "$SENSE_INTRON" \
| best_by_overlap "sense_intronic" "intron" | sort -u > "$tmp_dir/sense_intronic.tsv" || true

# Bidirectional: opposite-strand promoters without any gene/exon/intron overlap
cat "$tmp_dir/antisense.tsv" "$tmp_dir/sense_overlapping.tsv" "$tmp_dir/sense_intronic.tsv" 2>/dev/null \
| awk '{print $1}' | sort -u > "$tmp_dir/_has_geneoverlap.lnc" || true
awk 'NR==FNR{h[$1]=1;next} {if(!( $4 in h)) print $0}' "$tmp_dir/_has_geneoverlap.lnc" "$BIDI_RAW" \
| best_by_overlap "bidirectional" "promoter" | sort -u > "$tmp_dir/bidirectional.tsv" || true

# -------------------------------
# 8) Merge with precedence (no subshell; truly global marked[])
# -------------------------------
prec=("antisense" "sense_overlapping" "sense_intronic" "bidirectional")
echo -e "lncRNA\tclass\tchrom\tstart\tend\tstrand\tgene_id\tgene_name\ttranscript_id\ttarget_feature" > "$OUT"

declare -A marked
for cls in "${prec[@]}"; do
  f="$tmp_dir/${cls}.tsv"; [[ -s "$f" ]] || continue
  while IFS=$'\t' read -r lname c chr s e st gid gname tid tgt; do
    [[ -z "${marked[$lname]:-}" ]] || continue
    echo -e "${lname}\t${c}\t${chr}\t${s}\t${e}\t${st}\t${gid}\t${gname}\t${tid}\t${tgt}" >> "$OUT"
    marked["$lname"]=1
  done < "$f"
done

# Intergenic = not assigned yet
cut -f4 "$LNC_BED6" | sort -u > "$tmp_dir/all.lnc"
cut -f1 "$OUT" | tail -n +2 | sort -u > "$tmp_dir/assigned.lnc" || true
comm -23 "$tmp_dir/all.lnc" "$tmp_dir/assigned.lnc" \
| awk -v OFS="\t" '{print $1,"intergenic",".",".",".",".",".",".",".","."}' >> "$OUT"

# -------------------------------
# 9) Summary
# -------------------------------
echo "✅ Wrote: $OUT"
echo "🔢 Class counts:"
tail -n +2 "$OUT" | cut -f2 | sort | uniq -c | awk '{printf "  %-18s %d\n",$2,$1}'

echo "🧪 QC files saved in: $OUT_DIR"
ls -1 "$OUT_DIR" | grep -E 'lnc_vs_|genes\.pc|exons\.pc|introns\.pc|promoters\.pc|lnc\.bed6' || true
