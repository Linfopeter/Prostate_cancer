#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# --- Config ---
INPUT="${1:-Resultado_Final.TOP50_consensus.csv}"
OUTPUT="${2:-annotated_TOP50.csv}"   # <-- nombre por defecto solicitado
GTF="Homo_sapiens.GRCh38.115.gtf.gz"
URL="https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz"
DOWNLOADED=0

# --- Helpers (multi-core where possible) ---
downloader() {
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x16 -s16 -k1M -o "$GTF" "$URL"
  else
    wget -c "$URL" -O "$GTF"
  fi
}
decompress() {
  if [[ "$GTF" == *.gz ]]; then
    if command -v pigz >/dev/null 2>&1; then pigz -dc -- "$GTF"
    else zcat -f -- "$GTF"; fi
  else
    cat -- "$GTF"
  fi
}

# --- 0) Download GTF if missing ---
if [[ ! -f "$GTF" ]]; then
  echo "🌐 Downloading Ensembl GTF (r115)…"
  downloader
  DOWNLOADED=1
  echo "✅ Downloaded: $GTF"
else
  echo "ℹ️ Using existing $GTF"
fi

# --- 1) Check CSV ---
if [[ ! -f "$INPUT" ]]; then
  echo "❌ CSV not found: $INPUT" >&2; exit 1
fi
if grep -q ',"[^"]\+,[^"]\+"' "$INPUT"; then
  echo "⚠️  Detected quoted fields with commas in $INPUT. This simple parser assumes no embedded commas." >&2
fi

echo "🧩 Mapping ENST → (ENSG,gene_name) and writing $OUTPUT …"

# --- 2) Build annotated CSV ---
decompress | gawk -F'\t' -v OFS=',' -v CSV="$INPUT" '
  # PASS 1: read GTF -> maps
  {
    if ($0 ~ /gene_id "[^"]+"/ && $0 ~ /transcript_id "[^"]+"/) {
      match($0, /gene_id "([^"]+)"/, gid_m)
      match($0, /transcript_id "([^"]+)"/, tid_m)
      gnm_found = match($0, /gene_name "([^"]+)"/, gnm_m)
      if (tid_m[1] != "") {
        base = tid_m[1]; sub(/\.[0-9]+$/, "", base)
        if (!(base in GID)) { GID[base]=gid_m[1]; GNM[base]=(gnm_found?gnm_m[1]:"NA") }
      }
    }
  }
  END {
    # PASS 2: CSV -> insert gene_id,gene_name after col 2
    while ((getline line < CSV) > 0) {
      nf = split(line, A, ",")
      if (NRcsv == 0) {
        NRcsv = 1
        out = ""
        for (i=1; i<=nf; i++) {
          out = out ((i==1)?"":OFS) A[i]
          if (i==2) out = out OFS "gene_id" OFS "gene_name"
        }
        print out; continue
      }
      t=A[2]; sub(/\.[0-9]+$/, "", t)
      gid=(t in GID?GID[t]:"NA"); gnm=(t in GNM?GNM[t]:"NA")
      out = A[1] OFS A[2] OFS gid OFS gnm
      for (i=3; i<=nf; i++) out = out OFS A[i]
      print out
    }
    close(CSV)
  }
' > "$OUTPUT"

# --- 3) Quick validation ---
rows_in=$(($(wc -l < "$INPUT")))
rows_out=$(($(wc -l < "$OUTPUT")))
echo "📏 Rows in: $rows_in ; Rows out: $rows_out"
if [[ "$rows_in" -ne "$rows_out" ]]; then
  echo "⚠️  Row count changed; check for embedded commas/quoting in CSV." >&2
fi

# --- 4) Cleanup GTF if we downloaded it now ---
if [[ "$DOWNLOADED" -eq 1 ]]; then
  rm -f -- "$GTF"
  echo "🧹 Removed downloaded file: $GTF"
fi

echo "✅ Done: $OUTPUT"
