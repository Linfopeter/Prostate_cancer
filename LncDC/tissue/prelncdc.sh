#!/bin/bash
# =================================================================================
# LncDC Preprocessing Pipeline Script v3.6 (Portable + Idempotente + Robustez)
# =================================================================================

set -euo pipefail

# --- Variables Definidas por el Usuario ---
PROJECT_DIR=$(pwd)
RAW_DATA_TUMOR_DIR="$PROJECT_DIR/rawdata/cancer"
RAW_DATA_NORMAL_DIR="$PROJECT_DIR/rawdata/control"

PARENT_DIR=$(dirname "$PROJECT_DIR")
REF_GENOME="$PARENT_DIR/RefGen/GRCh38.primary_assembly.genome.fa"
REF_GTF="$PARENT_DIR/RefGen/gencode.v49.primary_assembly.annotation.gtf"
STAR_INDEX_DIR="$PARENT_DIR/RefGen/star_index"
HISAT2_INDEX_DIR="$PARENT_DIR/RefGen/hisat2_index"
NONCODE_DB_FASTA="$PARENT_DIR/RefGen/NONCODEv6_human.fa"   # archivo FASTA
NONCODE_DB_PREFIX="$PARENT_DIR/RefGen/NONCODEv6_human"     # prefijo DB si se crea

NUM_THREADS=30
FPKM_CUTOFF=0.2   # <<< Ajusta si es necesario (ej. 0.5, 0.7, 1.0)

# --- Directorios de Salida ---
CLEANED_TUMOR_DIR="$PROJECT_DIR/cleaned_reads/tumor"
CLEANED_NORMAL_DIR="$PROJECT_DIR/cleaned_reads/normal"
ROUTE1_ASSEMBLIES_TUMOR_DIR="$PROJECT_DIR/route1_assemblies/tumor"
ROUTE1_ASSEMBLIES_NORMAL_DIR="$PROJECT_DIR/route1_assemblies/normal"
ROUTE2_ASSEMBLIES_TUMOR_DIR="$PROJECT_DIR/route2_assemblies/tumor"
ROUTE2_ASSEMBLIES_NORMAL_DIR="$PROJECT_DIR/route2_assemblies/normal"
STAR_OUTPUT_DIR="$PROJECT_DIR/star_output"
HISAT2_OUTPUT_DIR="$PROJECT_DIR/hisat2_output"
CONSENSUS_DIR="$PROJECT_DIR/consensus_transcripts"

# --- Verificaciones Iniciales y Creación de Directorios ---
mkdir -p "$CLEANED_TUMOR_DIR" "$CLEANED_NORMAL_DIR"
mkdir -p "$ROUTE1_ASSEMBLIES_TUMOR_DIR" "$ROUTE1_ASSEMBLIES_NORMAL_DIR"
mkdir -p "$ROUTE2_ASSEMBLIES_TUMOR_DIR" "$ROUTE2_ASSEMBLIES_NORMAL_DIR"
mkdir -p "$STAR_OUTPUT_DIR" "$HISAT2_OUTPUT_DIR"
mkdir -p "$CONSENSUS_DIR"

cd "$PROJECT_DIR"

# --- Filtro FPKM (parametrizado) ---
fpkm_filter_logic=$(cat <<EOF
BEGIN{FS="\t"; OFS="\t"}
FNR==NR {
    if (\$3=="transcript") {
        match(\$9, /FPKM "([^"]+)"/, arr);
        if (arr[1] > $FPKM_CUTOFF) {
            match(\$9, /transcript_id "([^"]+)"/, tid);
            good_tids[tid[1]]=1;
        }
    }
    next;
}
{
    match(\$9, /transcript_id "([^"]+)"/, tid);
    if (tid[1] in good_tids) {
        print \$0;
    }
}
EOF
)

# =================================================================================
# Paso 1: Control de Calidad y Limpieza
# =================================================================================
echo "--- Paso 1: Realizando control de calidad y limpieza... 🧹 ---"
for sample_pair_1 in "$RAW_DATA_TUMOR_DIR"/*_1.fq.gz; do
    sample_name=$(basename "$sample_pair_1" "_1.fq.gz")
    CLEAN_FILE="$CLEANED_TUMOR_DIR/${sample_name}_1_val_1.fq.gz"
    if [ ! -f "$CLEAN_FILE" ]; then
        trim_galore --paired --fastqc --cores "$NUM_THREADS" \
                    --output_dir "$CLEANED_TUMOR_DIR" \
                    "$sample_pair_1" "$RAW_DATA_TUMOR_DIR/${sample_name}_2.fq.gz"
    fi
done
for sample_pair_1 in "$RAW_DATA_NORMAL_DIR"/*_1.fq.gz; do
    sample_name=$(basename "$sample_pair_1" "_1.fq.gz")
    CLEAN_FILE="$CLEANED_NORMAL_DIR/${sample_name}_1_val_1.fq.gz"
    if [ ! -f "$CLEAN_FILE" ]; then
        trim_galore --paired --fastqc --cores "$NUM_THREADS" \
                    --output_dir "$CLEANED_NORMAL_DIR" \
                    "$sample_pair_1" "$RAW_DATA_NORMAL_DIR/${sample_name}_2.fq.gz"
    fi
done

# =================================================================================
# Paso 2: Ruta 1 - STAR y Cufflinks + Merge (robusto, no toca otros GTF)
# =================================================================================
echo "--- Paso 2: Ejecutando Ruta 1 (STAR + Cufflinks)... 🚀 ---"
for sample_type in tumor normal; do
    CLEANED_DIR="$PROJECT_DIR/cleaned_reads/$sample_type"
    for sample_pair_1 in "$CLEANED_DIR"/*_1_val_1.fq.gz; do
        sample_name=$(basename "$sample_pair_1" "_1_val_1.fq.gz")
        BAM_FILE="$STAR_OUTPUT_DIR/${sample_type}_${sample_name}_Aligned.sortedByCoord.out.bam"
        CUFFLINKS_OUT_DIR="$PROJECT_DIR/route1_assemblies/${sample_type}/${sample_name}"
        GTF_FILTRADO="$CUFFLINKS_OUT_DIR/transcripts.filtered.gtf"
        if [ ! -f "$GTF_FILTRADO" ]; then
            mkdir -p "$CUFFLINKS_OUT_DIR"
            if [ ! -f "$BAM_FILE" ]; then
                STAR --runThreadN "$NUM_THREADS" --genomeDir "$STAR_INDEX_DIR" \
                     --readFilesIn "$sample_pair_1" "$CLEANED_DIR/${sample_name}_2_val_2.fq.gz" \
                     --readFilesCommand zcat \
                     --outFileNamePrefix "$STAR_OUTPUT_DIR/${sample_type}_${sample_name}_" \
                     --outSAMtype BAM SortedByCoordinate \
                     --outSAMstrandField intronMotif
            fi
            cufflinks -p "$NUM_THREADS" -o "$CUFFLINKS_OUT_DIR" "$BAM_FILE"
            gawk "$fpkm_filter_logic" \
                "$CUFFLINKS_OUT_DIR/transcripts.gtf" \
                "$CUFFLINKS_OUT_DIR/transcripts.gtf" > "$GTF_FILTRADO"
        fi
    done
done

# --- Merge Ruta 1: Tumor (conserva si ya existe) ---
if [ ! -f "$CONSENSUS_DIR/route1_tumor_gtf_list.txt" ]; then
    find "$ROUTE1_ASSEMBLIES_TUMOR_DIR" -name "transcripts.filtered.gtf" -size +1k \
        | sort > "$CONSENSUS_DIR/route1_tumor_gtf_list.txt"
fi
if [ ! -f "$CONSENSUS_DIR/merged_route1_tumor.gtf" ]; then
    cuffmerge -p "$NUM_THREADS" -o "$CONSENSUS_DIR" \
              -g "$REF_GTF" -s "$REF_GENOME" \
              "$CONSENSUS_DIR/route1_tumor_gtf_list.txt"
    mv "$CONSENSUS_DIR/merged.gtf" "$CONSENSUS_DIR/merged_route1_tumor.gtf"
fi

# --- Merge Ruta 1: Normal (robusto; NO toca los otros GTF) ---
if [ ! -s "$CONSENSUS_DIR/route1_normal_gtf_list.txt" ]; then
    find "$ROUTE1_ASSEMBLIES_NORMAL_DIR" -name "transcripts.filtered.gtf" -size +1k \
        | sort > "$CONSENSUS_DIR/route1_normal_gtf_list.txt"
fi

if [ ! -s "$CONSENSUS_DIR/merged_route1_normal.gtf" ]; then
    echo "⚠️  Generando merged_route1_normal.gtf (merge limpio y validado)..."
    TMP_MERGE_DIR="$CONSENSUS_DIR/tmp_normal_merge"
    rm -rf "$TMP_MERGE_DIR"
    cuffmerge -p "$NUM_THREADS" -o "$TMP_MERGE_DIR" \
              -g "$REF_GTF" -s "$REF_GENOME" \
              "$CONSENSUS_DIR/route1_normal_gtf_list.txt"
    if [ -s "$TMP_MERGE_DIR/merged.gtf" ]; then
        mv "$TMP_MERGE_DIR/merged.gtf" "$CONSENSUS_DIR/merged_route1_normal.gtf"
        echo "✅ merged_route1_normal.gtf generado correctamente."
    else
        echo "❌ Error: merged.gtf vacío para normales (Ruta 1). Revisa $TMP_MERGE_DIR"
        # No salimos; dejamos que el pipeline continúe (process_candidates saltará normal si falta)
    fi
    rm -rf "$TMP_MERGE_DIR"
else
    echo "✅ merged_route1_normal.gtf ya existe y no está vacío. Se conserva."
fi

# =================================================================================
# Paso 3: Ruta 2 - HISAT2 y StringTie + Merge (robusto)
# =================================================================================
echo "--- Paso 3: Ejecutando Ruta 2 (HISAT2 + StringTie)... 🚀 ---"
for sample_type in tumor normal; do
    CLEANED_DIR="$PROJECT_DIR/cleaned_reads/$sample_type"
    for sample_pair_1 in "$CLEANED_DIR"/*_1_val_1.fq.gz; do
        sample_name=$(basename "$sample_pair_1" "_1_val_1.fq.gz")
        SORTED_BAM="$HISAT2_OUTPUT_DIR/${sample_type}_${sample_name}_sorted.bam"
        GTF_FILTRADO="$PROJECT_DIR/route2_assemblies/${sample_type}/${sample_name}.filtered.gtf"
        if [ ! -f "$GTF_FILTRADO" ]; then
            mkdir -p "$PROJECT_DIR/route2_assemblies/${sample_type}"
            if [ ! -f "$SORTED_BAM" ]; then
                hisat2 -p "$NUM_THREADS" --dta \
                       --known-splicesite-infile "$HISAT2_INDEX_DIR/splicesites.txt" \
                       -x "$HISAT2_INDEX_DIR/genome_index" \
                       -1 "$sample_pair_1" -2 "$CLEANED_DIR/${sample_name}_2_val_2.fq.gz" | \
                       samtools sort -@ "$NUM_THREADS" -o "$SORTED_BAM"
                samtools index "$SORTED_BAM"
            fi
            STRINGTIE_OUT_GTF="$PROJECT_DIR/route2_assemblies/${sample_type}/${sample_name}.gtf"
            stringtie -p "$NUM_THREADS" -G "$REF_GTF" -o "$STRINGTIE_OUT_GTF" "$SORTED_BAM"
            gawk "$fpkm_filter_logic" \
                "$STRINGTIE_OUT_GTF" "$STRINGTIE_OUT_GTF" > "$GTF_FILTRADO"
        fi
    done
done

# --- Merge Ruta 2: Tumor (robusto) ---
if [ ! -f "$CONSENSUS_DIR/merged_route2_tumor.gtf" ]; then
    echo "Fusionando ensamblajes de tumor (Ruta 2) — excluyendo GTFs vacíos..."
    gtf_list_file="$CONSENSUS_DIR/route2_tumor_gtf_list.txt"
    rm -f "$gtf_list_file"
    for f in $(find "$ROUTE2_ASSEMBLIES_TUMOR_DIR" -name "*.filtered.gtf" -size +1k | sort); do
        if grep -qP '\ttranscript\t' "$f"; then
            echo "$f" >> "$gtf_list_file"
        else
            echo "⚠️ Muestra sin transcritos, excluida: $f"
        fi
    done
    if [ -s "$gtf_list_file" ]; then
        stringtie --merge -p "$NUM_THREADS" -G "$REF_GTF" \
            -o "$CONSENSUS_DIR/merged_route2_tumor.gtf" "$gtf_list_file"
    else
        echo "⚠️ Ningún ensamblaje válido encontrado para tumor (Ruta 2)."
    fi
fi

# --- Merge Ruta 2: Normal (robusto) ---
if [ ! -f "$CONSENSUS_DIR/merged_route2_normal.gtf" ]; then
    echo "Fusionando ensamblajes normales (Ruta 2) — excluyendo GTFs vacíos..."
    gtf_list_file="$CONSENSUS_DIR/route2_normal_gtf_list.txt"
    rm -f "$gtf_list_file"
    for f in $(find "$ROUTE2_ASSEMBLIES_NORMAL_DIR" -name "*.filtered.gtf" -size +1k | sort); do
        if grep -qP '\ttranscript\t' "$f"; then
            echo "$f" >> "$gtf_list_file"
        else
            echo "⚠️ Muestra sin transcritos, excluida: $f"
        fi
    done
    if [ -s "$gtf_list_file" ]; then
        stringtie --merge -p "$NUM_THREADS" -G "$REF_GTF" \
            -o "$CONSENSUS_DIR/merged_route2_normal.gtf" "$gtf_list_file"
    else
        echo "⚠️ Ningún ensamblaje válido encontrado para normal (Ruta 2)."
    fi
fi

# =================================================================================
# Paso 4–7: Generar candidatos (tumor y normal) con tolerancia a vacíos
# =================================================================================
process_candidates () {
    sample_type=$1  # tumor | normal
    echo "--- Paso 4–7: Procesando ensamblajes de $sample_type ---"

    # Comprobar que existan y no estén vacíos los merges
    if [ ! -s "$CONSENSUS_DIR/merged_route1_${sample_type}.gtf" ] || \
       [ ! -s "$CONSENSUS_DIR/merged_route2_${sample_type}.gtf" ]; then
        echo "⚠️  Faltan merges no vacíos para ${sample_type}. Se omite generación de candidatos."
        return
    fi

    # Paso 4: Comparación de ensamblajes con gffcompare (blindado)
    if ! gffcompare -r "$REF_GTF" \
         -o "$CONSENSUS_DIR/gffcmp_r1_vs_r2_${sample_type}" \
         "$CONSENSUS_DIR/merged_route1_${sample_type}.gtf" \
         "$CONSENSUS_DIR/merged_route2_${sample_type}.gtf"; then
        echo "⚠️  gffcompare falló para ${sample_type}. Se omite."
        return
    fi

    COMBINED="$CONSENSUS_DIR/gffcmp_r1_vs_r2_${sample_type}.combined.gtf"
    if [ ! -s "$COMBINED" ]; then
        echo "❌ gffcompare no generó combined.gtf para ${sample_type}. Se omite."
        return
    fi

    # Paso 5: Filtro de novedad (class_code u/i/o/x)
    awk 'BEGIN{FS=OFS="\t"}
         $3=="transcript" && $9 ~ /class_code "[uiox]"/ {print $0}
         $3=="exon"       && $9 ~ /transcript_id/       {print $0}' \
         "$COMBINED" > "$CONSENSUS_DIR/step_class_novel_${sample_type}.gtf"

    # Paso 6: Filtro de longitud + multi-exón
    awk 'BEGIN{FS=OFS="\t"}
         $3=="transcript" {
             match($9,/transcript_id "([^"]+)"/,a);
             len=$5-$4+1;
             if(a[1]!="" && len>=200 && len<=100000) keep[a[1]]=1;
         }
         $3=="exon" {
             match($9,/transcript_id "([^"]+)"/,a);
             if(a[1]!="") exon_count[a[1]]++;
         }
         END {
             for(t in exon_count)
                 if(exon_count[t] >= 2 && keep[t]==1)
                     print t;
         }' "$CONSENSUS_DIR/step_class_novel_${sample_type}.gtf" \
         > "$CONSENSUS_DIR/novel_multi_ids_${sample_type}.txt"

    echo "IDs de transcritos candidatos ($sample_type): $(wc -l < "$CONSENSUS_DIR/novel_multi_ids_${sample_type}.txt")"

    # Paso 7: Extraer GTF y FASTA
    awk 'BEGIN{FS=OFS="\t"}
         NR==FNR {ok[$1]=1; next}
         {match($9,/transcript_id "([^"]+)"/,a);
          if(a[1] in ok) print $0}' \
         "$CONSENSUS_DIR/novel_multi_ids_${sample_type}.txt" \
         "$CONSENSUS_DIR/step_class_novel_${sample_type}.gtf" \
         > "$CONSENSUS_DIR/novel_${sample_type}_candidates.gtf"

    gffread "$CONSENSUS_DIR/novel_${sample_type}_candidates.gtf" \
      -g "$REF_GENOME" \
      -w "$CONSENSUS_DIR/novel_${sample_type}_candidates.multi_exon.fasta"

    echo "Secuencias en FASTA ($sample_type): $(grep -c '^>' "$CONSENSUS_DIR/novel_${sample_type}_candidates.multi_exon.fasta")"

    # Limpieza intermedios
    rm -f "$CONSENSUS_DIR/step_class_novel_${sample_type}.gtf" \
          "$CONSENSUS_DIR/novel_multi_ids_${sample_type}.txt"
}

# ⚠️ IMPORTANTE: Invocar la función para generar los FASTA de tumor y normal
process_candidates "tumor"
process_candidates "normal"

# =================================================================================
# Paso 8: Filtrado contra NONCODE (optimizado si hay makeblastdb)
# =================================================================================
echo "🧬 Filtrando candidatos tumor contra NONCODE..."
NONCODE_DB_ARG=()
if command -v makeblastdb >/dev/null 2>&1; then
    # Crear DB si no existe (formato BLAST: .nhr/.nin/.nsq)
    if [ ! -f "${NONCODE_DB_PREFIX}.nhr" ]; then
        makeblastdb -in "$NONCODE_DB_FASTA" -dbtype nucl -out "$NONCODE_DB_PREFIX"
    fi
    NONCODE_DB_ARG=(-db "$NONCODE_DB_PREFIX" -num_threads "$NUM_THREADS")
else
    # Fallback a -subject (nota: -num_threads se ignora con -subject)
    NONCODE_DB_ARG=(-subject "$NONCODE_DB_FASTA")
fi

if [ -s "$CONSENSUS_DIR/novel_tumor_candidates.multi_exon.fasta" ]; then
    blastn -query "$CONSENSUS_DIR/novel_tumor_candidates.multi_exon.fasta" \
           "${NONCODE_DB_ARG[@]}" \
           -out "$CONSENSUS_DIR/noncode_hits.txt" \
           -evalue 1e-5 -perc_identity 95 -outfmt 6
    grep -v -F -f <(cut -f1 "$CONSENSUS_DIR/noncode_hits.txt") \
        "$CONSENSUS_DIR/novel_tumor_candidates.multi_exon.fasta" \
        > "$CONSENSUS_DIR/novel_tumor_noNONCODE.fasta"
    echo "Secuencias tumor tras NONCODE: $(grep -c '^>' "$CONSENSUS_DIR/novel_tumor_noNONCODE.fasta")"
else
    echo "⚠️  No existe FASTA de candidatos tumor. Se omite filtro NONCODE."
fi

# =================================================================================
# Paso 9: Filtrado contra controles normales
# =================================================================================
echo "🧬 Eliminando transcritos presentes en controles normales..."
if [ -s "$CONSENSUS_DIR/novel_tumor_noNONCODE.fasta" ] && \
   [ -s "$CONSENSUS_DIR/novel_normal_candidates.multi_exon.fasta" ]; then
    blastn -query "$CONSENSUS_DIR/novel_tumor_noNONCODE.fasta" \
           -subject "$CONSENSUS_DIR/novel_normal_candidates.multi_exon.fasta" \
           -out "$CONSENSUS_DIR/tumor_vs_normal_hits.txt" \
           -evalue 1e-5 -perc_identity 95 -outfmt 6
    grep -v -F -f <(cut -f1 "$CONSENSUS_DIR/tumor_vs_normal_hits.txt") \
        "$CONSENSUS_DIR/novel_tumor_noNONCODE.fasta" \
        > "$CONSENSUS_DIR/novel_tumor_unique.fasta"
    echo "Secuencias exclusivas tumor: $(grep -c '^>' "$CONSENSUS_DIR/novel_tumor_unique.fasta")"
else
    echo "⚠️  Faltan archivos para el filtro contra normales. Se omite Paso 9."
fi

# =================================================================================
# Paso 10: Redundancia con CD-HIT
# =================================================================================
echo "🧬 Eliminando redundancia con CD-HIT..."
if [ -s "$CONSENSUS_DIR/novel_tumor_unique.fasta" ]; then
    cd-hit-est -i "$CONSENSUS_DIR/novel_tumor_unique.fasta" \
               -o "$CONSENSUS_DIR/final_candidates_for_lncDC.fasta" \
               -c 0.95 -n 10 -M 16000 -T "$NUM_THREADS"
    echo "Secuencias finales (no redundantes): $(grep -c '^>' "$CONSENSUS_DIR/final_candidates_for_lncDC.fasta")"
else
    echo "⚠️  No hay FASTA de tumor único. Se omite CD-HIT."
fi

echo "✅ Pipeline finalizado."
