#!/bin/bash
# ============================================================
# validation.sh - Validación de FASTA candidatos lncRNA
# Ejecutar desde prepro/pbmcs
# ============================================================
set -e
set -o pipefail

CONS_DIR="consensus_transcripts"
INPUT="$CONS_DIR/final_candidates_for_lncDC.fasta"
REPORT="$CONS_DIR/validation_report.txt"

echo "Validando $INPUT ..."
{
    echo "Validando $INPUT ..."
    echo "Reporte generado en: $REPORT"
    echo "==================================================="

    # Contar número de secuencias
    num_seqs=$(grep -c "^>" "$INPUT")
    echo "🔢 Número de secuencias: $num_seqs"

    # Detectar líneas sin encabezado
    echo "⚠️ Encabezados inválidos encontrados (líneas sin '>'):"
    awk '($0 !~ /^>/ && $0 !~ /^[ACGTNacgtn]+$/) {print NR ":" $0}' "$INPUT"

    # Verificar que todas las secuencias contengan solo ACGTN
    echo
    echo "⚠️ Secuencias con caracteres no estándar:"
    grep -n -E -v "^[ACGTNacgtn>]" "$INPUT" || echo "Ninguna encontrada ✅"

} | tee "$REPORT"

echo "✅ Validación finalizada. Revisa $REPORT"
