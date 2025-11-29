# ==============================================================================
#      Script Independiente para Anotación, Reporte y Visualización
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Cargar las librerías necesarias
library(DESeq2)
library(biomaRt)
library(ReportingTools)
library(apeglm)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)

## Limpiar el entorno de trabajo
rm(list = ls())

## Configuración de Rutas y Carga de Datos
# Asume que el script se encuentra en la misma estructura de carpetas que la original.
dir_proyecto <- "."
dir_analisis <- file.path(dir_proyecto, "tissue2")
dir_resultados <- file.path(dir_analisis, "Resultados")
dir_deseq_out <- file.path(dir_resultados, "DESeq2_align")
dir.create(dir_deseq_out, showWarnings = FALSE, recursive = TRUE)

message("Cargando el objeto DESeqDataSet (dds) y la matriz de experimento...")

# La ruta al objeto 'dds' que fue guardado en el script de análisis inicial.
path_to_dds <- file.path(dir_deseq_out, "./Res_DS2_dds_Align.rds")
if (!file.exists(path_to_dds)) {
  stop("Error: El archivo 'Res_DS2_dds_Align.rds' no se encontró en la ruta esperada.
    Asegúrese de ejecutar el script de análisis inicial primero.")
}
dds <- readRDS(path_to_dds)

# La ruta a la matriz 'experimento.csv'
path_to_exp <- file.path(dir_analisis, 'experimento.csv')
if (!file.exists(path_to_exp)) {
  stop("Error: El archivo 'experimento.csv' no se encontró en la ruta esperada.")
}
experimento <- read.csv(path_to_exp, header = TRUE, row.names = 1)
experimento$Condition <- factor(experimento$Condition)

# ==============================================================================
#      RESULTADOS DEGs y Anotación de DEGs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Extrayendo resultados de DESeq2 y anotando genes...")

# Extraer los resultados de expresión diferencial.
# El contraste se define para "Prostate cancer" vs. "Healthy subjects".
res <- results(dds, contrast=c("Condition", "Case", "Control"))
mcols(res)$description
summary(res)

# Anotar resultados DESeq2 (ENSEMBL IDs -> Gene Names)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Cargar la base de datos de anotación de biomaRt para humanos.
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

Res_DS2_align <- data.frame(res)
Res_DS2_align$Ensembl <- rownames(Res_DS2_align)

message("Consultando biomaRt para la anotación...")
Lista_align <- getBM(filters = "ensembl_gene_id_version",
                     attributes = c("ensembl_gene_id_version",
                                    "external_gene_name",
                                    "description"),
                     values = Res_DS2_align$Ensembl,
                     mart = mart)

res_align <- merge(Res_DS2_align, Lista_align, by.x = "Ensembl", by.y = "ensembl_gene_id_version", all.x = TRUE)
write.csv(res_align, file.path(dir_deseq_out, "./Res_DS2_align_annot.csv"), row.names = TRUE)
message("Resultados anotados guardados en 'Res_DS2_align_annot.csv'.")

# ==============================================================================
#      Reporte con ReportingTools
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Generando reporte HTML...")
dir_report <- file.path(dir_deseq_out, "reporte_align")
dir.create(dir_report, showWarnings = FALSE)

add.anns <- function(df, mart) {
  nm <- rownames(df)
  anns <- getBM(attributes = c("ensembl_gene_id_version", "external_gene_name", "description"),
                filters = "ensembl_gene_id_version",
                values = nm,
                mart = mart)
  anns <- anns[match(nm, anns[, 1]), ]
  colnames(anns) <- c("Ensembl", "Gene_Name", "Gene_Description")
  df <- cbind(anns, df[, 2:ncol(df)])
  rownames(df) <- nm
  df
}

des2Report <- HTMLReport(shortName = "Reporte_align",
                         title = "Expresion diferencial DESeq2",
                         reportDirectory = dir_report)

publish(dds,
        des2Report,
        .modifyDF = list(add.anns),
        mart = mart,
        pvalueCutoff = 0.05,
        lfc = 1.5,
        factor = experimento$Condition,
        n = 10000,
        reportDir = dir_report)
finish(des2Report)
message("Reporte HTML generado en '", dir_report, "'.")

# =============================================================================
#      Compactación de datos para mejor visualización y gráficos
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

resultsNames(dds)

# Shrinkage de los LogFC para visualización mejorada
resApg <- lfcShrink(dds, coef=2, type="apeglm")
resAsh <- lfcShrink(dds, coef=2, type="ashr")
resNor <- lfcShrink(dds, coef=2, type="normal")

# Gráficos MA Plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plotMA(res, ylim = c(-8,8), main = "MA Plot (Raw)")
abline(h=c(-1,1), col="red", lwd=1)

plotMA(resApg, ylim = c(-10,10), main = "MA Plot (apeglm)")
abline(h=c(-1,1), col="red", lwd=1)

plotMA(resAsh, ylim = c(-6,6), main = "MA Plot (ashr)")
abline(h=c(-1,1), col="red", lwd=1)

plotMA(resNor, ylim = c(-6,6), main = "MA Plot (Normal)")
abline(h=c(-1,1), col="red", lwd=1)

# Gráficos de Conteos de Genes Específicos
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Generando gráficos de conteos para genes clave...")

idx <- identify(res$baseMean, res$log2FoldChange)
if (length(idx) > 0) {
  message("Se han identificado los siguientes genes en el MA Plot: ", paste(rownames(res)[idx], collapse = ", "))
} else {
  message("No se seleccionaron genes en el MA Plot.")
}
rownames(res)[idx]

plotCounts(dds, gene = which.max(res$log2FoldChange), intgroup = "Condition", main = "Gen con Máximo LogFC")

plotCounts(dds, gene = which.min(res$padj), intgroup = "Condition", main = "Gen con Mínimo p.adj")

plotCounts(dds, gene = "ENSG00000145794.17", intgroup = "Condition", main = "Gen 'ENSG00000145794.17'")
