# ==============================================================================
#      Análisis de Expresión Diferencial con DESeq2
#      (Flujo de Trabajo Unificado y Mejorado)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Cargar las librerías necesarias
library(DESeq2)
library(apeglm)
library(pheatmap)
library(vsn)
library(biomaRt)
library(ReportingTools)
library(lattice)
library(RColorBrewer)
library(readr)
library(ggplot2)

# Limpiar las variables del sistema
rm(list = ls())

# ==============================================================================
# 1. Configuración de Rutas y Carga de Datos
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 🚨 Ajuste de la ruta del proyecto.
# Asume que el script se ejecuta desde la raíz del proyecto.
dir_proyecto <- "./tissue1"

# Establece el directorio de trabajo una sola vez
setwd(dir_proyecto)

# Construye las rutas de manera más simple y directa
dir_analisis <- "." # Ya estamos en dir_proyecto
dir_resultados <- file.path(dir_analisis, "Resultados")
dir_deseq_out <- file.path(dir_resultados, "DESeq2_align")
dir_counts_out <- file.path(dir_resultados, "counts", "align")


message("Cargando la matriz de conteo generada por featureCounts...")

# La ruta completa al archivo se construye a partir de la raíz del proyecto.
path_to_counts <- file.path(dir_counts_out, 'Matriz_conteo_align.csv')
dsdata <- read.csv(path_to_counts, header = TRUE, row.names = 1)

# Cargar la información del experimento.
path_to_exp <- file.path(dir_analisis, 'experimento.csv')
experimento <- read.csv(path_to_exp, header = TRUE, row.names = 1)

# ==============================================================================
# 2. Preparación y Filtrado de Datos
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Eliminar la columna 'Length'
if ("Length" %in% colnames(dsdata)) {
  message("Columna 'Length' detectada. Se eliminará antes del análisis DESeq2.")
  dsdata_counts <- subset(dsdata, select = -c(Length))
} else {
  dsdata_counts <- dsdata
}

# Verificación de la coherencia entre las muestras.
message("Verificando la coherencia entre las muestras de la matriz de conteo y la de experimento...")
if (!all(colnames(dsdata_counts) %in% rownames(experimento)) || !all(colnames(dsdata_counts) == rownames(experimento))) {
  message("Se reordenará la matriz de conteo para que coincida con la matriz de experimento.")
  dsdata_counts <- dsdata_counts[, rownames(experimento)]
}
all(colnames(dsdata_counts) == rownames(experimento))

# Establecer la variable de condición como factor
experimento$Condition <- factor(experimento$Condition)

# Creación del objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(dsdata_counts),
                              colData = experimento,
                              design = ~ Condition)

# Especificar el grupo de referencia.
dds$Condition <- relevel(dds$Condition, ref = "Healthy subjects")
message(paste("El grupo de referencia para el análisis es:", levels(dds$Condition)[1]))

# Filtrado de features: Mantener solo genes con suficientes conteos.
Nmuestras = 30
keep <- rowSums(counts(dds)) >= (Nmuestras * 2)
dds <- dds[keep,]
message(paste("Se han filtrado genes con conteos bajos. El número de genes restantes es:", nrow(dds)))

# ==============================================================================
# 3. Análisis de Expresión Diferencial con DESeq2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Ejecutando DESeq2...")
dds <- DESeq(dds)
resultsNames(dds)

# Guardar el objeto 'dds' para análisis futuros.
saveRDS(dds, file = file.path(dir_deseq_out, "./Res_DS2_dds_Align.rds"))
message("Análisis DESeq2 completado y el objeto 'dds' se ha guardado.")

# ==============================================================================
# 4. Normalización, Transformación y Control de Calidad
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Iniciando normalización y transformación de datos para visualización...")

# Normalización y transformación de los datos para visualización
dds_lg2 <- normTransform(dds)
dds_vst <- vst(dds, blind = FALSE)
dds_rlog <- rlog(dds, blind = FALSE)

# Exportar la matriz de datos normalizada para uso futuro
norm_dds <- counts(dds, normalized=TRUE)
write.csv(norm_dds, file.path(dir_deseq_out, "./Data_DS2_align_norm.csv"), row.names = TRUE)

# ----------------- Visualización y Control de Calidad -----------------

message("Generando gráficos de control de calidad...")

# Gráfico de la distribución de distancias de Cook
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, main = "Distribución de Cook's distances")

# Gráficos de Varianza vs. Media (Mean-SD Plots)
meanSdPlot(assay(dds_lg2), main = "Mean-SD Plot (log2)")
meanSdPlot(assay(dds_vst), main = "Mean-SD Plot (VST)")
meanSdPlot(assay(dds_rlog), main = "Mean-SD Plot (rlog)")

# Gráficos de PCA (Análisis de Componentes Principales)
pca_rlog <- plotPCA(dds_rlog, intgroup = "Condition", ntop = 500, returnData = TRUE)
print(ggplot(pca_rlog, aes(PC1, PC2, color=Condition)) +
        geom_point(size=3) +
        ggtitle("PCA con rlog") +
        theme_bw())

pca_vst <- plotPCA(dds_vst, intgroup = "Condition", ntop = 500, returnData = TRUE)
print(ggplot(pca_vst, aes(PC1, PC2, color=Condition)) +
        geom_point(size=3) +
        ggtitle("PCA con VST") +
        theme_bw())

# Heatmap de distancia entre las muestras
sampleDists <- dist(t(assay(dds_vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(dds_vst), sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main = "Heatmap de Distancia de Muestras (VST)")

