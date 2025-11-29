# ==============================================================================
#      Script Independiente para Generar Volcano Plot
#      (Visualización en el panel de RStudio)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Cargar las librerías necesarias
library(dplyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

## Limpiar el entorno de trabajo
rm(list = ls())

## Configuración de Rutas y Carga de Datos
dir_proyecto <- "."
dir_analisis <- file.path(dir_proyecto, "tissue2")
dir_resultados <- file.path(dir_analisis, "Resultados")
dir_deseq_out <- file.path(dir_resultados, "DESeq2_align")

message("Cargando el archivo de resultados anotado...")
path_to_results <- file.path(dir_deseq_out, 'Res_DS2_align_annot.csv')
if (!file.exists(path_to_results)) {
  stop("Error: El archivo 'Res_DS2_align_annot.csv' no se encontró en la ruta esperada.
    Asegúrese de ejecutar el script de análisis y anotación primero.")
}
ds_volcano <- read.csv(path_to_results, header = TRUE, row.names = 1)

## Preparación de datos y selección de genes para etiquetar
# Reemplazar valores 'NA' en la columna de nombres de genes para evitar errores
ds_volcano$external_gene_name <- as.character(ds_volcano$external_gene_name)
ds_volcano$external_gene_name[is.na(ds_volcano$external_gene_name)] <- ""



ds_volcano[ds_volcano$external_gene_name == "CRISP3", ]


genes_a_etiquetar <- c(
  "CDH2","FIRRE","DAPL1","MKRN3","LRRTM3","REG1A",
  "LINC01235","ALOX12B","REG1A","PURPL","ADAM7",
  "PITX2","MUC6","MLC1","LINC01187","IGFN1",
  "CPB1","SEMG1","ADAMTS16","ATP6V0A4","LGSN",
  "LINC02571","LINC01612",
  #DOWNREGULATED
  "DQX1", "LCN2","ETV1","ATOH1","HBA2","GNG13","RLN1",
  "CST1","SLC9A3","PRIM2","HSPA6","RNA5SP202","MT2P1",
  "MT1P1",  "GC","ASIC5","GATA4","SPINK1","CEACAM20",
  "VSTM2A","SPON2","MT1G","PIGR","MT1E","MT1M","LCN15"
  #OVERLAPPED
  #"LRRTM3",
)

# Encontrar el valor máximo de -log10(pvalue) para ajustar dinámicamente el eje Y
max_y_value <- max(-log10(ds_volcano$pvalue), na.rm = TRUE) + 2

## Generar el gráfico de volcán y mostrarlo directamente
message("Generando el gráfico de volcán...")
x11()
EnhancedVolcano(ds_volcano,
                lab = ds_volcano$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = genes_a_etiquetar,
                xlim = c(min(ds_volcano$log2FoldChange, na.rm = TRUE) - 1.5, max(ds_volcano$log2FoldChange, na.rm = TRUE) + 1.5),
                ylim = c(0, 6),
                title = "Volcano Plot: Recurrence vs Non-Recurrence",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 50,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                legendLabels = c("NS", "Log2 FC", "p-value", "p-value & Log2 FC")
)

message("El Volcano Plot se ha visualizado en el panel de RStudio.")

#
#
#
#
#
#
#

















# Generar el gráfico de volcán y mostrarlo directamente
message("Generando el gráfico de volcán...")
x11()
EnhancedVolcano(ds_volcano,
                lab = ds_volcano$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                #selectLab = genes_a_etiquetar,
                # 🚨 Define the new limits here 🚨
                xlim = c(-4, 4),
                ylim = c(0, 4),
                title = "Volcano Plot: Prostate Cancer Tissue — Recurrence vs. Non-Recurrence ",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 25,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                legendLabels = c("NS", "Log2 FC", "p-value", "p-value & Log2 FC")
)

message("El Volcano Plot se ha visualizado en el panel de RStudio.")

#################################




#SPECIFIC GENE NUMBER
# 🔍 Filtrar los genes que cumplen con los criterios de significancia
genes_significativos <- ds_volcano[
  ds_volcano$pvalue < 0.01 & abs(ds_volcano$log2FoldChange) > 1,
]

genes_a_etiquetar <- head(genes_significativos[order(-abs(genes_significativos$log2FoldChange)), "external_gene_name"], 90)

# 📊 Graficar el Volcano Plot con todos los genes, pero solo etiquetar 150
EnhancedVolcano(ds_volcano,
                lab = ds_volcano$external_gene_name,
                selectLab = genes_a_etiquetar,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-4, 4),
                ylim = c(0, 4),
                title = "Volcano Plot: Etiquetando 150 genes más significativos",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 150,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                legendLabels = c("NS", "Log2 FC", "p-value", "p-value & Log2 FC")
)











































# ==============================================================================
#      Script Independiente para Generar Volcano Plot
#      (Visualización en el panel de RStudio)
#      Muestra SÓLO los genes subexpresados
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Cargar las librerías necesarias
library(dplyr)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

## Limpiar el entorno de trabajo
rm(list = ls())

## Configuración de Rutas y Carga de Datos
dir_proyecto <- "."
dir_analisis <- file.path(dir_proyecto, "pbmcs")
dir_resultados <- file.path(dir_analisis, "Resultados")
dir_deseq_out <- file.path(dir_resultados, "DESeq2_align")

message("Cargando el archivo de resultados anotado...")
path_to_results <- file.path(dir_deseq_out, 'Res_DS2_align_annot.csv')
if (!file.exists(path_to_results)) {
  stop("Error: El archivo 'Res_DS2_align_annot.csv' no se encontró en la ruta esperada.
    Asegúrese de ejecutar el script de análisis y anotación primero.")
}
ds_volcano <- read.csv(path_to_results, header = TRUE, row.names = 1)

## Filtrar el dataframe para mostrar solo genes subexpresados
ds_volcano_downregulated <- ds_volcano %>%
  dplyr::filter(log2FoldChange < 0)

## Preparación de datos
# Reemplazar valores 'NA' en la columna de nombres de genes
ds_volcano_downregulated$external_gene_name <- as.character(ds_volcano_downregulated$external_gene_name)
ds_volcano_downregulated$external_gene_name[is.na(ds_volcano_downregulated$external_gene_name)] <- ""

# Encontrar el valor máximo de -log10(pvalue) para ajustar dinámicamente el eje Y
max_y_value <- max(-log10(ds_volcano_downregulated$pvalue), na.rm = TRUE) + 2

## Generar el gráfico de volcán y mostrarlo directamente
message("Generando el gráfico de volcán...")
EnhancedVolcano(ds_volcano_downregulated,
                lab = ds_volcano_downregulated$external_gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NULL, # Esto asegura que no se etiquete ningún gen.
                xlim = c(min(ds_volcano_downregulated$log2FoldChange, na.rm = TRUE) - 1.5, 0.5),
                ylim = c(0, max_y_value),
                title = "Volcano Plot: Genes Subexpresados",
                pCutoff = 0.05,
                FCcutoff = 2.5,
                pointSize = 2,
                labSize = 3,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 160,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                legendLabels = c("NS", "Log2 FC", "p-value", "p-value & Log2 FC")
)

message("El Volcano Plot se ha visualizado en el panel de RStudio.")
