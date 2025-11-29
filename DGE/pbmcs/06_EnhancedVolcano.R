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

## Preparación de datos y selección de genes para etiquetar
# Reemplazar valores 'NA' en la columna de nombres de genes para evitar errores
ds_volcano$external_gene_name <- as.character(ds_volcano$external_gene_name)
ds_volcano$external_gene_name[is.na(ds_volcano$external_gene_name)] <- ""



ds_volcano[ds_volcano$external_gene_name == "LGALS9", ]
ds_volcano[ds_volcano$external_gene_name == "PDCD1", ]
ds_volcano[ds_volcano$external_gene_name == "IL12A", ]
ds_volcano[ds_volcano$external_gene_name == "IFNG", ]
ds_volcano[ds_volcano$external_gene_name == "CD177", ]
ds_volcano[ds_volcano$external_gene_name == "MTOR", ]

# Definir la lista de genes que deseas etiquetar en el plot
genes_a_etiquetar <- c(
  "PTGS2", "MMP28", "GADD45A", "MPO",  "ARG1", 
  "MMP9", "FOSB", "ALPL",  "PPP1R15A", "LUCAT1", 
  "EGR1", "SCIRT", "EREG",  "JUN", "JUP", "MMP17", "MMP25", "LINC00200", "CXCL2",
  "CXCL1", "CXCL8", "CCL2", "TNF", "IL1A", "IL1B", "CCL20", "CCL3",
  "IL6", "PTGS2", "CXCL3", "IFNB1",  "HLA-DOB",
  "EGFR",   "VEGFA", "MMP9","PDCD1", "CD177","CXCR2", "CCL4"
  #TENTATIVOS "CCL4L2", "DDN", "ASIP", "THBD", "BAMBI","S100B", "S100A5", "RRAD",
  #"ATM","MYC","KRAS",  "BRAF", "TP53", "RB1","BRCA1","BRCA2","PTEN","STAT3","TGFB1","PARP1", 
  #"PD-L1","AR","SIRPA","CD47","CTLA4","IL10", "TP53", "IL12A","LGALS9" ,"CD274",  "IFNG",
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
                xlim = c(-5,10),
                ylim = c(0, max_y_value),
                title = "Volcano Plot: Prostate Cancer vs. Healthy Subjects",
                pCutoff = 0.05,
                FCcutoff = 1,
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

#"MMP9","FOSB","ALPL","RRAD","PPP1R15A","LUCAT1","DDN","BAMBI","EGR1","SCIRT","EREG","THBD","JUN","JUP","MMP17","MMP25"
#"CXCL1","CXCL8","CCL2","TNF","IL1A","IL1B","CCL20","CCL3","CCL4L2","IL6","COX2","CXCL3","","","",""






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
