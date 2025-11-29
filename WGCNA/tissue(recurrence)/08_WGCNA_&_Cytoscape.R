# CYTOSCAPE



# ------------------------------------------------------------------------------
# Cargar las LIBRERIAS necesarias
library(fields)
library(impute)
library(dynamicTreeCut)
library(qvalue)
library(flashClust)
library(Hmisc)
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(data.table)
library(dplyr)
library(ggplot2)

# =============================================================================
# Cargar los datos para el analisis
# Se necesita la matriz de conteo con todas las muestras
# Se necesita adicionalmente un archivo con la descripción de los fenotipos
setwd("./tissue2")

dir.create("./Resultados/WGCNA")
dir.create("./Resultados/WGCNA/Modulos")
dir.create("./Resultados/WGCNA/Expresion_Modulos")
dir.create("./Resultados/WGCNA/Cytoscape")
dir_analisis <- "./Resultados/WGCNA"
setwd(dir_analisis)
# En la carpeta WGCNA debe estar el archivo: 
#     "Matriz_conteo_align.csv"
#     "experimento.csv"
#     "traits.csv"
# El archivo "traits.csv" en el siguiente formato:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SampleID     Healthy    Cancer 
# muestra_1      1         0
# muestar_2      1         0
# muestra_3      0         1
# muestra_4      0         1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#
# ==============================================================================
# Cargar los datos
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read.csv("./Matriz_conteo_align.csv", header = T, row.names = 1)
data$Length <- NULL
traits <- read.csv("traits.csv", header = T, row.names = 1)
traits1 <- read.csv("experimento.csv", header = T, row.names = 1)
sum(is.na(data))

# ------------------------------------------------------------------------------
# Detectar genes que se comportan como outliers
gsg <- goodSamplesGenes(t(data), verbose = 3)
summary(gsg)
gsg$allOK

table(gsg$goodGenes) # genes outliers
table(gsg$goodSamples) # muestras outliers

data <- data[gsg$goodGenes == TRUE,] # eliminar genes outliers

# Detectar muestras outliers con arbol de jerarquia
# *************************************************
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# Detectar muestras outliers con PCA
# **********************************
pca <- prcomp(t(data))
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# ------------------------------------------------------------------------------
# Eliminar muestras outliers
samples.to.be.excluded <- c("SRR7949418.BAM","SRR7949463.BAM","SRR7949473.BAM","SRR7949450.BAM")
data <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# **CRÍTICO: Eliminar las mismas muestras de la tabla de metadatos**
# Esto es necesario para que el número de columnas de 'data' coincida con el número de filas de 'traits'
traits <- traits[!(rownames(traits) %in% samples.to.be.excluded), ]
# ------------------------------------------------------------------------------
# Establecer los datos para el DeSeq2 y filtrar la informacion
# ************************************************************
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = traits,
                              design = ~ 1)

# Eliminar genes con bajos conteos
# Conteos menores a 15 en mas del 75% de las muestras (muestras*0.75)
muestras <- ncol(data)
cuartil <- muestras*0.75
dds75 <- dds[rowSums(counts(dds) >= 15) >= cuartil,]

# Realizar estabilización de varianzas vst 
dds_norm <- vst(dds75)

# ------------------------------------------------------------------------------
# Obtener conteos normalizados y escribir los datos filtrados y normalizados:
norm.counts <- assay(dds_norm)

# Escribir el archivo de conteos normalizados filtrado
write.csv(norm.counts, "./data_wgcna_vst_norm.csv" , row.names = TRUE)

# ==============================================================================
# Si se cuenta con la tabla normalizada y filtrada
# se puede iniciar desde este paso:

# 1.- cargar los datos o seguir en el paso 2.- 
data_wgcna<- read.csv('./data_wgcna_vst_norm.csv', header = TRUE, row.names = 1)
head(data_wgcna)
# 2.- los datos estan en forma de matrix, convertirlos a dataframe
data_wgcna = data.frame (norm.counts)

# ------------------------------------------------------------------------------
# Generar graficos de agrupamiento y distribucion de las muestras normalizadas
# y filtradas
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

htree <- hclust(dist(t(data_wgcna)), method = "average")
plot(htree)

# Checar comportamiento de las muestras con PCA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pca <- prcomp(t(data_wgcna))
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)

Condition <- traits1$Condition

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point(aes(fill = Condition), shape = 21, size=6, stroke=1) +
  geom_text(label = rownames(pca.dat), size = 3) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


#  Graficar la expresion normalizada de las muestras
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_wgcna_plot <- data_wgcna %>%
  mutate(Gene_id = row.names(data_wgcna)) %>%
  pivot_longer(-Gene_id)

data_wgcna_plot %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text( angle = 90)) +
  ylim(0, NA) +
  labs(
    title = "Normalized and Filtered WGCNA Expression Data",
    x = "Muestras",
    y = "Normalized Expression"
  )

# ------------------------------------------------------------------------------
# *******************      Analisis de WGCNA    ********************************
# Para el analisis de WGCNA los datos deben de estar transpuestos
# en filas las muestras y en las columnas los genes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# transponer los datos:
data_wgcna_t <- data_wgcna %>% t()

# Determinar el valor a usar de softthreshold
allowWGCNAThreads(20)
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

#funcion para llamar topologia de red
sft <- pickSoftThreshold(data_wgcna_t,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
sft.data <- sft$fitIndices

# visualizacion para elegir el valor de Beta (poder)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# asignar el poder beta de acuerdo a la gráfica
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
soft_power <- 9

# convierte la matriz de cuentas normalizadas a valores numéricos
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_wgcna_t[] <- sapply(data_wgcna_t, as.numeric)

# Modifica la funcion de cor (correlacion) por la de WGCNA
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cor <- WGCNA::cor

# Calcular la Matriz de topologia "topological overlap matrix" (TOM)
# unsigned -> nodes with positive & negative correlation are treated equally 
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero
bwnet <- blockwiseModules(data_wgcna_t,
                          maxBlockSize = 20000, # de acuerdo a la compu 4G ~ 5000
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

# escribir el resultado para no tener que realizarlo cada vez
# cargar el archivo de resultados en la variable bwnet
saveRDS(bwnet, file = "./bwnet.Data")
bwnet <- readRDS ("bwnet.Data")

# Devolver la funcion cor a su original
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cor<- stats::cor

# ------------------------------------------------------------------------------
# Evaluar los Modulos de eigengenes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MEs <- bwnet$MEs
head(MEs)
str(MEs)

# Para ver cuantos genes tiene cada módulo
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TablaModulos <- data.frame(table(bwnet$colors))
colnames(TablaModulos) <- c("Modulo","No.Genes")
write.csv(TablaModulos, "./Tabla_modulos.csv" , row.names = TRUE)

# Generar el dendograma
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(bwnet$unmergedColors)
length(bwnet$colors)
length(bwnet$dendrograms[[1]]$order)
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# Definir número de genes y muestras
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nSamples <- nrow(data_wgcna_t)
nGenes <- ncol(data_wgcna_t)

# ------------------------------------------------------------------------------
# Genera las tablas de correlacion de cada modulo por fenotipo
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Mod_Trait_corr <- as.data.frame (cor(MEs, traits, use = 'p'))
Mod_Trait_pvals <- as.data.frame(corPvalueStudent(as.matrix(Mod_Trait_corr), nSamples))

# visualizar asociacion modulo-fenotipo con heatmap
# Definir adecuadamente los valos de x y de y en las columnas de heatmap.data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
heatmap.data <- merge(MEs, traits, by = 'row.names')
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
str(heatmap.data)
# Modificar los valores de "x" y de "y" de acuerdo a la tabla heatmap.data  
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[21:22],  # Muestra los fenotipos (Healthy y Cancer)
             y = names(heatmap.data)[1:20],  # Muestra todos los módulos (del 1 al 22)
             col = c("blue1", "skyblue", "white", "pink", "red"))

# -----------------------------------------------------------------------------
# Extraer y anotar los genes de un Modulo determinado
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cargar los datos para anotación 
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# Extraer los datos y depositarlos en una data frame
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene_mod <- as.data.frame(bwnet$colors)
gene_mod <- data.frame(
  "ENSEMBL_ID" = gene_mod %>% rownames(),
  "Modulo" = gene_mod$`bwnet$colors`)

temp_ann <- getBM(filters= "ensembl_gene_id_version", 
                  attributes= c("ensembl_gene_id_version",
                                "external_gene_name",
                                "description"),
                  values=gene_mod$ENSEMBL_ID, mart= mart)

# Fusiona las tablas, el numero de genes disminuye ya que hay algunos que no se conocen
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene_mod_ann <-merge(gene_mod,temp_ann,by.x="ENSEMBL_ID",by.y="ensembl_gene_id_version")

# Escribir las tablas con los genes de cada módulo
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mod_colors <- substring(names(bwnet$MEs), 3)

for (val in mod_colors)
{
  nombre1 <-paste("./Modulos/Modulo_",val,".csv", sep = "", collapse = NULL)
  genes_por_modulo <-  gene_mod_ann %>%
    filter(Modulo == val)
  
  write.csv(genes_por_modulo, nombre1, row.names = TRUE)
}

# ------------------------------------------------------------------------------
# Identificar genes altamente relacionados con un Modulo o Fenotipo
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calcular la correlacion de cada gen con el módulo (pearson) y los valores p asociados
Gen_Mod_Membership_cor = as.data.frame (cor(data_wgcna_t, MEs, use = "p"))
Gen_Mod_Membership_pvals = as.data.frame (corPvalueStudent(as.matrix(Gen_Mod_Membership_cor), nSamples))
names(Gen_Mod_Membership_cor) = paste("c.GMM_", mod_colors, sep="")
names(Gen_Mod_Membership_pvals) = paste("p.GMM_", mod_colors, sep="")

# Calcular la correlación de cada gen con el fenotipo (pearson) y los valores p asociados
Gen_Trait_cor = as.data.frame(cor(data_wgcna_t, traits, use = "p"))
Gen_Trait_pvals = as.data.frame(corPvalueStudent(as.matrix(Gen_Trait_cor), nSamples))
names(Gen_Trait_cor) = paste("c.GT_", names(traits), sep="")
names(Gen_Trait_pvals) = paste("p.GTS_", names(traits), sep="")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Escribir y annotar las tablas con los valores de p de cada gen en relación con
# su pertenencia a cada módulo o fenotipo, las tablas se puede abrir en excel
# y ordenar de menor a mayor para seleccionar los genes con menor valor de "p"
# estos son los genes mas significativos para cada módulo o fenotipo.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Para Genes y Modulos
temp_ann2 <- merge(Gen_Mod_Membership_pvals, temp_ann, by.x=0, by.y="ensembl_gene_id_version")
write.csv (temp_ann2, "./Genes_Modules_Membership_pvals.csv", row.names = TRUE)

# Para Genes y Fenotipo
temp_ann2 <- merge(Gen_Trait_pvals, temp_ann, by.x=0, by.y="ensembl_gene_id_version")
write.csv (temp_ann2, "./Genes_Traits_Membership_pvals.csv", row.names = TRUE)

# Evaluar la relación de los fenotipos y los modulos de manera gráfica
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MET = orderMEs(cbind(MEs, traits))

# Plot the dendrogram
plotEigengeneNetworks(MET, "Modules and Traits Dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,1,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

# ------------------------------------------------------------------------------
# Calcular la significancia de cada gen en relación con un fenotipo
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Para Control(Healthy):
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene.signf.corr.Control <- cor(data_wgcna_t, traits$Control, use = 'p')
colnames(gene.signf.corr.Control) <- c('Per_cor')
gene.signf.corr.pvals.Control <- corPvalueStudent(gene.signf.corr.Control, nSamples)
colnames(gene.signf.corr.pvals.Control) <- c('pvals')
gene.signf.Control <- cbind(gene.signf.corr.Control, gene.signf.corr.pvals.Control)
gene.signf.Control <- data.frame(gene.signf.Control)
write.csv(gene.signf.Control, "./tissue2_Control_Assoc_genes_total.csv" , row.names = TRUE)

setDT(gene.signf.Control, keep.rownames = "ENSEMBL_ID")  
gene.signf.Control_anno <- merge(gene.signf.Control,temp_ann,by.x= "ENSEMBL_ID", by.y="ensembl_gene_id_version")
write.csv(gene.signf.Control_anno, "./tissue2_Control_Assoc_genes_anno.csv" , row.names = TRUE)

# Para Paciente:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gene.signf.corr.Case <- cor(data_wgcna_t, traits$Case, use = 'p')
colnames(gene.signf.corr.Case) <- c('Per_cor')
gene.signf.corr.pvals.Case <- corPvalueStudent(gene.signf.corr.Case, nSamples)
colnames(gene.signf.corr.pvals.Case) <- c('pvals')
gene.signf.Case <- cbind(gene.signf.corr.Case, gene.signf.corr.pvals.Case)
gene.signf.Case <- data.frame(gene.signf.Case)
write.csv(gene.signf.Case, "./tissue2_Case_Assoc_genes_total.csv", row.names = TRUE)

setDT(gene.signf.Case, keep.rownames = "ENSEMBL_ID")  
gene.signf.Case_anno <- merge(gene.signf.Case, temp_ann, by.x= "ENSEMBL_ID", by.y="ensembl_gene_id_version")
write.csv(gene.signf.Case_anno, "./tissue2_Case_Assoc_genes_anno.csv" , row.names = TRUE)

# ------------------------------------------------------------------------------
# Seleccionar modulos de interes y graficar su expresion
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modules_of_interest = c("magenta",
                        "pink")

# Generar una lista con los genes de los modulos seleccionados
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
submod = gene_mod %>% subset(Modulo %in% modules_of_interest)
row.names(submod) = submod$ENSEMBL_ID

# Extraer la expresion normalizada de los genes de los modulos seleccionados
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
expr_norm_sel_mod = data_wgcna[submod$ENSEMBL_ID,]

expresion_df = data.frame(expr_norm_sel_mod) %>%
  mutate(ENSEMBL_ID = row.names(.)) %>%
  pivot_longer(-ENSEMBL_ID) %>%
  mutate(module = submod[ENSEMBL_ID,]$Modulo)

# generar el grafico
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
expresion_df %>% ggplot(., aes(x=name, y=value, group=ENSEMBL_ID)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Fenotipo",
       y = "normalized expression")

# ==============================================================================
# Generar tablas con los datos de expresión normalizada para los genes de cada
# Modulo.
# ------------------------------------------------------------------------------

for (val in mod_colors)
{
  nombre1 <-paste("./Expresion_Modulos/Exp_Data_Modulo_",val,"_all_samples.csv", sep = "", collapse = NULL)
  genes_of_interest = gene_mod %>% subset(Modulo %in% val)
  expr_genes_of_interest = data_wgcna[genes_of_interest$ENSEMBL_ID,]
  expr_genes_of_interest <- data.frame (expr_genes_of_interest)
  write.csv(expr_genes_of_interest, nombre1, row.names = TRUE)
}

# ------------------------------------------------------------------------------
# Crear la lista de genes de los Modulos para creación de Redes en Cytoscape
# 1.- Crear la matriz de topología para todos los genes
# topological overlap matrix (TOM) 

TOM <- TOMsimilarityFromExpr(data_wgcna_t, power = soft_power)
row.names(TOM) = row.names(data_wgcna)
colnames(TOM) = row.names(data_wgcna)

dissTOM = 1-TOM

# Extraer los genes de los diferentes modulos y asignarles valores para crear
# las redes y las conecciones

for (i in 1:length(bwnet$MEs))
{
  modulo = c(substring(names(bwnet$MEs)[i], 3))
  nombre_modulo = paste("./Modulos/Modulo_",modulo,".csv", sep = "", collapse = NULL)
  print(nombre_modulo)
  modulo_temp = read.csv(nombre_modulo)
  genes = colnames(data_wgcna_t)
  modGenes = c(modulo_temp$ENSEMBL_ID)
  modGenes2 = c(modulo_temp$external_gene_name)
  modTOM=TOM[modGenes,modGenes]
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("./Cytoscape/Cytoscape_edges_", modulo, ".txt", sep=""),
                                 nodeFile = paste("./Cytoscape/Cytoscape_nodes_", modulo, ".txt", sep=""),
                                 weighted = TRUE, threshold = 0.1, nodeNames = modGenes2, nodeAttr = modulo);
}

# =====================================================================================
# Guardar todo el trabajo hasta aqui
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save.image(file="./wgcna_work_data.RData")
load(file="./wgcna_work_data.RDataTmp") # <- para cargar el trabajo
setwd("./tissue2/Resultados/WGCNA"
# =====================================================================================
# *******************************  Cytoscape   ****************************************
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#BiocManager::install("RCy3")
library(RCy3)
# Asegurarse de que Cytoscape esta corriendo
cytoscapePing ()
cytoscapeVersionInfo ()

# ------------------------------------------------------------------------------
# Cargar los archivos del modulo que se desea analizar
# ------------------------------------------------------------------------------
Modulo_Interes <- "cyan"
Cytofile1 <- paste("./Cytoscape/Cytoscape_edges_",Modulo_Interes,".txt", sep = "", collapse = NULL)
Cytofile2 <- paste("./Cytoscape/Cytoscape_nodes_",Modulo_Interes,".txt", sep = "", collapse = NULL)

edge <- read.delim(Cytofile1)
colnames(edge)
colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName")

node <- read.delim(Cytofile2)
colnames(node)  
colnames(node) <- c("id","altName","node_attributes") 

# Numero de interacciones a mostrar, cambiar el valor si es necesario
Vertices = round(length(edge$source)/5)

createNetworkFromDataFrames(node,edge[1:Vertices,], title=Modulo_Interes, collection="DataFrame Example")




# Calculate connectivity within the module
kWithin <- intramodularConnectivity(adjacency, moduleColors)$kWithin

# Find top hub gene in your module
topHub <- names(sort(kWithin[moduleColors == "yourModuleColor"], decreasing = TRUE))[1]
































#COMBINACION DE MODULOS RELEVANTES ****(IGNORAR)****

# 1. Crear la matriz de adyacencia si no está creada
adjacency <- adjacency(data_wgcna_t, power = soft_power, type = "signed")

# 2. Asegurar que la matriz tenga nombres de genes
rownames(adjacency) <- colnames(data_wgcna_t)
colnames(adjacency) <- colnames(data_wgcna_t)

# 3. Obtener los colores de los módulos
moduleColors <- bwnet$colors

# 4. Convertir moduleColors a data frame
module_df <- data.frame(ENSEMBL_ID = names(moduleColors),
                        color = moduleColors)

# 5. Unir con anotaciones para obtener ENSEMBL IDs
module_df_anno <- merge(module_df, temp_ann, by.x = "ENSEMBL_ID", by.y = "external_gene_name")

# 6. Filtrar genes de los módulos combinados
modulos_combinados <- c("greenyellow", "magenta")
genes_combinados <- module_df_anno$ensembl_gene_id_version[
  module_df_anno$color %in% modulos_combinados
]

# 7. Verificar que los genes existan en la matriz de adyacencia
genes_combinados <- intersect(genes_combinados, rownames(adjacency))

# 8. Crear la submatriz de adyacencia para esos genes
adjacency_comb <- adjacency[genes_combinados, genes_combinados]

# 9. Asegurar que la submatriz tenga nombres
rownames(adjacency_comb) <- genes_combinados
colnames(adjacency_comb) <- genes_combinados

# 10. Exportar la red combinada a Cytoscape
exportNetworkToCytoscape(adjacency_comb,
                         edgeFile = "./Cytoscape/Cytoscape_edges_combined.txt",
                         nodeFile = "./Cytoscape/Cytoscape_nodes_combined.txt",
                         weighted = TRUE,
                         threshold = 0.01,
                         nodeNames = genes_combinados,
                         altNodeNames = genes_combinados,
                         nodeAttr = moduleColors[match(genes_combinados, names(moduleColors))])

# 11. Visualizar en Cytoscape
Modulo_Interes <- "combined"
Cytofile1 <- "./Cytoscape/Cytoscape_edges_combined.txt"
Cytofile2 <- "./Cytoscape/Cytoscape_nodes_combined.txt"

edge <- read.delim(Cytofile1)
colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName")

node <- read.delim(Cytofile2)
colnames(node) <- c("id","altName","node_attributes") 

Vertices <- round(length(edge$source)/5)

createNetworkFromDataFrames(node, edge[1:Vertices,], title = Modulo_Interes, collection = "DataFrame Example")





# Check if the submatrix has any non-zero values
summary(as.vector(adjacency_comb))

# Optional: visualize the distribution of weights
hist(as.vector(adjacency_comb), breaks = 50, main = "Connectivity values", xlab = "Weight")
