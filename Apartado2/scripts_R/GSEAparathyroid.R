rm(list = ls()) # R version 4.3.1 (2023-06-16)
install.packages("tidyverse", version = "2.0.0")
library(tidyverse) # tidyverse_2.0.0
library(DESeq2)
# --- Data ---
getwd()

# GENERACIÓN DE .RNK DISEÑO 1
# Leo los datos de cuentas y metadatos
rawcounts <- read.table("~/Documentos/Transcriptómica/transcriptomic-final-exercise/Apartado2/input/rawcounts.tsv", header = TRUE, sep = "\t", row.names = 1)
metadata <- read.table("~/Documentos/Transcriptómica/transcriptomic-final-exercise/Apartado2/input/metadata.tsv", header = TRUE, sep = "\t", row.names = 1)
# Establecemos como factores las variables
metadata$patient <- factor(metadata$patient)
metadata$agent <- factor(metadata$agent)
metadata$time <- factor(metadata$time)
# Elimino los datos de la muestra GSM913896 (patient=4, agent=Control, time= 24h)
# y me quedo solo con datos a 24h para hacer el contraste que me piden.
metadata_fil= metadata[-24,]
rawcounts_fil= rawcounts[,-24]
metadata_filtered <- metadata_fil[metadata_fil$time == "24h", ]
rawcounts_filtered <- rawcounts[, rownames(metadata_filtered)]
# Creo un objeto de DESeq utilizando las cuentas y metadatos filtrados y un diseño complejo ~patient + agent
dds3 <- DESeqDataSetFromMatrix(countData = rawcounts_filtered,
                               colData = metadata_filtered,
                               design = ~ patient + agent)
# Vemos el contenido. Matriz de 53160 filas y 11 columnas
dds3
# Filtro para quedarme con genes de al menos 10 cuentas
keep <- rowSums(counts(dds3)) >= 10 
dds3 <- dds3[keep, ]
# Hago una transformación con vst() útil para representaciones gráficas.
vsd <- vst(dds3, blind = TRUE)
# Análisis de expresión diferencial. Normalización, control de la dispersión
dds4 <- DESeq(dds3, test = "Wald")
resultsNames(dds4)
design(dds4)
res_DPN_vs_Control <- results(dds4, alpha = 0.05, contrast = c("agent", "DPN", "Control"))
summary(res_DPN_vs_Control)
res_DPN_vs_Control
# Shrunken LFC: Necessary for ranking!!! Also useful to compare LFC across 
# experiments. It is recommended to use these shrunken results when doing DEA.
# Please note that the p-values might vary a little.

### Low count genes have a lot of variance and thus their estimated LFCs can be 
### very extreme and imprecise.
### lfcShrink allows to shrunk the LFC of genes with low counts while conserving
### the LFC of genes with a high number of counts.
### apeglm is the default method, but there are others available:
### https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#moreshrink
resultsNames(dds4)
res_DPN_vs_Control24h.ape <- lfcShrink(dds4, coef = "agent_DPN_vs_Control", type = "apeglm",
                     res = res_DPN_vs_Control)
summary(res_DPN_vs_Control24h.ape) # 9 genes up y 3 down con padj < 0.05

# Visual explanation of shrunken LFCs

png("plots_DPN_vs_Control_ShrinkGSEA.png")
par(mfrow = c(1, 2))
plotMA(res_DPN_vs_Control, ylim = c(-3, 3))
plotMA(res_DPN_vs_Control24h.ape, ylim = c(-3, 3))
dev.off()
# Create .rnk. crea data frame con 2 columnas.feature es el gen ,lfc es el fold change
rnk1 <- data.frame(Feature = rownames(res_DPN_vs_Control24h.ape), 
                   LFC = res_DPN_vs_Control24h.ape$log2FoldChange)
head(rnk1)
#elimino todo lo que haya desde un . hasta el final
rnk1$Feature <- str_remove(rnk1$Feature, "\\..*$")
head(rnk1)

# Save .rnk (without header and tab separated)
write.table(rnk1, file = "GSEAparathyroid.rnk", sep = "\t", quote = FALSE, 
            col.names = FALSE, row.names = FALSE)





