setwd("C:/Users/nmishra/SCB2005 Dropbox/Nitish Mishra/PC/Desktop/EMT-2022-August/Translational efficiency analysis/")

suppressMessages(suppressWarnings({
  # Global functions related code
  #source("../../jupyter_common.R")
  #source("../../function/enrichplot/gseaplot3.R")
  # Excel file related tools
  library("openxlsx"); library('readxl'); library('writexl'); library('data.table')
  # tidyverse related tools
  library('dplyr'); library('tidyverse'); library('stringr')
  # ggplot related tools
  library('ggplot2'); library('ggVennDiagram'); library('ggrepel'); library('ggpubr'); library('ggthemes')
  # Others tools
  
  # EMT/Ribosome/Translation related genes data
  #library('nichenetr')
  library('msigdbr'); library('gplots'); library('ComplexHeatmap')
  ##
  library('edgeR'); library('DESeq2')
  
}))

options(rstudio.help.showDataPreview = FALSE)

verb <- function(...) cat(sprintf(...), sep='', file=stdout())


#samples <- data.frame(
#  SampleID = c("JBF001", "JBF002", "JBF003", "JBF004", "JBF005", "JBF006", "JBF007", "JBF008", "JBF009"),
#  Condition = factor(c("unt", "unt", "unt", "tgfb", "tgfb", "tgfb", "tgfbCX", "tgfbCX", "tgfbCX"))
#)

# Create the design matrix
#design <- model.matrix(~ 0 + Condition, data=samples)
#colnames(design) <- levels(samples$Condition)
#print(design)
#

#lmat_pc <- lmat[rownames(lmat)%in%genes_symbol,]
#gene_lengths[gene_lengths$Gene %in% genes_symbol,]
#gene_lengths_pc <- gene_lengths[match(genes_symbol, gene_lengths$Gene), ]

################################################################################
### Globa; parameters
################################################################################
mtx_count = NULL 
min.expr.counts = 3
min.expr.log2cpm = -Inf
min.expr.num.samples = 2
ebayes.trend = TRUE 
ebayes.robust = TRUE 
MHadjust = "BH"
th_log2fc = log2(1.2)



gene_lengths <- read.table("Mus_musculus.GRCm39.104.rdna_rn18s_geneLength.txt", header=FALSE) |>
  mutate(Gene=V1, Transcript=V2, biotype=V3, Length=V4) |>
  select(c("Gene", "Transcript", "biotype", "Length")) |>
  #filter(biotype=="protein_coding") |>
  group_by(Gene) %>%
  filter(Length == max(Length)) %>%
  ungroup()


################################################################################
### Input RNAseq 
################################################################################

# Read the design matrix
design_matrix <- read.table("blancgrp_257956_RNAseq_polysome_input.design.txt", header=TRUE, row.names=1)
df <- design_matrix
df$SampleID <- rownames(df)
long_df <- pivot_longer(df, cols = -SampleID, names_to = "Condition", values_to = "value")
# Filter rows where value is 1
result_df <- long_df %>% filter(value == 1) %>% select(SampleID, Condition)
print(result_df)

### begin of addition by Nitish Mishra
# Exclused rDNA and mitochondial rRNA gene
total_RNAseq <-  fread("blancgrp_257956_RNAseq_polysome_input.counts.raw.txt", header = TRUE) |>
  filter(!grepl("^__|rDNA_promoter|rDNA|Rn18|mt-Rnr|Metazoa_SRP", V1)) |>
  column_to_rownames("V1") 


#### Filtering ####
mtx_count <- round(total_RNAseq)
mtx_count <- mtx_count[!grepl("^NA$", rownames(mtx_count)),]
mtx_count_inut <- mtx_count

# make DEGList for further analysis
myData.counts.dge <- DGEList( counts = mtx_count  ,  group = result_df$Condition  ,  genes = row.names(mtx_count) )

## remove genes with low expression
logi.good.count.expr  <-  rowSums(getCounts(myData.counts.dge) >= min.expr.counts)  >=  min.expr.num.samples
logi.good.log2cpm.expr  <-  rowSums(cpm(myData.counts.dge ) >= min.expr.log2cpm) >= min.expr.num.samples
logi.is.expr  <-  logi.good.log2cpm.expr  &  logi.good.count.expr


# remove genes with "__"
f_meta_tags <- grepl("^__", rownames(myData.counts.dge))
logi.is.expr <- logi.is.expr & !f_meta_tags
logi.is.nonexpr <- !logi.is.expr & !f_meta_tags
ngene <- length(logi.is.expr) - length(which(f_meta_tags))


myData.counts.filt <- myData.counts.dge[logi.is.expr,, keep.lib.sizes=FALSE]

dge <- calcNormFactors(myData.counts.filt)
cpm.total.RNAseq  <-  cpm( dge , normalized.lib.sizes=TRUE, log=FALSE )


## Gene length of genes
genes_symbol <- intersect(rownames(cpm.total.RNAseq), gene_lengths$Gene)
gene_lengths.selected <- gene_lengths[match(genes_symbol, gene_lengths$Gene), ]


myData.counts.filt <- myData.counts.filt[genes_symbol,]

## rpkm calculation
dge <- calcNormFactors(myData.counts.filt)
rpkm.total.RNAseq <- rpkm(dge, gene.length=gene_lengths.selected$Length) # rpkm for individual samples
rpkm_group.total.RNAseq <- rpkmByGroup(dge, gene.length=gene_lengths.selected$Length, group = design_matrix$Condition) # rpkm based on samples group

################################################################################
##### Read Polysome RNAseq data
################################################################################


# Read the design matrix
design_matrix <- read.table("blancgrp_257956_RNAseq_polysome_Fraction5_10.design.txt", header=TRUE, row.names=1)
df <- design_matrix
df$SampleID <- rownames(df)
long_df <- pivot_longer(df, cols = -SampleID, names_to = "Condition", values_to = "value")
# Filter rows where value is 1
result_df <- long_df %>% filter(value == 1) %>% select(SampleID, Condition)
print(result_df)

### begin of addition by Nitish Mishra
# Read epression matrix
# Exclused rDNA and mitochondial rRNA gene
total_Riboseq <-  fread("blancgrp_257956_RNAseq_polysome_Fraction5_10.counts.raw.txt", header = TRUE) |>
  filter(!grepl("^__|rDNA_promoter|rDNA|Rn18|mt-Rnr|Metazoa_SRP", V1)) |>
  column_to_rownames("V1") 


#### Filtering ####
mtx_count <- round(total_Riboseq)
mtx_count <- mtx_count[!grepl("^NA$", rownames(mtx_count)),]
mtx_count_polysome <- mtx_count

# make DEGList for further analysis
myData.counts.dge <- DGEList( counts = mtx_count  ,  group = result_df$Condition  ,  genes = row.names(mtx_count) )

## remove genes with low expression
logi.good.count.expr  <-  rowSums(getCounts(myData.counts.dge) >= min.expr.counts)  >=  min.expr.num.samples
logi.good.log2cpm.expr  <-  rowSums(cpm(myData.counts.dge ) >= min.expr.log2cpm) >= min.expr.num.samples
logi.is.expr  <-  logi.good.log2cpm.expr  &  logi.good.count.expr


# remove genes with "__"
f_meta_tags <- grepl("^__", rownames(myData.counts.dge))
logi.is.expr <- logi.is.expr & !f_meta_tags
logi.is.nonexpr <- !logi.is.expr & !f_meta_tags
ngene <- length(logi.is.expr) - length(which(f_meta_tags))


myData.counts.filt <- myData.counts.dge[logi.is.expr,, keep.lib.sizes=FALSE]

dge <- calcNormFactors(myData.counts.filt)
cpm.total.Riboseq  <-  cpm( dge , normalized.lib.sizes=TRUE, log=FALSE )


## Gene length of genes
genes_symbol <- intersect(rownames(cpm.total.Riboseq), gene_lengths$Gene)
gene_lengths.selected <- gene_lengths[match(genes_symbol, gene_lengths$Gene), ]


myData.counts.filt <- myData.counts.filt[genes_symbol,]

## rpkm calculation
dge <- calcNormFactors(myData.counts.filt)
rpkm.total.Riboseq <- rpkm(dge, gene.length=gene_lengths.selected$Length) # rpkm for individual samples
rpkm_group.total.Riboseq <- rpkmByGroup(dge, gene.length=gene_lengths.selected$Length, group = design_matrix$Condition) # rpkm based on samples group


################################################################################
### Common genes in Input RNAseq/Polysome RNAseq
################################################################################
gene_lengths_pc <- 
  gene_lengths |>
  filter(biotype=="protein_coding")

cpm.common <- intersect(rownames(cpm.total.RNAseq),rownames(cpm.total.Riboseq))
rpkm.common <- intersect(rownames(rpkm.total.RNAseq), rownames(rpkm_group.total.Riboseq))
rpkm.genes.pc <- sort(intersect(gene_lengths_pc$Gene,rpkm.common))

rpkm_group.total.RNAseq.pc <- rpkm_group.total.RNAseq[rpkm.genes.pc,]
rpkm.total.RNAseq.pc <- rpkm.total.RNAseq[rpkm.genes.pc,]
cpm.total.Riboseq.pc <- cpm.total.Riboseq[rpkm.genes.pc,]

rpkm_group.total.Riboseq.pc <- rpkm_group.total.Riboseq[rpkm.genes.pc,]
rpkm.total.Riboseq.pc <- rpkm.total.Riboseq[rpkm.genes.pc,]
cpm.total.Riboseq.pc <- cpm.total.Riboseq[rpkm.genes.pc,]

TE.matrix <- rpkm.total.Riboseq.pc/(rpkm.total.Riboseq.pc+rpkm.total.RNAseq.pc)
TE.matrix.1 <- rpkm.total.Riboseq.pc/rpkm.total.RNAseq.pc


TE.matrix <-   TE.matrix[complete.cases(TE.matrix) & !apply(TE.matrix, 1, function(row) any(is.infinite(row))), ]
TE.matrix.1 <-   TE.matrix.1[complete.cases(TE.matrix.1) & !apply(TE.matrix.1, 1, function(row) any(is.infinite(row))), ]

################################################################################
## DE analysis
################################################################################

# Create the design matrix

contrast_matrix <- makeContrasts(
  tgfb_vs_unt = tgfb - unt,
  tgfbCX_vs_tgfb = tgfbCX - tgfb,
  tgfbCX_vs_unt = tgfbCX - unt,
  levels = design_matrix
)

fit <- lmFit( object = TE.matrix , design = design_matrix)

### contrast fit
verb("\tcontrast fit.\n")

cfit <- contrasts.fit( fit = fit , contrasts = contrast_matrix )


### ebayes
verb("\teBayes\n")

cfit <- eBayes(fit = cfit , trend = ebayes.trend , robust = ebayes.robust )

#test.res  <-  limma::decideTests( cfit ,  adjust = MHadjust , method = "separate" , p.value = FDR.thresh )
# Extract the top table for each contrast
topTable_tgfb_vs_unt <- topTable(cfit, coef = "tgfb_vs_unt", number = Inf)
topTable_tgfbCX_vs_unt <- topTable(cfit, coef = "tgfbCX_vs_unt", number = Inf)
topTable_tgfbCX_vs_tgfb <- topTable(cfit, coef = "tgfbCX_vs_tgfb", number = Inf)



# Define the p-value and FDR thresholds
pvalue_threshold <- 0.05
fdr_threshold <- 0.05

# Extract the top table for each contrast with thresholds
topTable_tgfb_vs_unt <- topTable(cfit, coef = "tgfb_vs_unt", number = Inf, p.value = pvalue_threshold, adjust.method = "fdr")
topTable_tgfbCX_vs_unt <- topTable(cfit, coef = "tgfbCX_vs_unt", number = Inf, p.value = pvalue_threshold, adjust.method = "fdr")
topTable_tgfbCX_vs_tgfb <- topTable(cfit, coef = "tgfbCX_vs_tgfb", number = Inf, p.value = pvalue_threshold, adjust.method = "fdr")

# Filter results based on FDR threshold
topTable_tgfb_vs_unt <- topTable_tgfb_vs_unt[topTable_tgfb_vs_unt$adj.P.Val <= fdr_threshold & abs(topTable_tgfb_vs_unt$logFC) >= th_log2fc, ]
#topTable_tgfbCX_vs_unt <- topTable_tgfbCX_vs_unt[topTable_tgfbCX_vs_unt$adj.P.Val <= fdr_threshold & abs(topTable_tgfbCX_vs_unt$logFC) >= th_log2fc, ]
topTable_tgfbCX_vs_tgfb <- topTable_tgfbCX_vs_tgfb[topTable_tgfbCX_vs_tgfb$adj.P.Val <= fdr_threshold & abs(topTable_tgfbCX_vs_tgfb$logFC) >= th_log2fc, ]

# Print the filtered results
#print(topTable_tgfb_vs_unt)
#print(topTable_tgfbCX_vs_unt)
#print(topTable_tgfbCX_vs_tgfb)




# Remove genes with NA, NaN, Inf and -Inf
#df_clean <-   TE.matrix.1[complete.cases(TE.matrix.1) & !apply(TE.matrix.1, 1, function(row) any(is.infinite(row))), ]
fit <- lmFit( object = TE.matrix.1, design = design_matrix)
### contrast fit
verb("\tcontrast fit.\n")

cfit <- contrasts.fit( fit = fit , contrasts = contrast_matrix )

### ebayes
verb("\teBayes\n")
cfit <- eBayes(fit = cfit , trend = ebayes.trend , robust = ebayes.robust )

#test.res  <-  limma::decideTests( cfit ,  adjust = MHadjust , method = "separate" , p.value = FDR.thresh )
# Extract the top table for each contrast
topTable_tgfb_vs_unt_1 <- topTable(cfit, coef = "tgfb_vs_unt", number = Inf)
#topTable_tgfbCX_vs_unt <- topTable(cfit, coef = "tgfbCX_vs_unt", number = Inf)
topTable_tgfbCX_vs_tgfb_1 <- topTable(cfit, coef = "tgfbCX_vs_tgfb", number = Inf)

# Filter results based on FDR threshold
topTable_tgfb_vs_unt_1 <- topTable_tgfb_vs_unt_1[topTable_tgfb_vs_unt_1$adj.P.Val <= fdr_threshold & abs(topTable_tgfb_vs_unt_1$logFC) >= th_log2fc, ]
#topTable_tgfbCX_vs_unt <- topTable_tgfbCX_vs_unt[topTable_tgfbCX_vs_unt$adj.P.Val <= fdr_threshold & abs(topTable_tgfbCX_vs_unt$logFC) >= th_log2fc, ]
topTable_tgfbCX_vs_tgfb_1 <- topTable_tgfbCX_vs_tgfb_1[topTable_tgfbCX_vs_tgfb_1$adj.P.Val <= fdr_threshold & abs(topTable_tgfbCX_vs_tgfb_1$logFC) >= th_log2fc, ]






################################################################################
## Davide
################################################################################



exp.common.genes <- intersect(rownames(mtx_count_inut),rownames(mtx_count_polysome))
exp.common.genes.pc <- sort(intersect(gene_lengths_pc$Gene, exp.common.genes))


mtx_count_mono.pc <- mtx_count_inut[exp.common.genes.pc,]
mtx_count_poly.pc <- mtx_count_polysome[exp.common.genes.pc,]
#mtx_count_poly.pc <- log2(mtx_count_mono.pc)


contrast_matrix <- makeContrasts(
  tgfb_vs_unt = tgfb - unt,
  tgfbCX_vs_tgfb = tgfbCX - tgfb,
  tgfbCX_vs_unt = tgfbCX - unt,
  levels = design_matrix
)

matrix_davide <- mtx_count_poly.pc+1/(mtx_count_mono.pc + mtx_count_poly.pc)+1

matrix_davide <-   matrix_davide[complete.cases(matrix_davide) & !apply(matrix_davide, 1, function(row) any(is.infinite(row))), ]
matrix_davide <- log2(matrix_davide)

fit <- lmFit( object = matrix_davide , design = design_matrix)

### contrast fit
verb("\tcontrast fit.\n")

cfit <- contrasts.fit( fit = fit , contrasts = contrast_matrix )


### ebayes
verb("\teBayes\n")

cfit <- eBayes(fit = cfit , trend = ebayes.trend , robust = ebayes.robust )

#test.res  <-  limma::decideTests( cfit ,  adjust = MHadjust , method = "separate" , p.value = FDR.thresh )
# Extract the top table for each contrast
topTable_tgfb_vs_unt_2 <- topTable(cfit, coef = "tgfb_vs_unt", number = Inf)
topTable_tgfbCX_vs_unt_2 <- topTable(cfit, coef = "tgfbCX_vs_unt", number = Inf)
topTable_tgfbCX_vs_tgfb_2 <- topTable(cfit, coef = "tgfbCX_vs_tgfb", number = Inf)


topTable_tgfb_vs_unt_2 <- topTable_tgfb_vs_unt_2[topTable_tgfb_vs_unt_2$adj.P.Val <= 0.005 & abs(topTable_tgfb_vs_unt_2$logFC) >= 1, ]
#topTable_tgfbCX_vs_unt <- topTable_tgfbCX_vs_unt[topTable_tgfbCX_vs_unt$adj.P.Val <= fdr_threshold & abs(topTable_tgfbCX_vs_unt$logFC) >= th_log2fc, ]
topTable_tgfbCX_vs_tgfb_2 <- topTable_tgfbCX_vs_tgfb_2[topTable_tgfbCX_vs_tgfb_2$adj.P.Val <= 0.005 & abs(topTable_tgfbCX_vs_tgfb_2$logFC) >= 1, ]







################################################################################
## Save excel and Rdata image
################################################################################


#TGFb_Vs_Unt :: TE.matrix <- rpkm.total.Riboseq.pc/(rpkm.total.Riboseq.pc+rpkm.total.RNAseq.pc)
# TGFb_Vs_Unt_1 :: TE.matrix.1 <- rpkm.total.Riboseq.pc/rpkm.total.RNAseq.pc
# matrix_davide <-   matrix_davide[complete.cases(matrix_davide) & !apply(matrix_davide, 1, function(row) any(is.infinite(row))), ]
# matrix_davide <- log2(matrix_davide)

xl_list <- list('TGFb_Vs_Unt'= topTable_tgfb_vs_unt, 'CX5461_Vs_TGFb'= topTable_tgfbCX_vs_tgfb, 'TGFb_Vs_Unt_1'= topTable_tgfb_vs_unt_1, 'CX5461_Vs_TGFb_1'= topTable_tgfbCX_vs_tgfb_1,
                'TGFb_Vs_Unt_davide' = topTable_tgfb_vs_unt_2, 'CX541_Vs_TGFb_davide' = topTable_tgfbCX_vs_tgfb_2)
write.xlsx(xl_list, file = "Translational efficiency in Input RNAseq Vs Polysome RNAseq profiling.xlsx", rowNames=TRUE)



save.image("translational_eficiency_input_vs_polysome.Rdata")


