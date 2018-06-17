# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4387895/

arg <- commandArgs(trailingOnly = TRUE)

print(arg)

setwd(arg)

# load edgeR library
library(edgeR)

# read in csv containing counts
datain <- read.delim(file='featureCounts_matrix.csv', sep=',', row.names = 'Geneid')

# visualize the data
head(datain)


# defines grouping factors for the treaments
group <- factor(c('control', 'control', 'msp_ron', 'control', 'msp_ron', 'msp_ron'))

# generate DGElist object, genes that have all zero gene counts
dge <- DGEList(counts=datain, group=group)

# visualize the DGElist object
dge


# filter DGElist object to keep only the genes that have at least .5 count per million (cpm) in at least 2 samples

keep <- rowSums(cpm(dge) > 0.5) >= 2
dge <- dge[keep, keep.lib.sizes=FALSE]


# normalization of gene counts by the weighted trimmed mean of of M-values method
dge <- calcNormFactors(dge, method="TMM")

# MDS plot of DGE list object

par("mar")
par(mar=c(1,1,1,1))
plotMDS(dge)

# MD plots of all samples to detect skews in gene expression
for (column in c(1,2,3,4,5,6)){
  plotMD(dge, column=column)
  abline(h=0, col="red", lty=2, lwd=2)
}

# Generation fo the design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

# Dispersion estimate

dge <- estimateDisp(dge, design, robust=TRUE)

# Visualize dispersion estimate
plotBCV(dge)

# Estimation of QL dispersions

fit <- glmQLFit(dge, design, robust=TRUE)
head(fit$coefficients)

# Visualize QL dispersions

plotQLDisp(fit)

# summarize the df prior
summary(fit$df.prior)

# make differential gene expression comparison

controlvsmsp_ron <- makeContrasts(control-msp_ron,levels=design)
result <- glmQLFTest(fit, contrast=controlvsmsp_ron)

# observe DE genes

topTags(result)

# total number of DE genes identified at the default FDR of 5%

is.de <- decideTestsDGE(result)
summary(is.de)

# magnitude of differential expression plots visualized with fitted model MD plot

plotMD(result, status=is.de, values=c(1,-1), col=c('red', 'blue'), legend='topright')

# To reduce the number of genes to inlcude only those with log fold changes grater than 1.5

tr <- glmTreat(fit, contrast=controlvsmsp_ron, lfc=log2(1.5))
topTags(tr)

# total number of DE genes identified with logFC > 1.5 at the default FDR of 5%

is.de <- decideTestsDGE(tr)
summary(is.de)

# save the results table
write.table(is.de, file='key.csv', sep=',')
write.table(topTags(tr, n=Inf), file='trimmed_results.csv', sep=',')

# visualize more stringent list of DE genes 

plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")





