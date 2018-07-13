# Script is closely adapted from the following publications: doi 10.1007/978-1-4939-2444-8_24 and doi 10.12688/f1000research.8987.2

setwd(arg)

# Load the edgeR library
library(edgeR)

# Read in csv containing raw gene counts
datain <- read.delim(file='featureCounts_matrix.csv', sep=',', row.names = 'Geneid')

# Visualize the data by printing the head of the data frame
head(datain)

# Define grouping factors for the treaments, groups here are cells transfected with vector or and MSP/RON construct
group <- factor(c('control', 'control', 'msp_ron', 'control', 'msp_ron', 'msp_ron'))

# Instantiate the DGElist object (an EdgeR specific class)
dge <- DGEList(counts=datain, group=group)

# Visualize the DGElist object
dge

# Filter the DGElist object to keep only the genes that have at least .5 count per million (cpm) in at least 2 samples

keep <- rowSums(cpm(dge) > 0.5) >= 2
dge <- dge[keep, keep.lib.sizes=FALSE]

# Normalize the gene counts between groups by the weighted trimmed mean of of M-values method (doi 10.1186/gb-2010-11-3-r25)
dge <- calcNormFactors(dge, method="TMM")

# Generate an MDS plot of DGE list object to determine if there are outliers or batch effects
par("mar")
par(mar=c(1,1,1,1))
plotMDS(dge, col=as.numeric(dge$samples$group))

# Generate MD plots of all samples to detect skews in the gene expression profile
for (column in c(1,2,3,4,5,6)){
  plotMD(dge, column=column)
  abline(h=0, col="red", lty=2, lwd=2)
}

# Generation of the design matrix, which is required for linear modeling
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

# Dispersion estimate parameter is created to account for variability between biological replicates

dge <- estimateDisp(dge, design, robust=TRUE)

# Visualize the dispersion estimate
plotBCV(dge)

# Estimation of the Quasi-likelihood (QL) dispersions

fit <- glmQLFit(dge, design, robust=TRUE)
head(fit$coefficients)

# Visualize QL dispersions

plotQLDisp(fit)

# Summarize the degrees of freedom prior
summary(fit$df.prior)

# Calculate differential gene expression comparison between vector control and MSP/RON transfected cells

controlvsmsp_ron <- makeContrasts(control-msp_ron,levels=design)
result <- glmQLFTest(fit, contrast=controlvsmsp_ron)

# Observe a list of differentially expressed (DE) genes

topTags(result)

# Print the total number of DE genes identified at the default false discovery rate (FDR) cutoff of 5%

is.de <- decideTestsDGE(result)
summary(is.de)

# Visualize the magnitude of differential expression with a fitted model MD plot

plotMD(result, status=is.de, values=c(1,-1), col=c('red', 'blue'), legend='topright')

# Filter DE gene list to inlcude only those with log fold changes greater than 1.5

tr <- glmTreat(fit, contrast=controlvsmsp_ron, lfc=log2(1.5))
topTags(tr)

# Print the total number of DE genes identified with logFC > 1.5 at the default FDR of 5%

is.de <- decideTestsDGE(tr)
summary(is.de)

# Visualize the more stringently filtered list of DE genes 

plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")

# Save the filtered DE gene list results as a csv
write.table(topTags(tr, n=Inf), file='trimmed_results.csv', sep=',')







