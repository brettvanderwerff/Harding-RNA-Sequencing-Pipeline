# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4387895/
setwd('/media/vande060/101A04EA101A04EA/MSP_RNA_Seq/featureCounts_output')

# load edgeR library
library(edgeR)

# read in csv containing counts
datain <- read.delim(file='featureCounts_matrix.csv', sep=',', row.names = 'Geneid')

# visualize the data
head(datain)


# defines grouping factors for the treaments
group <- c('control', 'control', 'msp_ron', 'control', 'msp_ron', 'msp_ron')

# generate DGElist object
dge <- DGEList(counts=datain, group=group)

# visualize the DGElist object
dge

dge = calcNormFactors(dge)
dge
par("mar")
par(mar=c(1,1,1,1))
main='MDS Plot for Count Data'
plotMDS(dge,main=main,labels=colnames(dge$counts),col=as.numeric(dge$samples$group),las=1)
normalized.counts = cpm(dge)
transposed = t(normalized.counts)
distance = dist(transposed)
clusters=hclust(distance)
plot(clusters)

dge=estimateCommonDisp(dge)
dge=estimateTagwiseDisp(dge)
dge
dex=exactTest(dge,pair=c("control","PDGF"),dispersion="tagwise")

fdrvalues=p.adjust(dex$table$PValue, method='BH') 
dex$table$fdr=fdrvalues
summary(decideTestsDGE(dex,p=0.05))
summary(decideTestsDGE(dex,p=0.01))
summary(decideTestsDGE(dex,p=0.005))
cutoff=0.005 
de = decideTestsDGE(dex, p = cutoff, adjust = "BH")
detags = rownames(dex)[as.logical(de)]
plotSmear(dex, de.tags = detags)
abline(h = c(-1, 1), col = "blue")


# get normalized counts per million
cpms=cpm(dge$counts) 
# find out which columns have controls
control=grep('control',colnames(cpms))
# find out which columns have treatments
msp_ron=grep('msp_ron',colnames(cpms))
# calculate the average expression for control samples
ave.control=apply(cpms[,control],1,mean)
# calculate average expression for treatment samples
ave.msp_ron=apply(cpms[,msp_ron],1,mean)
# make a data frame 
res=data.frame(gene=row.names(dex$table),
               fdr=dex$table$fdr,
               logFC=dex$table$logFC,
               control=ave.control,
               msp_ron=ave.msp_ron)
# read gene information (originally from IGB quickload site)
# this is a "bed detail" file

# keep gene id and gene description columns
annots=read.delim(file='featureCounts_matrix.csv', sep=',')[1]
# name the columns
names(annots)=c('gene')
# combine gene expression and annotations
res=merge(res,annots,by.x='gene',by.y='gene')
# put the results in order of signficance (fdr)
res=res[order(res$fdr),]
print(res)
# put columns in an order easy to browse
res=res[,c('fdr','logFC','control','msp_ron','gene')]
# write DE genes to a file:
out_file='results.tsv'
# select just the rows that made the cutoff
de=res$fdr<=cutoff
# write DE genes only
write.table(res[de,],file=out_file,row.names=F,sep='\t',quote=F)
# Make a file we can load into LycoCyc for pathways visualization
out_file='results/forLycoCyc.tsv'
write.table(res[de,c('gene','logFC')],file=out_file,quote=FALSE,
            sep='\t',col.names=FALSE,row.names=FALSE)
# write all DE genes
out_file='results/tomatoAllDE.txt'
write.table(res,file=out_file,row.names=F,sep='\t',quote=F)