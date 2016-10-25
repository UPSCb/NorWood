###  Import data from raw read files
library(DESeq2)
library(scatterplot3d)
library(RColorBrewer)
library(vsn)
library(LSD)
library(genefilter)
library(gplots)
library(parallel)
library(pander)
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#### Data handeling

### Remember to run functions before running data.

### save and load latest data
## Latest saved 2015-10-14
#save.image("latest.R")
#load("latest.R")

## set workdir
setwd("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata")
indir <- "../HTSeq/"

res <- mclapply(dir(indir,pattern="*.txt",full.names=TRUE),
                read.delim,header=FALSE,stringsAsFactors=FALSE,mc.cores=8L)

names(res) <- sub("_sortmerna_trimmomatic_STAR\\.txt","",dir(indir,pattern="*.txt"))

addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])

count.table <- do.call(cbind,lapply(res,"[",3))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]

count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]
pander(apply(count.stats,2,function(co){round(co*100/sum(co))}))

mar <- par("mar")
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=c(0,16000000)
        #range(count.stats)
        ,cex.names=.6)
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

## Check gene order

stopifnot(all(sapply(res[1:length(res)],function(r,o){
    all(r[,1] == o)
},res[[1]][,1])))
count.table <- do.call(cbind,lapply(res,"[",2))
tail(count.table)
head(res)

#addInfo<-c("__no_feature","__ambigous","__too_low_aQual","__not_aligned","__alignment_not_unique")
#sel <- match(addInfo,res[[1]][,1])
#count.table <- res[-sel,]
## remove last lines!!  Statistics data from files Only needed for gene models
colnames(count.table) <- names(res)
rownames(count.table) <- sub("\\.0$","",res[[1]][,1])
nmodels = length(count.table[,1])-5
count.table2 <- count.table[1:nmodels,]
rm(nmodels)
tail(count.table)  #confirm correct rows are removed
count.table <- count.table2
rm(count.table2)
annot <- colnames(count.table)
## Store imported counttable
count.table.raw <- count.table

## 
#count.table = count.table.merge
#rownames(count.table) <- merge(rownames(count.table.ref),row.names(count.table.new),byrow=TRUE)
annot <- as.data.frame(annot)
colnames(annot) <- c("sample")
conditions <- factor(sub("-[0-9][0-9]","",annot[match(colnames(count.table),annot$sample),"sample"]))
annot$conditions <- conditions

count.table.full.raw <- count.table.raw
#count.table.raw <- count.table.full.raw
rm(count.table)
rm(res)
####
##  analysis with images and VST calculation.
###

rmphloem = TRUE
rmphloemall = TRUE
count.table.raw <- filterBadSamples(count.table.full.raw,rmphloem=rmphloem,rmphloemall=rmphloemall)
#count.table.raw <- filterBadSamples(count.table.raw,rmphloem=rmphloem)
head(count.table.raw)

#count.table.base <- count.table.raw
### Do variance stabilization transformation
## Save exon models
#save.image("spruce_wood_count_table.RData")
if(rmphloemall){
    #annot_filt <- annot[-c(38,37,20,21,1),]
    annot_filt <- annot[-c(37,20,21,1),]
}else if(rmphloem){
    annot_filt <- annot[-c(20,21,1),]
}else{
    annot_filt <- annot
}
###  Do the other dataset and then continue with both


#dds.raw <- DESeqDataSetFromMatrix(countData=count.table.raw, colData=annot, design=~conditions)
dds.raw <- DESeqDataSetFromMatrix(countData=count.table.raw, colData=annot_filt, design=~conditions)
#dds.raw.2 <- DESeqDataSetFromMatrix(countData=count.table.raw.2, colData=annot_filt, design=~conditions)

## Only needed if gene model otw continue with DEXSeq_new.R

dds.raw.2 <- estimateSizeFactors(dds.raw.2)
dds.raw.2 <- estimateDispersions(dds.raw.2)

#png(file="analysis/Dispersion-estimation-heatmapplotVST.png",height=600,width=600,pointsize=16)
#comparisonplot(cds,main="Dispersion Estimation")
#dev.off()
plotDispEsts(dds.raw.2, main="Dispersion estimation")
.lib

### VST
vst.raw <- varianceStabilizingTransformation(dds.raw)
vst.raw.exprs <- as.data.frame(assay(vst.raw))
colnames(vst.raw.exprs) <- colnames(count.table.raw)
vst.raw.exprs <- vst.raw.exprs-min(vst.raw.exprs)
## vst.raw.exprs.2 contains T3-02

vst.raw.2 <- varianceStabilizingTransformation(dds.raw.2)
vst.raw.exprs.2 <- as.data.frame(assay(vst.raw.2))
vst.raw.exprs.2 <- vst.raw.exprs.2-min(vst.raw.exprs.2)
colnames(vst.raw.exprs.2) <- colnames(count.table.raw.2)

## Rlog2
rlog.raw <- rlog(dds.raw)
save.image("2016_03_23_withrlog.raw.Rdata")
rlog2.df <- as.data.frame(assay(rlog.raw))

colnames(rlog2.df) <- colnames(vst.raw.exprs)
library("limma")
library("LSD")
par(mfrow=c(2,1))
plotMA(rlog2.df,array=5,ylab="this vs other",main=paste("rlog2 ",colnames(rlog2.df)[5], " vs other"))
abline(h=median(as.numeric(t(as.vector(rlog2.df[5])))),col="red")
plotMA(vst.raw.exprs,array=5,ylab="this vs other",main=paste("vst ",colnames(rlog2.df)[5], " vs other"))
abline(h=median(as.numeric(t(as.vector(vst.raw.exprs[5])))),col="red")
plotMA(rlog2.df,array=3,ylab="this vs other",main=paste("rlog2 ",colnames(rlog2.df)[3], " vs other"))
abline(h=median(as.numeric(t(as.vector(rlog2.df[3])))),col="red")
plotMA(vst.raw.exprs,array=3,ylab="this vs other",main=paste("vst ",colnames(rlog2.df)[3], " vs other"))
abline(h=median(as.numeric(t(as.vector(vst.raw.exprs[3])))),col="red")

par(mfrow=c(1,2))
samp = "T2-17"
samp = match(samp,colnames(vst.raw.exprs))
y <- rlog2.df
heatscatter(x=rowMeans(y), y=y[,samp]-rowMeans(y[-c(samp)]),main="rlog2",xlab="Average expression",ylab=paste(colnames(vst.raw.exprs)[samp]," vs all other"))
abline(h=median(as.numeric(t(as.vector(rlog2.df[samp])))),col="red")
abline(lm(y[,samp]-rowMeans(y[-c(samp)])~rowMeans(y)),col="blue")
y <- vst.raw.exprs
heatscatter(x=rowMeans(y), y=y[,samp]-rowMeans(y[-c(samp)]),main="vst",xlab="Average expression",ylab=paste(colnames(vst.raw.exprs)[samp]," vs all other"))
abline(h=median(as.numeric(t(as.vector(vst.raw.exprs[samp])))),col="red")
abline(lm(y[,samp]-rowMeans(y[-c(samp)])~rowMeans(y)),col="blue")




LSD::plotmatrix(as.matrix(rlog2.df))
LSD::plotmatrix(as.matrix(vst.raw.exprs))
detach(package:limma)

### filter dataset
rlog2.df.filtered <- rlog2.df[apply(rlog2.df,1,function(x) sum(x==0)==0),]
rlog2.df.filtered2 <- rlog2.df.filtered[apply(rlog2.df.filtered,1,function(x) sum(x<3)<length(colnames(rlog2.df.filtered))-6),]
head(rlog2.df.filtered2)

sum(head(rlog2.df)[2,]==0)

vignette("vsn")
system.file(package="DESeq2")
#####  Find matching genes between the two sets and check the difference in expression!
match(c("MA_34359g0010"),rownames(vst.raw.filtered.2))
plot(x=seq(1,52),vst.raw.filtered.2[8012,seq(1,52)],type="l",ylim=c(0,13),xlab="sample",ylab="vst exprs")
lines(x=seq(1,52),vst.raw.filtered.2[2636,seq(1,52)],col="red")
abline(v=18)
abline(v=33)

plot(x=seq(34,52),vst.raw.filtered.2[9143,seq(34,52)],type="l",ylim=c(0,13),xlab="sample",ylab="vst exprs")
lines(x=seq(34,52),vst.raw.filtered.2[11491,seq(34,52)],col="red")
abline(v=18)
abline(v=33)


#filter all zero genes

vst.raw.exprs.expressed <- vst.raw.exprs[apply(vst.raw.exprs[,-1], 1, function(x) all(mean(x)>0)),]

vst.raw.filtered <- filter_exprs(vst.raw.exprs,n=51,R=c("T1","T2","T3"),lim=5,k=2,r=2)
vst.raw.filtered.2 <- filter_exprs(vst.raw.exprs.2,n=52,R=c("T1","T2","T3"),lim=3,k=2,r=2)


# Write matrix
if(FALSE){
    write.table(x=vst.raw.filtered,quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t", file="spruce.wood.expression.table_2015-09-24.txt")
    write.table(x=rownames(vst.raw.filtered),quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t", file="spruce.wood.genes_2015-09-24.txt")
    write.table(x=t(colnames(vst.raw.filtered)),quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t", file="spruce.wood.samples_2015-09-24.txt")
}
if(TRUE){
    write.table(x=vst.raw.filtered.2,quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t", file="spruce.wood.expression.table_2016-06-13.txt")
    write.table(x=rownames(vst.raw.filtered.2),quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t", file="spruce.wood.genes_2016-06-13.txt")
    write.table(x=t(colnames(vst.raw.filtered.2)),quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t", file="spruce.wood.samples_2016-06-13.txt")
}

## 44 72
k=seq(71,80)
match("MA_10433465g0030",rownames(vst.raw.exprs.expressed))

k=c(9370,8470,7726,7489,7085,6951,3309,1521)
k=c(10685,7085,6951,1521)
plotfast(k,vst.raw.exprs.expressed, n=53,main="Novel expressed genes (mean > 3) 1-5")

### Plot gene profiles
#### PCA  VST

## Long non Coding RNAs
if(rmphloemall){
    new_conditions <- conditions[-c(38,37,20,21,1)]
} else if(rmphloem){
    new_conditions <- conditions[-c(20,21,1)]
}else{
    new_conditions <- conditions
}
tmp_conditions <- conditions[1:18]
tmp_conditions <- factor(tmp_conditions)

vst.raw.filtered <- filter_exprs(vst.raw.exprs,n=51,R=c("T1","T2","T3"),lim=3,k=2,r=2)
pc <- prcomp(t(vst.raw.exprs))
makePCA(pc,main="PCA raw T1 T2 phloem filtered",conditions = new_conditions,comp=c(1,2)) 
pc <- prcomp(t(vst.raw.filtered))
makePCA(pc,main="PCA with phloem filtered vst5 minS 2",legPos="top",conditions=new_conditions,comp=c(1,2),lbs=annot_filt$sample)
typeof(annot_filt$sample)
head(vst.raw.filtered[16])  

plot.multidensity(log10(count.table.raw),col=rep(pal,each=6),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts log10"
)

plot.multidensity(vst.raw.exprs[1:20], col=rep(pal,each=6),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)

plot.multidensity(vst.raw.exprs[14:18], col=rep(pal,each=1),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)

##Plot multidensity T3 Log 10
plot.multidensity(vst.raw.exprs[34:44], col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)

plot.multidensity(vst.raw.exprs[45:53], col=rep(pal,each=1),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)
    ## Plot multidensity T3 filtered
vst.raw.filtered <- filter_exprs(vst.raw.exprs,n=51,R=c("T1","T2","T3"),lim=3,k=2,r=2)
plot.multidensity(vst.raw.filtered, col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)
plot.multidensity(vst.raw.filtered[1:11], col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)
plot.multidensity(vst.raw.filtered[12:23], col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)
plot.multidensity(vst.raw.filtered[24:33], col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)
plot.multidensity(vst.raw.filtered[34:44], col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)


plot.multidensity(vst.raw.filtered[45:51], col=rep(pal,each=1),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)

## Plot multidensity filtered T3 log 10
plot.multidensity(log10(vst.raw.filtered[34:44]), col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)

plot.multidensity(log10(vst.raw.filtered[45:53]), col=rep(pal,each=1),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood",
                  xlab="per gene counts (vst)"
)

plot.multidensity(rlog2.df.filtered2, col=rep(pal,each=6),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood (filtered)",
                  xlab="per gene exprs (rlog2)",
                  xlim=c(-6,25)
)

plot.multidensity(test, col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood (filtered by rep similarity)",
                  xlab="per gene exprs (vst)",
                  xlim=c(-6,25)
)



## All samples
plot.multidensity(rlog2.df, col=rep(pal,each=6),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood (all)",
                  xlab="per gene exprs (rlog2)",
                  xlim=c(-6,25)
)

plot.multidensity(vst.raw.exprs, col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution spruce wood (filtered)",
                  xlab="per gene counts (vst)"
)

plot.multidensity(vst.raw.filtered, col=rep(pal,each=3),
                  legend.x="topright", legend.cex=0.5,
                  main="sample vst counts distribution novel genes",
                  xlab="per gene counts (vst)"
)


## Create vedd diagram of rlog vs vst filtered
library("VennDiagram")

length(intersect(rownames(vst.raw.exprs.expressed.Vfilt),rownames(test2)))
length(intersect(rownames(vst.raw.filtered),rownames(rlog2.df.filtered2)))
setdiff(union(rownames(vst.raw.filtered),rownames(rlog2.df.filtered2)), intersect(rownames(rlog2.df.filtered2),rownames(vst.raw.filtered)))
setdiff(union(rownames(rlog2.df.filtered2),rownames(vst.raw.filtered)), intersect(rownames(vst.raw.filtered),rownames(rlog2.df.filtered2)))


###
#   Select aligned subsamples for expression filtering
###

##subsetSamples <- c("T1-06","T1-07","T1-08","T1-09","T1-10","T1-12","T1-13","T1-14","T1-16","T1-18","T2-04","T2-05","T2-06","T2-07","T2-08","T2-09","T2-10","T2-11","T2-13","T2-15","T3-04","T3-05","T3-06","T3-07","T3-09","T3-10","T3-11","T3-12","T3-14","T3-16")
#subsetConditions <- rep(seq(1,9),3)

## Set 1
if(FALSE){
    subsetSamples <- c("T1-06","T1-07","T1-08","T1-09","T1-10","T1-12","T1-13","T1-14","T1-16","T2-03","T2-04","T2-05","T2-06","T2-07","T2-09","T2-10","T2-13","T2-14","T3-04","T3-05","T3-06","T3-07","T3-09","T3-10","T3-11","T3-12","T3-14")
    endsamp <- length(subsetSamples)
    nsamp <- endsamp/3
    subsetConditions <- c(seq(1,endsamp,nsamp),seq(2,endsamp,nsamp),seq(3,endsamp,nsamp),seq(4,endsamp,nsamp),seq(5,endsamp,nsamp),seq(6,endsamp,nsamp),seq(7,endsamp,nsamp),seq(8,endsamp,nsamp),seq(9,endsamp,nsamp))
}

## Set 2
if(FALSE){
    subsetSamples <- c("T1-06","T1-07","T1-08","T1-09","T1-10","T1-12","T1-13","T1-14","T2-04","T2-05","T2-06","T2-07","T2-08","T2-09","T2-10","T2-11","T3-04","T3-05","T3-06","T3-07","T3-09","T3-10","T3-11","T3-12")
    endsamp <- length(subsetSamples)
    nsamp <- endsamp/3
    subsetConditions <- c(seq(1,endsamp,nsamp),seq(2,endsamp,nsamp),seq(3,endsamp,nsamp),seq(4,endsamp,nsamp),seq(5,endsamp,nsamp),seq(6,endsamp,nsamp),seq(7,endsamp,nsamp),seq(8,endsamp,nsamp))
}

## Set 3
if(FALSE){
    subsetSamples <- c("T1-06","T1-07","T1-08","T1-09","T1-10","T1-12","T1-13","T1-14","T1-16","T2-03","T2-05","T2-06","T2-07","T2-08","T2-09","T2-10","T2-11","T2-14","T3-04","T3-05","T3-06","T3-07","T3-09","T3-10","T3-11","T3-12","T3-14")
    endsamp <- length(subsetSamples)
    nsamp <- endsamp/3
    subsetConditions <- c(seq(1,endsamp,nsamp),seq(2,endsamp,nsamp),seq(3,endsamp,nsamp),seq(4,endsamp,nsamp),seq(5,endsamp,nsamp),seq(6,endsamp,nsamp),seq(7,endsamp,nsamp),seq(8,endsamp,nsamp),seq(9,endsamp,nsamp))
}

## Set 4 as many as possible
if(FALSE){
    subsetSamples <- c("T1-05","T1-06","T1-07","T1-08","T1-09","T1-10","T1-11","T1-12","T1-13","T1-14","T1-15","T1-16","T1-17","T1-18","T1-19","T2-03","T2-04","T2-05","T2-06","T2-07","T2-08","T2-09","T2-10","T2-11","T2-12","T2-13","T2-14","T2-15","T2-16","T2-17","T3-03","T3-04","T3-05","T3-06","T3-07","T3-08","T3-09","T3-10","T3-11","T3-12","T3-13","T3-14","T3-15","T3-16","T3-17")
    endsamp <- length(subsetSamples)
    nsamp <- endsamp/3
    subsetConditions <- c(seq(1,endsamp,nsamp),seq(2,endsamp,nsamp),seq(3,endsamp,nsamp),seq(4,endsamp,nsamp),seq(5,endsamp,nsamp),seq(6,endsamp,nsamp),seq(7,endsamp,nsamp),seq(8,endsamp,nsamp),seq(9,endsamp,nsamp)#,seq(nsamp,endsamp,nsamp)
                          ,seq(10,endsamp,nsamp),seq(11,endsamp,nsamp),seq(12,endsamp,nsamp),seq(13,endsamp,nsamp),seq(14,endsamp,nsamp),seq(15,endsamp,nsamp))
}

##Set 5 Variance filtered
vst.raw.exprs.expressed.Vfilt <- filterByVariance(vst.raw.filtered.2,filt=1)

vst.raw.filtered <- filter_exprs(vst.raw.exprs,n=51,R=c("T1","T2","T3"),lim=3,k=2,r=2)
#setdiff(rownames(vst.raw.exprs.expressed.Vfilt),rownames(vst.raw.exprs.expressed))

vst.raw.exprs.sample.filt <- vst.raw.filtered[,!is.na(match(colnames(vst.raw.filtered),subsetSamples))]

lm.smp.filt.s <- apply(vst.raw.exprs.sample.filt[,
                                                 c(seq(1,endsamp,nsamp),seq(2,endsamp,nsamp),seq(3,endsamp,nsamp),seq(4,endsamp,nsamp),seq(5,endsamp,nsamp),
                                                   seq(6,endsamp,nsamp),seq(7,endsamp,nsamp),seq(8,endsamp,nsamp)
                                                   ,seq(9,endsamp,nsamp)#,seq(nsamp,endsamp,nsamp)
                                                   ,seq(10,endsamp,nsamp),seq(11,endsamp,nsamp),seq(12,endsamp,nsamp),seq(13,endsamp,nsamp),seq(14,endsamp,nsamp),seq(15,endsamp,nsamp)
                                                 )]
                       ,1,function(y) loess(y~subsetConditions)$s)
#test <- vst.raw.exprs.sample.filt[lm.smp.filt.s>1,]
#test1 <- vst.raw.filtered[lm.smp.filt.s>1,]
#test2 <- vst.raw.filtered[lm.smp.filt.s>1,]
#test3 <- vst.raw.filtered[lm.smp.filt.s>1,]
test4 <- vst.raw.filtered[lm.smp.filt.s>1,]

#pcorr <- cor(t(vst.raw.exprs.sample.filt),method="pearson")
pcorr <- cor(t(test),method="pearson")

if(FALSE){
    nsamp <- seq(1,nsamp)
    matT1 <- vst.raw.exprs.sample.filt[,nsamp]
    matT2 <- vst.raw.exprs.sample.filt[,nsamp+length(nsamp)]
    matT3 <- vst.raw.exprs.sample.filt[,nsamp+length(nsamp)*2]
}

if(TRUE){
    nsamp <- seq(1,nsamp)
    matT1 <- test[,nsamp]
    matT2 <- test[,nsamp+length(nsamp)]
    matT3 <- test[,nsamp+length(nsamp)*2]
}
## Filtered
#resSmp1 <- calCorr(matT1,matT2,matT3)
#resSmp2 <- calCorr(matT1,matT2,matT3)
#resSmp3 <- calCorr(matT1,matT2,matT3)
#resSmp4 <- calCorr(matT1,matT2,matT3)
#resSmp1vf <- calCorr(matT1,matT2,matT3)
#resSmp2vf <- calCorr(matT1,matT2,matT3)
#resSmp3vf <- calCorr(matT1,matT2,matT3)
#resSmp4vf <- calCorr(matT1,matT2,matT3)


#write.table(file="correlate_profiles_w9samp.txt", x=res, quote=FALSE,col.names=TRUE,row.names=TRUE, sep="\t")
res <- resSmp4vf

tmp1 <- res[!is.na(res$T1vsT2),]
tmp2 <- res[!is.na(res$T1vsT3),]
tmp3 <- res[!is.na(res$T2vsT3),]
#par(mfrow=c(2,2))
plot(density(tmp1$T1vsT2),col="pink",
     #ylim=c(0,1.4),
     ylim=c(0,3.2),
     main="Density plots pw corr, Set4")
lines(density(tmp2$T1vsT3),col="red")
lines(density(tmp3$T2vsT3),col="blue")
legend("topleft", c("T1vsT2","T1vsT3","T2vsT3"),lwd=c(2.5,2.5,2.5),lty=c(1,1,1),cex=0.7,col=c("pink","red","blue"))


rm(tmp1,tmp2,tmp3)


tmpres <- as.data.frame(res>.7)
tmpfinal <- tmpres[rowSums(tmpres)>2,]
match("NA",rownames(tmpfinal))
rownames(tmpfinal)
match("MA_10430902g0010",rownames(tmpres))

length(union(rownames(tmpfinal),rownames(test)))
length(setdiff(rownames(tmpfinal),rownames(test)))
length(setdiff(rownames(test),rownames(tmpfinal)))
length(intersect(rownames(tmpfinal),rownames(test)))
#vennDiagram()
png(file="SetVfiltExprsFilt_vst3_Var1.5_vst3.png")
draw.pairwise.venn(length(rownames(vst.raw.exprs.expressed.Vfilt)),length(rownames(vst.raw.filtered)),length(intersect(rownames(vst.raw.exprs.expressed.Vfilt),rownames(vst.raw.filtered)))
                   ,category = c("rowVars", "exprs"), lty = rep("blank",2), 
                   fill = c("light blue", "pink"), alpha = rep(0.5, 2), 
                   cat.pos = c(0,0), cat.dist = rep(0.025, 2))
dev.off()
test <- vst.raw.filtered[match(intersect(rownames(vst.raw.exprs.expressed.Vfilt),rownames(vst.raw.filtered)),rownames(vst.raw.filtered)),]

pcorrT1T2 <- cor(t(vst.raw.exprs.sample.filt[,seq(1,16)]),method="pearson")
pcorrT2T3 <- cor(t(vst.raw.exprs.sample.filt[,seq(9,24)]),method="pearson")
pcorrT1T3 <- cor(t(vst.raw.exprs.sample.filt[,c(seq(1,8),seq(17,24))]),method="pearson")
pcorr[2,9981]
pcorrT1T3[2,9981]
pcorrT1T2[2,9981]
pcorrT2T3[2,9981]
match("MA_204448g0010",rownames(vst.raw.exprs.sample.filt))

write.table(file="test_corr_correlate_profiles.txt", x=pcorr, quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t")
write.table(file="test_corr_correlate_profilesT1T2.txt", x=pcorrT1T2, quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t")
write.table(file="test_corr_correlate_profilesT1T3.txt", x=pcorrT1T3, quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t")
write.table(file="test_corr_correlate_profilesT2T3.txt", x=pcorrT2T3, quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t")
library(affy)
par(mfrow=c(1,1))
pcorrNONA <- pcorrT2T3[,colSums(is.na(pcorrT2T3))<19013]
colnames(pcorrNONA) <- NULL
rownames(pcorrNONA) <- NULL
pcorrNONA <- pcorrNONA[rowSums(is.na(pcorrNONA))==0,]
pcorrNONA[2,9981]
plot(density(pcorrNONA),main="T2-T3")

#write.table(file="genes_corr.txt", x=row.names(vst.raw.exprs.sample.filt), quote=FALSE,col.names=FALSE,row.names=FALSE, sep="\t")


cat(rownames(test))

#loess(as.formula(x~subsetConditions),vst.raw.exprs.sample.filt)
summary(lm.smp.filt[[2614]])
summary(lm.smp.filt[[13]])
lm.smp.filt[[2614]]$s

plot(lm.smp.filt[[2614]])
plot(lm.smp.filt[[13]])
vst.raw.exprs.sample.filt[1,]

plot(x=seq(1,endsamp),y=vst.raw.exprs.sample.filt[2614,
                                                  c(seq(1,endsamp,nsamp),seq(2,endsamp,nsamp),seq(3,endsamp,nsamp),seq(4,endsamp,nsamp),seq(5,endsamp,nsamp),
                                                    seq(6,endsamp,nsamp),seq(7,endsamp,nsamp),seq(8,endsamp,nsamp)
                                                    #,seq(9,endsamp,nsamp),seq(nsamp,endsamp,nsamp)
                                                  )],xlab=c(0,30),ylab=c(0,10),col=rep(pal,each=3))
text(x=seq(1,endsamp),y=vst.raw.exprs.sample.filt[2614,],col=rep(pal,each=nsamp))

match("MA_10429177g0010",rownames(vst.raw.exprs.sample.filt))
###
head(count.table.raw[,c(25,26,27)])


###  Functions

plotfast <- function(x,dataset,n=25,main="Novel expressed genes (mean > 3)"){
    plot(xlim=c(1,n),ylim=c(1,13),x=rep(1:n),y=dataset[x[1],1:n],xlab="Sample",ylab="gene expression (vst)",main=main,type="l")
    col1=rep(brewer.pal(8,"Dark2"),each=1)
    colX <- 1
    for (k in x){
        if (colX == 9){
            return("max nr of genes in the same plot is 8")
        }
        rownames(dataset[,])[k]
        lines(x=rep(1:n), y=dataset[k,1:n],col=col1[colX])
        colX = colX + 1
        Sys.sleep(1)
    }
}

filter_genes_by_exprs <- function(dataset,n=50,R=c("T1","T2","T3"),lim=3,k=2,r=2){
    ## n = number of samples, 
    ## r = n samples in each replicate
    ## lim = min vst value, 
    ## k = n sampeles that needs to pass criteria in each replicate
    ## r = n replicates that need to pass criteria
    
    test <- sapply(R,function(coln,dat){
        rowSums(dat[,grep(coln,colnames(dat))])
    },dataset>lim)
    test[rowSums(test>k) >= r,]
}

filter_exprs <- function(dat,n=106,R=c("T1","T2","T3"),lim=1,k=2,r=2){
    ## n = number of samples, 
    ## r = n samples in each replicate
    ## lim = min vst value, 
    ## k = n sampeles that needs to pass criteria in each replicate
    ## r = n replicates that need to pass criteria
    filt <- names(rowSums(filter_genes_by_exprs(dat,n=n,R=R,lim=lim,k=k,r=r))>0)
    dat[row.names(dat) %in% filt,]
}


makePCA <- function(pc,main="PCA lncRNA",legPos="bottomleft",conditions=conditions,comp=c(1,2),lbs=FALSE)
{
    if(typeof(lbs) != "integer"){
        lbs=unlist(sapply(sapply(table(conditions),":",1),rev))
    }
    percent <- round(summary(pc)$importance[2,]*100)
    smpls <- colnames(vst.raw.exprs)
    cols <- c(1,rep(brewer.pal(8,"Dark2"),2))
    
    plot(pc$x[,comp[1]],pc$x[,comp[2]],main= main, type="n",col=cols[as.integer(factor(conditions))],xlim=range(pc$x[,comp[1]]),ylim=range(pc$x[,comp[2]]),
         xlab=paste("Comp. ",comp[1],"(",percent[comp[1]],"%)",sep=""),
         ylab=paste("Comp. ",comp[2],"(",percent[comp[2]],"%)",sep=""))
    text(labels=lbs,x=pc$x[,comp[1]],y=pc$x[,comp[2]],col=cols[as.integer(factor(conditions))])
    legend(legPos,pch=c(15,15,15,15,15),col=c(NA,cols[1:5],box.lwd = 0,box.col = "white",bg = "white"),
           legend=c("Color:",levels(factor(conditions))),cex=0.8,pt.cex = 1.0,text.font=1,box.lwd = 0,box.col = "white",bg = "transparent",y.intersp=0.7)
}

filterBadSamples <- function(count.table,rmphloem=TRUE,rmphloemall=FALSE){
    #Filter bad samples in 
    #count.table <- count.table.raw
    count.table2 <- count.table
    if(rmphloem){
        count.table2$"T2-01" <- NULL
        count.table2$"T2-02" <- NULL
        count.table2$"T1-01" <- NULL
        if(rmphloemall){
            count.table2$"T3-01" <- NULL
            count.table2$"T3-02" <- NULL
        }
    }
    count.table2$"T2-10" <- rowMeans(subset(count.table, select = c("T2-09", "T2-11")), na.rm = TRUE)
    count.table2$"T2-10" <- as.integer(count.table2$"T2-10")
    count.table <-  count.table2
    return(count.table)
}

### Calculate correlations

calCorr <- function(T1,T2,T3){
    ndf <- T1[,c(1,2,3)]
    for (i in 1:length(T1[,1])){
        ndf[i,1] <- cor(t(matT1[i,]),t(matT2[i,]),method="pearson")
        ndf[i,2] <- cor(t(matT2[i,]),t(matT3[i,]),method="pearson")
        ndf[i,3] <- cor(t(matT1[i,]),t(matT3[i,]),method="pearson")
    }
    colnames(ndf) <- c("T1vsT2","T2vsT3","T1vsT3")
    return(ndf)
}

filterByVariance <- function (mat,filt=1){
    res <- mat[rowVars(mat,na.rm=TRUE)>filt,]
    return(res)
}

k=c(10685,7085,6951,1521)
plotfast(k,tmp, n=51,main="Novel expressed genes (mean > 3) 1-5")
tmp <- vst.raw.filtered.2
tmp$"T1-15" <- NA
tmp$"T3-02" <- NULL