#data <- read.table("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/expression/spruce.wood.expression.table_2016-06-23.txt", header=F, quote="\"")
#data <- read.table("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/expression/spruce.wood.expression.table_2016-06-27_no18.txt", header=F, quote="\"")

#rownames(data) <- as.vector(read.table("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/expression/spruce.wood.genes_2016-06-13.txt", header=F, quote="\"")$V1)
#colnames(data) <- as.vector(t(read.table("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/expression/spruce.wood.samples_2016-06-13.txt", header=F, quote="\""))[,1])
#data$Genes <- rownames(data)
#head(data)
# 26: K1
# 107: K1, K2, K3, K5

#data <- vst.raw.filtered.2
#data <- test
#data <- test1
#data <- vst.raw.exprs.expressed.Vfilt
#data <- vst.raw.filtered.variance
#data <- vst.raw.filtered.variance0.5
#data <- filterByVariance(tmp,filt=0.5)
#data <- filterByVariance(data,filt=1)

setwd("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/")
data <- read.table("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/expression/spruce.wood.expression.table_2016-06-27.txt", header=F, quote="\"")
dim(data)
#data <- read.table("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/expression/spruce.wood.expression.table_2016-06-27_no18.txt", header=F, quote="\"")
rownames(data) <- as.vector(read.table("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/expression/spruce.wood.genes_2016-06-27.txt", header=F, quote="\"")$V1)
colnames(data) <- as.vector(t(read.table("/mnt/picea/projects/spruce/htuominen/spruce-wood-cross-section/Rdata/expression/spruce.wood.samples_2016-06-27.txt", header=F, quote="\""))[,1])
#data$"T1-15" <- NULL
#data <- filterByVariance(data,filt=0.5)
data <- filterByVariance(data,filt=1)
n <- 51
data2 <- data

###############################################################################

# Plot heatmap with tree

#library(gplots)

# Cluster rows
head(data2)
#library(dendextend)
#rdist <- dist(data2, method="euclidean")h

rdist <- as.dist(1-cor(t(data2), method="pearson"))
rhc <- hclust(rdist,method="ward.D")
rdendro <- as.dendrogram(rhc)

if (TRUE){
    weights.Ord <- c(rep(0,length(rhc$order)))
    head(rhc$order)
    tail(rhc$order[ct==2])
    head(rhc$order[ct==1])
    weights.Ord[match(rhc$order[ct==1],rhc$order)] <- 300
    weights.Ord[match(rhc$order[ct==2],rhc$order)] <- 4000
    weights.Ord[match(rhc$order[ct==4],rhc$order)] <- 50000
    weights.Ord[match(rhc$order[ct==3],rhc$order)] <- 20
    weights.Ord[match(rhc$order[ct==5],rhc$order)] <- 1
    weights.Ord[match(rhc$order[ct==7],rhc$order)] <- 700000
    weights.Ord[match(rhc$order[ct==6],rhc$order)] <- 60000
    
    rhcOrder <- c(  rhc$order[ct==7],rhc$order,
                    rhc[3]$order[ct==6],
                    rhc[3]$order[ct==5],
                    rhc[3]$order[ct==4],
                    rhc[3]$order[ct==3],
                    rhc[3]$order[ct==2],
                    rhc[3]$order[ct==1]
                    )
    #save.image("ImageWithCorrectClusterOrder.R")
    rdendro <- reorder(rdendro,weights.Ord,agglo.FUN = mean)
}
#rdendro <- reorder(rdendro,rhcOrder, agglo.FUN = mean)
rhc <- as.hclust(rdendro)
# Cluster columns
cdist <- dist(t(data2), method="euclidean")
chc <- hclust(cdist,method="ward.D")
cdendro <- as.dendrogram(chc)

# Gene cluster color bar
nclust = 7
ct <- cutree(rhc, k=nclust)

#colgenes <- c("mediumpurple3","gray8","#CC5151","orange","#B26F2C","lightsteelblue3","#A3CC51","#51A3CC");
#colgenes <- c("gray8","#A3CC51","orange","#B26F2C","#51A3CC","lightsteelblue3","#CC5151","mediumpurple3","green","white","purple");
#colgenes <- c("gray8","#A3CC51","#A3CC51","orange","orange","#B26F2C","#B26F2C","#51A3CC","#51A3CC","lightsteelblue3","#CC5151","#CC5151","mediumpurple3","green","white","purple");
#colgenes <- terrain.colors(nclust, alpha = 1)
### spruce clusters
##,"#ad1f1f""#9f2d2d",
#colgenes <- c("gray8","#A3CC51","orange","#B26F2C","#51A3CC","lightsteelblue3","#CC5151","mediumpurple3","green","white","purple");
#colgenes <- c("#51A3CC","#708bdb","#d87474","mediumpurple3","#c63939","gray8","#B26F2C")
#colgenes <- c("#51A3CC","mediumpurple3","#c63939","mediumpurple3","#c63939","gray8","#B26F2C")
colgenes <- c("#51A3CC","lightsteelblue3","mediumpurple3","orange3","#c63939","#B26F2C","gray8")
colgenes <- c("#51A3CC","#8968CD","#c63939","#6B468B","#8B1A1A","gray8","#B26F2C")


c <- 0
prev <- 0
ct2 <- ct
for (j in nrow(data2):1) {
  
  g <- rhc$order[j]
  
  if (ct[g] != prev) {
    c <- c + 1
    prev <- ct[g]
  }
  
  ct2[g] <- c
  
}
ct <- ct2
vcol <- colgenes[ct]

order <- rep(0, (n-1))
h <- colnames(data2)
z <- strsplit(as.character(h), "[-]")
h <- unlist(lapply(z,"[",2))
h <- as.integer(h)
k <- 1
for(i in 1:20) {
  for(j in 1:(n-1)) {
     if (h[j] == i) {
       order[j] = k
      k <- k+1
    }
  }
}
cdendro <- reorder(cdendro,order, agglo.FUN = mean)
chc <- as.hclust(cdendro)
#chc$order

nclust2 = 3
ct2 <- cutree(chc, k=nclust2)
#col <- c("#A3CC51","#51A3CC","#CC5151","#B26F2C")
col <- c("#51A3CC","#CC5151","#B26F2C")
hcol <- col[ct2]
#save(hcol,file = "/mnt/picea/home/torgeir/wood/centrality_clustering/sample_colors.RData")

#pdf(file = "/mnt/picea/home/david/Rdata/SpruceWood/heatmap_var1/hclust_wood_samples_Kall_subclust.pdf")
#plot(cdendro, nodePar = list(pch = c(NA), lab.cex = 0.3), axes = FALSE)
#dev.off()

if (FALSE) { # Merge cluster a
  # 1-7
  ct <- gsub("\\<2\\>", 1, ct)
  ct <- gsub("\\<3\\>", 1, ct)
  ct <- gsub("\\<4\\>", 1, ct)
  ct <- gsub("\\<5\\>", 2, ct)
  ct <- gsub("\\<6\\>", 3, ct)
  ct <- gsub("\\<7\\>", 4, ct)
  ct <- gsub("\\<8\\>", 5, ct)
  ct <- gsub("\\<9\\>", 6, ct)
  ct <- gsub("\\<10\\>", 7, ct)
  ct <- gsub("\\<11\\>", 8, ct)
  ct <- gsub("\\<12\\>", 9, ct)
  ct <- gsub("\\<13\\>", 10, ct)
  ct <- gsub("\\<14\\>", 11, ct)
  ct <- gsub("\\<15\\>", 12, ct)
  ct <- gsub("\\<16\\>", 13, ct)
  
  ct <- as.integer(ct)
  vcol <- colgenes[ct]
  
  nclust <- 13
}

colorR <- colorRampPalette(c("blue","white","red"))(20)
#colorR[2:6]  <- colorR[1] 
#colorR[15:19] <- colorR[20]

#scaled <- scale(data2)
#min(scaled)
#max(scaled)
#hist(scaled)
#hist(as.matrix(data2))
#png(file = "/mnt/picea/home/david/Rdata/SpruceWood/hclust_wood_both_Kall_vst5_wname.png",width=15,height=15,units="cm",res=1000)
#pdf(file = "/mnt/picea/home/david/Rdata/SpruceWood/hclust_wood_both_3cl_Kall_vst3_52samp.pdf",width=15,height=15)
#pdf(file = "/mnt/picea/home/david/Rdata/SpruceWood/Thclust_wood_gclust_7cl_Kall_vst3_51samp_vfilt1.pdf",width=15,height=15)
#pdf(file = "/mnt/picea/home/david/Rdata/SpruceWood/hclust_wood_both_7cl_Kall_vst3_51samp_vfilt1.pdf",width=15,height=15)
#png(file = "/mnt/picea/home/david/Rdata/SpruceWood/hclust_wood_both_7cl_Kall_vst3_51samp_vfilt1.png",width=15,height=15,units="cm",res=1000)
png(file = "/mnt/picea/home/david/Rdata/SpruceWood/hclust_wood_gclust_7cl_Kall_vst3_51samp_vfilt_orgOrder_joined.png",width=15,height=15,units="cm",res=1000)


h <- heatmap.2(as.matrix(data2),
           Colv = cdendro, #cdendro, # FALSE
           dendrogram = "both", # both, row
           ColSideColors = hcol,
           #Colv = F,
           #dendrogram = "row",
           #colsep = c(18,33),
        Rowv = rdendro,
        scale = "row",
        margins = c(5, 5), 
        trace = "none",
        cexRow = 0.2,
        cexCol = 0.4, #0.4 0.75
        labRow = rep("", nrow(data2)), 
        #labRow = rownames(data2), 
        labCol = colnames(data2),  #as.matrix(header[t(1:n-1)]), 
        main = NULL,
        xlab = NULL, ylab = NULL,
        #col = redblue.colors(max(data2)),
        col = colorR,
        key=FALSE, keysize=1, density.info = "none",
        RowSideColors = vcol
)
dev.off()

####################
# PCA
shapes <- as.character(colnames(data2))
shapes <- gsub("\\.", "-", shapes)
shapes <- gsub("-0", "-", shapes)
#data2 <- vst.raw.exprs.expressed.Vfilt
data2 <- data
pc <- prcomp(t(data2), scale=TRUE)
white <- rep("white", ncol(data2))
pos <- rep(4, ncol(data2))
tmp <- summary(pc)
tmp$importance[2,1]*100
#pos[26] <- 1
#pos[52] <- 2
#pos[5] <- 2
#pos[30] <- 3
pdf(file = "/mnt/picea/home/david/Rdata/SpruceWood/pca_wood_both_Kall_pc12_VarianceFiltered_extra.pdf")
plot(pc$x[,1], pc$x[,2], 
        xlab=paste0("Principle component ", 1, " (", round(tmp$importance[2,1]*100,0), "%)")
        ,ylab=paste0("Principle component ", 2, " (", round(tmp$importance[2,2]*100,0), "%)"), 
     col=white, pch=1, cex=0.5, xlim=c(-150,150))
text(pc$x[,1], pc$x[,2], shapes, cex=0.55, pos=pos, offset=0, col=hcol)

dev.off()

pdf(file = "/mnt/picea/home/david/Rdata/SpruceWood/pca_wood_both_Kall_pc13_VarianceFiltered_extra.pdf")

plot(pc$x[,1], pc$x[,3],     
     xlab=paste0("Principle component ", 1, " (", round(tmp$importance[2,1]*100,0), "%)")
     ,ylab=paste0("Principle component ", 3, " (", round(tmp$importance[2,3]*100,0), "%)"),
     col=white, pch=1, cex=0.5, xlim=c(-150,150))
text(pc$x[,1], pc$x[,3], shapes, cex=0.55, pos=pos, offset=0, col=hcol)

dev.off()

pdf(file = "/mnt/picea/home/david/Rdata/SpruceWood/pca_wood_both_Kall_pc23_VarianceFiltered_extra.pdf")

plot(pc$x[,2], pc$x[,3], 
     xlab=paste0("Principle component ", 2, " (", round(tmp$importance[2,2]*100,0), "%)")
     ,ylab=paste0("Principle component ", 3, " (", round(tmp$importance[2,3]*100,0), "%)"), 
     col=white, pch=1, cex=0.5, xlim=c(-150,150))
text(pc$x[,2], pc$x[,3], shapes, cex=0.55, pos=pos, offset=0, col=hcol)

dev.off()

####

z <- strsplit(as.character(colnames(data2)), "\\.")
shapes <- unlist(lapply(z,"[",1))
names  <- unlist(lapply(z,"[",2))
              
shapes <- gsub("T1", "1", as.character(shapes))
shapes <- gsub("T2", "2", as.character(shapes))
shapes <- gsub("T3", "3", as.character(shapes))
shapes <- as.integer(shapes)

summary(pc <- prcomp(t(data2), scale=TRUE))
plot(pc$x[,1], pc$x[,2], xlab="Principle component 1", ylab="Principle component 2", col=hcol, pch=shapes, cex=0.5) #, xlim=c(-100,100))
text(pc$x[,1], pc$x[,2], names, cex=0.25, pos=4, col=hcol)


# Alternative
plot(pc$loadings[,1:2], pch=2, col=hcol, xlab="Principle component 1", ylab="Principle component 2")
text(pc$loadings[,1], pc$loadings[,2],colnames(data2), cex=0.6, pos=4, col=hcol)

####################
# Correct clustering
ct_org <- ct # ct <- ct_org

#clustername <- c("A","B1","B2","C1","C2","D1","D2","E1","E2","F","G1","G2","H")
#clustername <- c("A","B","C1","D1","C2","D2","E")
clustername <- c("A","B","C","D","E","F","G","A")
nclust <- length(clustername)

tableclust <- data.frame(data2,ct)

centroids <- data.frame()
for(j in 1:nclust) {
  
  group <- tableclust[tableclust$ct==j,]
  meangroup <- colMeans(group)
  
  if (j == 1) {
    centroids <- data.frame(meangroup[1:(n-1)])
  } else {
    centroids[,j] <- meangroup[1:(n-1)]
  }
}

# Novel genes
if (FALSE) {
  # novel_genes, lncRNA, fragments
  datan <- read.table("/mnt/picea/home/david/Rdata/SpruceWood/wood.filtered.expression_3_10_novel_genes.txt", header=T, quote="\"")
  data2n <- datan[,2:n]
  nrow(data2n)
  tableclust <- data.frame(data2n,ct=rep("a",nrow(data2n)))
  tableclust$ct
}

shift <- matrix(0, nclust, nclust)
clustcorr <- data.frame()
for(j in 1:nrow(tableclust)) {
#for(j in 1:10) {
  
  #gene <- data[data$Genes=="Potri.013G001000",2:n]
  gene <- tableclust[j,1:(n-1)]
  #data[rownames(gene),1]
  
  Rs <- cor(t(gene),centroids)
  R <- max(Rs)
  i_now <- tableclust[j,n]
  i <- which.max(Rs)
  #clustername[i]
  
  if (R < 0.50000) {
    i <- 1
  }
  
  shift[i_now,i] <- shift[i_now,i] + 1
  ct[j] <- i
  
  # OBS: data or datan
  g <- c(toString(data[j,1]),Rs)
  if (j == 1) {
    clustcorr <- data.frame(g)
  } else {
    clustcorr[,j] <- g
  }
}

if (TRUE) {
  clustcorr <- t(clustcorr)
  colnames(clustcorr) <- c("Gene",clustername)
  f <- paste0("/mnt/picea/home/david/Rdata/SpruceWood/T3heatmap_var1/cluster_correlations_novel_genes.txt")
  write.table(as.matrix(clustcorr), f, quote = FALSE, row.names = FALSE, col.names = TRUE)
}

shift
sum(diag(shift))
sum(shift)-sum(diag(shift))

####################
# Print clusters

filename <- "cluster"
#clustername <- c("A","B1","B2","C1","C2","D1","D2","E1","E2","F","G1","G2","H")
clustername <- c("A","B","C","D","E","F","G")
#clustername <- c("B","C1","D1","C2","D2","A","E")
clustername <- c("A","B","C","D","E","F","G","A")

#clustername <- c("A","B","C","D","E","F")
colnames(data)
m <- 19 # 26 107
#dim(data2)
data3 <- data2[,1:(m-1)] # only print profiles from K1
### Do K3 instead
#m <- 52 # 26 107
#dim(data2)
#data3 <- data2[,34:(m-1)] # only print profiles from K1
#m <- 19
#data2[,1:(m-1)]
#colnames(data3)
tableclust <- data.frame(data3,ct)
#tableclust <- tableclust[rhc$order,]

if (TRUE) { # pie diagram of clust dist
  
  n.genes <- rep(0,nclust)
  for (i in 1:nclust) {
    
    n.genes[i] <- nrow(tableclust[tableclust$ct==i,])
    
  }
  sum(n.genes)
  
  f <- paste0("/mnt/picea/home/david/Rdata/SpruceWood/T3heatmap_var1/pie_clusterdistr.png")
  png(file = f)
  pie(n.genes,labels = paste(clustername, " (",n.genes,")", sep = ""), col = colgenes, clockwise = TRUE, cex=0.7)
  dev.off()
}

if (FALSE) {
  tableclust2 <- data.frame(data,ct,rhc$order)
  write.table(tableclust2[rhc$order,], "/mnt/picea/home/david/Rdata/SpruceWood/T3heatmap_var1/wood.filtered.expression_3_10_clustered.txt");
}
centroids <- data.frame()
for(j in 1:nclust) {
  
  group <- tableclust[tableclust$ct==j,]
  
  # Plot correlation distributions
  if (TRUE) {
    Rm <- cor(t(group))
    R <- Rm[upper.tri(Rm, diag = FALSE)]
    f <- paste0("/mnt/picea/home/david/Rdata/SpruceWood/T3heatmap_var1/",filename,clustername[j],"-hist.png")
    png(file = f)
    hist(R, xlab = "Expression", ylab = "Gene pairs", labels = FALSE, main = "")
    dev.off()
  }
  
  # Plot
  if (TRUE) {
    f <- paste0("/mnt/picea/home/david/Rdata/SpruceWood/T3heatmap_var1/",filename,clustername[j],".png")
    png(file = f)
    
    max <- max(group[,1:(m-1)])
    min <- min(group[,1:(m-1)])
    
    meangroup <- colMeans(group)
    meangroup <- scale(meangroup) # SCALE!
    if (j == 1) {
      centroids <- data.frame(meangroup[1:(m-1)])
    } else {
      centroids[,j] <- meangroup[1:(m-1)]
    }
    
    max <- max(meangroup[1:(m-1)])
    min <- min(meangroup[1:(m-1)])
    nog <- nrow(group)
    
    par(bg = colgenes[j])
    op <- par(mar = rep(0, 4))
    plot(t(1:(m-1)),  t(meangroup[1:(m-1)]),type="l", pch="", lwd=4, xaxt="n", yaxt="n", ann=FALSE, ylim=c(min,max), xlim=c(1,m-1), col="white"
         , main = paste0("Cluster ",clustername[j],": ",nog," genes"))
    
    #axis(1, at=1:(n-1), lab=colnames(data2), las=2, cex.axis=1)
    #axis(1, at=1:(n-1), las=2)
    #axis(2, cex.axis=1);
    #box(col = colgenes[j], lwd=4)
    
    #title(ylab="expression", mgp=c(2.3,1,0))
    
    abline(v = 4.5, lty = 3,lwd=2)
    abline(v = 12.5, lty = 3,lwd=2)
    #abline(v = 19.5, lty = 3)
    
    par(op)
    
    if (TRUE)  { # Add genes
      
      for(i in 1:nrow(group)) {
        
        gene <- scale(t(group[i,1:(m-1)]))
        
        lines(1:(m-1), gene, type="l", lwd=0.1, col="snow3")
      }
      
      lines(t(1:(m-1)),  t(meangroup[1:(m-1)]),type="l", pch="", lwd=4, col="white")
      
    }
    
    dev.off()
  } 
  
  # Print genes
  if (TRUE) {
    genes <- rownames(group)
    f <- paste0("/mnt/picea/home/david/Rdata/SpruceWood/T3heatmap_var1/",filename,clustername[j],".txt")
    write.table(as.matrix(genes), f, quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}












#########################
# Sub-clusters in the same plot!
# B 2 3
# C 4 5 
# D 6 7
# E 8 9
# F 10
# G 11 12
# H 13

i <- 5
j <- 5
clustername <- "H-SUB"

f <- paste0("/mnt/picea/home/david/Rdata/SpruceWood/heatmap_var1/",filename,clustername,".png")
png(file = f)

group <- tableclust[tableclust$ct==i | tableclust$ct==j,]
meangroup <- colMeans(group)

group <- tableclust[tableclust$ct==i,]
meangroupi <- colMeans(group)

group <- tableclust[tableclust$ct==j,]
meangroupj <- colMeans(group)

max <- max(meangroup[1:(m-1)],meangroupi[1:(m-1)],meangroupj[1:(m-1)])
min <- min(meangroup[1:(m-1)],meangroupi[1:(m-1)],meangroupj[1:(m-1)])
    
par(bg = colgenes[i])
op <- par(mar = rep(0, 4))

plot(t(1:(m-1)),  t(meangroup[1:(m-1)]),type="l", pch="", lwd=4, xaxt="n", yaxt="n", ann=FALSE, ylim=c(min,max), xlim=c(1,m), col="black")

lines(t(1:(m-1)),  t(meangroupi[1:(m-1)]),type="l", lty=2, pch="", lwd=4, col="black")

lines(t(1:(m-1)),  t(meangroupj[1:(m-1)]),type="l", lty=2, pch="", lwd=4, col="black")

abline(v = 5.5, lty = 3)
abline(v = 12.5, lty = 3)
abline(v = 19.5, lty = 3)
par(op)
dev.off()

########################
# Clusters in the same plot

f <- paste0("/mnt/picea/home/david/Rdata/SpruceWood/heatmap_var1/cluster_profiles.pdf")
pdf(file = f)

centroids <- scale(centroids)
max <- max(centroids)
min <- min(centroids)

plot(1:(m-1), centroids[,1], type="l", pch="", lwd=3, col = colgenes[1], xaxt="n", yaxt="n", ann=FALSE, ylim=c(min,max), xlim=c(1,m))

abline(v = 5.5, lty = 3)
abline(v = 12.5, lty = 3)
abline(v = 19.5, lty = 3)

axis(1, at=1:(m-1), lab=colnames(data3), las=2, cex.axis=1)
#axis(2, cex.axis=0)
box()

#for(i in 1:ncol(centroids)) {
for(i in c(1,3,4,5,6,7,8)) {
    lines(1:(n-1), centroids[,i], type="l", lwd=3, col = colgenes[i])
}

#cluster_names <- c("h","g","c","d","f","b","e")
#colgenes <- c("mediumpurple3","#CC5151","orange","#B26F2C","lightsteelblue3","#A3CC51","#51A3CC");
cluster_names <- c("b","c","d","e","f","g","h")
cluster_col <- c("#A3CC51","orange","#B26F2C","#51A3CC","lightsteelblue3","#CC5151","mediumpurple3");

#legend(1, 1, c("BLA","BLA2"), lty=c(1,1), lwd=c(2.5,2.5), col = c("red","blue"))
legend(9, 2.2, cluster_names, lty=c(1,1,1,1,1,1,1,1), lwd=c(2.5,2.5), col = cluster_col, ncol=2, bty = "n")

dev.off()

write.table(t(centroids), "/mnt/picea/home/david/Rdata/SpruceWood/cluster_profiles.txt");

# Plot genes
############

n <- 19 # 26 107

#gene_names <- c("Potri.018G101400","Potri.011G114500")
gene_names <- c("Potri.006G216700","Potri.006G002000","Potri.001G401100","Potri.013G014200","Potri.003G088800")

#gene_col <- c("#A3CC51","orange","#B26F2C","#51A3CC","lightsteelblue3","#CC5151","mediumpurple3");
gene_col <- c("#A3CC51","#51A3CC","mediumpurple3","orange","#B26F2C");

#,"","lightsteelblue3","",

profiles <- data.frame()
for(j in 1:length(gene_names)) {
  
  #group <- tableclust[tableclust$ct==j,]
  #meangroup <- colMeans(group)
  #if (j == 1) {
  #  centroids <- data.frame(meangroup[1:(n-1)])
  #} else {
  #  centroids[,j] <- meangroup[1:(n-1)]
  #}
  
  p <- data[data$Genes==gene_names[j],2:n]
  p <- colMeans(p)
  
  if (j == 1) {
    profiles <- data.frame(p[1:(n-1)])
  } else {
    profiles[,j] <- p[1:(n-1)]
  }
}

f <- paste0("/mnt/picea/home/david/Rdata/SpruceWood/gene_profiles_all.pdf")
pdf(file = f)

#profiles <- scale(profiles)
max <- max(profiles)
min <- min(profiles)

plot(1:(n-1), profiles[,1], type="l", pch="", lwd=3, col = gene_col[1], 
     xaxt="n", yaxt="n", ann=FALSE, ylim=c(min,max), xlim=c(1,n), pin = c(20,5))

abline(v = 5.5, lty = 3)
abline(v = 12.5, lty = 3)
abline(v = 19.5, lty = 3)

axis(1, at=1:(n-1), lab=colnames(data3), las=2, cex.axis=1)
axis(2, cex.axis=1)
box()

for(i in 2:length(gene_names)) {
  lines(1:(n-1), profiles[,i], type="l", lwd=3, col = gene_col[i])
}

legend(0.2, 14, gene_names, lty=1, lwd=c(2.5,2.5), col = gene_col, ncol=2, bty = "n")

dev.off()

###############################################################################

# Plot sample tree

library(hyperSpec)

z <- strsplit(as.character(header), "-")
header2 <- unlist(lapply(z,"[",3))

d <- dist(as.matrix(t(data3)), method = "euclidean")
#d <- pearson.dist (t(data3))
hc <- hclust(d , method = "ward")

hc$height <- hc$height / sum(hc$height) # Trick to make mark.dendrogram work
#Ugly hack!
#hc$order[21:49] <- rev(hc$order[21:49])
#hc$order[25:26] <- rev(hc$order[25:26])
#hc$order[42:49] <- rev(hc$order[42:49])

hc <- as.dendrogram(hc)
order <- rep(0, (n-1))
h <- colnames(data3)
z <- strsplit(as.character(h), "[.]")
h <- unlist(lapply(z,"[",2))
h <- gsub("00", "01", as.character(h))
h <- as.integer(h)
k <- 1
for(i in 1:30) {
  for(j in 1:(n-1)) {
    if (h[j] == i) {
      order[j] = k
      k <- k+1
    }
  }
}
hc <- reorder(hc,order, agglo.FUN = mean)
hc <- as.hclust(hc)

pdf(file = "/mnt/picea/home/david/Rdata/SpruceWood/hclust_wood_samples_Kall_Euclidian.pdf")

par (xpd = TRUE, mar = c(2, 2, 2, 2)) # allows plotting into the margin

plot(hc, xlab = NA, ylab = NA, hang=-1, ann = FALSE, cex=0.25, labels = colnames(data3), axes = FALSE)
col <- c("#848484", "#9F5050","#A6BC6A","#864E11")
clusters <- as.factor (as.matrix(header2))
levels (clusters) <- seq_along(unique(clusters))
mark.dendrogram (hc, clusters, border = TRUE, text.col = NA, label=NA, pos.marker = -0.013,  col = col)





















###############################################################################
###############################################################################
###############################################################################

# Plot dendrogram

#d <- dist(as.matrix(t(data3)), method = "euclidean")
d <- as.dist(1-cor(data3))
ord=order.dendrogram(as.dendrogram(hclust(d, method = "ward")))

hc <- hclust(d, method = "ward")

plot(hc, xlab = NA, ylab = NA, hang=-1, ann = FALSE, cex=0.25, labels = as.matrix(header[t(1:n-1)]))

###############################################################################

# Plot level plot

library(lattice)

# order genes
dg=dist(as.matrix(data3), method = "euclidean")

ordg=order.dendrogram(as.dendrogram(hclust(dg, method = "ward")))

# order samples
ds=dist(as.matrix(t(data3)), method = "euclidean")

ords=order.dendrogram(as.dendrogram(hclust(ds, method = "ward")))

# plot heatmap
plot(data3,ylab="Genes",xlab="Samples",pretty=T,aspect=2,
     panel = function(...){
       panel.levelplot(...)
       panel.abline(v = 2.5,lty = 2)
       panel.abline(v = 10.5,lty = 2)
       panel.abline(v = 14.5,lty = 2)
     }, 
     col.regions=redblue.colors(18),
     scales=list(x=list(labels=as.matrix(data3[0,ords]), rot=90,tck=0,cex=0.5)
                 ,y=list(labels=data[ordg,1],tck=0,cex=0.5) # at=seq(0,0,0)
     )
)




