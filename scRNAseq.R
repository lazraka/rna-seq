install.packages("Seurat")
install.packages("tidyverse")
install.packages("cowplot")
install.packages("RColorBrewer")

library(Seurat)
library(tidyverse)
library(cowplot)
library(RColorBrewer)

pbmc.data <- Read10X(data.dir = "")

install.packages("Matrix")
install.packages("ISLR")
install.packages("factoextra")
install.packages("dendextend")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")
install.packages("raster")
BiocManager::install("edgeR")
install.packages("irlba")


library(Matrix)
library(ISLR)
library(factoextra)
library(dendextend)
library(SingleCellExperiment)
library(scater)
library(raster)
library(edgeR)
library(irlba)

matrix.path = "~/GSE132044_cortex_mm10_count_matrix.mtx"

cortex_mat <- readMM(file = matrix.path)

wells.path <- "~/GSE132044_cortex_mm10_cell.tsv"
features.path <- "~/GSE132044_cortex_mm10_gene.tsv"

feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
wells.names = read.delim(wells.path, header = FALSE, stringsAsFactors = FALSE)
colnames(cortex_mat) = wells.names$V1
rownames(cortex_mat) = feature.names$V1
rownames(cortex_mat) = substr(rownames(cortex_mat),start=20, stop=100)

cortex_mat<-as.matrix(cortex_mat)

smartseq_subset<-cortex_mat[,grep("Cortex1.Smart",colnames(cortex_mat))]

tenx_subset<-cortex_mat[,grep("Cortex1.10",colnames(cortex_mat))]
droncseq_subset<-cortex_mat[,grep("Cortex1.DroNc",colnames(cortex_mat))]
sciseq_subset<-cortex_mat[,grep("Cortex1.sci-RNA",colnames(cortex_mat))]

sample_counts_smartseq<-colSums(smartseq_subset)
sample_counts_tenx<-colSums(tenx_subset)
sample_counts_droncseq<-colSums(droncseq_subset)
sample_counts_sciseq<-colSums(sciseq_subset)

par(mfrow=c(2,2))
hist(sample_counts_smartseq, main="Total Counts/Cell - Smartseq", xlab="Counts", breaks="FD")
hist(sample_counts_tenx, main="Total Counts/Cell - 10X", xlab="Counts", breaks="FD")
hist(sample_counts_droncseq, main="Total Counts/Cell - DroNcSeq", xlab="Counts", breaks="FD")
hist(sample_counts_sciseq, main="Total Counts/Cell - SciSeq", xlab="Counts", breaks="FD")

nonzero_smartseq <- rowSums(smartseq_subset) > 0
smartseq_keep<-smartseq_subset[nonzero_smartseq,]
dge_smartseq <- DGEList(smartseq_keep)
dge_smartseq$genes <- rownames(smartseq_keep)
smartseq_norm<-calcNormFactors(dge_smartseq)

nonzero_tenx <- rowSums(tenx_subset) > 0
tenx_keep<-tenx_subset[nonzero_tenx,]
dge_tenx <- DGEList(tenx_keep)
dge_tenx$genes <- rownames(tenx_keep)
tenx_norm<-calcNormFactors(dge_tenx, method="TMM")

nonzero_droncseq <- rowSums(droncseq_subset) > 0
droncseq_keep<-droncseq_subset[nonzero_droncseq,]
dge_droncseq <- DGEList(droncseq_keep)
dge_droncseq$genes <- rownames(droncseq_keep)
droncseq_norm<-calcNormFactors(dge_droncseq)

nonzero_sciseq <- rowSums(sciseq_subset) > 0
sciseq_keep<-sciseq_subset[nonzero_sciseq,]
dge_sciseq <- DGEList(sciseq_keep)
dge_sciseq$genes <- rownames(sciseq_keep)
sciseq_norm<-calcNormFactors(dge_sciseq)

#plot normalization factors
par(mfrow=c(2,2))
hist(smartseq_norm$samples$norm.factors, main="Smart-Seq Normalization", xlab="Normalization Factors", breaks="FD")
abline(v=1, col="green", lwd=3, lty=2)
hist(tenx_norm$samples$norm.factors, main="10X Normalization", xlab="Normalization Factors", breaks="FD")
abline(v=1, col="blue", lwd=3, lty=2)
hist(droncseq_norm$samples$norm.factors, main="DroNcSeq Normalization", xlab="Normalization Factors", breaks="FD")
abline(v=1, col="red", lwd=3, lty=2)
hist(sciseq_norm$samples$norm.factors, main="Sci-Seq Normalization", xlab="Normalization Factors", breaks="FD")
abline(v=1, col="yellow", lwd=3, lty=2)

smartseq_a <- edgeR::aveLogCPM(smartseq_norm)
tenx_a <- edgeR::aveLogCPM(tenx_norm)
tenx_a_norm<-edgeR::aveLogCPM(dge_tenx_norm)

droncseq_a <- edgeR::aveLogCPM(droncseq_norm)
sciseq_a <- edgeR::aveLogCPM(sciseq_norm)

smartseq_avelogcpm_thresh <- 2

hist_smartseq <- list(
  Histogram=ggplot(data.frame(logCPM=smartseq_a)) + aes(x=logCPM) +
    geom_histogram(aes(y=100*(..count..)/sum(..count..)), binwidth=0.25, boundary=0) +
    geom_vline(xintercept=smartseq_avelogcpm_thresh, color="red", linetype="dashed") +
    xlab("Average logCPM") + ylab("Percent of genes in bin") + coord_cartesian(xlim=quantile(smartseq_a, c(0, 0.995)), ylim=c(0,15)) +
    labs(title="Average gene LogCPM distribution - Smartseq")
  )

hist_tenx <- list(
  Histogram=ggplot(data.frame(logCPM=tenx_a_norm)) + aes(x=logCPM) +
    geom_histogram(aes(y=100*(..count..)/sum(..count..)), binwidth=0.05, boundary=0) +
    geom_vline(xintercept=avelogcpm.presence.threshold, color="red", linetype="dashed") +
    xlab("Average logCPM") + ylab("Percent of genes in bin") + coord_cartesian(xlim=quantile(tenx_a, c(0, 0.995)), ylim=c(0,15)) +
    labs(title="Average gene LogCPM distribution - 10X")
)

hist_droncseq <- list(
  Histogram=ggplot(data.frame(logCPM=droncseq_a)) + aes(x=logCPM) +
    geom_histogram(aes(y=100*(..count..)/sum(..count..)), binwidth=0.1, boundary=0) +
    geom_vline(xintercept=avelogcpm.presence.threshold, color="red", linetype="dashed") +
    xlab("Average logCPM") + ylab("Percent of genes in bin") + coord_cartesian(xlim=quantile(droncseq_a, c(0, 0.995)), ylim=c(0,15)) +
    labs(title="Average gene LogCPM distribution - DroNc-Seq")
)

hist_sciseq <- list(
  Histogram=ggplot(data.frame(logCPM=sciseq_a)) + aes(x=logCPM) +
    geom_histogram(aes(y=100*(..count..)/sum(..count..)), binwidth=0.1, boundary=0) +
    geom_vline(xintercept=avelogcpm.presence.threshold, color="red", linetype="dashed") +
    xlab("Average logCPM") + ylab("Percent of genes in bin") + coord_cartesian(xlim=c(8, 12), ylim=c(0,15)) +
    labs(title="Average gene LogCPM distribution - Smartseq")
)

par(mfrow=c(2,2))


gene_counts_smartseq<-rowSums(smartseq_subset)
gene_counts_tenx<-rowSums(tenx_subset)
gene_counts_droncseq<-rowSums(droncseq_subset)
gene_counts_sciseq<-rowSums(sciseq_subset)

#non-normalize
smartseq_cpm<- edgeR::cpm(smartseq_subset)
tenx_cpm<- edgeR::cpm(tenx_subset)
droncseq_cpm<- edgeR::cpm(droncseq_subset)
sciseq_cpm<- edgeR::cpm(sciseq_subset)

#normalized
smartseq_cpm_norm<-edgeR::cpm(smartseq_norm)
tenx_cpm_norm<- edgeR::cpm(tenx_norm)
droncseq_cpm_norm<- edgeR::cpm(droncseq_norm)
sciseq_cpm_norm<- edgeR::cpm(sciseq_norm, log=F)


plot(smartseq_cpm[,1],smartseq_subset[,1])
plot(tenx_cpm[,1],tenx_subset[,1])
plot(droncseq_cpm[,1],droncseq_subset[,1])
plot(sciseq_cpm_norm[,1],sciseq_keep[,1])

par(mfrow=c(1,1))
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(smartseq_cpm[,1],smartseq_subset[,1],ylim=c(0,50),xlim=c(0,20), main="Smartseq CPM Threshold", ylab ="Raw Counts", xlab="CPM") #threshold is 20
plot(tenx_cpm[,3],tenx_subset[,3],ylim=c(0,50),xlim=c(0,500)) #threshold is 500
plot(droncseq_cpm[,10],droncseq_subset[,10],ylim=c(0,50),xlim=c(0,15000), main="DroNcSeq CPM Threshold", ylab ="Raw Counts", xlab="CPM") #threshold is 10,000
plot(sciseq_cpm[,1],sciseq_subset[,1],ylim=c(0,50),xlim=c(0,5000))#threshold is 4000

abline(h=5)
abline(v=20)

thresh_smartseq <- smartseq_cpm > 2
table(rowSums(thresh_smartseq))

thresh_tenx <- tenx_cpm > 500
table(rowSums(thresh_tenx))

thresh_droncseq <- droncseq_cpm > 6000
table(rowSums(thresh_droncseq))

thresh_sciseq <- sciseq_cpm > 2000
table(rowSums(thresh_sciseq))

keep_smartseq <- rowSums(thresh_smartseq) >= 2 # 2 replicates per sample
keep_tenx <- rowSums(thresh_tenx) >= 2 # 2 replicates per sample
keep_droncseq <- rowSums(thresh_droncseq) >= 2 # 2 replicates per sample
keep_sciseq <- rowSums(thresh_sciseq) >= 2 # 2 replicates per sample

counts.keep.smartseq <- smartseq_subset[keep_smartseq,]
counts.keep.tenx<-tenx_subset[keep_tenx,]
counts.keep.droncseq<-droncseq_subset[keep_droncseq,]
counts.keep.sciseq<-sciseq_subset[keep_sciseq,]

summary(keep)
dim(counts.keep.smartseq)
dim(counts.keep.tenx)
dim(counts.keep.droncseq)
dim(counts.keep.sciseq)

#this is for the 10th sample because the plot before was only for one sample
plot(smartseq_cpm[,10],smartseq_subset[,10],ylim=c(0,50),xlim=c(0,20))

dge_smartseq <- DGEList(counts.keep.smartseq)
dge_tenx <- DGEList(counts.keep.tenx)
dge_droncseq <- DGEList(counts.keep.droncseq)
dge_sciseq <- DGEList(counts.keep.sciseq)

dge_tenx_norm<-DGEList(t(tenx_normalized))

barplot(dge_smartseq$samples$lib.size,names=colnames(dge_smartseq),las=2)
title("Barplot of library sizes")

logcounts_smartseq <- edgeR::cpm(dge_smartseq,log=TRUE)
logcounts_tenx <- edgeR::cpm(dge_tenx,log=TRUE)
logcounts_droncseq <- edgeR::cpm(dge_droncseq,log=TRUE)
logcounts_sciseq <- edgeR::cpm(dge_sciseq,log=TRUE)

logcounts_tenx_norm<-edgeR::cpm(dge_tenx_norm,log=TRUE)

boxplot(logcounts_smartseq, xlab="", ylab="Log2 counts per million",las=2)
barplot(logcounts_tenx_norm[,1:50], xlab="", ylab="Log2 counts per million",las=2)

# We estimate the variance for each row in the logcounts matrix
var_genes_smartseq <- apply(logcounts_smartseq, 1, var)
var_genes_tenx <- apply(logcounts_tenx, 1, var)
var_genes_droncseq <- apply(logcounts_droncseq, 1, var)
var_genes_sciseq <- apply(logcounts_sciseq, 1, var)

var_genes_tenx_norm<- apply(dge_tenx_norm, 1, var)
head(var_genes_smartseq)

# Get the gene names for the top 500 most variable genes
select_var_smartseq <- names(sort(var_genes_smartseq, decreasing=TRUE))[1:10]
select_var_tenx <- names(sort(var_genes_tenx, decreasing=TRUE))[1:10]
select_var_droncseq <- names(sort(var_genes_droncseq, decreasing=TRUE))[1:10]
select_var_sciseq <- names(sort(var_genes_sciseq, decreasing=TRUE))[1:10]

select_var_tenx_norm <- names(sort(var_genes_tenx_norm, decreasing=TRUE))[1:10]
head(select_var)

highly_variable_lcpm_smartseq <- logcounts_smartseq[select_var_smartseq,]
highly_variable_lcpm_tenx <- logcounts_tenx[select_var_tenx,]
highly_variable_lcpm_droncseq <- logcounts_droncseq[select_var_droncseq,]
highly_variable_lcpm_sciseq <- logcounts_sciseq[select_var_sciseq,]

highly_variable_tenx_norm <-dge_tenx_norm$counts[select_var_tenx_norm,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm_smartseq)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]

par(mar = c(1, 1, 1, 3))
heatmap.2(highly_variable_lcpm_smartseq,col=rev(morecols(50)),trace="none", main="Smartseq - Top 10 Genes",scale="row", key=TRUE, keysize=1.5)
heatmap.2(highly_variable_lcpm_tenx,col=rev(morecols(50)),trace="none", main="10X - Top 10 Genes",scale="row")
heatmap.2(highly_variable_lcpm_droncseq,col=rev(morecols(50)),trace="none", main="DroNcSeq - Top 10 Genes",scale="row")
heatmap.2(highly_variable_lcpm_sciseq,col=rev(morecols(50)),trace="none", main="Sciseq - Top 10 Genes",scale="row")

heatmap.2(highly_variable_tenx_norm,col=rev(morecols(50)),trace="none", main="10X - Top 10 DEG",scale="row")


geo_lib_size <- colSums(log(smartseq_subset + 1))
barplot(geo_lib_size, xlab = "Cell", ylab = "Geometric Lib Size", las = 2)

lcounts <- log(smartseq_subset + 1)

# PC1 and Geometric library size correlation
pc1 <- irlba(lcounts - rowMeans(lcounts), 1)$v[, 1]
smartseq_cor<-cor(geo_lib_size, pc1)

# create a summarizedexperiment class
summexp<-SummarizedExperiment(assays = smartseq_subset)
seq_type<-c(rep("Smartseq",length(colnames(smartseq_subset))),rep("10X",length(colnames(tenx_subset))),rep("DroNcseq",length(colnames(droncseq_subset))),rep("Sciseq",length(colnames(sciseq_subset))))
seq_type<-as.data.frame(seq_type)


cortex_normalized<-cbind()
tenx_transpose<-t(tenx_subset)
tenx_normalized<-library.size.normalize(tenx_transpose, verbose = FALSE)
