#Wls Knockout Experiment
#Chulan Kwon Laboratory
#Primary code authors: Suraj Kannan, Matthew Miyamoto
#January 26, 2021
#Document: Code + Figures 3.1

library(Matrix)
library(Seurat)
library(DropletUtils)
library(ggplot2)
library(ggthemes)
library(scales)
library(dplyr)
library(deMULTIplex)
library(patchwork)
library(clipr)
library(reshape2)
library(stringr)
library(monocle)
library(tradeSeq)
library(pheatmap)

#Set the working directory
#Please note that the uploaded files on Synapse contain the counts tables from Kallisto|Bustools, but not the fastqs for the MultiSeq/deMULTIplex barcode assignment steps! If you want to do that step, you will have to download the appropriate fastqs (off of GEO or wherever you get your fastqs) and put them in the working directory.
setwd("~/Documents/Research/Wls/FinalFiles")

#####Prepare 03/04/20 Dataset#####

###Load in data

data = readMM("MMiy030420/cells_x_genes.mtx")
rownames(data) = read.table("MMiy030420/cells_x_genes.barcodes.txt", as.is = TRUE)$V1
colnames(data) = read.table("MMiy030420/cells_x_genes.genes.txt", as.is = TRUE)$V1
data = t(data)
pheno = data.frame(row.names = colnames(data), depth = colSums(data), genes = colSums(data > 0))
bc_rank = barcodeRanks(data)
qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
           label = c("knee", "inflection"), color = c("blue", "green")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count")
cell_droplets = rownames(pheno)[pheno$depth >= metadata(bc_rank)$inflection]
data = data[, cell_droplets]
pheno = pheno[cell_droplets, ]
ttg = read.table("transcripts_to_genes.txt", as.is = TRUE)
genes = ttg[match(rownames(data), ttg$V2), ]$V3
data = data[!duplicated(genes), ]
genes = genes[!duplicated(genes)]
rownames(data) = genes

###Demultiplex barcodes
#Note: This is a computationally intensive step. What I did was clear out everything except for cell_droplets, and then run the script for demultiplexing the barcodes. I then kept the final calls, cleared out everything else, and reloaded the data back in. Please note that this code was taken from Chris McGinnis, except that I have vindictively replaced all arrows with equals signs.

cell.id.vec = cell_droplets
bar.ref = read.table("~/NGSdir/MMiy102219-MultiTag/multi_barcodes.txt", as.is = TRUE)$V1
readTable = MULTIseq.preProcess(R1 = 'MMiy-789-BC_S10_L001_R1_001.fastq.gz', R2 = 'MMiy-789-BC_S10_L001_R2_001_trimmed.fastq.gz', cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8))
bar.table = MULTIseq.align(readTable, cell.id.vec, bar.ref)
bar.table.full = bar.table
bar.table = bar.table.full[, seq(19, 30)]
bar.table_sweep.list = list()
n = 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n = n + 1
  bar.table_sweep.list[[n]] = classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] = paste("q=",q,sep="")
}

# Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 = findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
round1.calls = classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells = names(round1.calls)[which(round1.calls == "Negative")]
bar.table = bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list = list()
n = 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n = n + 1
  bar.table_sweep.list[[n]] = classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] = paste("q=",q,sep="")
}

threshold.results2 = findThresh(call.list=bar.table_sweep.list)
round2.calls = classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells = c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
bar.table = bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list = list()
n = 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n = n + 1
  bar.table_sweep.list[[n]] = classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] = paste("q=",q,sep="")
}

threshold.results2 = findThresh(call.list=bar.table_sweep.list)
round2.calls = classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells = c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
final.calls = c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls) = c(names(round2.calls),neg.cells)

rm(bar.table, bar.table_sweep.list, readTable, threshold.results1, threshold.results2, n, neg.cells, q, round1.calls, round2.calls, bar.table.full)

#Here would be where to load back in the rest of the data
final.calls.order = final.calls[colnames(data)]
pheno$barcode = final.calls.order

pheno$timepoint = "e9.5"
pheno$condition = "wildtype"
pheno[pheno$barcode %in% c("Bar24", "Bar25", "Bar26", "Bar27", "Bar28", "Bar29", "Bar30"), ]$timepoint = "e8.0"
pheno[pheno$barcode %in% c("Bar19", "Bar20", "Bar21", "Bar25", "Bar27", "Bar28", "Bar30"), ]$condition = "knockout"
pheno$run = "030420"

#####Prepare 10/21/19 Dataset 1#####

###Load in data

data2 = readMM("MMiy102219-10xscRNA/cells_x_genes.mtx")
rownames(data2) = read.table("MMiy102219-10xscRNA/cells_x_genes.barcodes.txt", as.is = TRUE)$V1
colnames(data2) = read.table("MMiy102219-10xscRNA/cells_x_genes.genes.txt", as.is = TRUE)$V1
data2 = t(data2)
pheno2 = data.frame(row.names = colnames(data2), depth = colSums(data2), genes = colSums(data2 > 0))
bc_rank = barcodeRanks(data2)
qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
           label = c("knee", "inflection"), color = c("blue", "green")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count")
cell_droplets2 = rownames(pheno2)[pheno2$depth >= metadata(bc_rank)$inflection]
data2 = data2[, cell_droplets2]
pheno2 = pheno2[cell_droplets2, ]
ttg = read.table("transcripts_to_genes.txt", as.is = TRUE)
genes = ttg[match(rownames(data2), ttg$V2), ]$V3
data2 = data2[!duplicated(genes), ]
genes = genes[!duplicated(genes)]
rownames(data2) = genes

###Demultiplex barcodes
#As above. Note that for this run, we had to make a minor tweak to the classifyCells code, since one barcode dominated. Please use the attached custom classifyCells.

cell.id.vec = cell_droplets2
bar.ref = read.table("multi_barcodes.txt", as.is = TRUE)$V1
readTable = MULTIseq.preProcess(R1 = 'e9_5_S1_R1_001.fastq.gz', R2 = 'e9_5_S1_R2_001.fastq.gz', cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8))
bar.table = MULTIseq.align(readTable, cell.id.vec, bar.ref)
bar.table.full = bar.table
bar.table = as.data.frame(bar.table.full[, seq(1,6)])
bar.table_sweep.list = list()
n = 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n = n + 1
  bar.table_sweep.list[[n]] = classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] = paste("q=",q,sep="")
}

# Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 = findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
round1.calls = classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells = names(round1.calls)[which(round1.calls == "Negative")]
bar.table = bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list = list()
n = 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n = n + 1
  bar.table_sweep.list[[n]] = classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] = paste("q=",q,sep="")
}

threshold.results2 = findThresh(call.list=bar.table_sweep.list)
round2.calls = classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells = c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
bar.table = bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list = list()
n = 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n = n + 1
  bar.table_sweep.list[[n]] = classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] = paste("q=",q,sep="")
}

threshold.results2 = findThresh(call.list=bar.table_sweep.list)
round2.calls = classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells = c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
final.calls2 = c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls2) = c(names(round2.calls),neg.cells)

rm(bar.table, bar.table_sweep.list, readTable, threshold.results1, threshold.results2, n, neg.cells, q, round1.calls, round2.calls, bar.table.full)

#Here would be where to load back in the rest of the data
final.calls.order = final.calls2[colnames(data2)]
pheno2$barcode = final.calls.order

pheno2$timepoint = "e9.5"
pheno2$condition = "wildtype"
pheno2[pheno2$barcode %in% c("Bar1", "Bar2", "Bar3"), ]$timepoint = "e8.5"
pheno2[pheno2$barcode %in% c("Bar1", "Bar5"), ]$condition = "knockout"
pheno2$run = "102219-1"

#####Prepare 10/21/19 Dataset 2#####

###Load in data

data3 = readMM("MMiy102219-10xscRNA2/cells_x_genes.mtx")
rownames(data3) = read.table("MMiy102219-10xscRNA2/cells_x_genes.barcodes.txt", as.is = TRUE)$V1
colnames(data3) = read.table("MMiy102219-10xscRNA2/cells_x_genes.genes.txt", as.is = TRUE)$V1
data3 = t(data3)
pheno3 = data.frame(row.names = colnames(data3), depth = colSums(data3), genes = colSums(data3 > 0))
bc_rank = barcodeRanks(data3[, colSums(data3) > 105]) #the shape of this curve was weird without filtering and thus it seemed like the estimation of the inflection point wasn't correct; I think this inflection point is closer to what we were expecting
qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
           label = c("knee", "inflection"), color = c("blue", "green")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count")
cell_droplets3 = rownames(pheno3)[pheno3$depth >= metadata(bc_rank)$inflection]
data3 = data3[, cell_droplets3]
pheno3 = pheno3[cell_droplets3, ]
ttg = read.table("transcripts_to_genes.txt", as.is = TRUE)
genes = ttg[match(rownames(data3), ttg$V2), ]$V3
data3 = data3[!duplicated(genes), ]
genes = genes[!duplicated(genes)]
rownames(data3) = genes

###Demultiplex barcodes
#As above.

cell.id.vec = cell_droplets3
bar.ref = read.table("multi_barcodes.txt", as.is = TRUE)$V1
readTable = MULTIseq.preProcess(R1 = 'e9_5-2_S2_R1_001.fastq.gz', R2 = 'e9_5-2_S2_R2_001.fastq.gz', cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8))
bar.table = MULTIseq.align(readTable, cell.id.vec, bar.ref)
bar.table.full = bar.table
bar.table = as.data.frame(bar.table.full[, seq(7,27)])
bar.table = bar.table[, -c(5, 6, 7, 8, 9, 10, 11, 16)] #These are low quality barcodes that could be readily removed from the data (in fact, they need to be removed for good classification)
bar.table_sweep.list = list()
n = 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n = n + 1
  bar.table_sweep.list[[n]] = classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] = paste("q=",q,sep="")
}

# Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 = findThresh(call.list=bar.table_sweep.list)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
round1.calls = classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells = names(round1.calls)[which(round1.calls == "Negative")]
bar.table = bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list = list()
n = 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n = n + 1
  bar.table_sweep.list[[n]] = classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] = paste("q=",q,sep="")
}

threshold.results2 = findThresh(call.list=bar.table_sweep.list)
round2.calls = classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells = c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
bar.table = bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list = list()
n = 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n = n + 1
  bar.table_sweep.list[[n]] = classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] = paste("q=",q,sep="")
}

threshold.results2 = findThresh(call.list=bar.table_sweep.list)
round2.calls = classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
neg.cells = c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
final.calls3 = c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls3) = c(names(round2.calls),neg.cells)
final.calls3 = final.calls3[!duplicated(names(final.calls3))]

#Here would be where to load back in the rest of the data
final.calls.order = final.calls3[colnames(data3)]
pheno3$barcode = final.calls.order

pheno3$timepoint = "e9.5"
pheno3$condition = "wildtype"
pheno3[pheno3$barcode %in% c("Bar15", "Bar16", "Bar18"), ]$timepoint = "e7.5"
pheno3[pheno3$barcode %in% c("Bar8", "Bar10"), ]$condition = "knockout"
pheno3[pheno3$barcode %in% c("Bar19", "Bar20", "Bar21", "Bar23", "Bar24", "Bar25", "Bar26", "Bar27"), ]$condition = "+/+"
pheno3$run = "102219-2"

#####Workspace Cleanup and Data Combination#####
rm(bc_rank, ttg, bar.ref, cell_droplets, cell_droplets2, cell_droplets3, cell.id.vec, final.calls.order, genes, classifyCells)
colnames(data) = paste(colnames(data), "_1", sep = "")
colnames(data2) = paste(colnames(data2), "_2", sep = "")
colnames(data3) = paste(colnames(data3), "_3", sep = "")
rownames(pheno) = colnames(data)
rownames(pheno2) = colnames(data2)
rownames(pheno3) = colnames(data3)
full_data = cbind(data, data2, data3)
full_pheno = rbind(pheno, pheno2, pheno3)
full_data = full_data[, !full_pheno$barcode %in% c("Doublet", "Negative") & !full_pheno$condition == "+/+"]
full_pheno = full_pheno[colnames(full_data), ]
full_pheno$full_condition = paste(full_pheno$timepoint, full_pheno$condition)

#####Create Full Seurat Object and QC#####
full_seurat = CreateSeuratObject(counts = full_data, meta.data = full_pheno, min.cells = 3, min.features = 200)
full_seurat[["percent.mt"]] = PercentageFeatureSet(full_seurat, pattern = "^mt-")
VlnPlot(full_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
full_seurat = subset(full_seurat, subset = nFeature_RNA > 2500 & nFeature_RNA < 9000 & percent.mt < 22 & nCount_RNA < 62500)
good_cells = colnames(full_seurat)

#####Integrated Seurat Object#####
#We integrate the three datasets using the Seurat SCTransform workflow. Note that this code is taken more or less directly from the Seurat vignette, except again I have replaced arrows with equals signs because I am a bad person.
options(future.globals.maxSize = 4000 * 1024^2)
seurat.list = SplitObject(full_seurat, split.by = "run")
for (i in 1:length(seurat.list)) {
  seurat.list[[i]] = NormalizeData(seurat.list[[i]], verbose = FALSE)
  seurat.list[[i]] = FindVariableFeatures(seurat.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}
for (i in 1:length(seurat.list)) {
  seurat.list[[i]] = SCTransform(seurat.list[[i]], verbose = FALSE)
}
seurat.features = SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
seurat.list = PrepSCTIntegration(object.list = seurat.list, anchor.features = seurat.features, 
                                    verbose = FALSE)
seurat.anchors = FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", 
                                           anchor.features = seurat.features, verbose = FALSE)
seurat.integrated = IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)
seurat.integrated = RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated = RunUMAP(seurat.integrated, dims = 1:30)
DimPlot(seurat.integrated, group.by = "run")
seurat.integrated = FindNeighbors(seurat.integrated, dims = 1:30)
seurat.integrated = FindClusters(seurat.integrated, resolution = 1.2)
DimPlot(seurat.integrated, label = TRUE)
rm(seurat.anchors, full_seurat) #clear some memory

#For sake of visualization, we remove the RBCs. We then recluster and identify clusters here. Our strategy was to effectively "overcluster" - e.g. generate more clusters than we expected to be biologically relevant, and subsequently merge based on marker expression. Our rationale was that "unbiased" clustering isn't often fully biological, particularly for small numbers of clusters that might separate based on genes that aren't as biologically relevant. By generating more clusters, we can gain the benefits of unbiased clustering but then use biological intuition to merge clusters.

seurat.integrated.clean = seurat.integrated[, !seurat.integrated$seurat_clusters %in% c(3, 9, 23, 24, 30, 31, 33)] #RBC clusters
seurat.integrated.clean = RunPCA(seurat.integrated.clean, verbose = FALSE)
seurat.integrated.clean = RunUMAP(seurat.integrated.clean, dims = 1:30)
#We used the default integration vignette parameters above, and in particular on the UMAP plot; visually, using anywhere from 10-30 PCs didn't change the broad results of the UMAP plot but 30 appeared visually pleasing. However, probably only the first 10-15 PCs contribute effectively to the biology (e.g. from the elbow plot), and including more PCs in the clustering seemed to produce some abiological clusters. So we used 15 dims below.
seurat.integrated.clean = FindNeighbors(seurat.integrated.clean, dims = 1:15)
seurat.integrated.clean = FindClusters(seurat.integrated.clean, resolution = 2.0) #Resolution set very high to produce large numbers of clusters
DimPlot(seurat.integrated.clean, label = TRUE, pt.size = 1)
label = c("aSHF", "PhM", "Forelimb", "EC", "PhM", "pSHF", "Epicardium", "pSHF", "Dermamyotome", "Sclerotome", "Hindlimb", "Forelimb", "Sclerotome", "Epicardium", "Forelimb", "OT CM", "Extraembryonic", "CM", "Forelimb", "PhM", "PhM", "Sclerotome", "Mesenchyme", "Dermamyotome", "FHF", "Hindlimb", "PhM", "pSHF", "Forelimb", "aSHF", "Extraembryonic", "Forelimb", "Hindlimb", "pSHF", "Extraembryonic")
names(label) = seq(0, 34)
seurat.integrated.clean$base_label = label[seurat.integrated.clean$seurat_clusters]
DimPlot(seurat.integrated.clean, group.by = "base_label", label = TRUE, pt.size = 1)
DimPlot(seurat.integrated.clean, group.by = "condition", pt.size = 1)
DefaultAssay(seurat.integrated.clean) = "SCT"

#####Differential Gene Expression Testing#####
#While we have a candidate heatmap, we might also want a heatmap that is based on unbiased markers. To this end, we compute top genes as per the code in the Seurat vignette.
idents = seurat.integrated.clean$base_label
seurat.integrated.clean@active.ident = factor(idents)
seurat.markers <- FindAllMarkers(seurat.integrated.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top5 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

#We now use Seurat's diff gene testing to identify differentially regulated genes in each tissue and at each timepoint. We only identify diff genes if there are at least 10 cells in the category.
seurat.integrated.clean$condition = factor(seurat.integrated.clean$condition, levels = c("wildtype", "knockout"))
diff_genes = list()
for(i in unique(seurat.integrated.clean$base_label)){
  print(i)
  sub = seurat.integrated.clean[, seurat.integrated.clean$base_label == i]
  DefaultAssay(sub) = "RNA"
  test = character()
  if(sum(sub$full_condition == "e9.5 wildtype") >= 10 & sum(sub$full_condition == "e9.5 knockout") >= 10){
    print("e9.5")
    test = FindMarkers(sub[, sub$timepoint == "e9.5"], group.by = "condition", ident.1 = "wildtype", ident.2 = "knockout", min.pct = 0.25, only.pos = FALSE) #e9.5 diff genes
    test = test[test$p_val_adj < 0.05, ]
  }
  test2 = character()
  if(sum(sub$full_condition == "e8.5 wildtype") >= 10 & sum(sub$full_condition == "e8.5 knockout") >= 10){
    print("e8.5")
    test2 = FindMarkers(sub[, sub$timepoint == "e8.5"], group.by = "condition", ident.1 = "wildtype", ident.2 = "knockout", min.pct = 0.25, only.pos = FALSE) #e8.5 diff genes
    test2 = test2[test2$p_val_adj < 0.05, ]
  }
  test3 = character()
  if(sum(sub$full_condition == "e8.0 wildtype") >= 10 & sum(sub$full_condition == "e8.0 knockout") >= 10){
    print("e8.0")
    test3 = FindMarkers(sub[, sub$timepoint == "e8.0"], group.by = "condition", ident.1 = "wildtype", ident.2 = "knockout", min.pct = 0.25, only.pos = FALSE) #e8.0 diff genes
    test3 = test3[test3$p_val_adj < 0.05, ]
  }
  diff_genes[[i]] = list(e9.5 = test, e8.5 = test2, e8.0 = test3)
}
rm(i, test, test2, test3, sub)

#We now combine to make one dataframe with the number of differentially expressed genes for each celltype at each timepoint
diff_no = data.frame(row.names = unique(seurat.integrated.clean$base_label))
diff_no$e8_0 = ""
diff_no$e8_5 = ""
diff_no$e9_5 = ""
diff_no[names(unlist(lapply(diff_genes, function(l) nrow(l[[3]])))), ]$e8_0 = unlist(lapply(diff_genes, function(l) nrow(l[[3]])))
diff_no[names(unlist(lapply(diff_genes, function(l) nrow(l[[2]])))), ]$e8_5 = unlist(lapply(diff_genes, function(l) nrow(l[[2]])))
diff_no[names(unlist(lapply(diff_genes, function(l) nrow(l[[1]])))), ]$e9_5 = unlist(lapply(diff_genes, function(l) nrow(l[[1]])))
diff_no$celltype = rownames(diff_no)
diff_no = melt(diff_no, id.vars = "celltype")
diff_no$variable = str_replace(diff_no$variable, "_", ".")
diff_no$value = as.numeric(diff_no$value)

#This allows us to check how many differentially expressed genes are commonly expressed across all celltypes
diff_count = table(unlist(lapply(diff_genes, function(l) rownames(l[[1]]))))
diff_count = sort(diff_count, decreasing = TRUE)
conserved_diff = names(diff_count)[diff_count >= 7]

#This allows us to check how many differentially expressed genes are unique to one celltype. singles gets a list with all of the genes (this was for exploration), while singles_counts gets the count. We add this diff_no for plotting later.
singles = list()
singles_count = vector()
for(i in names(diff_genes)){
  singles[i] = list(names(diff_count)[diff_count <= 1 & names(diff_count) %in% rownames(diff_genes[[i]][["e9.5"]])])
  singles_count[i] = sum(names(diff_count)[diff_count <= 1] %in% rownames(diff_genes[[i]][["e9.5"]]))
}
diff_no$singles = singles_count[diff_no$celltype]

####Celltype Percentage Analysis#####
#Our study is likely underpowered to detect anything other than very large differences in cell number per celltype, but we tried anyway. To ensure some semblence of reasonableness, we only took animals that had >300 cells and only used e8.0/e9.5 (since we only had n =1, 2 for the e8.5 conditions).

cell_pct = vector()
seurat.integrated.clean$base_label = factor(seurat.integrated.clean$base_label)
for(i in names(table(seurat.integrated.clean$barcode))[table(seurat.integrated.clean$barcode) > 300]){
  sub = seurat.integrated.clean[, seurat.integrated.clean$barcode == i]
  cell_pct = rbind(cell_pct, table(sub$base_label))
}
cell_pct = as.data.frame(cell_pct)
cell_pct = sweep(cell_pct, 1, rowSums(cell_pct), "/")
cell_pct$barcode = names(table(seurat.integrated.clean$barcode))[table(seurat.integrated.clean$barcode) > 300]
barcode_condition = seurat.integrated.clean@meta.data[!duplicated(seurat.integrated.clean@meta.data$barcode), ]$full_condition
names(barcode_condition) = seurat.integrated.clean@meta.data[!duplicated(seurat.integrated.clean@meta.data$barcode), ]$barcode
cell_pct = melt(cell_pct)
cell_pct$full_condition = barcode_condition[cell_pct$barcode]
cell_pct$timepoint = sapply(strsplit(cell_pct$full_condition, " "), "[[", 1)
cell_pct$condition = sapply(strsplit(cell_pct$full_condition, " "), "[[", 2)
cell_pct$condition = factor(cell_pct$condition, levels = c("wildtype", "knockout"))

#####Trajectory Reconstruction#####
#We perform trajectory reconstruction in Monocle 2, followed by differential gene expression analysis with Trade-seq. We selected Monocle 2 because it appeared to give the most biologically sensible trajectories as compared to a couple of others we tested (Monocle 3, Slingshot + PCA). With Monocle 2, there is a decision point in what genes to give as input for ordering. A totally unbiased way would be to use dpFeature, which uses genes differentially expressed across clusters on tSNE. Here, however, we can admit to being interested in a particular biological comparison - wildtype vs. knockout; therefore, we used differentially expressed genes between wildtype and knockout at e9.5. Note, however, that the trajectories were often similar to those generated with dpFeature, leading us to be confident that our approach wasn't overtly biased.

#We focused on aSHF and pSHF for the sake of our analysis.

###pSHF
trajectory_pshf = newCellDataSet(full_data[rownames(seurat.integrated.clean), colnames(seurat.integrated.clean)[seurat.integrated.clean$base_label %in% c("pSHF")]], phenoData = new("AnnotatedDataFrame", data = seurat.integrated.clean@meta.data[colnames(seurat.integrated.clean)[seurat.integrated.clean$base_label %in% c("pSHF")], ]), featureData = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(seurat.integrated.clean), gene1 = rownames(seurat.integrated.clean), gene2 = rownames(seurat.integrated.clean))), lowerDetectionLimit = 0.1, expressionFamily = negbinomial.size())
trajectory_pshf = estimateSizeFactors(trajectory_pshf)
trajectory_pshf = estimateDispersions(trajectory_pshf)
#It would take forever to run diff expression on all genes, and a lot of them are likely not informative. So we used two steps - first, filter out only reasonably highly expressed genes, then selected the top 8000 most dispersed genes (arbitrary cutoff).
trajectory_pshf = detectGenes(trajectory_pshf, min_expr = 0.1)
disp_table = dispersionTable(trajectory_pshf)
unsup_clustering_genes = subset(disp_table, mean_expression > 0.1)
unsup_clustering_genes = unsup_clustering_genes[order(unsup_clustering_genes$dispersion_empirical, decreasing = TRUE), ]
clustering_DEG_genes = differentialGeneTest(trajectory_pshf[as.character(head(unsup_clustering_genes$gene_id, 8000)), trajectory_pshf$timepoint == "e9.5"], fullModelFormulaStr = "~condition + run", reducedModelFormulaStr = "~run", cores = 1)
ordering_genes = row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:3000] #We then used 3000 most differentially expressed genes to build trajectory (again, arbitrary cutoff)
trajectory_pshf = setOrderingFilter(trajectory_pshf, ordering_genes = ordering_genes)
trajectory_pshf = reduceDimension(trajectory_pshf, reduction_method = 'DDRTree', residualModelFormulaStr = "~run")
trajectory_pshf = orderCells(trajectory_pshf) #by default this should give the correct ordering
#So one thing to observe is that, by default, Monocle assigns different path lengths and thus different max pseudotimes to both branches. Monocle understands this, and stretches the pseudotimes when plotting genes. However, Trade-Seq doesn't do this, which in a couple of cases results in what are spurious diff gene calls; also, I rather like plotting the smoothers from Trade-Seq. We corrected the pseudotimes by stretching the short branch to match the long branch.
trajectory_pshf$fixed_pseudotime = trajectory_pshf$Pseudotime
pData(trajectory_pshf)[trajectory_pshf$State == 2, ]$fixed_pseudotime = pData(trajectory_pshf)[trajectory_pshf$State == 2, ]$Pseudotime * max(pData(trajectory_pshf)[trajectory_pshf$State == 3, ]$Pseudotime) / max(pData(trajectory_pshf)[trajectory_pshf$State == 2, ]$Pseudotime)
trajectory_pshf$old_pseudotime = trajectory_pshf$Pseudotime
trajectory_pshf$Pseudotime = trajectory_pshf$fixed_pseudotime
#Differential gene expression analysis
sce_pshf <- fitGAM(trajectory_pshf[fData(trajectory_pshf)$num_cells_expressed > (0.2 * ncol(trajectory_pshf)), ], verbose = TRUE) #Limited here to genes expressed in at least 20% of cells; again, arbitrary threshold, this was mostly just to reduce my annoyance at plotting "differentially expressed genes" that were barely expressed
#Trade-seq is a really neat package with plenty of options for differential gene expression testing. However, at least in this case, end testing was basically all we needed.
endRes_pshf = diffEndTest(sce_pshf)
endRes_pshf = endRes_pshf[order(endRes_pshf$waldStat, decreasing = TRUE), ]
endRes_pshf$padj = p.adjust(endRes_pshf$pvalue, method = "fdr")

###aSHF
trajectory_ashf = newCellDataSet(full_data[rownames(seurat.integrated.clean), colnames(seurat.integrated.clean)[seurat.integrated.clean$base_label %in% c("aSHF")]], phenoData = new("AnnotatedDataFrame", data = seurat.integrated.clean@meta.data[colnames(seurat.integrated.clean)[seurat.integrated.clean$base_label %in% c("aSHF")], ]), featureData = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(seurat.integrated.clean), gene1 = rownames(seurat.integrated.clean), gene2 = rownames(seurat.integrated.clean))), lowerDetectionLimit = 0.1, expressionFamily = negbinomial.size())
trajectory_ashf = estimateSizeFactors(trajectory_ashf)
trajectory_ashf = estimateDispersions(trajectory_ashf)
trajectory_ashf = detectGenes(trajectory_ashf, min_expr = 0.1)
disp_table = dispersionTable(trajectory_ashf)
unsup_clustering_genes = subset(disp_table, mean_expression > 0.1)
unsup_clustering_genes = unsup_clustering_genes[order(unsup_clustering_genes$dispersion_empirical, decreasing = TRUE), ]
clustering_DEG_genes = differentialGeneTest(trajectory_ashf[as.character(head(unsup_clustering_genes$gene_id, 8000)), trajectory_ashf$timepoint == "e9.5"], fullModelFormulaStr = "~condition + run", reducedModelFormulaStr = "~run", cores = 1)
ordering_genes = row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:3000]
trajectory_ashf = setOrderingFilter(trajectory_ashf, ordering_genes = ordering_genes)
trajectory_ashf = reduceDimension(trajectory_ashf, reduction_method = 'DDRTree', residualModelFormulaStr = "~run")
trajectory_ashf = orderCells(trajectory_ashf, root_state = 6)
pData(trajectory_ashf)[trajectory_ashf$State %in% c(1, 2, 7), ]$fixed_pseudotime = pData(trajectory_ashf)[trajectory_ashf$State %in% c(1, 2, 7), ]$Pseudotime * max(pData(trajectory_ashf)[trajectory_ashf$State == 4, ]$Pseudotime) / max(pData(trajectory_ashf)[trajectory_ashf$State %in% c(1, 2, 7), ]$Pseudotime)
trajectory_ashf$old_pseudotime = trajectory_ashf$Pseudotime
trajectory_ashf$Pseudotime = trajectory_ashf$fixed_pseudotime
#Differential gene expression analysis
#Unlike the pSHF, the aSHF has a couple of spurrious branches that really aren't of much interest to us. We could in theory keep them for Trade-seq, but truthfully they are a bit annoying. So I manually pruned and removed here.
info <- extract_monocle_info(trajectory_ashf)
info$cellWeights = info$cellWeights[, c("Y_63", "Y_103")]
info$pseudotime = info$pseudotime[, c("Y_63", "Y_103") ]
info$cellWeights = info$cellWeights[rowSums(info$cellWeights) > 0, ]
info$pseudotime = info$pseudotime[rownames(info$cellWeights), ]
sce_ashf <- fitGAM(counts = Biobase::exprs(trajectory_ashf[fData(trajectory_ashf)$num_cells_expressed > (0.2 * ncol(trajectory_ashf)), rownames(info$cellWeights)]), cellWeights = info$cellWeights, pseudotime = info$pseudotime, verbose = TRUE)
rm(info)
endRes_ashf = diffEndTest(sce_ashf)
endRes_ashf = endRes_ashf[order(endRes_ashf$waldStat, decreasing = TRUE), ]
endRes_ashf$padj = p.adjust(endRes_ashf$pvalue, method = "fdr")
rm(disp_table, unsup_clustering_genes, clustering_DEG_genes, ordering_genes)

#Yes, I am aware of the arguments in terms of cautious use of the p-value with Trade-Seq. For our interests, however, they will suffice, at least to gear us towards putatively interesting results. I haven't thought too deeply about the utility of testing for the log2FC or filtering after testing against a threshold of 1 (which is what we did) - I suppose the former better models fold change uncertainty but also may be too conservative.

pshf_diff = rownames(endRes_pshf)[endRes_pshf$padj < 0.05 & abs(endRes_pshf$logFC1_2) >= 0.8]
ashf_diff = rownames(endRes_ashf)[endRes_ashf$padj < 0.05 & abs(endRes_ashf$logFC1_2) >= 0.8]

###Clustering genes and making heatmaps

#Trade-seq offers an option for clustering using RSEC but personally, after trying it, I didn't find it particularly neat. I took a simpler approach, which was to extract the smoothers from Trade-seq and use the standard clustering in pheatmap to generate my clusters. The end heatmaps look rather like the Monocle 2 BEAM heatmaps, just using the smoothers from Trade-seq instead. There is also code here for generating kinetic maps from the clusters.

#For the sake of space, I didn't recopy the code for ashf and pshf; the two are more or less identical, just with changing the variable names and one tweak noted below to get the order of the clusters correct.

ysmooth <- predictSmooth(models = sce_ashf, gene = ashf_diff, nPoints = 100, tidy = FALSE)
ysmooth2 = ysmooth[, c(seq(100, 1), seq(101, 200))]
breaks = c(seq(-6, -3), seq(-2, 2, 0.1), seq(3, 6))
cols = c(rep("lightgrey", 4), colorRampPalette(c("lightgrey", "blue"))(41), rep("blue", 4))
test = pheatmap(log(ysmooth2 + 1), scale = "row", cluster_cols = FALSE, show_rownames = FALSE, gaps_col = 100, color = cols, breaks = breaks)
tree = cutree(test$tree_row, k = 5)
rm(test)
#The next couple of lines are done to set the orders for the clusters. Uncomment the first three lines for pSHF, or use the last four for aSHF.
#tree[tree == 2] = 5 #for pshf
#tree[tree == 1] = 2 #for pshf
#tree[tree == 5] = 1 #for pshf
tree[tree == 3] = "a"
tree[tree == 2] = "b"
tree[tree == 4] = "c"
tree[tree %in% c(1, 5)] = "d"
#Generating kinetic maps
ysmooth_tidy <- predictSmooth(models = sce_ashf, gene = ashf_diff, nPoints = 100)
plots = list()
ysmooth_tidy$yhat_norm = 0
for(g in unique(ysmooth_tidy$gene)){
  ysmooth_tidy[ysmooth_tidy$gene == g, ]$yhat_norm = as.numeric(scale(log(ysmooth_tidy[ysmooth_tidy$gene == g, ]$yhat + 1)))
}
for (xx in sort(unique(tree))) {
  p = ggplot()
  p <- p + geom_density2d(data = ysmooth_tidy[ysmooth_tidy$gene %in% names(tree)[tree == xx], ], aes(x = time, yhat_norm, col = as.character(lineage), group = lineage), lwd = 0.7/fig_factor)
  p <- p + guides(color = FALSE) + scale_color_manual(values = c("royalblue", "indianred1"), breaks = c("0", "1")) + labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") + theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 0), axis.text = element_text(size = 0), title = element_text(size = 10), axis.ticks = element_line(size= 0))
  plots[[xx]] = p
}
do.call("grid.arrange", c(plots, ncol=2)) #The list plots() should have a list of each individiual kinetic plot; we simply extracted these for the final figures.

#####Code for Making Pretty Figures#####
fig_factor = 1 #1 for plotting for computer; 2.5 for making figures

###Integrated Cleaned Plot by Celltype
DimPlot(seurat.integrated.clean, group.by = "base_label", label = TRUE, pt.size = 1/fig_factor, label.size = 4.5/fig_factor) + xlab("UMAP Dim. 1") + ylab("UMAP Dim. 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + guides(color = guide_legend(title = "Celltype", override.aes = list(size = 4/fig_factor)))

###Integrated Cleaned Plot by Timepoint
#Note - I took out the e7.5 ones here because there was barely any of them and the labeling could confuse
DimPlot(seurat.integrated.clean[, seurat.integrated.clean$timepoint != "e7.5"], group.by = "timepoint", label = FALSE, pt.size = 1/fig_factor, label.size = 4.5/fig_factor) + xlab("UMAP Dim. 1") + ylab("UMAP Dim. 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + guides(color = guide_legend(title = "Timepoint", override.aes = list(size = 4))) + scale_color_manual(values = c("green3", "darkorange", "darkorchid3"))

###Integrated Cleaned Plot by Condition
seurat.integrated.clean$condition = factor(seurat.integrated.clean$condition, levels = c("wildtype", "knockout"))
labels = c("Control", "Knockout")
names(labels) = c("wildtype", "knockout")
DimPlot(seurat.integrated.clean, group.by = "condition", label = FALSE, pt.size = 1/fig_factor, label.size = 4.5/fig_factor) + xlab("UMAP Dim. 1") + ylab("UMAP Dim. 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + guides(color = guide_legend(title = "Condition", override.aes = list(size = 4))) + scale_color_manual(label = labels, values = c("royalblue", "indianred1"))

###Integrated Cleaned Plot by Run
DimPlot(seurat.integrated.clean, group.by = "run", label = FALSE, pt.size = 1/fig_factor, label.size = 4.5/fig_factor) + xlab("UMAP Dim. 1") + ylab("UMAP Dim. 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 15, l = 0)), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + guides(color = guide_legend(title = "Run", override.aes = list(size = 4)))

###Integrated Heatmap or Dotplot using Candidate Markers
seurat.integrated.clean$order_label = factor(seurat.integrated.clean$base_label, levels = c("Dermamyotome", "Sclerotome", "PhM", "aSHF", "pSHF", "FHF", "CM", "OT CM", "Epicardium", "Mesenchyme", "Forelimb", "Hindlimb", "EC", "Extraembryonic"))
DoHeatmap(seurat.integrated.clean, group.by = "order_label", features = c("Fst", "Meox1", "Pax3", "Pax1", "Nr2f2", "Ebf1", "Tbx1", "Fgf8", "Fgf10", "Isl1", "Osr1", "Wnt2", "Sfrp5", "Tbx5", "Tnnt2", "Nkx2-5", "Myh6", "Wt1", "Tbx18", "Postn", "Vim", "Pitx1", "Tshz1", "Tbx3", "Twist1", "Tbx4", "Msx1", "Msx2", "Hoxb8" ,"Hoxc6", "Kdr", "Tek"), lines.width = 50, draw.lines = TRUE, group.colors = hue_pal()(14)[match(levels(seurat.integrated.clean$order_label), levels(factor(seurat.integrated.clean$base_label)))])

seurat.integrated.clean$order_label = factor(seurat.integrated.clean$base_label, levels = rev(c("Dermamyotome", "Sclerotome", "PhM", "aSHF", "pSHF", "FHF", "CM", "OT CM", "Epicardium", "Mesenchyme", "Forelimb", "Hindlimb", "EC", "Extraembryonic")))
DotPlot(seurat.integrated.clean, group.by = "order_label", features = rev(c("Fst", "Meox1", "Pax3", "Pax1", "Nr2f2", "Ebf1", "Tbx1", "Fgf8", "Fgf10", "Isl1", "Osr1", "Wnt2", "Sfrp5", "Tbx5", "Tnnt2", "Nkx2-5", "Myh6", "Wt1", "Tbx18", "Postn", "Vim", "Pitx1", "Tshz1", "Tbx3", "Twist1", "Tbx4", "Msx1", "Msx2", "Hoxb8" ,"Hoxc6", "Kdr", "Tek")), dot.scale = 6/fig_factor) + RotatedAxis() + theme(axis.text.x = element_text(size = 18/fig_factor), axis.text.y = element_text(size = 16/fig_factor), legend.position = "none")

###Integrated Dotplot for Highly Expressed Wnts
DotPlot(seurat.integrated.clean[, !seurat.integrated.clean$timepoint == "e7.5"], group.by = "order_label", features = rev(c("Wnt2", "Wnt4", "Wnt5a", "Wnt5b", "Wnt6", "Wnt11")), dot.scale = 6/fig_factor) + RotatedAxis() + theme(axis.text.x = element_text(size = 18/fig_factor), axis.text.y = element_text(size = 16/fig_factor), legend.position = "none")

###Integrated Heatmap using Unbiased Markers
DoHeatmap(seurat.integrated.clean, features = top5$gene) + NoLegend()

###VlnPlots
VlnPlot(seurat.integrated.clean, features = c("Wls"), group.by = "order_label", split.by = "condition", multi.group = TRUE, pt.size = 0) + scale_fill_manual(labels = labels, values = c("royalblue", "indianred1")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(), axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), title = element_text(size = 28/fig_factor), legend.text = element_text(size = 18/fig_factor))

VlnPlot(seurat.integrated.clean, features = c("Isl1"), group.by = "order_label", pt.size = 0) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(), axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), title = element_text(size = 28/fig_factor), legend.text = element_text(size = 18/fig_factor))

VlnPlot(seurat.integrated.clean[, seurat.integrated.clean$timepoint == "e9.5"], features = c("Ddit4"), group.by = "order_label", split.by = "condition", multi.group = TRUE, pt.size = 0) + scale_fill_manual(labels = labels, values = c("royalblue", "indianred1")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line = element_blank(), axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), title = element_text(size = 28/fig_factor), legend.text = element_text(size = 18/fig_factor))

###Differentially Expressed Genes by Celltype and Timepoint; this code can also be adapted to plot singles
ggplot(diff_no, aes(x = celltype, y = value, fill = variable)) + geom_bar(stat = "identity", position = position_dodge()) + xlab("Celltype") + ylab("Number of DEGs") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_fill_manual(name = "Timepoint", values = c("green3", "darkorange", "darkorchid3"))

###Celltype Percentages
ggplot(cell_pct[!cell_pct$timepoint == "e8.5", ], aes(x = timepoint, y = value, fill = condition)) + geom_boxplot() + geom_point(position = position_jitterdodge()) + facet_wrap(~variable, scale = "free_y") + theme_bw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 22/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), strip.text = element_text(size = 24/fig_factor)) + guides(fill = guide_legend(title = "Condition")) + scale_fill_manual(label = labels, values = c("royalblue", "indianred1")) + ylab("Fraction") + xlab("Timepoint") + expand_limits(y = 0)

###Trajectory plots for pSHF (first by timepoint, then by condition)
plot_cell_trajectory(trajectory_pshf, color_by = "timepoint", cell_size = 1.5/fig_factor, show_branch_points = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Timepoint") + scale_color_manual(values = c("green3", "darkorange", "darkorchid3"))
plot_cell_trajectory(trajectory_pshf, color_by = "condition", cell_size = 1.5/fig_factor, show_branch_points = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Condition") + scale_color_manual(values = c("royalblue", "indianred1"))

###Trajectory plots for aSHF (first by timepoint, then by condition)
plot_cell_trajectory(trajectory_ashf, color_by = "timepoint", cell_size = 1.5/fig_factor, show_branch_points = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Timepoint") + scale_color_manual(values = c("green3", "darkorange", "darkorchid3"))
plot_cell_trajectory(trajectory_ashf, color_by = "condition", cell_size = 1.5/fig_factor, show_branch_points = FALSE) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Condition") + scale_color_manual(values = c("royalblue", "indianred1"))
trajectory_ashf$fixed_pseudotime = trajectory_ashf$Pseudotime

###Individual gene trajectory plots (using smoothers from Trade-seq) - note, to make this code work in this manner, where the colours could be adjusted, I had to fork a branch of Trade-seq and make a tiny tweak. The authors may correct this issue soon, however.
colData(sce_pshf)$condition = pData(trajectory_pshf)[colnames(sce_pshf), ]$condition
colData(sce_pshf)$condition = factor(colData(sce_pshf)$condition, levels = c("wildtype", "knockout"))
colData(sce_ashf)$condition = pData(trajectory_ashf)[colnames(sce_ashf), ]$condition
colData(sce_ashf)$condition = factor(colData(sce_ashf)$condition, levels = c("wildtype", "knockout"))
plotSmoothers(sce_ashf, counts(sce_ashf), "Ddit4", pointCol = "condition", size = 1/fig_factor, lwd = 2/fig_factor) +xlab("Pseudotime") + ylab("Scaled Expression") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_manual(values = c("royalblue", "indianred1", "indianred1", "royalblue"))
plotSmoothers(sce_pshf, counts(sce_pshf), "Ddit4", pointCol = "condition", size = 1/fig_factor, lwd = 2/fig_factor) +xlab("Pseudotime") + ylab("Scaled Expression") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_manual(values = c("royalblue", "indianred1", "indianred1", "royalblue"))

###Two-sided pseudotime heatmaps for aSHF (can be modified to use for pSHF)
pheatmap(log(ysmooth2[names(sort(tree)), ] + 1), scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, gaps_col = 100, gaps_row = cumsum(as.numeric(table(tree))), color = cols, breaks = breaks)

#####Helper Functions#####
#I quickly defined some helper functions during exploratory analysis; they have been left here in case they are useful later.

#all_genes: Returns all of the differentially expressed genes for a given celltype at a certain timepoint
all_genes = function(celltype, timepoint = "e9.5"){
  return(rownames(diff_genes[[celltype]][[timepoint]]))
}

#good_genes: Returns all of the differentially expressed genes for a given celltype at a given timepoint outside of those in the conserved list (e.g. diff exp. in at least 7 celltypes)
good_genes = function(celltype, timepoint = "e9.5"){
  genes = all_genes(celltype, timepoint)
  return(genes[!genes %in% conserved_diff])
}

#find_gene: Finds which tissues a given gene may be differentially expressed in at e9.5
find_gene = function(gene){
  return(names(diff_genes)[unlist(lapply(diff_genes, function(l) (gene %in% rownames(l[[1]]))))])
}
