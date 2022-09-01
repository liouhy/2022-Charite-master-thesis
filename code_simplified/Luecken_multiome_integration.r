# this session follows the instruction from: https://satijalab.org/seurat/articles/integration_large_datasets.html

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)

# create GEX object
expression_matrix <- ReadMtx(
  mtx = '../processed/Luecken_GEX_counts.mtx', features = '../processed/Luecken_GEX_var.tsv',
  cells = '../processed/Luecken_GEX_obs.tsv',feature.column=1,skip.cell=1,skip.feature=1
)
GEX <- CreateSeuratObject(counts = expression_matrix)

obs = read.csv('../processed/Luecken_GEX_obs.tsv',sep='\t')
Idents(GEX) <- obs$batch


# create ATAC assay
atac_counts <- Matrix::readMM('../processed/Luecken_ATAC_counts.mtx')
colnames(atac_counts) = obs[,1]

genome <- read.csv('../processed/Luecken_ATAC.bed',sep='\t',header=FALSE)
vec_g<-paste(genome$V1,genome$V2,genome$V3,sep='-')
grange.counts <- StringToGRanges(vec_g, sep = c(":", "-"))

rownames(atac_counts) = vec_g

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   genome = 'hg38',
   ranges = grange.counts,
   annotation = annotations
 )

# save assay
GEX[['ATAC']] <- chrom_assay

# find integration anchors to correct for batch effect based on GEX
GEX.list <- SplitObject(GEX, split.by = "ident")
GEX.list <- lapply(X = GEX.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = GEX.list)
GEX.list <- lapply(X = GEX.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = GEX.list, reference = c(1, 2), reduction = "rpca",
    dims = 1:50)

# integrate GEX (batch correction)
GEX.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

GEX.integrated <- ScaleData(GEX.integrated, verbose = FALSE)
GEX.integrated <- RunPCA(GEX.integrated, verbose = FALSE)
GEX.integrated <- RunUMAP(GEX.integrated, dims = 1:50, n.neighbors = 55L,min.dist = 0.45)

saveRDS(GEX.integrated,'../processed/Luecken_GEX_int.rds')

GEX.integrated <- readRDS('../processed/Luecken_GEX_int.rds')

# integrate ATAC
# LSI
DefaultAssay(GEX.integrated) <- "ATAC"
GEX.integrated <- RunTFIDF(GEX.integrated)
GEX.integrated <- FindTopFeatures(GEX.integrated, min.cutoff = 'q0')
GEX.integrated <- RunSVD(GEX.integrated)

# integrate ATAC at LSI space using the anchors obtained from GEX
multi.integrated <- IntegrateEmbeddings(
  anchorset = anchors,
  reductions = GEX.integrated[['lsi']],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50)

GEX.integrated[['integrated_lsi']]<-multi.integrated[['integrated_lsi']]

GEX.integrated <- RunUMAP(GEX.integrated,
 reduction = "integrated_lsi",
 dims = 2:50,
 reduction.name="umap_atac",
 reduction.key = "atacUMAP_",
 n.neighbors = 55L,
 min.dist = 0.45)

# create a joint WNN embedding based on both GEX and ATAC
# find WNN and run umap
GEX.integrated <- FindMultiModalNeighbors(GEX.integrated,
 reduction.list = list("pca", "integrated_lsi"),
 dims.list = list(1:50, 2:50))

GEX.integrated <- RunUMAP(GEX.integrated,
 nn.name = "weighted.nn",
 reduction.name = "wnn_umap",
 reduction.key = "wnnUMAP_",
 n.neighbors = 55L,
 min.dist = 0.45)

# 25D UMAP for SEMITONES reference cell selection
GEX.integrated <- RunUMAP(GEX.integrated,
 nn.name = "weighted.nn",
 reduction.name = "wnn_umap25",
 reduction.key = "wnnUMAP25_",
 n.components =25L,
 n.neighbors = 55L,
 min.dist = 0.45)


# also run umap25 for GEX and ATAC for SEMITONES enrichment score calculation
GEX.integrated <- RunUMAP(GEX.integrated,
 verbose = FALSE,
 reduction.name = "umap25",
 reduction.key = "UMAP25_",
 n.components =25L,
 dims=1:50,
 n.neighbors = 55L,
 min.dist = 0.45)

GEX.integrated <- RunUMAP(GEX.integrated,
 reduction = "integrated_lsi",
 dims = 2:50,
 reduction.name="atac_umap25",
 reduction.key = "atacUMAP25_",
 n.components =25L,
 n.neighbors = 55L,
 min.dist = 0.45)

saveRDS(GEX.integrated, '../processed/Luecken_integrated_multiome.rds')

# saving data for python
multi_integrated <- readRDS('../processed/Luecken_integrated_multiome.rds')

write.table(Embeddings(multi_integrated[['wnn_umap25']]),'../processed/Luecken_wnn_umap25.csv',sep=',',
            col.names=FALSE,row.names=FALSE)

write.table(Embeddings(multi_integrated[['umap25']]),'../processed/Luecken_GEX_int_umap25.csv',sep=',',
            col.names=FALSE,row.names=FALSE)

write.table(Embeddings(multi_integrated[['atac_umap25']]),'../processed/Luecken_ATAC_int_umap25.csv',sep=',',
            col.names=FALSE,row.names=FALSE)


write.table(Embeddings(multi_integrated[['wnn_umap']]),'../processed/Luecken_wnn_umap.csv',sep=',',
            col.names=FALSE,row.names=FALSE)

write.table(Embeddings(multi_integrated[['umap']]),'../processed/Luecken_GEX_int_umap.csv',sep=',',
            col.names=FALSE,row.names=FALSE)

write.table(Embeddings(multi_integrated[['umap_atac']]),'../processed/Luecken_ATAC_int_umap.csv',sep=',',
            col.names=FALSE,row.names=FALSE)
