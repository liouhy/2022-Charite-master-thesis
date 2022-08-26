# the data preprocessing follows the instruction from: https://satijalab.org/signac/articles/pbmc_multiomic.html

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SeuratDisk)

# create a Seurat object containing GEX and ATAC data
counts = Read10X_h5('../raw/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5')
fragpath = '../raw/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'

annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) = 'UCSC'

pbmc = CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  min.cells = 20
)

pbmc[["ATAC"]] = CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation,
  min.cells = 15
)

# QC
DefaultAssay(pbmc) = "ATAC"

pbmc = NucleosomeSignal(pbmc)
pbmc = TSSEnrichment(pbmc)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


# dimensionality reduction for GEX
DefaultAssay(pbmc) = "RNA"
pbmc = SCTransform(pbmc)
pbmc = RunPCA(pbmc)
pbmc = RunUMAP(pbmc, dims = 1:50, n.neighbors = 55L,min.dist = 0.45)
pbmc = RunUMAP(pbmc, reduction.name = "umap25", reduction.key = "UMAP25_",
               n.components =25L, dims=1:50, n.neighbors = 55L, min.dist = 0.45)

# dimensionality reduction for ATAC
DefaultAssay(pbmc) = "ATAC"
pbmc = FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc = RunTFIDF(pbmc)
pbmc = RunSVD(pbmc)
pbmc = RunUMAP(pbmc, reduction = "lsi", dims = 2:50, reduction.name="umap_atac",
               reduction.key = "atacUMAP_", n.neighbors = 55L, min.dist = 0.45)
pbmc = RunUMAP(pbmc, reduction = "lsi", dims = 2:50, reduction.name="atac_umap25",
               reduction.key = "atacUMAP25_", n.components =25L, n.neighbors = 55L,
               min.dist = 0.45)

# wnn
pbmc = FindMultiModalNeighbors(pbmc,
                               reduction.list = list("pca", "lsi"),
                               dims.list = list(1:50, 2:50))
pbmc = RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn_umap",
                reduction.key = "wnnUMAP_", n.neighbors = 55L, min.dist = 0.45)

pbmc = RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn_umap25",
               reduction.key = "wnnUMAP25_", n.components =25L, n.neighbors = 55L,
               min.dist = 0.45)

# transfer cell type labels from Hao et al. (2020) since the data set didn't contain the information of cell types
# load PBMC reference
reference = LoadH5Seurat("../raw/pbmc_multimodal.h5seurat")

DefaultAssay(pbmc) = "SCT"

# find anchors between reference and query
transfer_anchors = FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

# transfer cell type labels
predictions = TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc = AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the predicted cell types
Idents(pbmc) = "predicted.id"

# reorder cell type labels
levels(pbmc) = c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                  "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                  "NK Proliferating", "gdT",
                  "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                  "CD14 Mono", "CD16 Mono",
                  "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")


# visualization
VlnPlot(pbmc, features = "SCT.weight", group.by = 'predicted.id', sort = TRUE, pt.size = 0.1) +
  NoLegend()
VlnPlot(pbmc, features = "ATAC.weight", group.by = 'predicted.id', sort = TRUE, pt.size = 0.1) +
  NoLegend()

library(ggplot2)
p1 = DimPlot(pbmc, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 = DimPlot(pbmc, reduction = "umap_atac", group.by = "predicted.id", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 = DimPlot(pbmc, reduction = "wnn_umap", group.by = "predicted.id", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


# save data

saveRDS(pbmc,'../processed/10X_multiome/seurat_obj.rds')

pbmc = readRDS('../processed/10X_multiome/seurat_obj.rds')

write.table(Embeddings(pbmc[['wnn_umap25']]),'../processed/10X_multiome/wnn_umap25.csv',sep=',',
            col.names=FALSE,row.names=FALSE)

write.table(Embeddings(pbmc[['umap25']]),'../processed/10X_multiome/RNA_umap25.csv',sep=',',
            col.names=FALSE,row.names=FALSE)

write.table(Embeddings(pbmc[['atac_umap25']]),'../processed/10X_multiome/ATAC_umap25.csv',sep=',',
            col.names=FALSE,row.names=FALSE)


write.table(Embeddings(pbmc[['wnn_umap']]),'../processed/10X_multiome/wnn_umap.csv',sep=',',
            col.names=FALSE,row.names=FALSE)

write.table(Embeddings(pbmc[['umap']]),'../processed/10X_multiome/RNA_umap.csv',sep=',',
            col.names=FALSE,row.names=FALSE)

write.table(Embeddings(pbmc[['umap_atac']]),'../processed/10X_multiome/ATAC_umap.csv',sep=',',
            col.names=FALSE,row.names=FALSE)

write.csv(pbmc[[]], '../processed/10X_multiome/metadata.csv')

# change seurat object into h5ad format so it can be read in python
DefaultAssay(pbmc) = "RNA"
SaveH5Seurat(pbmc, filename = "../processed/10X_multiome/pbmc_rna.h5Seurat")
Convert("../processed/10X_multiome/pbmc_rna.h5Seurat", dest = "h5ad")


DefaultAssay(pbmc) = "ATAC"
SaveH5Seurat(pbmc, filename = "../processed/10X_multiome/pbmc_atac.h5Seurat")
Convert("../processed/10X_multiome/pbmc_atac.h5Seurat", dest = "h5ad")


