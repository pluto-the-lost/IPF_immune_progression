require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(Hmisc)
library(sctransform)
library(monocle3)
library(monocle)
library(harmony)
library(RColorBrewer)

grid.col<-colorRampPalette(brewer.pal(8, "Set3"))(16)[1:16]
rds<-readRDS('/Volumes/LENOVO_USB_HDD/IPF/B_TET_074/Figure4/Macrophage_new.rds')
sub<-rds
DefaultAssay(object = sub) <- "RNA"
sub<- FindVariableFeatures(object = sub)
all.genes <- rownames(sub)
sub <- ScaleData(sub, features = all.genes, verbose = FALSE)
sub <- RunPCA(sub, verbose = FALSE)
sub <- RunHarmony(sub, group.by.vars = "sampleID")
sub <- RunUMAP(sub, dims = 1:20, verbose = FALSE,reduction = "harmony",reduction.name = 'umap', reduction.key='umap')
sub <- FindNeighbors(sub, dims = 1:20, verbose = FALSE,reduction = "harmony")
sub <- FindClusters(sub, resolution = 1, verbose = FALSE)
current.cluster.ids <- levels(Idents(sub))
new.cluster.ids <- as.numeric(current.cluster.ids)+1
Idents(sub) <- plyr::mapvalues(x = Idents(sub), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(sub))
sub@meta.data$seurat_clusters<-Idents(sub)[rownames(sub@meta.data)]
saveRDS(sub,'Mono.umap.rds')

rds<-sub
DefaultAssay(object = rds) <- "RNA"
Idents(rds)<-rds@meta.data$seurat_clusters
DimPlot(rds, reduction = "umap", group.by = "seurat_clusters",pt.size=2,cols=grid.col,raster=TRUE)
ggsave('Mono_dimplot.pdf', dpi=900,width=5, height=4, units="in")
DimPlot(rds, reduction = "umap", group.by = "seurat_clusters",split.by='sampleID',pt.size =4,cols=grid.col,raster=TRUE,ncol=4)
ggsave('Mono_dimplot_group.pdf', dpi=900,width=16, height=4, units="in")

markers <- FindAllMarkers(object = rds, group.by = 'seurat_clusters',logfc.threshold = 0.25,test.use = "wilcox",min.pct = 0.1,p_val_adj=0.05)
top5 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.table(top5,'Mono.marker.xls',sep='\t',quote=F,col.names=T,row.names=F)
