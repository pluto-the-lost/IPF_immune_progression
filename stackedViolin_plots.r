require(Seurat)
require(dplyr)
require(Matrix)
require(magrittr)
library(scales)
library(ggplot2)
library(configr)
library(cowplot)
library(Hmisc)
library(RColorBrewer)
library(tidyverse)
library(GSA
library(GSA)


rds<-readRDS('/Volumes/LENOVO_USB_HDD/IPF/B_TET_074/AJR_hs_immuneCells_withAnno.rds')
genes<-c('Dpp4','C5ar1')
FeaturePlot(rds,features=genes,ncol =2,cols = c("lightgrey", "#1C00FF"),reduction='umap',raster=TRUE)
ggsave('Fig1_feature_sup.pdf', dpi=900,width=10, height=5, units="in")

rds<-readRDS('/Volumes/LENOVO_USB_HDD/IPF/B_TET_074/Figure4/Macrophage_new.rds')
genes<-c('Cd68'，'Itgam', 'Itgax', 'Ly6c2', 'Siglecf', 'Pparg', 'C1qa', 'C1qb', 'C1qc','Fn1','Cx3cr1')

my36colors<-c('#E31F1B','#498FB3','#6F9D6B','#C86F68')
modify_vlnplot1<- function(obj,
                          features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          cols=my36colors,
                          slot="data",
                          assay="RNA",
                          ...) {
  p<- VlnPlot(obj, features = features,pt.size = pt.size, cols=cols, slot=slot,
              assay=assay,... )  +
    xlab("") + ylab(features) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(0.7), angle = 0, vjust = 0.5,color="black"),
          plot.margin = plot.margin )
  return(p)
}
## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          slot="data",
                          assay="RNA",
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot1(obj = obj, features = x, slot=slot, assay=assay,...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(color="black",size=10,angle=90),axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

rds@meta.data$cell.type <- factor(rds@meta.data$cell.type, levels = c("AM","MonMac","Mono","IntMacs"))
DefaultAssay(object = rds) <- "RNA"
Idents(rds)<-rds@meta.data$cell.type

p<-StackedVlnPlot(obj = rds, features=genes,slot="data",assay="RNA")
pdf('Macrophage.pdf',w=4,h=11)
print(p)
dev.off()

#rds<-readRDS('/Volumes/LENOVO_USB_HDD/IPF/B_TET_074/Figure4/AM.anno.rds')
#tab<-table(rds@meta.data$cell_type,rds@meta.data$sampleID)
#write.table(tab,'AM.number.xls',sep='\t',quote=F)

#rds<-readRDS('/Volumes/LENOVO_USB_HDD/IPF/B_TET_074/Figure4/MonMac.anno.rds')
#table(rds@meta.data$cell_type)
#rds<-readRDS('/Volumes/LENOVO_USB_HDD/IPF/B_TET_074/Figure4/AM.anno.rds')
#table(rds@meta.data$cell.type)

#AM      MonMac  Mono    IntMacs
#3105    1604    489     121
#B    DC   Mac    Neu   T.NKT
#413  1146 14851  3055  1513

#Fabp5+Gpnmd+       Lung resident Ams         CCL17-producing
#  1106                    1071                     401
#Nr4a1+ Fabp1-Car4-Krt79+Cidec+ Fabp1+Car4+Krt79+Cidec+
#370        79                      49
#Isg15+Ly6e+
#29

#Cx3cr1+       CD36+      C5ar1+      F13a1+ Isg15+Ly6a+
#  868         397         389         388          51
