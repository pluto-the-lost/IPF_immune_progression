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
library('ggalluvial')

DC_colors <- c("#66C2A5","#AB98C8","#BDA37D","#A89BB0","#DF8BC3","#E99073")
plot_sank <- function(dada_per, condition, groups, count, colors){
    p <- ggplot(dada_per,
    aes_string(x = condition, stratum = groups, alluvium = groups, y=count,
    fill = groups, label = groups)) +
    scale_x_discrete(expand = c(0, 0)) +
    geom_flow(width = 1/8) + #线跟方块间空隙的宽窄
    geom_stratum(alpha = .9,width = 1/10) + #方块的透明度、宽度
    #geom_text(stat = "stratum", size = 3, color="black") +  #文字大小、颜色
    scale_fill_manual(values = colors) +
    
    xlab("") + ylab("") +
    theme_bw() + #去除背景色
    theme(panel.grid =element_blank()) + #去除网格线
    theme(panel.border = element_blank()) + #去除外层边框
    theme(legend.title = element_blank()) +
    theme(axis.line = element_blank(),axis.ticks = element_blank()) + #去掉坐标轴
    ggtitle("")
    return(p)
}


rds<-readRDS('IPF_dbfree.rds')
Idents(rds)<-rds@meta.data$cell.type
rds@meta.data$orig.ident<-as.vector(rds@meta.data$orig.ident)
options(repr.plot.width=5, repr.plot.height=5)
DC_per <- round(prop.table(table(rds@meta.data[,c("cell.type","orig.ident")]), margin=1),3)
DC_per <- data.frame(DC_per)
#plot_sank(DC_per,"cell.type","orig.ident",,"Freq", DC_colors)
#ggsave("Macrophage_prop1.pdf",w=5,h=5)
Idents(rds)<-rds@meta.data$orig.ident
options(repr.plot.width=5, repr.plot.height=5)
DC_per <- round(prop.table(table(rds@meta.data[,c("orig.ident","cell.type")]), margin=1),3)
DC_per <- data.frame(DC_per)
print (unique(DC_per$orig.ident))
DC_per$orig.ident <- factor(DC_per$orig.ident, levels = c("D0", "D7", "D14","D21"))
plot_sank(DC_per, "orig.ident", "cell.type","Freq", DC_colors)
ggsave("Macrophage_prop2.pdf",w=5,h=5)
