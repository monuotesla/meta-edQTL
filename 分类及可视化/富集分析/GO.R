rm(list=ls())
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
install.packages("devtools")
BiocManager::install("clusterProfiler")
install.packages("openxlsx")
install.packages("readxl")
library(openxlsx)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
setwd("C:\\Users\\dell\\Desktop")
#####
data_hESCs_KO <- read.xlsx("hESCs-WT2-vs-hESCs-KO-Hom.all.annot2.xlsx")
data_hESCs_Q <- read.xlsx("hESCs-WT2-vs-hESCs-Q268R-Hom.all.annot2.xlsx")
data_hESCs_Y <- read.xlsx("hESCs-WT2-vs-hESCs-Y240X-Hom.all.annot2.xlsx")
data_NPC_KO <- read.xlsx("NPC-WT2-vs-NPC-KO-Hom.all.annot.xlsx")
data_NPC_Y <- read.xlsx("NPC-WT2-vs-NPC-Y240X-Hom.all.annot.xlsx")
#####
data_hESCs_KO_filter <- data_hESCs_KO[data_hESCs_KO$FDR < 0.05, ]
data_hESCs_KO_filter <- data_hESCs_KO_filter[data_hESCs_KO_filter$`log2(fc)` > 1, ]
data_hESCs_Q_filter <- data_hESCs_Q[data_hESCs_Q$FDR < 0.05, ]
data_hESCs_Q_filter <- data_hESCs_Q_filter[data_hESCs_Q_filter$`log2(fc)` > 1, ]
data_hESCs_Y_filter <- data_hESCs_Y[data_hESCs_Y$FDR < 0.05, ]
data_hESCs_Y_filter <- data_hESCs_Y_filter[data_hESCs_Y_filter$`log2(fc)` > 1, ]
data_NPC_KO_filter  <- data_NPC_KO[data_NPC_KO$FDR < 0.05, ]
data_NPC_KO_filter <- data_NPC_KO_filter[data_NPC_KO_filter$`log2(fc)` > 1, ]
data_NPC_Y_filter  <- data_NPC_Y[data_NPC_Y$FDR < 0.05, ]
data_NPC_Y_filter <- data_NPC_Y_filter[data_NPC_Y_filter$`log2(fc)` > 1, ]
#####
entrezIDs1=mget(data_hESCs_KO_filter$Symbol,org.Hs.egSYMBOL2EG,ifnotfound=NA)
entrezIDs2=mget(data_hESCs_Q_filter$Symbol,org.Hs.egSYMBOL2EG,ifnotfound=NA)
entrezIDs3=mget(data_hESCs_Y_filter$Symbol,org.Hs.egSYMBOL2EG,ifnotfound=NA)
entrezIDs4=mget(data_NPC_KO_filter$Symbol,org.Hs.egSYMBOL2EG,ifnotfound=NA)
entrezIDs5=mget(data_NPC_Y_filter$Symbol,org.Hs.egSYMBOL2EG,ifnotfound=NA)
entrezIDs1 = as.character(entrezIDs1)
entrezIDs1 <- gsub("NA", NA, entrezIDs1)
entrezIDs1 <- na.omit(entrezIDs1)
entrezIDs2 = as.character(entrezIDs2)
entrezIDs2 <- gsub("NA", NA, entrezIDs2)
entrezIDs2 <- na.omit(entrezIDs2)
entrezIDs3 = as.character(entrezIDs3)
entrezIDs3 <- gsub("NA", NA, entrezIDs3)
entrezIDs3 <- na.omit(entrezIDs3)
entrezIDs4 = as.character(entrezIDs4)
entrezIDs4 <- gsub("NA", NA, entrezIDs4)
entrezIDs4 <- na.omit(entrezIDs4)
entrezIDs5 = as.character(entrezIDs5)
entrezIDs5 <- gsub("NA", NA, entrezIDs5)
entrezIDs5 <- na.omit(entrezIDs5)
#####################data_hESCs_KO_filter#####
pdf("data_hESCs_KO_filter.pdf", width = 20, height = 15)
common_NPC_cc=enrichGO(gene = entrezIDs1,
                       OrgDb = org.Hs.eg.db,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       ont="ALL",
                       readable =T)
dotplot(common_NPC_cc, split="ONTOLOGY",showCategory = 10,label_format=50) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free", space = "free") +  # 将ONTOLOGY放在行上
  theme(panel.spacing = unit(0.2, "lines"),  # 调整面板间距
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 10))  # 调整面板标
dev.off()
write.xlsx(common_NPC_cc, file = "data_hESCs_KO_filter_GO.xlsx")
#kg
pdf("data_hESCs_KO_filter_kg.pdf", width = 10, height = 8)
gene=`gene_list_all(1)`
kk <- enrichKEGG(gene$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff =1,
                 qvalueCutoff =1)
dotplot(kk)
write.xlsx(kk, file = "data_hESCs_KO_filter_KEGG.xlsx")
dev.off()
#####################data_hESCs_Q_filter#####
pdf("data_hESCs_Q_filter.pdf", width = 20, height = 15)
common_NPC_cc=enrichGO(gene = entrezIDs2,
                       OrgDb = org.Hs.eg.db,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       ont="ALL",
                       readable =T)
dotplot(common_NPC_cc, split="ONTOLOGY",showCategory = 10,label_format=50) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free", space = "free") +  # 将ONTOLOGY放在行上
  theme(panel.spacing = unit(0.2, "lines"),  # 调整面板间距
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 10))  # 调整面板标
dev.off()
write.xlsx(common_NPC_cc, file = "data_hESCs_Q_filter_GO.xlsx")
#kg
pdf("data_hESCs_Q_filter_kg.pdf", width = 10, height = 8)
kk <- enrichKEGG(entrezIDs2,
                 organism = "hsa",
                 pvalueCutoff =1,
                 qvalueCutoff =1)
dotplot(kk)
dev.off()
write.xlsx(kk, file = "data_hESCs_Q_filter_KEGG.xlsx")
############data_hESCs_Y_filter
pdf("data_hESCs_Y_filter.pdf", width = 20, height = 15)
common_NPC_cc=enrichGO(gene = entrezIDs3,
                       OrgDb = org.Hs.eg.db,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       ont="ALL",
                       readable =T)
dotplot(common_NPC_cc, split="ONTOLOGY",showCategory = 10,label_format=50) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free", space = "free") +  # 将ONTOLOGY放在行上
  theme(panel.spacing = unit(0.2, "lines"),  # 调整面板间距
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 10))  # 调整面板标
dev.off()
write.xlsx(common_NPC_cc, file = "data_hESCs_Y_filter_GO.xlsx")
#kg
pdf("data_hESCs_Y_filter_kg.pdf", width = 10, height = 8)
kk <- enrichKEGG(entrezIDs3,
                 organism = "hsa",
                 pvalueCutoff =1,
                 qvalueCutoff =1)
dotplot(kk)
dev.off()
write.xlsx(kk, file = "data_hESCs_Y_filter_KEGG.xlsx")
###########data_NPC_KO_filter###############
pdf("data_NPC_KO_filter.pdf", width = 20, height = 15)
common_NPC_cc=enrichGO(gene = entrezIDs4,
                       OrgDb = org.Hs.eg.db,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       ont="ALL",
                       readable =T)
dotplot(common_NPC_cc, split="ONTOLOGY",showCategory = 10,label_format=50) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free", space = "free") +  # 将ONTOLOGY放在行上
  theme(panel.spacing = unit(0.2, "lines"),  # 调整面板间距
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 10))  # 调整面板标
dev.off()
write.xlsx(common_NPC_cc, file = "data_NPC_KO_filter_GO.xlsx")
#kg
pdf("data_NPC_KO_filter_kg.pdf", width = 10, height = 8)
kk <- enrichKEGG(entrezIDs4,
                 organism = "hsa",
                 pvalueCutoff =1,
                 qvalueCutoff =1)
dotplot(kk)
dev.off()
write.xlsx(kk, file = "data_NPC_KO_filter_KEGG.xlsx")
#################data_NPC_Y_filter
pdf("data_NPC_Y_filter.pdf", width = 20, height = 15)
common_NPC_cc=enrichGO(gene = entrezIDs5,
                       OrgDb = org.Hs.eg.db,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       ont="ALL",
                       readable =T)
dotplot(common_NPC_cc, split="ONTOLOGY",showCategory = 10,label_format=70) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free", space = "free") +  # 将ONTOLOGY放在行上
  theme(panel.spacing = unit(0.2, "lines"),  # 调整面板间距
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 10))  # 调整面板标
dev.off()
write.xlsx(common_NPC_cc, file = "data_NPC_Y_filter_GO.xlsx")
#kg
pdf("data_NPC_Y_filter_kg.pdf", width = 10, height = 8)
kk <- enrichKEGG(entrezIDs5,
                 organism = "hsa",
                 pvalueCutoff =1,
                 qvalueCutoff =1)
dotplot(kk)
dev.off()
write.xlsx(kk, file = "data_NPC_Y_filter_KEGG.xlsx")
#########################################




gene1=as.data.frame(cbind(data_hESCs_KO_filter$Symbol,entrezID=entrezIDs1))
gene2=as.data.frame(cbind(data_hESCs_Q_filter$Symbol,entrezID=entrezIDs2))
gene3=as.data.frame(cbind(data_hESCs_Y_filter$Symbol,entrezID=entrezIDs3))
gene4=as.data.frame(cbind(data_NPC_KO_filter$Symbol,entrezID=entrezIDs4))
gene5=as.data.frame(cbind(data_NPC_Y_filter$Symbol,entrezID=entrezIDs5))
gene1 = gene1[entrezIDs!="NA",]
gene2 = gene2[entrezIDs!="NA",]
gene3 = gene3[entrezIDs!="NA",]
gene4 = gene4[entrezIDs!="NA",]
gene5 = gene5[entrezIDs!="NA",]
gene1$entrezIDs <- gsub('c\\("([0-9]+)", "[0-9]+"\\)', '\\1', gene1$entrezID)
gene2$entrezIDs <- gsub('c\\("([0-9]+)", "[0-9]+"\\)', '\\1', gene2$entrezID)
gene3$entrezIDs <- gsub('c\\("([0-9]+)", "[0-9]+"\\)', '\\1', gene3$entrezID)
gene4$entrezIDs <- gsub('c\\("([0-9]+)", "[0-9]+"\\)', '\\1', gene4$entrezID)
gene5$entrezIDs <- gsub('c\\("([0-9]+)", "[0-9]+"\\)', '\\1', gene5$entrezID)

pdf("data_hESCs_KO_filter", width = 20, height = 15)
common_NPC_cc=enrichGO(gene = entrezIDs5,
                       OrgDb = org.Hs.eg.db,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       ont="ALL",
                       readable =T)
dotplot(common_NPC_cc, split="ONTOLOGY",showCategory = 10,label_format=50) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free", space = "free") +  # 将ONTOLOGY放在行上
  theme(panel.spacing = unit(0.2, "lines"),  # 调整面板间距
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 10))  # 调整面板标
dev.off()














entrezIDs = as.character(entrezIDs)
gene=as.data.frame(cbind(gene,entrezID=entrezIDs))
gene = gene[entrezIDs!="NA",]
gene$entrezIDs <- gsub('c\\("([0-9]+)", "[0-9]+"\\)', '\\1', gene$entrezID)
done

gene=read.table(file="edsite_gene")
gene=read.table(file="s2e")
entrezIDs=mget(data_hESCs_KO_filter$Symbol,org.Hs.egSYMBOL2EG,ifnotfound=NA)
entrezIDs = as.character(entrezIDs)
gene=as.data.frame(cbind(gene,entrezID=entrezIDs))
gene = gene[entrezIDs!="NA",]
# write.table(gene,"gene.txt",quote = F,sep = "\t") 
# gene = read.table("gene.txt",sep = "\t",header = T,check.names = 1,row.names = 1)
gene$entrezIDs <- gsub('c\\("([0-9]+)", "[0-9]+"\\)', '\\1', gene$entrezID)
pdf("path/to/your/directory/dotplot_output.pdf", width = 10, height = 8)
common_NPC_cc=enrichGO(gene = entrezIDs1,
                    OrgDb = org.Hs.eg.db,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    ont="ALL",
                    readable =T)
write.xlsx(go_results, file = "kegg_results.xlsx")

dotplot(common_NPC_cc, split="ONTOLOGY",showCategory = 10,label_format=50) +
  facet_grid(rows = vars(ONTOLOGY), scales = "free", space = "free") +  # 将ONTOLOGY放在行上
  theme(panel.spacing = unit(0.2, "lines"),  # 调整面板间距
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 10))  # 调整面板标题字体大小



kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism = "hsa",
                 pvalueCutoff =1,
                 qvalueCutoff =1)
dotplot(kk)
go_results <- as.data.frame(kk)


devtools::install_github("YuLab-SMU/clusterProfiler")




dotplot(common_NPC_cc)
cnetplot(common_NPC,categorySize="pvalue", foldChange=gene$entrezID,colorEdge = TRUE,shadowtext="none")
write.table(overlap2,"overlap2.txt",quote = F,row.names = F,col.names = F)
p=phyper(1, 155487, 364624, 137680, lower.tail=F)
phyper(10-1, 16, 13, 11, lower.tail=F)
520111 155487 137680  34817
1-phyper(34817, 155487, 364624, 137680)



ego_BP <- enrichGO(gene = s2e$ENTREZID, 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 1, 
                   qvalueCutoff = 1, 
                   readable = TRUE)

ego_MF <- enrichGO(gene = s2e$ENTREZID, 
                   OrgDb = org.Hs.eg.db, 
                   ont = "MF", 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 1, 
                   qvalueCutoff = 1, 
                   readable = TRUE)

ego_CC <- enrichGO(gene = s2e$ENTREZID, 
                   OrgDb = org.Hs.eg.db,
                   ont = "CC", 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 1, 
                   qvalueCutoff = 1, 
                   readable = TRUE)

# 将三个富集结果合并为一个数据框
combined_results <- rbind(
  cbind(ego_BP@result, ontology = "BP"),
  cbind(ego_MF@result, ontology = "MF"),
  cbind(ego_CC@result, ontology = "CC")
)
combined_results <- combined_results[, c("Description", "p.adjust", "ontology")]
colnames(combined_results) <- c("TERM", "p.adjust", "ontology")
# 使用ggplot2绘制综合点图
ggplot(combined_results, aes(x = -log10(p.adjust), y = reorder(TERM, -log10(p.adjust)), color = ontology)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("BP" = "blue", "MF" = "red", "CC" = "green")) +
  labs(title = "GO Enrichment Analysis",
       x = "-log10(Adjusted P-value)",
       y = "GO Terms",
       color = "Ontology") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 16, face = "bold"))
