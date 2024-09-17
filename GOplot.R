library(org.Hs.eg.db)
library(clusterProfiler) 
library(ggplot2)
library(RColorBrewer) 
library(patchwork)
library(dplyr)
library(stringr)

ori = data.frame(matrix(nrow = 410, ncol = 1))
colnames(ori)[1]<-"SYMBOL"
a<-fread("~/res")

ori$SYMBOL<-a$a
up_entrez <- bitr(ori$SYMBOL,
                  fromType = "SYMBOL", 
                  toType = "ENTREZID",
                  OrgDb = "org.Hs.eg.db")

GO_MF_diff <- enrichGO(gene = up_entrez$ENTREZID,  
                       OrgDb = org.Hs.eg.db,  
                       ont = "MF",  
                       pAdjustMethod = "BH",  
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)  
GO_MF_result <- GO_MF_diff@result

GO_CC_diff <- enrichGO(gene = up_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
GO_CC_result <- GO_CC_diff@result

GO_BP_diff <- enrichGO(gene = up_entrez$ENTREZI,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
GO_BP_result <- GO_BP_diff@result

write.table(GO_MF_result,file = "~/GO_MF_result_all.txt",col.names = TRUE,row.names = FALSE,sep = "\t",quote = FALSE)
write.table(GO_CC_result,file = "~/GO_CC_result_all.txt",col.names = TRUE,row.names = FALSE,sep = "\t",quote = FALSE)
write.table(GO_BP_result,file = "~/GO_BP_result_all.txt",col.names = TRUE,row.names = FALSE,sep = "\t",quote = FALSE)

MF <- GO_MF_result[1:20,]
CC <- GO_CC_result[1:20,]
BP <- GO_BP_result[1:20,]

MF$Description <- str_trunc(MF$Description,width = 50,side = "right")
MF$Description
MF$term <- factor(MF$Description,levels = rev(MF$Description))
CC$term <- factor(CC$Description,levels = rev(CC$Description))
BP$term <- factor(BP$Description,levels = rev(BP$Description))

rf<- apply(MF,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF<- round(GeneRatio/BgRatio,2)
  RF
})
MF$Rich_Factor <- rf

rf<- apply(CC,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF<- round(GeneRatio/BgRatio,2)
  RF
})
CC$Rich_Factor <- rf

rf<- apply(BP,1,function(x){
  GeneRatio <- eval(parse(text = x["GeneRatio"]))
  BgRatio <- eval(parse(text = x["BgRatio"]))
  RF<- round(GeneRatio/BgRatio,2)
  RF
})
BP$Rich_Factor <- rf

GO_dot <- function(x){
  y = get(x)
  ggplot(data = y,
         aes(x = Rich_Factor,
             y = term)) +
    geom_point(aes(size = Count,
                   color = -log10(pvalue))) + 
    scale_y_discrete(labels = function(y) str_wrap(y, width = 50) ) +
    labs(x = "Rich Factor",
         y = "Description",
         title = paste0(x," of GO enrichment Dotplot"), 
         size = "Count") +  
    theme_bw()+
    mytheme
}
mytheme <- theme_minimal() 
pp1 <- GO_dot("MF") + scale_color_distiller(palette = "YlOrRd",direction = 1)
pp2 <- GO_dot("CC") + scale_color_distiller(palette = "YlGnBu",direction = 1)
pp3 <- GO_dot("BP") + scale_color_distiller(palette = "YlOrBr",direction = 1)
setwd("~/GO_KEGG")
pdf("go_all.pdf",width = 30,height = 10)
print(pp1+pp2+pp3)
dev.off()
