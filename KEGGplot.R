library(dplyr)
library(org.Hs.eg.db) 
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer) 
a<-fread("~/res")
a<-unique(a$exposure)
a<-as.data.frame(a)

data = data.frame(matrix(nrow = 410, ncol = 1))
colnames(data)[1]<-"SYMBOL"
data$SYMBOL<-a$a
up_entrez <- bitr(data$SYMBOL,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

KEGG_diff <- enrichKEGG(gene = up_entrez$ENTREZID,organism = "hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
KEGG_result <- KEGG_diff@result
View(KEGG_result)

KEGG_top <- KEGG_result[1:20,]

KEGG_top$pathway <- factor(KEGG_top$Description,levels = rev(KEGG_top$Description))
KEGG_top$pathway <- str_wrap(KEGG_top$pathway, width = 40)

mytheme <- theme(axis.title = element_text(size = 12),
                 axis.text = element_text(size = 10),
                 plot.title = element_text(size = 13,
                                           hjust = 0.5,
                                           face = "bold"),
                 legend.title = element_text(size = 12),
                 legend.text = element_text(size = 10))
KEGG_plot <- ggplot(data = KEGG_top,
            aes(x = Count,
                y = pathway,
                fill = -log10(pvalue)))+
  scale_fill_distiller(palette = "RdPu",direction = 1) + 
  geom_bar(stat = "identity",width = 0.8) +
  theme_bw() +
  labs(x = "Count",
       y = "pathway",
       title = "KEGG enrichment barplot") + mytheme
pdf("kegg.pdf")
print(KEGG_plot)
dev.off()

