install.packages("devtools")
devtools::install_github("BioSenior/ggvolcano")
library(ggVolcano)
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)

deg_data<-fread("~/all_res")
deg_data <- deg_data[deg_data$method == "Inverse variance weighted" | deg_data$method == "Wald ratio", ]
#                row   baseMean log2FoldChange     lfcSE       stat       pvalue         padj
#GCR1           GCR1  7201.5782       2.244064 0.2004959  11.192564 4.434241e-29 2.153711e-25
data <- add_regulate(deg_data, log2FC_name = "b",
                     fdr_name = "pval",log2FC = 0, fdr = 0.05)
#                row   baseMean log2FoldChange     lfcSE       stat       pvalue         padj regulate
#GCR1           GCR1  7201.5782       2.244064 0.2004959  11.192564 4.434241e-29 2.153711e-25       Up
colnames(data)[4]<-"geneName"
colnames(data)[7]<-"Beta"
colnames(data)[9]<-"P-value"

EnhancedVolcano(data,
                lab = data$geneName,
                x = 'Beta',y = 'P-value',xlab = 'Beta',ylab = '-Log10(P-value)',
                xlim = c(-2, 2),
                ylim = c(0,9),
                title=NULL,
                subtitle = NULL,
                pCutoff = 0.001172275,
                FCcutoff = 0,
                pointSize = 1.5,
                labSize = 3.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                legendLabels = NULL,
                selectLab = c('NPL','RET','SERPINA10','UROD','TGOLN2','PRRT3',
                              'IL12RB1','LAG3','POSTN','CDH15','QPCTL','CPM',
                              'FBN2','DBI','CEACAM16','CXCL11','RCN1','LRRC37A2',
                              'SH2B3','VNN2','SBSN','CSF2','NUCB2','ADH7','ORM1',
                              'VWA1','CCN2','PLA2G2A','NCR3LG1','PER3','MRI1',
                              'TNFAIP8','MAPK9','SHANK3','INPP5B','OSTN','IFI30',
                              'CTRB2','CHMP6','FUT10','FABP2','PCDH9','PCOLCE2',
                              'CD14','TFRC','ELANE','DNAJB6','TIMD4','EGLN1',
                              'NAGLU','HABP4','NME3','ENTPD6','CLEC7A','SIGLEC12',
                              'SCG3','SNX15','CRELD1','C1QTNF9','PSG5','F3','IDI2',
                              'MAPRE1','SIGLEC14','FSHB'), 
                drawConnectors = TRUE,
                vlineWidth = 0.4,
                col = c('royalblue', 'grey30', 'royalblue', 'red2'))