library(gwasglue)
library(coloc)
library(data.table)
library(openxlsx)
library(tidyverse)

outcome <- fread("～/pcos.tsv")  
colnames(outcome) <- c('SNP','chrom',"pos",'effect_allele','other_allele',"eaf","beta","se","P","samplesize","number_cases")
outcome$MAF <- ifelse(outcome$eaf<0.5,outcome$eaf,1-outcome$eaf)
outcome$varbeta <- outcome$se^2
outcome <- subset(outcome, !duplicated(SNP))
outcome$z = outcome$beta/outcome$se
outcome$s <- outcome$number_cases/outcome$samplesize
outcome_GWASdata <- outcome %>% na.omit() %>% filter(P>0)

exposure <- fread("~/protein.txt.gz") 
exposure$end <- exposure$Pos
exposure <- exposure %>% dplyr::select("rsids","Chrom","Pos","end","effectAllele","otherAllele","ImpMAF","Beta","SE","Pval","N")
colnames(exposure) <- c('SNP','chrom',"start","end",'effect_allele','other_allele',"eaf","beta","se","P","samplesize")
exposure<- exposure %>% na.omit() %>% filter(P>0)
exposure$z = exposure$beta/exposure$se
exposure$varbeta <- exposure$se^2
exposure$MAF <- ifelse(exposure$eaf<0.5,exposure$eaf,1-exposure$eaf)

lead <- exposure %>% dplyr::arrange(P)
leadchr <- lead$chrom[1]
leadstart <- as.numeric(lead$start[1])
leadend <- as.numeric(lead$end[1])
pQTLdata <- exposure[exposure$chrom==leadchr,]
pQTLdata <- pQTLdata[pQTLdata$start>leadstart-500000 & pQTLdata$end<leadend+500000,]
pQTLdata <- subset(pQTLdata, !duplicated(SNP))
pQTLdata <- na.omit(pQTLdata)

sameSNP <- intersect(pQTLdata$SNP,outcome_GWASdata$SNP)
pQTLdata <- pQTLdata[pQTLdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
pQTLdata$P<-as.numeric(pQTLdata$P)
GWASdata <- outcome_GWASdata[outcome_GWASdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()

result <- coloc.abf(dataset1=list(pvalues=GWASdata$P, snp=GWASdata$SNP, type="cc", s=GWASdata$s[1], N=GWASdata$samplesize[1]), 
                    dataset2=list(pvalues=pQTLdata$P, snp=pQTLdata$SNP, type="quant", N=pQTLdata$samplesize[1]), MAF=pQTLdata$MAF)

coloc_result <- result$results %>% dplyr::arrange(desc(SNP.PP.H4))
write.table(coloc_result,file = "coloc_result.txt",sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
write.csv(result$summary,file = 'coloc_sum_result.csv')
