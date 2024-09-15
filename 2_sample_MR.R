library(data.table)
library(dplyr)
library(TwoSampleMR) 
library(MRInstruments)

trait1.trait2.outcome <- read_outcome_data(
  filename = "~/outcome.tsv",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "freq",
  pval_col = "p",
  samplesize_col = "n",
)

trait1.exposure_data <- read_exposure_data(
  "~/protein.txt.",
  sep = '\t',
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p",
  samplesize_col = "n"
)
trait1.exposure_data$exposure = "protein"
trait1.exposure_data = trait1.exposure_data[trait1.exposure_data$pval.exposure < 5e-8,]
if(length(trait1.exposure_data$SNP) == 0){
  next;
}
dat <- trait1.exposure_data[,c("SNP","pval.exposure")]
colnames(dat)[1:2]<-c("rsid","pval")
dat<-mutate(dat,pval = as.numeric(dat$pval)) 
expo_clump <- ld_clump(dat = dat, 
                       clump_r2 = 0.1, 
                       clump_kb = 5000,
                       plink_bin = "～/R/x86_64-pc-linux-gnu-library/4.2/plinkbinr/bin/plink_Linux",
                       bfile = "～/1000G/EUR")

trait1.exposure_data_clumped <- trait1.exposure_data[trait1.exposure_data$SNP %in% expo_clump$rsid,]
trait1.exposure_data_clumped <- mutate(trait1.exposure_data_clumped,F= beta.exposure^2/se.exposure^2) 
trait1.exposure_data_clumped <- subset(trait1.exposure_data_clumped,F>=10) 

trait1.trait2.outcome_data <- trait1.trait2.outcome[trait1.trait2.outcome$SNP %in% trait1.exposure_data_clumped$SNP,]
trait1.trait2.dat <- harmonise_data(exposure_dat = trait1.exposure_data_clumped, outcome_dat = trait1.trait2.outcome_data)
trait1.trait2.results <- mr(trait1.trait2.dat)  
trait1.trait2.dat <-steiger_filtering(trait1.trait2.dat)  
trait1.trait2.dat <- subset(trait1.trait2.dat,steiger_dir=="TRUE")
trait1.trait2.results.withOR = generate_odds_ratios(trait1.trait2.results)
trait1.trait2.dat.run_mr_presso = run_mr_presso(trait1.trait2.dat)
write.table(trait1.trait2.results,"～/protein_res",sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)

#fdr
all_res<-fread("~/all_res")
res_w<-subset(all_res,nsnp==1)
res_ivw<-subset(all_res,nsnp>1)
res_ivw2<-subset(res_ivw,method=="Inverse variance weighted")
res_rbind<-rbind(res_w,res_ivw2)
res_rbind<-arrange(res_rbind,res_rbind$pval)
res<-mutate(res_rbind,fdr=p.adjust(res_rbind$pval,method = "fdr",n=length(res_rbind$pval)))
resfdr<-subset(res,fdr<0.05)
write.table(resfdr,"~/result_fdr",sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
