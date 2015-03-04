library(dplyr)
library(data.table)
library(foreach)
library(iterators)
add_table=function(x){write.table(x, file="~/Documents/urate_variance/snps.csv",sep="\t", quote=FALSE, append=TRUE)}

tfam <- read.table("~/MurrayXsan/vGWAS/scratch/ARIC/Manual_lm/aric_var_set.tfam", header=FALSE, stringsAsFactors = FALSE)
names(tfam) <- c("SUBJECT", "FID", "V3", "V4", "V5","V6")

allaricfhsforvariance <- read.delim("~/Documents/urate_variance/allaricfhsdataforvariance.txt", stringsAsFactors=FALSE)

pheno <- left_join(x = tfam, y = allaricfhsforvariance, by="SUBJECT")
pheno$SUBJECT <- as.character(pheno$SUBJECT)

#tped1 = read.table("~/MurrayXsan/vGWAS/scratch/ARIC/Manual_lm/aric_var_set.tped", header=FALSE, nrows = 1, stringsAsFactors = FALSE)
tped_file <- ireadLines(file.path("~/MurrayXsan/vGWAS/scratch/ARIC/Manual_lm/aric_var_set.tped"))
i=0
foreach ( line = tped_file) %do% {
  i <- i+1
  tped <- line
  #tped=readLines("~/MurrayXsan/vGWAS/scratch/ARIC/Manual_lm/rs6449173.tped")
  tped <- unlist(strsplit( tped, split = " " ))

  a <- matrix(nrow = length(tfam[,1]), ncol = 2, byrow = TRUE, data = tped[5:length(tped)])
  gt <- paste0(a[,1], a[,2])
  gt <- cbind(tfam[,1], gt)
  gt[gt[,2] == "00",2] <- NA
  gt <- na.omit(gt)
  gt <- as.data.frame(gt)
  b<-table(gt[,2])

  colnames(gt) <- c("SUBJECT", "SNP")
  gt$addsnp = -1
  #assign major/minor allele for additive model
  if((b[1]*2 + b[2])/ (2*sum(b)) >= 0.5){
    gt[gt[,2] == names(b[1]),"addsnp"] <- 0
    gt[gt[,2] == names(b[2]),"addsnp"] <- 1
    gt[gt[,2] == names(b[3]),"addsnp"] <- 2
    major=names(b[1])
  } else{
    gt[gt[,2] == names(b[1]),"addsnp"] <- 2
    gt[gt[,2] == names(b[2]),"addsnp"] <- 1
    gt[gt[,2] == names(b[3]),"addsnp"] <- 0
    major = names(b[3])
  }
  snpid <-tped[2]
  pheno_gt <- left_join(x = pheno, y = as.data.frame(gt), by="SUBJECT")
  pheno_gt <- pheno_gt %>% select(SUBJECT, URICACID, AGE, BMI, PCA1, SEX, PCA2, SNP, addsnp)
  pheno_gt <- na.omit(pheno_gt)
  pheno_gt$rownames <- rownames(pheno_gt)
  test <- lm(BMI~AGE+as.factor(SEX)+URICACID+PCA1+PCA2+addsnp, data=pheno_gt)
  test2 <- lm(residuals.lm(test)^2 ~ pheno_gt$addsnp)



  #test <- lm(URICACID~AGE+as.factor(SEX)+BMI+PCA1+PCA2, data=pheno_gt)
  #test3 <- lm( rstandard(test)^2 ~ as.factor(pheno_gt$SNP) )


  #pheno_gt %>% group_by(SNP) %>% summarise(mean(URICACID), var(URICACID), IQR(URICACID))
  #summary(test)
  #residuals <- data.frame(test$residuals)
  #rownames(residuals) <- pheno_gt$SUBJECT #[(!is.na(pheno_gt$URICACID & pheno_gt$AGE & pheno_gt$BMI & pheno_gt$PCA1 & pheno_gt$PCA2 & pheno_gt$SNP))]
  #transformed <- qnorm((rank(residuals$test.residuals,na.last="keep")-0.5)/sum(!is.na(residuals$test.residuals)))
  #residuals$transformed <- (transformed)^2
  #residuals$SUBJECT <- row.names(residuals)
  print(paste0(i," ",snpid))
  add_table(rbind(c(snpid,major,snpid,snpid),cbind((cbind(beta= coef(test), confint(test))),Pvalue= summary(test)$coefficients[,4]),c("test2","test2","test2","test2"), cbind((cbind(beta= coef(test2), confint(test2))),Pvalue= summary(test2)$coefficients[,4])))

  #pheno_gt <- left_join(x = pheno_gt, y = residuals, by="SUBJECT")
  rm(a,tped,gt,snpid,pheno_gt, test, test2,b,major)
  #test2 <- lm(pheno_gt$transformed ~ pheno_gt$SNP)

}


#
# grep -v 'AGE*\|(Intercept*\|PC*\|BMI\|as.factor*\|test*' all_snps.csv >all_snps_filtered.csv
# cat all_snps_filtered.csv | cut -f5 > snp_p.txt
#

lines <- read.table("snp_p.txt", stringsAsFactors = FALSE)

df <- data.frame()
df <- cbind(lines[seq(from=1,to = length(lines[,1]),by = 3), 1 ], lines[seq(2,length(lines[,1]),3), 1 ] , lines[seq(3,length(lines[,1]),3), 1 ])
colnames(df) <- c("snp", "snp_p", "res_snp_p")

map <- read.table("~/MurrayXsan/vGWAS/archive/variance datasets/ARIC/ARICforvariance.bim", stringsAsFactors = FALSE)
colnames(map) <- c("chr", "snp", "gen","pos","a1", "a2")

df2 <- merge(map, df, by="snp")
df2$snp_p <- as.numeric(as.character(df2$snp_p))
df2$res_snp_p <- as.numeric(as.character(df2$res_snp_p))
df2$chr <- as.numeric(df2$chr)
df2$pos <- as.numeric(df2$pos)

library(qqman)
manhattan(x=df2, chr = 'chr', bp ='pos', p = 'snp_p', snp = 'snp', main= "lm(URICACID ~ SNP)")
manhattan(x=df2, chr = 'chr', bp ='pos', p = 'res_snp_p', snp = 'snp', main= "lm(residuals^2 ~ SNP")
