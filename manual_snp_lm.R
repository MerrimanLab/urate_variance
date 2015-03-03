library(dplyr)
library(data.table)
library(foreach)

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
  colnames(gt) <- c("SUBJECT", "SNP")
  snpid <-tped[2]
  pheno_gt <- left_join(x = pheno, y = as.data.frame(gt), by="SUBJECT")
  pheno_gt <- pheno_gt %>% select(SUBJECT, URICACID, AGE, BMI, PCA1, SEX, PCA2, SNP)
  pheno_gt <- na.omit(pheno_gt)
  pheno_gt$rownames <- rownames(pheno_gt)
  test <- lm(URICACID~AGE+as.factor(SEX)+BMI+PCA1+PCA2+as.factor(SNP), data=pheno_gt)
  test2 <- lm(residuals.lm(test)^2 ~ as.factor(pheno_gt$SNP))



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
  add_table(rbind(c(snpid,snpid,snpid,snpid),cbind((cbind(beta= coef(test), confint(test))),Pvalue= summary(test)$coefficients[,4]),c("test2","test2","test2","test2"), cbind((cbind(beta= coef(test2), confint(test2))),Pvalue= summary(test2)$coefficients[,4])))

  #pheno_gt <- left_join(x = pheno_gt, y = residuals, by="SUBJECT")
  rm(a,tped,gt,snpid,pheno_gt, test, test2)
  #test2 <- lm(pheno_gt$transformed ~ pheno_gt$SNP)

}

