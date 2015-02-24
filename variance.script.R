#subset data into males and females for each cohort
aric<-subset(allaricfhsforvarinace, COHORT=="ARIC")
fhs<-subset(allaricfhsforvarinace, COHORT=="FHS")
aricfemales<-subset(aric, SEX==2)
aricmales<-subset(aric, SEX==1)
fhsmales<-subset(fhs, SEX==1)
fhsfemales<-subset(fhs, SEX==2)

#initial regression and transformation of residuals
test <- lm(URICACID~AGE+BMI+PCA1+PCA2, data=FHSmale)
summary(test)
residuals <- data.frame(test$residuals)
rownames(residuals)=FHSmale$SUBJECT[(!is.na(FHSmale$URICACID&FHSmale$AGE&FHSmale$BMI&FHSmale$PCA1&FHSmale$PCA2))]
transformed <- qnorm((rank(residuals$test.residuals,na.last="keep")-0.5)/sum(!is.na(residuals$test.residuals)))
residuals$transformed <- (transformed)^2
residuals$ID <- row.names(residuals)
write.table(residuals,file="residuals_FHSmale.txt",quote=FALSE,row.names=FALSE)

#create a pheotype file with the z2 scores from both males and females compiled for use in PLINK
#regress z2 to SNPs in PLINK
#import plink results into R
#merge datasets to remove SNPs not common to both datasets
#compare reference alleles and convert Beta values for those that don't match between the datasets
results <- merge (plink.ARIC.results,plink.FHS.results, by="MARKER")
results$A1.x <- as.character(results$A1.x)
results$A1.y <- as.character(results$A1.y)
results$fix =results$BETA.y
results[results$A1.x != results$A1.y,]$fix=results[results$A1.x != results$A1.y,]$fix * -1
write.table(results[,c("MARKER","A1.x","NMISS.x","BETA.x","SE.x","L95.x","U95.x","STAT.x","P.x")], file="ARICresultsassoc.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(results[,c("MARKER","A1.x","NMISS.y","fix","SE.y","L95.y","U95.y","STAT.y","P.y")], file="FHSresultsassoc.txt",quote=FALSE,row.names=FALSE,sep="\t")

#feed this data into METAL. Run METAL with STDERR scheme and estimate of genomic control.
#renaming columns for running manhattan plots
colnames(metalresults)
[1] "MARKER"  "CHROM"   "POS"     "ALLELEA" "ALLELEB" "WEIGHT"  "ZSCORE"  "PVALUE" 
colnames(metalresults)[1]="SNP"
colnames(metalresults)[2]="CHR"
colnames(metalresults)[3]="BP"
colnames(metalresults)[8]="P"

#run manhattan plot
library(qqman)
metalresults$CHR <- as.character(metalresults$CHR)
metalresults$CHR <- as.numeric(metalresults$CHR)
metalresults2<-subset(metalresults, is.na(metalresults$CHR)==FALSE)
manhattan(metalresults2, ymax=10, suggestiveline = F, genomewideline = F)


