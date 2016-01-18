## Normal urate association
# URATE = SNP + Predictors + residuals

## Unadjusted urate
# URATE = Predictors + residuals
# then(Inverse normal transformed residuals)^2 correlated to SNP genotype

## Adjusted urate
# Adjusted urate = Predictors + residuals.
# then (Inverse normal transformed residuals)^2 correlated to SNP genotype
# Where Adjusted urate = individual urate measurements with the mean urate subtracted for each genotype group

#### Change DATASET = xxxx and xxxx_result = result to desired dataset before starting

library(dplyr)
print(paste0("start", date()))

load("Data.RData")
data = list(aric, cardia, chs, fhs)
names(data) <- c("aric", "cardia", "chs", "fhs")
results_list = list()
d = "aric"
DATASET = data[[d]]
print(d)

print(paste0("loaded", date()))

result = data.frame(matrix(nrow=(length(DATASET[1,])-9),ncol=10))
names(result) = c("SNP","Mean_Beta","Mean_SE","Mean_P","Unadj_Beta","Unadj_SE","Unadj_P","Adj_Beta","Adj_SE","Adj_P")
result$SNP = names(DATASET)[10:length(DATASET[1,])]

MALES<-subset(DATASET, SEX==1)
FEMALES<-subset(DATASET, SEX==2)

MALESunadj<-lm(URICACID~AGE+BMI+PCA1+PCA2, data=MALES)
MALES$unadjres=MALESunadj$residuals
MALES$unadjZ <- qnorm((rank(MALES$unadjres,na.last="keep")-0.5)/sum(!is.na(MALES$unadjres)))
MALES$unadjZ2 <- (MALES$unadjZ)^2

FEMALESunadj<-lm(URICACID~AGE+BMI+PCA1+PCA2, data=FEMALES)
FEMALES$unadjres=FEMALESunadj$residuals
FEMALES$unadjZ <- qnorm((rank(FEMALES$unadjres,na.last="keep")-0.5)/sum(!is.na(FEMALES$unadjres)))
FEMALES$unadjZ2 <- (FEMALES$unadjZ)^2
result2 <- result
a<-rbind(MALES, FEMALES)
print(date())
for(i in 10:length(DATASET[1,])){
#for(i in 10:12) {
  tmp <- data.frame()
for(geno in levels(factor(a[,i]))){
  tmp1 <- a[a[,i] == geno & a$SEX == 1 & !is.na(a$URICACID) & 
              !is.na(a[,i] == geno), 
            c("SUBJECT","SEX","AGE","BMI","PCA1","PCA2","URICACID", "unadjZ", 
              "unadjZ2", "unadjres", "MENOPAUSE", names(a)[i])]
  names(tmp1) <- 
    c("SUBJECT","SEX","AGE","BMI","PCA1","PCA2","URICACID", "unadjZ", 
      "unadjZ2", "unadjres", "MENOPAUSE", "SNP")
  if(length(tmp1[,1] )> 0){
    tmp1$meanuratebysnp <- mean(tmp1$URICACID, na.rm=TRUE)
    tmp1$adjuratebysnp=(tmp1$URICACID-tmp1$meanuratebysnp)
  }
  tmp2 <- a[a[,i] == geno & a$SEX == 2 & !is.na(a$URICACID) & 
              !is.na(a[,i] == geno), 
            c("SUBJECT","SEX","AGE","BMI","PCA1","PCA2","URICACID", "unadjZ", 
              "unadjZ2", "unadjres", "MENOPAUSE",names(a)[i])]
  names(tmp2) <- 
    c("SUBJECT","SEX","AGE","BMI","PCA1","PCA2","URICACID", "unadjZ", 
      "unadjZ2", "unadjres", "MENOPAUSE", "SNP")
  if(length(tmp2[,1]) > 0){
    tmp2$meanuratebysnp <- mean(tmp2$URICACID, na.rm=TRUE)
    tmp2$adjuratebysnp=(tmp2$URICACID-tmp2$meanuratebysnp)
  }
  tmp <- rbind(tmp, rbind(tmp1,tmp2))
}
  m <- lm(adjuratebysnp~AGE+BMI+PCA1+PCA2, data=tmp[tmp$SEX ==1,])
  f <- lm(adjuratebysnp~AGE+BMI+PCA1+PCA2, data=tmp[tmp$SEX ==2,])
  tmp$adjres = NA
  tmp$adjZ = NA
  tmp$adjZ2 = NA
  tmp[tmp$SEX ==1,]$adjres <- m$residuals
  tmp[tmp$SEX ==1,]$adjZ <- qnorm((rank(tmp[tmp$SEX ==1,]$adjres,na.last="keep")-0.5)/sum(!is.na(tmp[tmp$SEX ==1,]$adjres)))
  tmp[tmp$SEX ==1,]$adjZ2 <- (tmp[tmp$SEX ==1,]$adjZ)^2

  tmp[tmp$SEX ==2,]$adjres <- f$residuals
  tmp[tmp$SEX ==2,]$adjZ <- qnorm((rank(tmp[tmp$SEX ==2,]$adjres,na.last="keep")-0.5)/sum(!is.na(tmp[tmp$SEX ==2,]$adjres)))
  tmp[tmp$SEX ==2,]$adjZ2 <- (tmp[tmp$SEX ==2,]$adjZ)^2


  # Urate to RS10000104 association
  tryCatch(testCOHORT<-lm(URICACID~as.numeric(SNP)+SEX+AGE+BMI+PCA1+PCA2, data=tmp, subset=SEX==1),error=function(e) NULL)
  tryCatch(result[i-9,2:4] <- coef(summary(testCOHORT))[2,c(1,2,4)],error=function(e) NULL)

  # Unadjusted urate variance association
  tryCatch(COHORTunadj<-lm(unadjZ2~as.numeric(SNP), data=tmp, subset=SEX==1),error=function(e) NULL)
  tryCatch(result[i-9,5:7] <- coef(summary(COHORTunadj))[2,c(1,2,4)],error=function(e) NULL)

  # Adjusted for mean urate variance association
  tryCatch(COHORTadj<-lm(adjZ2~as.numeric(SNP), data=tmp, subset=SEX==1),error=function(e) NULL)
  tryCatch(result[i-9,8:10] <- coef(summary(COHORTadj))[2,c(1,2,4)],error=function(e) NULL)


  rm(tmp,testCOHORT,COHORTunadj,COHORTadj,COHORTunadj2,COHORTadj2,COHORT)
  if(i %% 100 == 0){
    print(paste(i, d, sep=" "))
  }
  if(i %% 1000 == 0){
    write.table(result,file=paste0(d,'2.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
  }
}

print(date())
write.table(result, file=paste0(d,'2.txt'), row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

rm(DATASET,i,result)