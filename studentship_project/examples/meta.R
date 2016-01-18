library(meta)

## name according to cohort
names(aric) = c("SNP","aric_Mean_Beta","aric_Mean_SE","aric_Mean_P","aric_Unadj_Beta","aric_Unadj_SE","aric_Unadj_P","aric_Adj_Beta","aric_Adj_SE","aric_Adj_P")
names(cardia) = c("SNP","cardia_Mean_Beta","cardia_Mean_SE","cardia_Mean_P","cardia_Unadj_Beta","cardia_Unadj_SE","cardia_Unadj_P","cardia_Adj_Beta","cardia_Adj_SE","cardia_Adj_P")
names(chs) = c("SNP","chs_Mean_Beta","chs_Mean_SE","chs_Mean_P","chs_Unadj_Beta","chs_Unadj_SE","chs_Unadj_P","chs_Adj_Beta","chs_Adj_SE","chs_Adj_P")
names(fhs) = c("SNP","fhs_Mean_Beta","fhs_Mean_SE","fhs_Mean_P","fhs_Unadj_Beta","fhs_Unadj_SE","fhs_Unadj_P","fhs_Adj_Beta","fhs_Adj_SE","fhs_Adj_P")

## merge into one data.frame
results = merge(aric,cardia,by="SNP")
results = merge(results,chs,by="SNP")
results = merge(results,fhs,by="SNP")

## mean
result_Mean <- data.frame(matrix(nrow=length(results[,1]),ncol=6)) ##make a dataframe with number of rows = number of SNPs in object results, and 6 columns
names(result_Mean) <- c("SNP","TE","seTE","P","Het.P","K") ## name the 6 columns
result_Mean$SNP <- results$SNP #add SNPs from results to result_mean
## fill in TE, seTE, P, Het.P, K from combining the 4 cohorts ()
for(i in 1:length(results[,1])){ 
  meta <- metagen(TE=c(results$aric_Mean_Beta[i],results$cardia_Mean_Beta[i],results$chs_Mean_Beta[i],results$fhs_Mean_Beta[i]),seTE=c(results$aric_Mean_SE[i],results$cardia_Mean_SE[i],results$chs_Mean_SE[i],results$fhs_Mean_SE[i]),studlab=c("ARIC","CARDIA","CHS","FHS"))
  result_Mean[i,2] <- ifelse((1-pchisq(meta$Q,(meta$k-1)))<0.05,meta$TE.random,meta$TE.fixed)
  result_Mean[i,3] <- ifelse((1-pchisq(meta$Q,(meta$k-1)))<0.05,meta$seTE.random,meta$seTE.fixed)
  result_Mean[i,4] <- ifelse((1-pchisq(meta$Q,(meta$k-1)))<0.05,(pnorm(abs(meta$zval.random),lower.tail=FALSE)*2),(pnorm(abs(meta$zval.fixed),lower.tail=FALSE)*2))
  result_Mean[i,5] <- 1-pchisq(meta$Q,(meta$k-1))
  result_Mean[i,6] <- meta$k
  if(i %% 100==0){
    print(i)
  }
  rm(meta)
}
rm(i)

result_Unadj <- data.frame(matrix(nrow=length(results[,1]),ncol=6))
names(result_Unadj) <- c("SNP","TE","seTE","P","Het.P","K")
result_Unadj$SNP <- results$SNP
for(i in 1:length(results[,1])){
  meta <- metagen(TE=c(results$aric_Unadj_Beta[i],results$cardia_Unadj_Beta[i],results$chs_Unadj_Beta[i],results$fhs_Unadj_Beta[i]),seTE=c(results$aric_Unadj_SE[i],results$cardia_Unadj_SE[i],results$chs_Unadj_SE[i],results$fhs_Unadj_SE[i]),studlab=c("ARIC","CARDIA","CHS","FHS"))
  result_Unadj[i,2] <- ifelse((1-pchisq(meta$Q,(meta$k-1)))<0.05,meta$TE.random,meta$TE.fixed)
  result_Unadj[i,3] <- ifelse((1-pchisq(meta$Q,(meta$k-1)))<0.05,meta$seTE.random,meta$seTE.fixed)
  result_Unadj[i,4] <- ifelse((1-pchisq(meta$Q,(meta$k-1)))<0.05,(pnorm(abs(meta$zval.random),lower.tail=FALSE)*2),(pnorm(abs(meta$zval.fixed),lower.tail=FALSE)*2))
  result_Unadj[i,5] <- 1-pchisq(meta$Q,(meta$k-1))
  result_Unadj[i,6] <- meta$k
  if(i %% 100==0){
    print(i)
  }
  rm(meta)
}
rm(i)

result_Adj <- data.frame(matrix(nrow=length(results[,1]),ncol=6))
names(result_Adj) <- c("SNP","TE","seTE","P","Het.P","K")
result_Adj$SNP <- results$SNP
for(i in 1:length(results[,1])){
  meta <- metagen(TE=c(results$aric_Adj_Beta[i],results$cardia_Adj_Beta[i],results$chs_Adj_Beta[i],results$fhs_Adj_Beta[i]),seTE=c(results$aric_Adj_SE[i],results$cardia_Adj_SE[i],results$chs_Adj_SE[i],results$fhs_Adj_SE[i]),studlab=c("ARIC","CARDIA","CHS","FHS"))
  result_Adj[i,2] <- ifelse((1-pchisq(meta$Q,(meta$k-1)))<0.05,meta$TE.random,meta$TE.fixed)
  result_Adj[i,3] <- ifelse((1-pchisq(meta$Q,(meta$k-1)))<0.05,meta$seTE.random,meta$seTE.fixed)
  result_Adj[i,4] <- ifelse((1-pchisq(meta$Q,(meta$k-1)))<0.05,(pnorm(abs(meta$zval.random),lower.tail=FALSE)*2),(pnorm(abs(meta$zval.fixed),lower.tail=FALSE)*2))
  result_Adj[i,5] <- 1-pchisq(meta$Q,(meta$k-1))
  result_Adj[i,6] <- meta$k
  if(i %% 100==0){
    print(i)
  }
  rm(meta)
}
rm(i)

## replace p=0 with p=1e-300 for MeanEffect
result_Mean$P[result_Mean$P==0] <- 1e-300
## remove SNPs with K less than 4
result_Mean <- result_Mean[result_Adj$K==4,]
result_Adj <- result_Adj[result_Adj$K==4,]
result_Unadj <- result_Unadj[result_Unadj$K==4,]

name1 <- sapply(strsplit(x=result_Mean$SNP, split="_", fixed=TRUE), "[[", 1)
result_Mean$SNP <- name1

name2 <- sapply(strsplit(x=result_Unadj$SNP, split="_", fixed=TRUE), "[[", 1)
result_Unadj$SNP <- name2

name3 <- sapply(strsplit(x=result_Adj$SNP, split="_", fixed=TRUE), "[[", 1)
result_Adj$SNP <- name3

setwd("~/GCKR/results")

write.table(result_Mean,file="meta_MeanEffect_GCKR.txt",quote=F,row.names=F,sep="\t")
write.table(result_Unadj,file="meta_vUnadjEffect_GCKR.txt",quote=F,row.names=F,sep="\t")
write.table(result_Adj,file="meta_vAdjEffect_GCKR.txt",quote=F,row.names=F,sep="\t")


# extra columns/details
#result_Adj3 = merge(SNP_Positions,result_Adj,by.x="MARKER",by.y="SNP",all.y=TRUE)
#result_Unadj3 = merge(SNP_Positions,result_Unadj,by.x="MARKER",by.y="SNP",all.y=TRUE)
#result_Mean3 = merge(SNP_Positions,result_Mean,by.x="MARKER",by.y="SNP",all.y=TRUE)
#result_Mean2 = merge(result_Mean3,Keep_SNPs,by.x="MARKER",by.y="SNP",all=TRUE)
#result_Unadj2 = merge(result_Unadj3,Keep_SNPs,by.x="MARKER",by.y="SNP",all=TRUE)
#result_Adj2 = merge(result_Adj3,Keep_SNPs,by.x="MARKER",by.y="SNP",all=TRUE)
#result_Adj_Keep = result_Adj2[result_Adj2$KeepOrNot=="Worth_Keeping",]
#result_Unadj_Keep = result_Unadj2[result_Adj2$KeepOrNot=="Worth_Keeping",]
#result_Mean_Keep = result_Mean2[result_Mean2$KeepOrNot=="Worth_Keeping",]
#result_Unadj_Keep = result_Unadj2[result_Unadj2$KeepOrNot=="Worth_Keeping",]

#write.table(result_Mean_Keep,file="meta_MeanEffect_KeepSNPs.txt",quote=F,row.names=F,sep="\t")
#write.table(result_Unadj_Keep,file="meta_vUnadjEffect_KeepSNPs.txt",quote=F,row.names=F,sep="\t")
#write.table(result_Adj_Keep,file="meta_vAdjEffect_KeepSNPs.txt",quote=F,row.names=F,sep="\t")
#write.table(result_Mean3,file="meta_MeanEffect2.txt",quote=F,row.names=F,sep="\t")
#write.table(result_Unadj3,file="meta_vUnadjEffect2.txt",quote=F,row.names=F,sep="\t")
#write.table(result_Adj3,file="meta_vAdjEffect2.txt",quote=F,row.names=F,sep="\t")