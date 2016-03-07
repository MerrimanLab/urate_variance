load("Non-additive interaction for diuretics.RData")
dataset <- rbind(ARIC_snps, CARDIA_snps, CHS, CHSDIU_snps, FHS_snps)
dataset$DIU <- NA
dataset$DIU[dataset$COHORT=="CHSDIU" | dataset$COHORT=="ARICDIU"] <- 1
dataset$DIU[dataset$COHORT=="CHS" | dataset$COHORT=="ARIC" | dataset$COHORT=="FHS" | dataset$COHORT=="CARDIA"] <- 0
##RS3775935 AA > CA > CC (AA urate raising)
##6449173 TT > TG > GG (TT urate raising)
##7657853 AA > AG > GG (AA urate raising)

## Non-additive interaction testing - interaction term is the coefficient extracted. 
dataset$RS3775935 <- as.character(dataset$RS3775935)
dataset$RS3775935_code[dataset$RS3775935=="CC" & !is.na(dataset$RS3775935)] <- 0 ## Change to RS3775935=="AA" for alternative allelic grouping
dataset$RS3775935_code[dataset$RS3775935=="CA" | dataset$RS3775935=="AA" & !is.na(dataset$RS3775935)] <- 1 ## Change to RS3775935=="CC" for alternative allelic grouping (AA vs. CA/AA)
table(dataset$RS3775935_code)
summary(lm(URICACID~(RS3775935_code)*DIU+SEX+AGE+BMI+PCA1+PCA2, data=dataset))

dataset$RS6449173 <- as.character(dataset$RS6449173)
dataset$RS6449173_code[dataset$RS6449173=="TT" & !is.na(dataset$RS6449173)] <- 0
dataset$RS6449173_code[dataset$RS6449173=="TG" | dataset$RS6449173=="GG" & !is.na(dataset$RS6449173)] <- 1
table(dataset$RS6449173_code)
summary(lm(URICACID~(RS6449173_code)*DIU+SEX+AGE+BMI+PCA1+PCA2, data=dataset))

dataset$RS7657853 <- as.character(dataset$RS7657853)
dataset$RS7657853_code[dataset$RS7657853=="AA" & !is.na(dataset$RS7657853)] <- 0 ## Change to RS7657853=="GG" for alternative
dataset$RS7657853_code[dataset$RS7657853=="GA" | dataset$RS7657853=="GG" & !is.na(dataset$RS7657853)] <- 1 ## Change to RS7657853=="AA" for alternative
table(dataset$RS7657853_code)
summary(lm(URICACID~(RS7657853_code)*DIU+SEX+AGE+BMI+PCA1+PCA2, data=dataset))

## GxE groups, with "grouP" as the variable regressed against. 
dataset$RS6449173_group <- NA
dataset$RS6449173_group[dataset$RS6449173_code==0 & dataset$DIU==0 & !is.na(dataset$RS6449173_code) & !is.na(dataset$DIU)] <- 0
dataset$RS6449173_group[dataset$RS6449173_code==1 & dataset$DIU==0 & !is.na(dataset$RS6449173_code) & !is.na(dataset$DIU)] <- 1
dataset$RS6449173_group[dataset$RS6449173_code==0 & dataset$DIU==1 & !is.na(dataset$RS6449173_code) & !is.na(dataset$DIU)] <- 2
dataset$RS6449173_group[dataset$RS6449173_code==1 & dataset$DIU==1 & !is.na(dataset$RS6449173_code) & !is.na(dataset$DIU)] <- 3
dataset$RS6449173_group <- as.factor(dataset$RS6449173_group)
table(dataset$RS6449173_group)

summary(lm(URICACID~RS6449173_group+SEX+AGE+BMI+PCA1+PCA2, data=dataset))

dataset$RS7657853_group <- NA
dataset$RS7657853_group[dataset$RS7657853_code==0 & dataset$DIU==0 & !is.na(dataset$RS7657853_code) & !is.na(dataset$DIU)] <- 0
dataset$RS7657853_group[dataset$RS7657853_code==1 & dataset$DIU==0 & !is.na(dataset$RS7657853_code) & !is.na(dataset$DIU)] <- 1
dataset$RS7657853_group[dataset$RS7657853_code==0 & dataset$DIU==1 & !is.na(dataset$RS7657853_code) & !is.na(dataset$DIU)] <- 2
dataset$RS7657853_group[dataset$RS7657853_code==1 & dataset$DIU==1 & !is.na(dataset$RS7657853_code) & !is.na(dataset$DIU)] <- 3
dataset$RS7657853_group <- as.factor(dataset$RS7657853_group)
table(dataset$RS7657853_group)

summary(lm(URICACID~RS7657853_group+SEX+AGE+BMI+PCA1+PCA2, data=dataset))

dataset$RS3775935_group <- NA
dataset$RS3775935_group[dataset$RS3775935_code==0 & dataset$DIU==0 & !is.na(dataset$RS3775935_code) & !is.na(dataset$DIU)] <- 0
dataset$RS3775935_group[dataset$RS3775935_code==1 & dataset$DIU==0 & !is.na(dataset$RS3775935_code) & !is.na(dataset$DIU)] <- 1
dataset$RS3775935_group[dataset$RS3775935_code==0 & dataset$DIU==1 & !is.na(dataset$RS3775935_code) & !is.na(dataset$DIU)] <- 2
dataset$RS3775935_group[dataset$RS3775935_code==1 & dataset$DIU==1 & !is.na(dataset$RS3775935_code) & !is.na(dataset$DIU)] <- 3
dataset$RS3775935_group <- as.factor(dataset$RS3775935_group)
table(dataset$RS3775935_group)

summary(lm(URICACID~RS3775935_group+SEX+AGE+BMI+PCA1+PCA2, data=dataset))
