setwd("~/Summer Research/Jaron/SLC16A9")
## ARIC
ARIC <- read.delim("ARIC.303382")
ARIC1_df <- read.delim("~/Summer Research/Jaron/vGWAS/ARIC1_df.txt")

aric <- merge(ARIC1_df[,c(1:9)], ARIC, by="SUBJECT")
table(aric$SUBJECT==ARIC1_df$SUBJECT) ## checking that there is no issue with duplicated participants

write.table(aric, "vGWAS_ARIC_SLC16A9.txt", row.names = F, col.names = T, sep = "\t", quote = F)

## CARDIA
CARDIA <- read.delim("CARDIA.303383")
CARDIA1_df <- read.delim("~/Summer Research/Jaron/vGWAS/CARDIA1_df.txt")

cardia <- merge(CARDIA1_df[,c(1:9)], CARDIA, by="SUBJECT")
table(cardia[,c(1:9)]==CARDIA1_df[,c(1:9)])

write.table(cardia, "vGWAS_CARDIA_SLC16A9.txt", row.names = F, col.names = T, sep = "\t", quote = F)

## CHS
CHS <- read.delim("CHS.303384")
CHS1_df <- read.delim("~/Summer Research/Jaron/vGWAS/CHS1_df.txt")

chs <- merge(CHS1_df[,c(1:9)], CHS, by="SUBJECT")
table(chs[,c(1:9)]==CHS1_df[,c(1:9)])

write.table(chs, "vGWAS_CHS_SLC16A9.txt", row.names = F, col.names = T, sep = "\t", quote = F)


## FHS
FHS <- read.delim("FHS.303394")
FHS1_df <- read.delim("~/Summer Research/Jaron/vGWAS/FHS1_df.txt")

fhs <- merge(FHS1_df[,c(1:9)], FHS, by="SUBJECT", all.x=TRUE)
table(fhs[,c(1:9)]==FHS1_df[,c(1:9)])

write.table(fhs, "vGWAS_FHS_SLC16A9.txt", row.names = F, col.names = T, sep = "\t", quote = F)

rm(ARIC)
rm(ARIC1_df)
rm(CARDIA)
rm(CARDIA1_df)
rm(CHS)
rm(CHS1_df)
rm(FHS)
rm(FHS1_df)

save.image("~/Summer Research/Jaron/SLC16A9/SLC16A9_Data.RData")
