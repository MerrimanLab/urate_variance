ARIC <- read.delim("ARIC_snpscov.txt", header=TRUE)
CARDIA <- read.delim("CARDIA_snpscov.txt", header=TRUE)
dataset <- rbind(ARIC, CARDIA)

dataset$RS7657853 <- as.character(dataset$RS7657853)
dataset$RS7657853_code[dataset$RS7657853=="AA" & !is.na(dataset$RS7657853)] <- 0
dataset$RS7657853_code[dataset$RS7657853=="GA" | dataset$RS7657853=="GG" & !is.na(dataset$RS7657853)] <- 1
table(dataset$RS7657853_code)

POST <- dataset[dataset$MENOPAUSE=="POST" & !is.na(dataset$MENOPAUSE),]
                
summary(lm(URICACID~PLATELET_COUNT+AGE+BMI+PCA1+PCA2, data=POST))
summary(lm(URICACID~PLATELET_COUNT+AGE+BMI+PCA1+PCA2, data=POST, subset=RS7657853_code==1))
