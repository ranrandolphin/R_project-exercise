
F1 <- as.numeric(as.character(public$V12)) + 0.5 * as.numeric(as.character(public$V13))
F2 <- as.numeric(as.character(public$V32)) + 0.5 * as.numeric(as.character(public$V33))
Y <- F2 - F1
indicator = as.factor(ifelse(Y > 0, 1, 0))
w1 <- as.numeric(as.character(public$V15))
w2 <- as.numeric(as.character(public$V35))
D <- w1 - w2

H <- as.factor(public$V2) # CHAIN
C <- as.factor(public$V3) #CO_OWNED
Z <- as.factor(public$V4) #STATE
EMPFT <- as.numeric(as.character(public$V12))
EMPPT <- as.numeric(as.character(public$V13))
NMGRS <- as.numeric(as.character(public$V14))
w1 <- as.numeric(as.character(public$V15))
INCTIME <- as.numeric(as.character(public$V16))
FIRSTINC <- as.numeric(as.character(public$V17))
BONUS <- as.factor(public$V18)
PCTAFF <- as.numeric(as.character(public$V18))
MEALS <- as.factor(public$V20)
OPEN <- as.numeric(public$V21)
HRSOPEN <- as.numeric(public$V22)
PSODA <- as.numeric(as.character(public$V23))
PFRY <- as.numeric(as.character(public$V24))
PENTREE <- as.numeric(as.character(public$V25))
NRGES <- as.numeric(as.character(public$V26))
NREGS11 <- as.numeric(as.character(public$V27))


wage <- na.omit(data.frame(Y,D,Z,C,H,EMPFT,EMPPT,NMGRS,w1,INCTIME,
                           FIRSTINC,BONUS,MEALS,OPEN,HRSOPEN,PSODA,
                           PFRY,PENTREE,NRGES,NREGS11))
row.names(wage) <- NULL
```{r eval =FALSE}
wage <- na.omit(data.frame(Y,D,Z,C,H,EMPFT,EMPPT,NMGRS,w1,INCTIME,
                           FIRSTINC,BONUS,MEALS,OPEN,HRSOPEN,PSODA,
                           PFRY,PENTREE,NRGES,NREGS11))
                           
###############################
####1. Random Forest Method####
library(randomForest)
library(ROCR)
wage.x <- wage[, 2:20]
wage.y <- wage[, 1]

for (i in c(50 ,200, 500, 1000)){
  for (j in c(3,6,9,12)){
    set.seed(2048)
    print(randomForest(x = wage.x, y = wage.y, ntree = i, mtry=j, nodesize = 15))
  }
}
```
```{r include = FALSE}
library(randomForest)
library(ROCR)
wage.x <- wage[, 2:20]
wage.y <- wage[, 1]

for (i in c(50 ,200, 500, 1000)){
  for (j in c(3,6,9,12)){
    set.seed(2048)
    print(randomForest(x = wage.x, y = wage.y, ntree = i, mtry=j, nodesize = 15))
  }
}

###############################
####2. Any other settings that characterize implementation of the method####

set.seed(13456)
bestmtry <-tuneRF(wage.x, wage.y, ntreTry = 500, 
                  stepFactor = 2, improve = 0.05, trace = TRUE, plot = TRUE)
plot(bestmtry, type="b",xlab = "optimal mtry parameter")

set.seed(13456)
bestmtry <-tuneRF(wage.x, wage.y, ntreTry = 500, stepFactor = 2, improve = 0.05, trace = TRUE, plot = TRUE)
plot(bestmtry, type="b",xlab = "optimal mtry parameter")

### From the mtry vs. OOB error plot, we find that there is tunning point at mtry = 12, becase we get the smallest oob error here. Thus, mtry = 12 is the best number with set.seed.

###############################
####3. A representation of how well the model predicts the outcome####

wage.index <- sample(1:nrow(wage), round(0.7*nrow(wage)))
train.x <- wage.x[wage.index, ]
train.y <- wage.y[wage.index]
test.x <- wage.x[-wage.index,]
test.y <- wage.y[-wage.index]

rf1 <- randomForest(x = train.x, y = train.y, xtest = test.x, 
                    ytest = test.y, ntree=500, mtry=3, importance=TRUE)
rf2 <- randomForest(x = train.x, y = train.y, xtest = test.x, 
                    ytest = test.y, ntree=500, mtry=6, importance=TRUE)
rf3 <- randomForest(x = train.x, y = train.y, xtest = test.x, 
                    ytest = test.y, ntree=500, mtry=9, importance=TRUE)
rf4 <- randomForest(x = train.x, y = train.y, xtest = test.x, 
                    ytest = test.y, ntree=500, mtry=12, importance=TRUE)

###############################
####4. A list of which elements of X are ‘important’ for the prediction####
par(mfrow=c(2,2))
plot(importance(rf1, type=1), main="Importance, mtry=3", 
     ylab="Mean Decrease in Accuracy", xaxt="n")
axis(1, at=1:14, las=2)
plot(importance(rf2, type=1), main="Importance, mtry=6", 
     ylab="Mean Decrease in Accuracy", xaxt="n")
axis(1, at=1:14, las=2)
plot(importance(rf3, type=1), main="Importance, mtry=9", 
     ylab="Mean Decrease in Accuracy", xaxt="n")
axis(1, at=1:14, las=2)
plot(importance(rf4, type=1), main="Importance, mtry=12", 
     ylab="Mean Decrease in Accuracy", xaxt="n")
axis(1, at=1:14, las=2)
