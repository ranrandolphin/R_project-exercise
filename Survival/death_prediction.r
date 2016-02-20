library(survival)
data(colon)
attach(colon)
#summary(colon)
# data for Ran's group Patient ids assigned 101-910

data <- colon[201:1820,] 
rownames(data) <- NULL
data[is.na(data$differ),"differ"] <- 2# replace 26 NA in differ with median value 2

# replace Age >= 60 to 1, and Age < 60 to 0
data[data$age >= 60, "age"] <- 1
data[data$age < 60 & data$age > 1, "age"] <- 0
new.data <- data
#cor(new.data[-c(1,2,3,10,15,16)]) 

# correlation between nodes(without NA) and node4 is high, and we decided to remove nodes because it has the NA missing data

#rownames(new.data) <- NULL
#new.data[is.na(new.data$differ),"differ"] <- 2
new.data$nodes <- NULL #node4 (>4), has relatively greater correlation with other variables compared with nodes


death <- subset(new.data, new.data$etype == 2)
rownames(death) <- NULL

death.obs <- subset(death, death$rx == "Obs")
rownames(death.obs) <- NULL
death.obs$rx <- droplevels(death.obs$rx)
death.lev <- death[which (death$rx == "Lev"),]
rownames(death.lev) <- NULL
death.lev$rx <- droplevels(death.lev$rx)
death.5fu <- subset(death, death$rx == "Lev+5FU")
rownames(death.5fu) <- NULL
death.5fu$rx <- droplevels(death.5fu$rx)

############################## Prediction: Death (Mortality) #############

########## Training set and Validation Set ##########
set.seed(5)
train.row.death <- sample(seq_len(nrow(death)), size = 540) #set randomly selected 183 rows in training set

train.death <- death[train.row.death,]
rownames(train.death) <- NULL
valid.death <- death[-train.row.death,]
rownames(valid.death) <- NULL

#######Mortality: Lev+5FU set in training model###########

train.obslev5.death <- subset(train.death, train.death$rx == "Obs" | train.death$rx == "Lev+5FU")
valid.obslev5.death <- subset(valid.death, valid.death$rx == "Obs" | valid.death$rx == "Lev+5FU")
rownames(train.obslev5.death) <- NULL
train.obslev5.death$rx <- droplevels(train.obslev5.death$rx)
rownames(valid.obslev5.death) <- NULL
valid.obslev5.death$rx <- droplevels(valid.obslev5.death$rx)
full.pre.lev5.death <- coxph(Surv(train.obslev5.death$time, train.obslev5.death $status)~rx + sex + age + obstruct + perfor + adhere + node4 + factor(differ) + factor(extent.new) + surg, data = train.obslev5.death) 
step(full.pre.lev5.death, direction = c("backward"))

# Thus, our final prediction model for lev5FU is
pre.lev5.death <- coxph(formula = Surv(train.obslev5.death$time, train.obslev5.death$status) ~ node4 + surg, data = train.obslev5.death)
summary(pre.lev5.death) #significant


### change value of rx (lev+5FU = 1, obs = 0) for calculation
X.train.death <- model.matrix(pre.lev5.death)
X.valid.death <- model.matrix(~adhere + node4, data = valid.obslev5.death)
X.valid.death <- X.valid.death[,-1]

#extract the coefficient from previous selected model without NA
coef.lev5.death <- coef(pre.lev5.death)
coef.lev5.death <- as.numeric(coef.lev5.death)
coef.lev5.death <- coef.lev5.death[!is.na(coef.lev5.death)]

risk.train.death <- X.train.death %*% coef.lev5.death
risk.valid.death <- X.valid.death %*% coef.lev5.death

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev5.death <- as.numeric(quantile(risk.train.death, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data


risk.valid.death <- as.factor(as.factor(cut(risk.valid.death, breaks =c(-Inf, quant.lev5.death, Inf), labels = FALSE)))
valid.obslev5.death$level <- risk.valid.death

plot(survfit(Surv(valid.obslev5.death$time, valid.obslev5.death$status) ~ level, data=valid.obslev5.death), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev5.death, data=valid.obslev5.death), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)


################## Obs+Lev5FU - different seed######################
########## Training set and Validation Set ##########
set.seed(818)
train.row.death <- sample(seq_len(nrow(death)), size = 540) #set randomly selected 183 rows in training set

train.death <- death[train.row.death,]
rownames(train.death) <- NULL
valid.death <- death[-train.row.death,]
rownames(valid.death) <- NULL

#######Mortality: Lev+5FU set in training model###########

train.obslev5.death <- subset(train.death, train.death$rx == "Obs" | train.death$rx == "Lev+5FU")
valid.obslev5.death <- subset(valid.death, valid.death$rx == "Obs" | valid.death$rx == "Lev+5FU")
rownames(train.obslev5.death) <- NULL
train.obslev5.death$rx <- droplevels(train.obslev5.death$rx)
rownames(valid.obslev5.death) <- NULL
valid.obslev5.death$rx <- droplevels(valid.obslev5.death$rx)
full.pre.lev5.death <- coxph(Surv(train.obslev5.death $time, train.obslev5.death $status)~rx + sex + age + obstruct + perfor + adhere + node4 + differ + extent + surg, data = train.obslev5.death) 
step(full.pre.lev5.death, direction = c("both"))

# Thus, our final prediction model for lev5FU is
pre.lev5.death <- coxph(formula = Surv(train.obslev5.death$time, train.obslev5.death$status) ~ rx + age + node4 + differ + extent + surg, data = train.obslev5.death)
summary(pre.lev5.death) #significant


### change value of rx (lev+5FU = 1, obs = 0) for calculation
X.train.death <- cbind(train.obslev5.death$rx, train.obslev5.death$age, train.obslev5.death$node4, train.obslev5.death$differ, train.obslev5.death$extent,train.obslev5.death$surg)
X.valid.death <- cbind(valid.obslev5.death$rx, valid.obslev5.death$age, valid.obslev5.death$node4, valid.obslev5.death$differ, valid.obslev5.death$extent, valid.obslev5.death$surg)

#extract the coefficient from previous selected model without NA
coef.lev5.death <- coef(pre.lev5.death)
coef.lev5.death <- as.numeric(coef.lev5.death)
coef.lev5.death <- coef.lev5.death[!is.na(coef.lev5.death)]

risk.train.death <- X.train.death %*% coef.lev5.death
risk.valid.death <- X.valid.death %*% coef.lev5.death

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev5.death <- as.numeric(quantile(risk.train.death, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data

risk.valid.death <- as.factor(cut(risk.valid.death, breaks =c(-Inf, quant.lev5.death, Inf), labels = FALSE))

valid.obslev5.death$level <- risk.valid.death

plot(survfit(Surv(valid.obslev5.death$time, valid.obslev5.death$status) ~ level, data=valid.obslev5.death), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability")
     
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev5.death, data=valid.obslev5.death), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

###################### Obs + Lev: Mortality ########################

set.seed(10)
train.row.death <- sample(seq_len(nrow(death)), size = 540) #set randomly selected 183 rows in training set

train.death <- death[train.row.death,]
rownames(train.death) <- NULL
valid.death <- death[-train.row.death,]
rownames(valid.death) <- NULL

#####################
train.obslev.death <- subset(train.death, train.death$rx == "Obs" | train.death$rx == "Lev")
valid.obslev.death <- subset(valid.death, valid.death$rx == "Obs" | valid.death$rx == "Lev")
rownames(train.obslev.death) <- NULL
train.obslev.death$rx <- droplevels(train.obslev.death$rx)
rownames(valid.obslev.death) <- NULL
valid.obslev.death$rx <- droplevels(valid.obslev.death$rx)

full.pre.lev.death <- coxph(Surv(train.obslev.death $time, train.obslev.death $status)~rx + sex + age + obstruct + perfor + adhere + node4 + factor(differ) + factor(extent) + surg, data = train.obslev.death) 
step(full.pre.lev.death, direction = c("both"))

# Thus, our final prediction model for lev5FU is
pre.lev.death <- coxph(formula = Surv(train.obslev.death$time, train.obslev.death$status) ~ age + obstruct + node4 + factor(differ) + factor(extent), data = train.obslev.death)
summary(pre.lev.death) #significant


### change value of rx (lev = 1, obs = 0) for calculation
X.train.death.lev <- model.matrix(pre.lev.death)
X.valid.death.lev <- model.matrix(~age + obstruct + node4 + factor(differ) + factor(extent), data = valid.obslev.death)
X.valid.death.lev <- X.valid.death.lev[,-1]
#extract the coefficient from previous selected model without NA
coef.lev.death <- coef(pre.lev.death)
coef.lev.death <- as.numeric(coef.lev.death)
coef.lev.death <- coef.lev.death[!is.na(coef.lev.death)]

risk.train.death.lev <- X.train.death.lev %*% coef.lev.death
risk.valid.death.lev <- X.valid.death.lev %*% coef.lev.death

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev.death <- as.numeric(quantile(risk.train.death.lev, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data


risk.valid.death.lev <- as.factor(cut(risk.valid.death.lev, breaks =c(-Inf, quant.lev.death, Inf), labels = FALSE))
valid.obslev.death$level <- risk.valid.death.lev

plot(survfit(Surv(valid.obslev.death$time, valid.obslev.death$status) ~ level, data=valid.obslev.death), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev.death, data=valid.obslev.death), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

########################################################
set.seed(1023)
train.row.death <- sample(seq_len(nrow(death)), size = 540) #set randomly selected 183 rows in training set

train.death <- death[train.row.death,]
rownames(train.death) <- NULL
valid.death <- death[-train.row.death,]
rownames(valid.death) <- NULL

#####################
train.obslev.death <- subset(train.death, train.death$rx == "Obs" | train.death$rx == "Lev")
valid.obslev.death <- subset(valid.death, valid.death$rx == "Obs" | valid.death$rx == "Lev")
rownames(train.obslev.death) <- NULL
train.obslev.death$rx <- droplevels(train.obslev.death$rx)
rownames(valid.obslev.death) <- NULL
valid.obslev.death$rx <- droplevels(valid.obslev.death$rx)

full.pre.lev.death <- coxph(Surv(train.obslev.death $time, train.obslev.death $status)~rx + sex + age + obstruct + perfor + adhere + node4 + differ + extent + surg, data = train.obslev.death) 
step(full.pre.lev.death, direction = c("both"))

# Thus, our final prediction model for lev5FU is
pre.lev.death <- coxph(formula = Surv(train.obslev.death$time, train.obslev.death$status) ~ age + obstruct + node4 + differ + extent + surg, data = train.obslev.death)
summary(pre.lev.death) #significant


### change value of rx (lev = 1, obs = 0) for calculation
X.train.death.lev <- cbind(train.obslev.death$age, train.obslev.death$obstruct, train.obslev.death$node4, train.obslev.death$differ, train.obslev.death$extent, train.obslev.death$surg)
X.valid.death.lev <- cbind(valid.obslev.death$age, valid.obslev.death$obstruct, valid.obslev.death$node4, valid.obslev.death$differ, valid.obslev.death$extent, valid.obslev.death$surg)

#extract the coefficient from previous selected model without NA
coef.lev.death <- coef(pre.lev.death)
coef.lev.death <- as.numeric(coef.lev.death)
coef.lev.death <- coef.lev.death[!is.na(coef.lev.death)]

risk.train.death.lev <- X.train.death.lev %*% coef.lev.death
risk.valid.death.lev <- X.valid.death.lev %*% coef.lev.death

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev.death <- as.numeric(quantile(risk.train.death.lev, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data


risk.valid.death.lev <- as.factor(cut(risk.valid.death.lev, breaks =c(-Inf, quant.lev.death, Inf), labels = FALSE))
valid.obslev.death$level <- risk.valid.death.lev

plot(survfit(Surv(valid.obslev.death$time, valid.obslev.death$status) ~ level, data=valid.obslev.death), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev.death, data=valid.obslev.death), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)
