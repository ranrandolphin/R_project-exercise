library(survival)
data(colon)
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

recurrence <- subset(new.data, new.data$etype == 1)
rownames(recurrence) <- NULL
death <- subset(new.data, new.data$etype == 2)
rownames(death) <- NULL

recur.obs <- subset(recurrence, recurrence$rx == "Obs")
rownames(recur.obs) <- NULL
recur.obs$rx <- droplevels(recur.obs$rx)
recur.lev <- recurrence[which (recurrence$rx == "Lev"),]
rownames(recur.lev) <- NULL
recur.lev$rx <- droplevels(recur.lev$rx)
recur.5fu <- subset(recurrence, recurrence$rx == "Lev+5FU")
rownames(recur.5fu) <- NULL
recur.5fu$rx <- droplevels(recur.5fu$rx)


######################### Prediction : Recurrence ####################################
########## Training set and Validation Set: Option 1 ##########
# randomly select 2/3 (540) of the patients from the observational group as training data, and use the rest 1/3(264) as the validation data
# We do not want to involve the treatment effect here to affect the proportional hazard model
set.seed(10)
train.row <- sample(seq_len(nrow(recurrence)), size = 540) #set randomly selected 183 rows in training set

train.rec <- recurrence[train.row,]
rownames(train.rec) <- NULL
valid.rec <- recurrence[-train.row,]
rownames(valid.rec) <- NULL

#######recurrence + Lev+5FU set in training model###########

train.obslev5 <- subset(train.rec, train.rec$rx == "Obs" | train.rec$rx == "Lev+5FU")
valid.obslev5 <- subset(valid.rec, valid.rec$rx == "Obs" | valid.rec$rx == "Lev+5FU")
rownames(train.obslev5) <- NULL
train.obslev5$rx <- droplevels(train.obslev5$rx)
rownames(valid.obslev5) <- NULL
valid.obslev5$rx <- droplevels(valid.obslev5$rx)
full.pre.lev5 <- coxph(Surv(train.obslev5$time, train.obslev5$status)~rx + sex + age + obstruct + perfor + adhere + node4 + factor(differ) + factor(extent) + surg, data = train.obslev5) 
step(full.pre.lev5, direction = c("both"))

# Thus, our final prediction model for lev5FU is
pre.lev5 <- coxph(formula = Surv(train.obslev5$time, train.obslev5$status) ~ rx + node4 + factor(differ) + factor(extent) + surg, data = train.obslev5)
summary(pre.lev5) #significant

# Diagnosis####


### change value of rx (lev+5FU = 1, obs = 0) for calculation
X.train <- model.matrix(pre.lev5)
X.valid <- model.matrix(~rx + node4 + factor(differ) + factor(extent) + surg , data = valid.obslev5)
X.valid <- X.valid[,-1]

#extract the coefficient from previous selected model without NA
coef.lev5 <- coef(pre.lev5)
coef.lev5 <- as.numeric(coef.lev5)
coef.lev5 <- coef.lev5[!is.na(coef.lev5)]

risk.train <- X.train %*% coef.lev5
risk.valid <- X.valid %*% coef.lev5

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev5.rec <- as.numeric(quantile(risk.train, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data

risk.valid <- as.factor(cut(risk.valid, breaks =c(-Inf, quant.lev5.rec, Inf), labels = FALSE))

valid.obslev5$level <- risk.valid

plot(survfit(Surv(valid.obslev5$time, valid.obslev5$status) ~ level, data=valid.obslev5), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev5, data=valid.obslev5), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

###########################Prediction: Obs + Lev+5FU############
######### different seed###################
set.seed(818)
train.row <- sample(seq_len(nrow(recurrence)), size = 540) #set randomly selected 183 rows in training set

train.rec <- recurrence[train.row,]
rownames(train.rec) <- NULL
valid.rec <- recurrence[-train.row,]
rownames(valid.rec) <- NULL

#######recurrence + Lev+5FU set in training model###########

train.obslev5 <- subset(train.rec, train.rec$rx == "Obs" | train.rec$rx == "Lev+5FU")
valid.obslev5 <- subset(valid.rec, valid.rec$rx == "Obs" | valid.rec$rx == "Lev+5FU")
rownames(train.obslev5) <- NULL
train.obslev5$rx <- droplevels(train.obslev5$rx)
rownames(valid.obslev5) <- NULL
valid.obslev5$rx <- droplevels(valid.obslev5$rx)
full.pre.lev5 <- coxph(Surv(train.obslev5$time, train.obslev5$status)~rx + sex + age + obstruct + perfor + adhere + node4 + factor(differ) + factor(extent) + surg, data = train.obslev5) 
step(full.pre.lev5, direction = c("both"))


# Thus, our final prediction model for lev5FU is
pre.lev5 <- coxph(formula = Surv(train.obslev5$time, train.obslev5$status) ~ rx + perfor + node4 + factor(differ) + factor(extent) + surg, data = train.obslev5)
summary(pre.lev5) #significant


### change value of rx (lev+5FU = 1, obs = 0) for calculation
X.train <- model.matrix(pre.lev5)
X.valid <- model.matrix(~rx + perfor + node4 + factor(differ) + factor(extent) + surg , data = valid.obslev5)
X.valid <- X.valid[,-1]

#extract the coefficient from previous selected model without NA
coef.lev5 <- coef(pre.lev5)
coef.lev5 <- as.numeric(coef.lev5)
coef.lev5 <- coef.lev5[!is.na(coef.lev5)]

risk.train <- X.train %*% coef.lev5
risk.valid <- X.valid %*% coef.lev5

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev5.rec <- as.numeric(quantile(risk.train, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data
risk.valid <- as.factor(cut(risk.valid, breaks =c(-Inf, quant.lev5.rec, Inf), labels = FALSE))

valid.obslev5$level <- risk.valid

plot(survfit(Surv(valid.obslev5$time, valid.obslev5$status) ~ level, data=valid.obslev5), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev5, data=valid.obslev5), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

##############################################
set.seed(1)
train.row <- sample(seq_len(nrow(recurrence)), size = 540) #set randomly selected 183 rows in training set

train.rec <- recurrence[train.row,]
rownames(train.rec) <- NULL
valid.rec <- recurrence[-train.row,]
rownames(valid.rec) <- NULL

#######recurrence + Lev set in training model###########

train.obslev <- subset(train.rec, train.rec$rx == "Obs" | train.rec$rx == "Lev")
valid.obslev <- subset(valid.rec, valid.rec$rx == "Obs" | valid.rec$rx == "Lev")
rownames(train.obslev) <- NULL
train.obslev$rx <- droplevels(train.obslev$rx)
rownames(valid.obslev) <- NULL
valid.obslev$rx <- droplevels(valid.obslev$rx)
full.pre.lev <- coxph(Surv(train.obslev$time, train.obslev$status)~rx + sex + age + obstruct + perfor + adhere + node4 + factor(differ) + factor(extent) + surg, data = train.obslev) 
step(full.pre.lev, direction = c("both"))

# Thus, our final prediction model for lev5FU is
pre.lev <- coxph(formula = Surv(train.obslev$time, train.obslev$status) ~ obstruct + adhere + node4 + factor(differ) + factor(extent) + surg, data = train.obslev)
summary(pre.lev) #significant


### change value of rx (lev+5FU = 1, obs = 0) for calculation
X.train.lev <- model.matrix(pre.lev)
X.valid.lev <- model.matrix(~obstruct + adhere + node4 + factor(differ) + factor(extent) + surg, data = valid.obslev)
X.valid.lev <- X.valid.lev[,-1]

#extract the coefficient from previous selected model without NA
coef.lev <- coef(pre.lev)
coef.lev <- as.numeric(coef.lev)
coef.lev <- coef.lev[!is.na(coef.lev)]

risk.train.lev <- X.train.lev %*% coef.lev
risk.valid.lev <- X.valid.lev %*% coef.lev

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev.rec <- as.numeric(quantile(risk.train.lev, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data
risk.valid.lev <- as.factor(cut(risk.valid.lev, breaks =c(-Inf, quant.lev.rec, Inf), labels = FALSE))
valid.obslev$level <- risk.valid.lev

plot(survfit(Surv(valid.obslev$time, valid.obslev$status) ~ level, data=valid.obslev), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev, data=valid.obslev), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

#################### Different set.seed for randomly selecting training set###############

set.seed(500)
train.row <- sample(seq_len(nrow(recurrence)), size = 540) #set randomly selected 183 rows in training set

train.rec <- recurrence[train.row,]
rownames(train.rec) <- NULL
valid.rec <- recurrence[-train.row,]
rownames(valid.rec) <- NULL

#######recurrence + Lev set in training model###########

train.obslev <- subset(train.rec, train.rec$rx == "Obs" | train.rec$rx == "Lev")
valid.obslev <- subset(valid.rec, valid.rec$rx == "Obs" | valid.rec$rx == "Lev")
rownames(train.obslev) <- NULL
train.obslev$rx <- droplevels(train.obslev$rx)
rownames(valid.obslev) <- NULL
valid.obslev$rx <- droplevels(valid.obslev$rx)
full.pre.lev <- coxph(Surv(train.obslev$time, train.obslev$status)~rx + sex + age + obstruct + perfor + adhere + node4 + differ + extent + surg, data = train.obslev) 
step(full.pre.lev, direction = c("both"))

# Thus, our final prediction model for lev is
pre.lev <- coxph(formula = Surv(train.obslev$time, train.obslev$status) ~ age + obstruct + perfor+ node4 + differ + surg, data = train.obslev)
summary(pre.lev) #significant


### change value of rx (lev+5FU = 1, obs = 0) for calculation
X.train.lev <- cbind(train.obslev$age, train.obslev$obstruct, train.obslev$perfor, train.obslev$node4, train.obslev$differ, train.obslev$extent)
X.valid.lev <- cbind(valid.obslev$age, valid.obslev$obstruct, valid.obslev$perfor, valid.obslev$node4, valid.obslev$differ, valid.obslev$extent)

#extract the coefficient from previous selected model without NA
coef.lev <- coef(pre.lev)
coef.lev <- as.numeric(coef.lev)
coef.lev <- coef.lev[!is.na(coef.lev)]

risk.train.lev <- X.train.lev %*% coef.lev
risk.valid.lev <- X.valid.lev %*% coef.lev

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev.rec <- as.numeric(quantile(risk.train.lev, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data
risk.valid.lev <- as.factor(cut(risk.valid.lev, breaks =c(-Inf, quant.lev.rec, Inf), labels = FALSE))

valid.obslev$level <- risk.valid.lev

plot(survfit(Surv(valid.obslev$time, valid.obslev$status) ~ level, data=valid.obslev), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev, data=valid.obslev), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

