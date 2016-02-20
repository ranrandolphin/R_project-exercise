##########Data Process################################################
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

########################################################################

#table(lev5$surg)
#data.5year <- subset(data, data$time < 1826.25)
#rownames(data.5year) <- NULL

#For male in 5 year data
#male5 <- subset(data.5year, data.5year$sex == 1)
#female5 <- subset(data.5year, data.5year$sex == 0)

#data.age <- subset(data.5year, data.5year$age > 61)
#rownames(data.age) <- NULL 

#2) Unadjusted and adjusted treatment effect
# unadjusted Group, the univariate Cox model including only treatment variable(Lev or Lev+5-Fu)

############################# recurrence group ########################
########################Observation v.s. Lev ##########################
# a. For observation and Lev group
recur.obs.lev <- rbind(recur.obs, recur.lev)

# a). Unadjusted effect
Surv.obs.lev <- Surv(recur.obs.lev$time, recur.obs.lev$status)
recur.olev.fit <- survfit(Surv.obs.lev ~ rx, conf.int = 0.95, se.fit = TRUE, type = "kaplan-meier", data = recur.obs.lev)
p.recur.olev <- plot(recur.olev.fit, col = 1:2, xlab = "Days from randomization", 
                ylab = "Recurrence", 
                main = "Figure 1: Survival curves over Lev Treatments & Observation Groups")
legend("topright", legend = c("Obs", "Lev"), lty = c(1,2), col = 1:2)
model.olev.un <- coxph(Surv.obs.lev ~ rx , data = recur.obs.lev)
summary(model.olev.un) #not significant
survdiff(Surv.obs.lev ~ rx, data = recur.obs.lev)

#b). Adjusted effect???????????
model.olev.ad <- coxph(Surv.obs.lev ~ rx + extent + time + node4, data = recur.obs.lev)
summary(model.olev.ad)
#c). secondary analysis
model.olev.sed <- coxph(Surv.obs.lev ~ rx + age + obstruct+ perfor + adhere + differ, data = recur.obs.lev)
summary(model.olev.sed)

###############Observation v.s. Lev+5-Fu#######################################
#b. For observation and Lev+5-FU
recur.obs.5fu <- rbind(recur.obs, recur.5fu)


#a). Unadjusted effect
Surv.obs.lev5 <- Surv(recur.obs.5fu$time, recur.obs.5fu$status)
recur.olev5.fit <- survfit(Surv.obs.lev5 ~ rx, conf.int = 0.95, se.fit = TRUE, type = "kaplan-meier", data = recur.obs.5fu)
p.recur.olev <- plot(recur.olev5.fit, col = c(1,3), xlab = "Days from randomization", 
                     ylab = "Recurrence", 
                     main = "Figure 1: Survival curves over Lev+5FU Treatments & Observation Groups")
legend("topright", legend = c("Obs", "Lev+5-FU"), lty = c(1,3), col = c(1,3))
model.olev5.un <- coxph(Surv.obs.lev5 ~ rx , data = recur.obs.5fu)
summary(model.olev5.un) # significant under LR, Wald, and Score test
survdiff(Surv.obs.lev5 ~ rx, data = recur.obs.5fu)

#b). Adjusted effect
model.olev5.ad <- coxph(Surv.obs.lev5 ~ rx + extent + time + node4, data = recur.obs.5fu)
summary(model.olev5.ad)
#c). secondary analysis
summary(coxph(Surv.obs.lev5 ~ rx + sex , data = recur.obs.5fu))
summary(coxph(Surv.obs.lev5 ~ rx + age , data = recur.obs.5fu))
summary(coxph(Surv.obs.lev5 ~ rx + obstruct, data = recur.obs.5fu))
summary(coxph(Surv.obs.lev5 ~ rx + perfor , data = recur.obs.5fu))
summary(coxph(Surv.obs.lev5 ~ rx + adhere , data = recur.obs.5fu))
summary(coxph(Surv.obs.lev5 ~ rx + differ, data = recur.obs.5fu))
summary(coxph(Surv.obs.lev5 ~ rx + extent, data = recur.obs.5fu))
summary(coxph(Surv.obs.lev5 ~ rx + node4, data = recur.obs.5fu))
summary(coxph(Surv.obs.lev5 ~ rx + surg, data = recur.obs.5fu))
#full.lev5 <- coxph(Surv.obs.lev5 ~ rx + age + obstruct+ perfor + adhere + differ, data = recur.obs.5fu)
#summary(model.olev5.sed)


###########################   Death Group    ##########################
death.obs <- subset(death, death$rx == "Obs")
rownames(death.obs) <- NULL
death.obs$rx <- droplevels(death.obs$rx)
death.lev <- death[which (death$rx == "Lev"),]
rownames(death.lev) <- NULL
death.lev$rx <- droplevels(death.lev$rx)
death.5fu <- subset(death, death$rx == "Lev+5FU")
rownames(death.5fu) <- NULL
death.5fu$rx <- droplevels(death.5fu$rx)
##################Death: obs v.s. Lev ################################
# a. For observation and Lev group
death.obs.lev <- rbind(death.obs, death.lev)
# a). Unadjusted effect
Surv.death.lev<- Surv(death.obs.lev$time, death.obs.lev$status)
death.olev.fit <- survfit(Surv.death.lev ~ rx, conf.int = 0.95, se.fit = TRUE, type = "kaplan-meier", data = death.obs.lev)
plot(death.olev.fit, col = 1:2, xlab = "Days from randomization", 
                     ylab = "Mortality ", 
                     main = "Figure 1: Survival curves over Lev Treatments & Observation Groups")
legend("topright", legend = c("Obs", "Lev"), lty = c(1,2), col = 1:2)
death.olev.un <- coxph(Surv.death.lev ~ rx , data = death.obs.lev)
summary(death.olev.un) #not significant - do not need to keep going

###################Death: obs v.s. Lev5FU#############################
#b. For observation and Lev+5-FU
death.obs.5fu <- rbind(death.obs, death.5fu)

#a). Unadjusted effect
death.obs.lev5 <- Surv(death.obs.5fu$time, death.obs.5fu$status)
death.olev5.fit <- survfit(death.obs.lev5 ~ rx, conf.int = 0.95, se.fit = TRUE, type = "kaplan-meier", data = death.obs.5fu)
plot(death.olev5.fit, col = c(1,3), xlab = "Days from randomization", 
                     ylab = "Mortality", 
                     main = "Figure 1: Survival curves over Lev+5FU Treatments & Observation Groups")
legend("topright", legend = c("Obs", "Lev+5-FU"), lty = c(1,3), col = c(1,3))
death.olev5.un <- coxph(death.obs.lev5 ~ rx , data = death.obs.5fu)
summary(death.olev5.un) # significant under LR, Wald, and Score test by 0.05
survdiff(death.obs.lev5 ~ rx, data = death.obs.5fu)

#b). Adjusted effect
death.olev5.ad <- coxph(death.obs.lev5 ~ rx + extent + time + node4, data = death.obs.5fu)
summary(model.olev5.ad)
#c). secondary analysis
summary(coxph(death.obs.lev5 ~ rx + sex , data = death.obs.5fu))
summary(coxph(death.obs.lev5 ~ rx + sex + age , data = death.obs.5fu))
summary(coxph(death.obs.lev5 ~ rx + sex + age + obstruct, data = death.obs.5fu))
summary(coxph(death.obs.lev5 ~ rx + sex + age + obstruct+ perfor , data = death.obs.5fu))
summary(coxph(death.obs.lev5 ~ rx + sex + age + obstruct+ perfor + adhere , data = death.obs.5fu))
summary(coxph(death.obs.lev5 ~ rx + sex + age + obstruct+ perfor + adhere + differ, data = death.obs.5fu))
summary(coxph(death.obs.lev5 ~ rx + sex + age + obstruct+ perfor + adhere + differ + extent, data = death.obs.5fu))
summary(coxph(death.obs.lev5 ~ rx + sex + age + obstruct+ perfor + adhere + differ + extent + node4, data = death.obs.5fu))
summary(coxph(death.obs.lev5 ~ rx + sex + age + obstruct+ perfor + adhere + differ + extent + node4 + surg, data = death.obs.5fu))


##############Feature Important : 
########1) Obsevation and lev Groups
# Feature Explanined: importance

obs <- subset(new.data, new.data$rx == "Obs" | new.data$rx == "Lev")
obs$rx <- droplevels(obs$rx)
rownames(obs) <- NULL
obs.feature<- obs[-c(1,2,3,9,14,15)]
# feature selection
apply(obs,2,var)
pr.out.obs = prcomp(obs.feature, scale = TRUE)
plot(pr.out.obs) # We found that first principle component has highest variace corresponding to each features, thus PC1 is the best choice compared with other principle components.
barplot(pr.out.obs$sdev/pr.out.obs$sdev[1])
pr.out.obs$scale
pr.out.obs$rotation

# Thus, the time, status, node4, nodes are relative importnat variables from PC1
biplot(pr.out.obs, scale = 0)
# We found that each variable are correlated with some of the rest variables. Although nodes and node4 almost overlap each other, we have checked the correlation matrix, and then decide to remain them for the remaining prediction.


pr.var = pr.out.obs$sdev^2 # the variance of each principal component
# to compute the propotion of variace explained by each principal component
pve = pr.var/sum(pr.var)
par(mfrow = c(1,2))
# plot the PVE explanied by each component
plot(pve, xlab = "Principal Component", ylab = "Proportion of Variance Explained",main = "Obs v.s. Lev Groups", ylim = c(0.05,0.2), type = "b")
# plot cucmulative PVE
plot(cumsum(pve), xlab = "Principal Component", ylab = "Cumulative Proportion of Variance Explained", main = "Obs v.s. Lev Groups", ylim = c(0,1), type = "b")

###########2) Feature Importantce (Obseravation & Lev+5FU groups)
obs2 <- subset(new.data, new.data$rx == "Obs" | new.data$rx == "Lev+5FU")
rownames(obs2) <- NULL
obs2$rx <- droplevels(obs2$rx)
obs2.feature<- obs2[-c(1,2,3,9,14,15)]
# feature selection
apply(obs2,2,var)
pr.out.obs2 = prcomp(obs2.feature, scale = TRUE)
plot(pr.out.obs2) # We found that first principle component has highest variace corresponding to each features, thus PC1 is the best choice compared with other principle components.
barplot(pr.out.obs2$sdev/pr.out.obs2$sdev[1])
pr.out.obs2$scale
pr.out.obs2$rotation

# Thus, the time, status, node4, nodes are relative importnat variables from PC1
par(mfrow = c(1,1))
biplot(pr.out.obs2, scale = 0)
# We found that each variable are correlated with some of the rest variables. Although nodes and node4 almost overlap each other, we have checked the correlation matrix, and then decide to remain them for the remaining prediction.


pr.var2 = pr.out.obs2$sdev^2 # the variance of each principal component
# to compute the propotion of variace explained by each principal component
pve2 = pr.var2/sum(pr.var2)
par(mfrow = c(1,2))
# plot the PVE explanied by each component
plot(pve2, xlab = "Principal Component", ylab = "Proportion of Variance Explained", main = "Obs v.s. Lev+5FU Groups", ylim = c(0.05,0.2), type = "b")
# plot cucmulative PVE
plot(cumsum(pve2), xlab = "Principal Component", ylab = "Cumulative Proportion of Variance Explained",main = "Obs v.s. Lev+5FU Groups", ylim = c(0,1), type = "b")
par(mfrow = c(1,1)) # back to original set of graphs

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
full.pre.lev5 <- coxph(Surv(train.obslev5$time, train.obslev5$status)~rx + sex + age + obstruct + perfor + adhere + node4 + differ + extent + surg, data = train.obslev5) 
step(full.pre.lev5, direction = c("both"))

# Thus, our final prediction model for lev5FU is
pre.lev5 <- coxph(formula = Surv(train.obslev5$time, train.obslev5$status) ~ rx + node4 + differ + extent + surg, data = train.obslev5)
summary(pre.lev5) #significant


### change value of rx (lev+5FU = 1, obs = 0) for calculation
X.train <- cbind(train.obslev5$rx, train.obslev5$node4, train.obslev5$differ, train.obslev5$extent, train.obslev5$surg)
X.valid <- cbind(valid.obslev5$rx,  valid.obslev5$node4, valid.obslev5$differ, valid.obslev5$extent, valid.obslev5$surg)

#extract the coefficient from previous selected model without NA
coef.lev5 <- coef(pre.lev5)
coef.lev5 <- as.numeric(coef.lev5)
coef.lev5 <- coef.lev5[!is.na(coef.lev5)]

risk.train <- X.train %*% coef.lev5
risk.valid <- X.valid %*% coef.lev5

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev5.rec <- as.numeric(quantile(risk.train, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data
q1.lev5.rec <- quant.lev5.rec[1]
q3.lev5.rec <- quant.lev5.rec[2]

risk.valid[risk.valid > q3.lev5.rec, ] <- 3
risk.valid[risk.valid <= q3.lev5.rec & risk.valid > q1.lev5.rec, ] <- 2
risk.valid[risk.valid <= q1.lev5.rec, ] <- 1
risk.valid <- as.factor(risk.valid)
valid.obslev5$level <- risk.valid

plot(survfit(Surv(valid.obslev5$time, valid.obslev5$status) ~ level, data=valid.obslev5), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability",
     main = "Observed and expected survival curves for test colon cancer patients in validation data")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev5, data=valid.obslev5), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

#######recurrence + Lev set in training model###########

train.obslev <- subset(train.rec, train.rec$rx == "Obs" | train.rec$rx == "Lev")
valid.obslev <- subset(valid.rec, valid.rec$rx == "Obs" | valid.rec$rx == "Lev")
rownames(train.obslev) <- NULL
rownames(valid.obslev) <- NULL
full.pre.lev <- coxph(Surv(train.obslev$time, train.obslev$status)~rx + sex + age + obstruct + perfor + adhere + node4 + differ + extent + surg, data = train.obslev) 
step(full.pre.lev, direction = c("both"))

# Thus, our final prediction model for lev5FU is
pre.lev <- coxph(formula = Surv(train.obslev$time, train.obslev$status) ~ obstruct + adhere + node4 + extent + surg, data = train.obslev)
summary(pre.lev) #significant


### change value of rx (lev+5FU = 1, obs = 0) for calculation
X.train.lev <- cbind(train.obslev$obstruct, train.obslev$adhere, train.obslev$node4, train.obslev$extent, train.obslev$surg)
X.valid.lev <- cbind(valid.obslev$obstruct, valid.obslev$adhere, valid.obslev$node4, valid.obslev$extent, valid.obslev$surg)

#extract the coefficient from previous selected model without NA
coef.lev <- coef(pre.lev)
coef.lev <- as.numeric(coef.lev)
coef.lev <- coef.lev[!is.na(coef.lev)]

risk.train.lev <- X.train.lev %*% coef.lev
risk.valid.lev <- X.valid.lev %*% coef.lev

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev.rec <- as.numeric(quantile(risk.train.lev, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data
q1.lev.rec <- quant.lev.rec[1]
q3.lev.rec <- quant.lev.rec[2]

#step by step
risk.valid.lev[risk.valid.lev > q3.lev.rec,] <- 3
risk.valid.lev[risk.valid.lev <= q3.lev.rec & risk.valid.lev > q1.lev.rec, ] <- 2
risk.valid.lev[risk.valid.lev <= q1.lev.rec, ] <- 1
risk.valid.lev <- as.factor(risk.valid.lev)
valid.obslev$level <- risk.valid.lev

plot(survfit(Surv(valid.obslev$time, valid.obslev$status) ~ level, data=valid.obslev), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability",
     main = "Observed and expected survival curves for test colon cancer patients in validation data")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev, data=valid.obslev), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

############################### 我是华丽丽的分界线#########################################################

############################## Prediction: Death (Mortality) #############

########## Training set and Validation Set ##########
set.seed(1023)
train.row.death <- sample(seq_len(nrow(death)), size = 540) #set randomly selected 183 rows in training set

train.death <- death[train.row.death,]
rownames(train.death) <- NULL
valid.death <- death[-train.row.death,]
rownames(valid.death) <- NULL

#######Mortality: Lev+5FU set in training model###########

train.obslev5.death <- subset(train.death, train.death$rx == "Obs" | train.death$rx == "Lev+5FU")
valid.obslev5.death <- subset(valid.death, valid.death$rx == "Obs" | valid.death$rx == "Lev+5FU")
rownames(train.obslev5.death) <- NULL
rownames(valid.obslev5.death) <- NULL
full.pre.lev5.death <- coxph(Surv(train.obslev5.death $time, train.obslev5.death $status)~rx + sex + age + obstruct + perfor + adhere + node4 + differ + extent + surg, data = train.obslev5.death) 
step(full.pre.lev5.death, direction = c("both"))

# Thus, our final prediction model for lev5FU is
pre.lev5.death <- coxph(formula = Surv(train.obslev5.death$time, train.obslev5.death$status) ~ age + node4 + differ + extent + surg, data = train.obslev5.death)
summary(pre.lev5.death) #significant


### change value of rx (lev+5FU = 1, obs = 0) for calculation
X.train.death <- cbind(train.obslev5.death$age, train.obslev5.death$node4, train.obslev5.death$differ, train.obslev5.death$extent, train.obslev5.death$surg)
X.valid.death <- cbind(valid.obslev5.death$age, valid.obslev5.death$node4, valid.obslev5.death$differ, valid.obslev5.death$extent, valid.obslev5.death$surg)

#extract the coefficient from previous selected model without NA
coef.lev5.death <- coef(pre.lev5.death)
coef.lev5.death <- as.numeric(coef.lev5.death)
coef.lev5.death <- coef.lev5.death[!is.na(coef.lev5.death)]

risk.train.death <- X.train.death %*% coef.lev5.death
risk.valid.death <- X.valid.death %*% coef.lev5.death

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev5.death <- as.numeric(quantile(risk.train.death, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data
q1.lev5.death <- quant.lev5.death[1]
q3.lev5.death <- quant.lev5.death[2]

risk.valid.death[risk.valid.death > q3.lev5.death, ] <- 3
risk.valid.death[risk.valid.death <= q3.lev5.death & risk.valid.death > q1.lev5.death, ] <- 2
risk.valid.death[risk.valid.death <= q1.lev5.death, ] <- 1
risk.valid.death <- as.factor(risk.valid.death)
valid.obslev5.death$level <- risk.valid.death

plot(survfit(Surv(valid.obslev5.death$time, valid.obslev5.death$status) ~ level, data=valid.obslev5.death), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability",
     main = "Observed and expected survival curves for test colon cancer patients in validation data")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev5.death, data=valid.obslev5.death), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

################ Mortality: Obs + Lev ########################
train.obslev.death <- subset(train.death, train.death$rx == "Obs" | train.death$rx == "Lev")
valid.obslev.death <- subset(valid.death, valid.death$rx == "Obs" | valid.death$rx == "Lev")
rownames(train.obslev.death) <- NULL
rownames(valid.obslev.death) <- NULL
full.pre.lev.death <- coxph(Surv(train.obslev.death $time, train.obslev.death $status)~rx + sex + age + obstruct + perfor + adhere + node4 + differ + extent + surg, data = train.obslev.death) 
step(full.pre.lev.death, direction = c("both"))

# Thus, our final prediction model for lev5FU is
pre.lev.death <- coxph(formula = Surv(train.obslev.death$time, train.obslev.death$status) ~ age + obstruct + node4 + extent + surg, data = train.obslev.death)
summary(pre.lev.death) #significant


### change value of rx (lev = 1, obs = 0) for calculation
X.train.death.lev <- cbind(train.obslev.death$age, train.obslev.death$obstruct, train.obslev.death$node4, train.obslev.death$extent, train.obslev.death$surg)
X.valid.death.lev <- cbind(valid.obslev.death$age, valid.obslev.death$obstruct, valid.obslev.death$node4, valid.obslev.death$extent, valid.obslev.death$surg)

#extract the coefficient from previous selected model without NA
coef.lev.death <- coef(pre.lev.death)
coef.lev.death <- as.numeric(coef.lev.death)
coef.lev.death <- coef.lev.death[!is.na(coef.lev.death)]

risk.train.death.lev <- X.train.death.lev %*% coef.lev.death
risk.valid.death.lev <- X.valid.death.lev %*% coef.lev.death

#levels data with 1, 2, 3 levels corresponding to 33th and 67th percentile rank
quant.lev.death <- as.numeric(quantile(risk.train.death.lev, c(0.33, 0.67))) #obtain the 33th and 67th percentile from training data
q1.lev.death <- quant.lev.death[1]
q3.lev.death <- quant.lev.death[2]

risk.valid.death.lev[risk.valid.death.lev > q3.lev.death, ] <- 3
risk.valid.death.lev[risk.valid.death.lev <= q3.lev.death & risk.valid.death.lev > q1.lev.death, ] <- 2
risk.valid.death.lev[risk.valid.death.lev <= q1.lev.death, ] <- 1
risk.valid.death.lev <- as.factor(risk.valid.death.lev)
valid.obslev.death$level <- risk.valid.death.lev

plot(survfit(Surv(valid.obslev.death$time, valid.obslev.death$status) ~ level, data=valid.obslev.death), mark.time=FALSE,
     lty = 1, col = 1:3, xlab = "Observation Times(days)", ylab = "Survival Probability",
     main = "Observed and expected survival curves for test colon cancer patients in validation data")
legend("topright", legend = c("Low Risk", "Moderate Risk", "High Risk"), lty = 1, col = 1:3)
lines(survexp( ~ level, ratetable=pre.lev.death, data=valid.obslev.death), lty = 2, col=1:3)
legend("bottomright", legend = c("Observed", "Expected"), lty = 1:2, col = 1)

#####################################################################################
