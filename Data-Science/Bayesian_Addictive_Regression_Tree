#### Clean Data ####
public <- read.table("C:/Users/Ran/Desktop/2690G/public.dat", header=F)

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


public <- read.table("C:/Users/Ran/Desktop/2690G/public.dat", header=F)

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

wage <- na.omit(data.frame(Y,D,Z,C,H,EMPFT,EMPPT,NMGRS,w1,INCTIME,
                           FIRSTINC,BONUS,MEALS,OPEN,HRSOPEN,PSODA,
                           PFRY,PENTREE,NRGES,NREGS11))


hist(wage$Y)
fit <- lm(Y ~  D+Z+C+H+EMPFT+EMPPT+NMGRS+w1+INCTIME+
            FIRSTINC+BONUS+MEALS+OPEN+HRSOPEN+PSODA+
            PFRY+PENTREE+NRGES+NREGS11, data = wage)
summary(fit)


####1. The method you used (BART Method)

library(BayesTree)
new_wage <- data.frame(wage[2:20])
new_wage$indicator <- NULL
Xind <- makeind(new_wage, all=TRUE)  
out <- bart(x.train=new_wage, y.train=wage$Y)
plot(out)
summary(out$sigma) 


par(mfrow=c(1,1), bty="n")
dim(out$varcount)
boxplot(out$varcount, names=dimnames(Xind)[[2]], cex.axis=0.5)



# Following code here is to estimate the effect of a single value
# change the number in xind = c(i), we can find the different plot for each variable
par(mfrow=c(1,2))
out <- pdbart(x.train=Xind, y.train=wage$Y, xind=c(6)) 
# EMPFT  Effect of a single variable
plot(out2)
names(out2)
out2$xlbs
out2$levs
dim(out2$fd[[1]])

out3 <- pd2bart(x.train=Xind, y.train=wage$Y, xind=c(6,10),
                levquants=c(.05,.1,.25,.5,.75,.9,.95),pl=FALSE) # EMPFT and HRSOPEN
plot(out3)

####From below, we use the bartMachine to analysis BART regression model:

library(rJava)
library(bartMachine)

# BART regression model 
bart_machine <- bartMachine(new_wage, wage$Y)
summary(bart_machine)

rmse_by_num_trees(bart_machine, num_replicates = 20)

# BART regression model with cross-validation
bart_machine_cv <- bartMachineCV(new_wage, wage$Y)
summary(bart_machine_cv)

#### A representation of how well the model predicts the outcome ####
calc_credible_intervals(bart_machine_cv, new_wage[1:19, ])
# Actual vs. Fitted for bart_machine
plot_y_vs_yhat(bart_machine_cv,credible_intervals = TRUE)
plot_y_vs_yhat(bart_machine_cv,  prediction_intervals = TRUE)
plot(resid(bart_machine_cv))
abline(0,0)

#### A list of which elements of X are ‘important’ for the prediction ####
#change the number of trees from 200 to 50 
#1) ntree = 200
vc_200 = out$varcount
vpercent_200 = vc_200/apply(vc_200,1,sum)
vppostmean_200 = apply(vpercent_200,2,mean)
plot(vppostmean_200)

new.out <- bart(x.train=new_wage, y.train=wage$Y,ntree=50,ndpost = 1000)
vc_50 = new.out$varcount
vpercent_50 = vc_50/apply(vc_50,1,sum)
vppostmean_50 = apply(vpercent_50,2,mean)
plot(vppostmean_50)

vs <- var_selection_by_permute(bart_machine, bottom_margin = 5, num_permute_samples = 5)
vs$important_vars_local_names
vs$important_vars_global_se_names

investigate_var_importance(bart_machine_cv, num_replicates_for_avg = 10)
vs_cv <- var_selection_by_permute_cv(bart_machine, k_folds = 3, num_reps_for_avg = 5, num_permute_samples = 100, num_trees_for_permute = 20, alpha = 0.05, num_trees_pred_cv = 50)
print(vs_cv$best_method)
print(vs_cv$important_vars_cv)
