load(file = "./dataset/COVID/WIMR_decon_cohort.rds")
load(file = "./dataset/COVID/WIMR_decon2_cohort.rds")

health <- wimr[grep("health", wimr$V3),]
healthy.avg <- apply(health[,342:ncol(health)], 2, median)
healthy.max <- apply(health[,342:ncol(health)], 2, max)
healthy.min <- apply(health[,342:ncol(health)], 2, min)
wimr[grep("Yes", wimr$bacteria_2nd),"Neutrophils"]
wimr[grep("Yes", wimr$bacteria_2nd),"Mono_C"]
wimr <- wimr[grep("health", wimr$V3, invert = T),]
library(randomForest)
library(caret)
library(e1071)
library(sigFeature)

##SVM####
wimr$diabetes
wimr$bacteria_2nd
wimr.df <- wimr[,grep("age|sex|bacteria_2nd|WHO_severity", colnames(wimr))]
wimr.df$bacteria_2nd <- as.factor(wimr.df$bacteria_2nd)
wimr.df <- cbind(wimr.df,  wimr[,342:ncol(wimr)])
wimr.df <- wimr.df[complete.cases(wimr.df), ]
wimr.df <- wimr.df[,-5]
smp_size <- floor(0.75 * nrow(wimr.df))
set.seed(123)
train_idx <- sample(seq_len(nrow(wimr.df)), size = smp_size)
train.df <- wimr.df[train_idx, ]
test.df <-  wimr.df[-train_idx, ]

x <- train.df[,-4]
x <- x[,-1]
y <- rep(1, nrow(x))
y[grep("No", train.df$bacteria_2nd)] = -1
pvals <- sigFeaturePvalue(x,y)
hist(unlist(pvals),col="skyblue", #breaks=seq(0,0.08,0.0015),
     xlab="p value",ylab="Frequency",main="")
system.time(sigFeature(x, y))
sigfeatureRankedList <- sigFeature(x[,1:30], y)

sigfeatureRankedList<- colnames(x)[sigfeatureRankedList]

#clinical data
set.seed(304)
rf <- svm(bacteria_2nd ~. , data=train.df[,1:4],  type = "C", kernel = "linear", probability=TRUE) 
print(rf)
summary(rf)
pred <- predict(rf , newdata = test.df[,1:3], probability = T)
table(pred, test.df[,4])



#cell type data
set.seed(304)
rf <- svm(bacteria_2nd ~. , data=train.df[,4:ncol(train.df)],  type = "C", kernel = "linear", probability=TRUE) 
print(rf)
summary(rf)
pred <- predict(rf , newdata = test.df[,5:ncol(train.df)], probability = T)
table(pred, test.df[,4])



#all data
set.seed(304)
rf <- svm(bacteria_2nd ~. , data=train.df[,c(sigfeatureRankedList[1:19],"bacteria_2nd")],  type = "C", kernel = "linear", probability=TRUE) 
print(rf)
summary(rf)
pred <- predict(rf , newdata = test.df[,c(1:3,5:ncol(train.df))], probability = T)
table(pred, test.df[,4])
(table(pred, test.df[,4])[1,1] + table(pred, test.df[,4])[2,2])/nrow(test.df)

accuracy <- c()
for (i in 2:length(sigfeatureRankedList)) {
  set.seed(304)
  rf <- svm(bacteria_2nd ~. , data=train.df[,c(sigfeatureRankedList[1:i],"bacteria_2nd")],  type = "C", kernel = "linear", probability=TRUE) 
  print(rf)
  summary(rf)
  pred <- predict(rf , newdata = test.df[,c(1:3,5:ncol(train.df))], probability = T)
  table(pred, test.df[,4])
  accuracy <- c(accuracy , (table(pred, test.df[,4])[1,1] + table(pred, test.df[,4])[2,2])/nrow(test.df))
}
accuracy
#Accuracy is 81% , not good as the published one?

set.seed(304)
rf <- svm(bacteria_2nd ~ Neutrophils + Mono_C+ WHO_severity+ Bex + CD8Tcm + age, data=train.df,  type = "C", kernel = "linear", probability=TRUE) 
print(rf)
summary(rf)
pred <- predict(rf , newdata = test.df, probability = T)
table(pred, test.df[,4])





##Lasso ####
library(dplyr)
library(tidyr)
wimr.df <- wimr[,grep("age|sex|bacteria_2nd|WHO_severity", colnames(wimr))]
wimr.df$bacteria_2nd <- as.factor(wimr.df$bacteria_2nd)
wimr.df <- cbind(wimr.df,  wimr[,342:ncol(wimr)])
wimr.df <- wimr.df[complete.cases(wimr.df), ]
wimr.df <- wimr.df[,-5]


set.seed(1)
train = wimr.df %>%
  sample_frac(0.8)
test = wimr.df %>%
  setdiff(train)
xfactors = model.matrix(bacteria_2nd ~ sex + WHO_severity , data = train)[, -1]
x_train  <- as.matrix(cbind(train[,c(2,5:ncol(train))], xfactors))
xfactors = model.matrix(bacteria_2nd ~ sex + WHO_severity , data = test)[, -1]
x_test  <- as.matrix(cbind(test[,c(2,5:ncol(test))], xfactors))
y_train = train %>%
  select(bacteria_2nd) %>%
  unlist() 
y_test = test %>%
  select(bacteria_2nd) %>%
  unlist() 

library(glmnet)
grid = 10^seq(10, -2, length = 100)
set.seed(2)
lasso_mod = glmnet(x_train, 
                   y_train, 
                   #lambda = grid,
                   family = "binomial",
                   alpha = 1 
                   ) # Fit lasso model on training data
plot(lasso_mod, xvar = "lambda", label = T)    # Draw plot of coefficients
print(lasso_mod)
coef(lasso_mod)[, 10]
#remove the variables with coef of 0 
lasso_mod = glmnet(x_train[,grep("mature|WHO_severity", colnames(x_train))], 
                   y_train, 
                   #lambda = grid,
                   family = "binomial",
                   alpha = 1 
) # Fit lasso model on training data
plot(lasso_mod, xvar = "lambda", label = T)    # Draw plot of coefficients
print(lasso_mod)
coef(lasso_mod)[, 10]

set.seed(1)
cv.out = cv.glmnet(x_train, y_train, alpha = 1, family = "binomial") # Fit lasso model on training data
plot(cv.out) # Draw plot of training MSE as a function of lambda
bestlam = cv.out$lambda.min # Select lamda that minimizes training MSE

lasso_pred = predict(lasso_mod, s = bestlam, newx = x_test[,grep("mature|WHO_severity", colnames(x_train))], type = "class") # Use best lambda to predict test data
#mean((lasso_pred - y_test)^2) # Calculate test MSE
table(lasso_pred, y_test)
(table(lasso_pred, y_test)[1,1] + table(lasso_pred, y_test)[2,2])/nrow(x_test)


lasso_train <- predict(lasso_mod, s = bestlam, newx = x_train[,grep("mature|WHO_severity", colnames(x_train))], type = "class") # Use best lambda to predict test data
table(lasso_train, y_train)
(table(lasso_train, y_train)[1,1] + table(lasso_train, y_train)[2,2])/nrow(x_train)
