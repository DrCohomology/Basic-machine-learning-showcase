library(readr)
library(party)
library(rpart)
library(rpart.plot)
library(magrittr)
library(glmnet)
library(boot)
library(PerformanceAnalytics)

# ADJUSTED rsq
rsq <- function(predictor, model, X, y){
	n <- nrow(X)
	p <- ncol(X)
	
	pred <- predictor(model, X)
	
	RSS <- sum((y-pred)^2)
	TSS <- sum((y-mean(y))^2)
	
	out <- 1 - RSS/(n-p-1) * (n-1)/TSS         
	
	return(out)
}

#---- Dataset, Multicollinearity ----

concrete_data <- read.csv("C:/Users/Io/Desktop/SML Supervised/Concrete_Data.csv", sep=';') #%>% standardize

X <- model.matrix(compressive_strength ~.-compressive_strength, data=concrete_data) %>% subset(select=c(-1))
y <- concrete_data$compressive_strength
Xpoly <- poly(X, 2, raw=TRUE) 

chart.Correlation(concrete_data, labels=names(concrete_data))	
boxplot(concrete_data)
#---- Model Predictivity ----

# Linear
linear_model <- glm(compressive_strength ~ .-compressive_strength, data=concrete_data)
summary(linear_model)
cv.glm(concrete_data, linear_model, K=10)$delta[1] 
# 109.3
rsq(predict, linear_model, data.frame(X), y)
# 0.612

# Lasso (has some stochastic component in the design)
test_lasso <- glmnet(X, y)
plot(test_lasso, xvar='lambda', label=TRUE)

cv_lasso <- cv.glmnet(X, y)
plot(cv_lasso)
min(cv_lasso$cvm) 
# 109.7

cv_lasso$lambda.1se
# 0.560

lasso_model <- glmnet(X, y, lambda=cv_lasso$lambda.1se)
coef(lasso_model)

rsq(predict, lasso_model, X, y)
# 0.598


# naive error for poly and tree - NOT USED 
train <- sample(nrow(concrete_data),nrow(concrete_data)*0.75,replace=FALSE)
poly_train <- glm(concrete_data$compressive_strength[train] ~ Xpoly[train,])

sum((y[-train]-predict(poly_train, Xpoly[-train,] %>% data.frame))^2)/length(y[-train])
#...

# Polynomial
poly_model <- glm(compressive_strength ~ Xpoly, data=concrete_data)
summary(poly_model)
cv.glm(data.frame(Xpoly, y), poly_model, K=10)$delta[1] 
#524???
sum((y-predict(poly_model, Xpoly))^2)/length(y)
#52.8

rsq(predict.glm, poly_model, Xpoly%>%data.frame, y)
# 0.802

# Tree
complete_tree <- rpart(compressive_strength ~ .-compressive_strength, data=concrete_data)
plotcp(complete_tree)
prp(complete_tree, roundint=FALSE, digits=4)

rsq(rpart.predict, complete_tree, data.frame(X), y)
# 0.738
sum((y-rpart.predict(complete_tree, data.frame(X)))^2)/length(y)
# 72.4

#---- Case Study ----
# Tree
# Let's assume we have at most 200 kg of cement per m^3 AND only 60 days to make it stiff
# What is the best combination of other materials?

target <- 200
time <- 60
case_data <- concrete_data %>% subset(cement < target) %>% subset(age < time)

boxplot(case_data)

Xc <- model.matrix(compressive_strength ~.-compressive_strength, data=case_data) 
yc <- case_data$compressive_strength

# does anything change in the coefficients? YES it does, a lot
linear_model <- glm(compressive_strength ~ .-compressive_strength, data=case_data)

rsq(predict, linear_model, data.frame(Xc), yc)
# 0.752

#---- CASE study tree ----
max(case_data$compressive_strength)

case_tree <- rpart(compressive_strength ~ .-compressive_strength, data=case_data)
plotcp(case_tree)
prp(case_tree, roundint=FALSE, digits=4)

rsq(rpart.predict, case_tree, data.frame(Xc), yc)
# 0.807

# best value
case_data[case_data$compressive_strength == max(case_data$compressive_strength), ]
