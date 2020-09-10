library(dplyr)

preproc.pacientes <- function(){
  pacientes <- readxl::read_xls(list.files(pattern = c(".xls")))
  
  colnames(pacientes)[1] <- "tel_p";
  colnames(pacientes)[14] <- "CCI"
  colnames(pacientes)[15] <- "CCI_cat"
  drops <- c("AP5")
  
  recod <- function(x) (max(x)+1)-x
  
  pacientes$AP1 <- recod(pacientes$AP1)
  pacientes$AP2 <- recod(pacientes$AP2)
  pacientes$AP3 <- recod(pacientes$AP3)
  
  pacientes <- pacientes %>% dplyr::select(-drops)
  pacientes <- pacientes %>% dplyr::select(-tel_p)
  pacientes$total_lcadl_baseline[116] <- mean(pacientes[[19]], na.rm = T)
  
  pacientes <- pacientes %>% 
    mutate(Sex = factor(Sex)) %>% 
    mutate(Sitl = factor(Sitl)) %>% 
    mutate(AP1 = factor(AP1)) %>% 
    mutate(AP2 = factor(AP2)) %>% 
    mutate(AP3 = factor(AP3)) %>% 
    mutate(AP4 = factor(AP4)) %>% 
    mutate(CCI_cat = factor(CCI_cat)) %>% 
    mutate(bmicat = factor(bmicat)) %>% 
    mutate(mortalidad_2years = factor(mortalidad_2years))
  
  return(pacientes)
}
y= pacientes$mortalidad_2years
pacientes <- preproc.pacientes()

for(i in 1:10){
  name <- paste0("csv",i,".csv")
}



fac.vars <- as.data.frame(pacientes[, sapply(pacientes, is.factor)])
pval <- c()
for(i in 1:8){
  tt <- table(fac.vars[,i], fac.vars[,9])
  if(any(tt==0)) tt = tt+0.5
  m1 <- chisq.test(tt)$expected
  G2.1<- 2*sum(tt*log(tt/m1))
  pval <-c(pval, 1-pchisq(G2.1,df=(2-1)*(2-1)) )
}

data.frame(variable =colnames(fac.vars)[-9],p_value=pval)


####################################### Logistic Regression ######################################################

full_model <- glm(mortalidad_2years ~ (AP1 + AP2 + AP3 + PF2 + PF4 + bmi + CCI + bode + w22_baseline)^2,
                  pacientes, family = "binomial")

# step(full_model, direction = "both", k=log(119)) #Step selection uncomment to perform

model.1 <- glm(formula = mortalidad_2years ~ AP1 + AP2 + AP3 + PF2 + PF4 + bmi + 
                 CCI + bode + w22_baseline + AP3:w22_baseline , 
               family = "binomial", data = pacientes)

model.2 <- glm(formula = mortalidad_2years ~ AP1+ AP2 + AP3 + bmi + PF4  + w22_baseline + AP3:w22_baseline , 
               family = "binomial", data = pacientes)

model.3 <- glm(formula = mortalidad_2years ~ AP1+ AP2 + AP3 + bmi + PF4, 
               family = "binomial", data = pacientes)

model.4 <- glm(formula = mortalidad_2years ~ AP2 + AP3 + bmi + PF4  + w22_baseline + AP3:w22_baseline , 
               family = "binomial", data = pacientes)

sumary <- summary(model.2)
sumary <- cbind(sumary$coefficients,exp = exp(sumary$coefficients[,1]))
sumary <- sumary[,c(1,5,2,3,4)]


library(boot)
glm.diag.plots(model.2)
outliers <- which(glm.diag(model.2)$cook>8/(119-2*7))



MKmisc::HLgof.test(fitted(model.2), as.numeric(as.character(y)))

library(pROC)
?plot.roc
plot(roc(fitted(model.2)>0.5, as.numeric(as.character(y))), reuse.auc=T)


library(caret)
confusionMatrix(table(as.numeric(fitted(model.2)>0.5), as.numeric(as.character(y))))
confusionMatrix(table(as.numeric(fitted(model.2)>0.5), as.numeric(as.character(y))))

library(Epi)


par(mfrow=c(1,2))
ROC(form = mortalidad_2years ~ AP1+ AP2 + AP3 + bmi + PF4  + w22_baseline + AP3:w22_baseline,
    data = pacientes, plot = "ROC", MI=T)
ROC(form = mortalidad_2years ~ AP1+ AP2 + AP3 + bmi + PF4  + w22_baseline + AP3:w22_baseline,
    data = pacientes[-outliers,], plot = "ROC", MI=T)


library(cvAUC)

set.seed(1)
datos <- pacientes
n <- nrow(datos)
d <- sample(1:n, size = n*0.6)
b=300
predictions <- matrix(NA, nrow = n-length(d), ncol = b)
observations <- predictions
i=1
for(i in 1:b){
  d <- sample(1:n, size = n*0.6)
  train <- datos[d,]
  test <- datos[-d,]
  
  model <- tryCatch(glm(formula = mortalidad_2years ~ AP1+ AP2 + AP3 + bmi + PF4  + w22_baseline + AP3:w22_baseline , 
                        family = "binomial", data = train),
                    error = function(e) break)
  predictions[,i] <- predict(model, test, type="response")
  observations[,i] <- test$mortalidad_2years
  i=i+1
}



out.1 <- cvAUC(predictions,
               observations)

plot(out.1$perf, col="grey82", lty=3)
plot(out.1$perf, col="red", avg="vertical", add=TRUE)
out.1$cvAUC

library(cvAUC)

set.seed(1)
datos <- pacientes[-outliers,]
n <- nrow(datos)
d <- sample(1:n, size = n*0.6)
b=300
predictions <- matrix(NA, nrow = n-length(d), ncol = b)
observations <- predictions
i=1
for(i in 1:b){
  d <- sample(1:n, size = n*0.6)
  train <- datos[d,]
  test <- datos[-d,]
  
  model <- glm(formula = mortalidad_2years ~ AP1+ AP2 + AP3 + bmi + PF4  + w22_baseline + AP3:w22_baseline , 
               family = "binomial", data = train)
  predictions[,i] <- tryCatch(predict(model, test, type="response"), error = function(e) rep(NA, length = n-length(d)))
  observations[,i] <- test$mortalidad_2years
}

observations <- observations[,colSums(is.na(predictions))==0]
predictions <- predictions[,colSums(is.na(predictions))==0]

out.2 <- cvAUC(predictions,
               observations)

plot(out.2$perf, col="grey82", lty=3)
plot(out.2$perf, col="red", avg="vertical", add=TRUE)

out.1$cvAUC
out.1$fold.AUC
out.2$cvAUC
out.2$fold.AUC

############################################# Poisson Regression #################################################

poisson.model <-glm(formula = prog_adm_COPD ~ bode , data = pacientes, offset = t_seguim)
poisson.model$fitted.values

plot(poisson.model)

plot(poisson.model)

summary(poisson.model)

anova(poisson.model, test = "Chisq")