######################### Figure 3.1 ##############################################
library(fda); library(MASS)

set.seed(45)
x <- c(10+cumsum(rnorm(150)), rep(NA,150), (c(20+cumsum(rnorm(150)))), rep(NA,150))

argvals = c(1:150, c(301:450))
breaks= c(seq(0,150, length.out = 10), seq(300,600,length.out = 20))
basis <- create.bspline.basis(rangeval = c(0,600),breaks = breaks, norder = 3)


phi <- eval.basis(na.omit(argvals),basis)
plot(phi[,1], type = "l", xlim = c(0,450),ylim = c(0,1))
for(i in 2:40) lines(phi[,i])
phi.2 <- eval.basis(1:600, basis)

phi.2 <- eval.basis(1:600, basis)

plot(x, ylim = c(-12,30), xlab = "")
lines(phi.2%*%ginv(t(phi)%*%phi)%*%t(phi)%*%na.omit(x), col = "blue", lwd=2)
abline(v= breaks, col ="grey", lty=3)

########################  Figure 3.2. ####################################

library(data.table);library(dplyr);library(plyr);library(stringr);library(fda)
archivos <- list.files(pattern = ".csv")
myfiles <- lapply(archivos, fread) # fread siempre es más rápido
pacientes <- readxl::read_xls(list.files(pattern = c(".xls")))
colnames(pacientes)[1] <- "tel_p";
podometro <- myfiles[[1]];podometro <- podometro %>% mutate(valor = as.numeric(as.character(valor)))
podometro <- na.omit(podometro)
podometro$fecha <- as.Date(podometro$fecha, format = "%d/%m/%Y")
podometro <- podometro[podometro$fecha > "2009-01-01",] ## 61
fechas <- podometro$fecha
podometro <- podometro[podometro$valor<=40000,]
podometro <- podometro %>% 
  dplyr::group_by(tel_p, fecha) %>% 
  dplyr::summarise(valor = sum(valor))
y <- pacientes$mortalidad_2years[pacientes$tel_p %in% unique(podometro$tel_p)]
podometro <- as.data.frame(podometro)

splited_podometro <- split(podometro, podometro$tel_p)
splited_podometro$`41` <- splited_podometro$`41`[-1,] # 2nd value three months later!
splited_podometro$`46` <- splited_podometro$`46`[-1,] # 2nd value one months later!
splited_podometro$`109` <- splited_podometro$`109`[-c(1:7),]
splited_podometro$`110` <- splited_podometro$`110`[-c(1:5),]
splited_podometro$`70` <- splited_podometro$`70`[-c(1:2),]

unique(podometro$tel_p)[64]

a <- as.numeric(min(podometro$fecha));
b <- as.numeric(max(podometro$fecha));
A <- podometro %>% dplyr::select(fecha, valor, tel_p)
splited <- split(A, podometro$tel_p)

splited$`41` <- splited$`41`[-1,] # 2nd value three months later!
splited$`46` <- splited$`46`[-1,] # 2nd value one months later!
splited$`109` <- splited$`109`[-c(1:7),]
splited$`110` <- splited$`110`[-c(1:5),]
splited$`70` <- splited$`70`[-c(1:2),]


splited <- lapply(splited, "[", -c(3))

rm(splited_podometro,A)

par(mfrow=c(1,2))

plot(rep(-6,length(splited$`1`$fecha))~splited$`1`$fecha, type="l", ylim= c(-6,5),
     ylab = "", xlab= "Date", yaxt="n")
pos.y <- -6
for(i in 2:112){
  pos.y <- pos.y +0.1
  lines(rep(pos.y,length(splited[[i]]$fecha))~splited[[i]]$fecha)
}

for(i in 1:112){ splited[[i]]$fecha <- splited[[i]]$fecha - min(splited[[i]]$fecha)}


max.fecha <- sapply(splited, function(x) max(x$fecha))
order(max.fecha)

lfechas <- lapply(splited, function(x) x = zoo::as.Date.numeric(x$fecha))

plot(rep(-6,length(splited$`1`$fecha))~splited$`1`$fecha, type="l", ylim= c(-6,5),
     ylab = "", xlab= expression("t"["j"]), yaxt="n")
pos.y = -6.1
for(i in order(max.fecha, decreasing = T)){
  pos.y <- pos.y +0.1
  lines(rep(pos.y,length(splited[[i]]$fecha))~splited[[i]]$fecha)
}


########################  Figure 3.3. ####################################

M <- fda::sparse.mat(splited)

## Outlier detection was performed looking graphically at all curves and using Boxplots for all curves.
# Many are obvious problems of data imputation (errors of the monitoring)

M[which(M[,4]>9000),4] <- 317
M[which(M[,6]>5000),6] <- c(596, 710, 978) ## Curve 6 change
M[which(M[,7]>30000),7] <- 3010 ## Curve 7 change
M[which(M[,12] >15000),12] <- c(2191, 2183, 2150) ## Curve 12
M[which(M[,15] > 25000),15] <- 8935 ## Curve 15
M[which(M[,16]>25000),16] <- 9512 ## Curve 16
M[which(M[,19]>15000),19] <- mean(M[,19], na.rm=T) ## Curve 19
M[which(M[,20] < 1000),20] <- mean(M[,20], na.rm=T)
M[which(M[,23] <1000)[-c(1,2)],23] <- mean(M[,23], na.rm = T)
M[which(M[,27]>15000),27] <- c(5208,3545,2692)
M[which(M[,33]>15000),33] <- 5303
M[which(M[,34]>30000),34] <- 3521
M[which(M[,65]>10000),65] <- mean(M[,65], na.rm = T)
M[which(M[,67]>9000),67] <- mean(M[,67], na.rm = T)
M[which(M[,70]>6000),70] <- mean(M[,70], na.rm = T)
M[which(M[,85]==540),85] <- 5400
M[which(M[,90] > 25000),90] <- 2754
M[which(M[,94] > 30000),94] <- 7008
M[which(M[,96]>40000),96] <- 8793
M[which(M[,110]>40000),110] <- mean(M[,110], na.rm=T)
M[1,112] <- 3480

y <- pacientes$mortalidad_2years[pacientes$tel_p %in% unique(podometro$tel_p)]


par(mfrow=c(1,2))
plot(M[,51], type = "l",xlab= expression("T"["j"]), ylab = "Number of steps", main = "Curve 51")
plot(M[,6], type = "l", main = "Curve 6", ylab = "", xlab= expression("T"["j"]))


######################### Figure 1.4, 1.5 and 1.6 ############################
############################### PEDOMETER ####################################
medias <- apply(M, 1, function(x) mean(x, na.rm=T))
medias.1 <- apply(M[,which(y==1)], 1, function(x) mean(x, na.rm=T))
medias.1 <- na.omit(cbind(1:1289, medias.1))
medias.0 <- apply(M[,which(y!=1)], 1, function(x) mean(x, na.rm=T))
sddata <- apply(M, 1, function(x) sd(x, na.rm=T))
sddata.1 <- apply(M[,which(y==1)], 1, function(x) sd(x, na.rm=T))
sddata.1 <- na.omit(cbind(1:1289, sddata.1))
sddata.0 <- apply(M[,which(y!=1)], 1, function(x) sd(x, na.rm=T))

par(mfrow=c(1,2))
plot(medias, type="l", ylim = c(0,9000), col = rgb(0, 0,0, alpha = 0.5), xlab= expression("t"["j"]),
     ylab = bquote(bar(x)(t)))
lines(smooth.basis(argvals = 1:length(medias),y= medias, 
                   create.bspline.basis(rangeval = c(0,length(medias)), nbasis = 10, norder = 3)
),
lwd = 2)
plot(sddata, type="l", ylim = c(0,9000), col = rgb(0, 0,0, alpha = 0.5), xlab= expression("t"["j"]),
     ylab = bquote(hat(sigma)(t)))
lines(smooth.basis(argvals = 1:length(sddata),y= sddata, 
                   create.bspline.basis(rangeval = c(0,length(sddata)), nbasis = 10, norder = 3)
),
lwd = 2)


par(mfrow=c(1,2))
plot(medias.1[,2], type="l",ylim=c(0,9000) ,xlim = c(1,1289),col = rgb(1, 0,0, alpha = 0.5), xlab =expression("t"["j"]), ylab = bquote(bar(x)(t)))
lines(smooth.basis(argvals = medias.1[,1],y= medias.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(medias.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(medias.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(medias.0),y= medias.0, 
                   create.bspline.basis(rangeval = c(0,length(medias.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### standard deviation



plot(sddata.1[,2], type="l", col = rgb(1, 0,0, alpha = 0.5), xlim = c(1,1289),xlab =expression("t"["j"]), ylab = bquote(hat(sigma)(t)))
lines(smooth.basis(argvals = sddata.1[,1],y= sddata.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(sddata.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(sddata.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(sddata.0),y= sddata.0, 
                   create.bspline.basis(rangeval = c(0,length(sddata.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### Covariance function

par(mfrow=c(1,1))

eval.grid <- as.integer(seq(1,1289,length.out = 100))
grid <- expand.grid(eval.grid, eval.grid)
expand.grid(eval.grid)
grid <-cbind(grid,Var3=NA)

M.2 <- M - rowMeans(M, na.rm = T)

for(i in 1:nrow(grid)){
  grid[i,3] <- sum(M.2[grid[i,1],]*M.2[grid[i,2],], na.rm = T)/112
  
}


grid
grid <- xtabs(Var3 ~ Var1+Var2, data = grid)

grid
image(grid, axes =F)
contour(grid, add=T, axes=T)
axis(side = 1, at = seq(0,1,length.out = 5), labels = seq(1,1289,length.out = 5))
axis(side = 2, at = seq(0,1,length.out = 5), labels = seq(1,1289,length.out = 5))

######################### Figure 1.7, 1.8 and 1.9 ############################
############################### SPO2 #########################################

SPO2 <- myfiles[[10]]

SPO2 <- SPO2 %>% mutate(valor = as.numeric(as.character(valor)))
SPO2 <- na.omit(SPO2)


SPO2$fecha <- as.Date(SPO2$fecha, format = "%d/%m/%Y")


SPO2 <- SPO2[SPO2$fecha > "2009-01-01",] ## 61


SPO2 <- SPO2 %>% 
  dplyr::group_by(tel_p, fecha) %>% 
  dplyr::summarise(valor = mean(valor)) %>% 
  dplyr::mutate(fecha = as.numeric(fecha))


SPO2 <- as.data.frame(SPO2)
SPO2$valor[10689] <- 97

splited_SPO2 <- split(SPO2, SPO2$tel_p)

a <- as.numeric(min(SPO2$fecha));
b <- as.numeric(max(SPO2$fecha));
A <- SPO2 %>% dplyr::select(fecha, valor, tel_p)
splited <- split(A, SPO2$tel_p)
splited <- lapply(splited, "[", -c(3))
for(i in 1:112){ splited[[i]]$fecha <- splited[[i]]$fecha - min(splited[[i]]$fecha) +a}


M.spo2 <- fda::sparse.mat(splited)
M <- M.spo2
M
plot(apply(M, 1, function(x) mean(x, na.rm=T)), type="l")
medias <- apply(M, 1, function(x) mean(x, na.rm=T))
medias.1 <- apply(M[,which(y==1)], 1, function(x) mean(x, na.rm=T))
medias.1 <- na.omit(cbind(1:1289, medias.1))
medias.0 <- apply(M[,which(y!=1)], 1, function(x) mean(x, na.rm=T))


par(mfrow=c(1,2))
plot(medias, type="l",ylim =c(90,98), col = rgb(0, 0,0, alpha = 0.5), xlab =expression("t"["j"]), ylab = bquote(bar(x)(t)))
lines(smooth.basis(argvals = 1:length(medias),y= medias, 
                   create.bspline.basis(rangeval = c(0,length(medias)), nbasis = 10, norder = 3)
),
lwd = 2)

lines(medias.1[,2], type="l", col = rgb(1, 0,0, alpha = 0.5))
lines(smooth.basis(argvals = medias.1[,1],y= medias.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(medias.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(medias.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(medias.0),y= medias.0, 
                   create.bspline.basis(rangeval = c(0,length(medias.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### standard deviation


sddata <- apply(M, 1, function(x) sd(x, na.rm=T))
sddata.1 <- apply(M[,which(y==1)], 1, function(x) sd(x, na.rm=T))
sddata.1 <- na.omit(cbind(1:1289, sddata.1))
sddata.0 <- apply(M[,which(y!=1)], 1, function(x) sd(x, na.rm=T))



plot(sddata, type="l", col = rgb(0, 0,0, alpha = 0.5),xlab =expression("t"["j"]), ylab = bquote(hat(sigma)(t)))
lines(smooth.basis(argvals = 1:length(sddata),y= sddata, 
                   create.bspline.basis(rangeval = c(0,length(sddata)), nbasis = 10, norder = 3)
),
lwd = 2)

lines(sddata.1[,2], type="l", col = rgb(1, 0,0, alpha = 0.5))
lines(smooth.basis(argvals = sddata.1[,1],y= sddata.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(sddata.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(sddata.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(sddata.0),y= sddata.0, 
                   create.bspline.basis(rangeval = c(0,length(sddata.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### Covariance function

eval.grid <- as.integer(seq(1,1289,length.out = 100))
grid <- expand.grid(eval.grid, eval.grid)
grid <-cbind(grid,Var3 =NA)

M.2 <- M
M.2 <- M.2 - rowMeans(M, na.rm = T)
rowMeans(M.2, na.rm = T)
grid

for(i in 1:nrow(grid)){
  grid[i,3] <- sum(M.2[grid[i,1],]*M.2[grid[i,2],], na.rm = T)/112
  
}
grid
grid <- xtabs(Var3 ~ Var1+Var2, data = grid)


persp(grid, theta = 20, phi = 20, col="deepskyblue", box = F)
image(grid)
contour(grid, add=T)


######################### Figure 1.7, 1.8 and 1.9 ############################
##############################################################################

myfiles[[9]]
RitmoResp <- myfiles[[9]]

RitmoResp <- RitmoResp %>% mutate(valor = as.numeric(as.character(valor)))
RitmoResp <- na.omit(RitmoResp)

sum(is.na(RitmoResp$valor))

RitmoResp$fecha <- as.Date(RitmoResp$fecha, format = "%d/%m/%Y")


RitmoResp <- RitmoResp[RitmoResp$fecha > "2009-01-01",] ## 61


RitmoResp <- RitmoResp %>% 
  dplyr::group_by(tel_p, fecha) %>% 
  dplyr::summarise(valor = mean(valor)) %>% 
  dplyr::mutate(fecha = as.numeric(fecha))


RitmoResp <- as.data.frame(RitmoResp)


plot(RitmoResp$valor)

splited_RitmoResp <- split(RitmoResp, RitmoResp$tel_p)

a <- as.numeric(min(RitmoResp$fecha));
b <- as.numeric(max(RitmoResp$fecha));
A <- RitmoResp %>% dplyr::select(fecha, valor, tel_p)
splited <- split(A, RitmoResp$tel_p)
splited <- lapply(splited, "[", -c(3))
for(i in 1:112){ splited[[i]]$fecha <- splited[[i]]$fecha - min(splited[[i]]$fecha) +a}


M.RitmoResp <- fda::sparse.mat(splited)
M <- M.RitmoResp
M
plot(apply(M, 1, function(x) mean(x, na.rm=T)), type="l")
medias <- apply(M, 1, function(x) mean(x, na.rm=T))
medias.1 <- apply(M[,which(y==1)], 1, function(x) mean(x, na.rm=T))
medias.1 <- na.omit(cbind(1:1289, medias.1))
medias.1
medias.0 <- apply(M[,which(y!=1)], 1, function(x) mean(x, na.rm=T))


par(mfrow=c(1,2))
plot(medias, type="l", col = rgb(0, 0,0, alpha = 0.5),ylim =c(14,26), xlab =expression("t"["j"]), ylab = bquote(bar(x)(t)))
lines(smooth.basis(argvals = 1:length(medias),y= medias, 
                   create.bspline.basis(rangeval = c(0,length(medias)), nbasis = 10, norder = 3)
),
lwd = 2)

lines(medias.1[,2], type="l", col = rgb(1, 0,0, alpha = 0.5))
lines(smooth.basis(argvals = medias.1[,1],y= medias.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(medias.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(medias.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(medias.0),y= medias.0, 
                   create.bspline.basis(rangeval = c(0,length(medias.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### standard deviation


sddata <- apply(M, 1, function(x) sd(x, na.rm=T))
sddata.1 <- apply(M[,which(y==1)], 1, function(x) sd(x, na.rm=T))
sddata.1 <- na.omit(cbind(1:1289, sddata.1))
sddata.0 <- apply(M[,which(y!=1)], 1, function(x) sd(x, na.rm=T))



plot(sddata, type="l", col = rgb(0, 0,0, alpha = 0.5), xlab =expression("t"["j"]), ylab = bquote(hat(sigma)(t)))
lines(smooth.basis(argvals = 1:length(sddata),y= sddata, 
                   create.bspline.basis(rangeval = c(0,length(sddata)), nbasis = 10, norder = 3)
),
lwd = 2)

lines(sddata.1[,2], type="l", col = rgb(1, 0,0, alpha = 0.5))
lines(smooth.basis(argvals = sddata.1[,1],y= sddata.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(sddata.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(sddata.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(sddata.0),y= sddata.0, 
                   create.bspline.basis(rangeval = c(0,length(sddata.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### Covariance function

eval.grid <- as.integer(seq(1,1289,length.out = 100))
grid <- expand.grid(eval.grid, eval.grid)
grid <-cbind(grid,Var3 =NA)

M.2 <- M
M.2 <- M.2 - rowMeans(M, na.rm = T)
rowMeans(M.2, na.rm = T)
grid

for(i in 1:nrow(grid)){
  grid[i,3] <- sum(M.2[grid[i,1],]*M.2[grid[i,2],], na.rm = T)/112
  
}
grid
grid <- xtabs(Var3 ~ Var1+Var2, data = grid)


persp(grid, theta = 20, phi = 20, col="deepskyblue", box = F)
image(grid)
contour(grid, add=T)

###################### Pulse #####################################
#####################################################################

myfiles[[8]]
Pulso <- myfiles[[8]]

Pulso <- Pulso %>% mutate(valor = as.numeric(as.character(valor)))
Pulso <- na.omit(Pulso)

sum(is.na(Pulso$valor))

Pulso$fecha <- as.Date(Pulso$fecha, format = "%d/%m/%Y")


Pulso <- Pulso[Pulso$fecha > "2009-01-01",] ## 61


Pulso <- Pulso %>% 
  dplyr::group_by(tel_p, fecha) %>% 
  dplyr::summarise(valor = mean(valor)) %>% 
  dplyr::mutate(fecha = as.numeric(fecha))


Pulso <- as.data.frame(Pulso)


splited_Pulso <- split(Pulso, Pulso$tel_p)

a <- as.numeric(min(Pulso$fecha));
b <- as.numeric(max(Pulso$fecha));
A <- Pulso %>% dplyr::select(fecha, valor, tel_p)
splited <- split(A, Pulso$tel_p)
splited <- lapply(splited, "[", -c(3))
for(i in 1:112){ splited[[i]]$fecha <- splited[[i]]$fecha - min(splited[[i]]$fecha) +a}


M.Pulso <- fda::sparse.mat(splited)
M <- M.Pulso
M
medias <- apply(M, 1, function(x) mean(x, na.rm=T))
medias.1 <- apply(M[,which(y==1)], 1, function(x) mean(x, na.rm=T))
medias.1 <- na.omit(cbind(1:1289, medias.1))
medias.0 <- apply(M[,which(y!=1)], 1, function(x) mean(x, na.rm=T))



plot(medias, type="l", col = rgb(0, 0,0, alpha = 0.5),xlab =expression("t"["j"]), ylab = bquote(bar(x)(t)))
lines(smooth.basis(argvals = 1:length(medias),y= medias, 
                   create.bspline.basis(rangeval = c(0,length(medias)), nbasis = 10, norder = 3)
),
lwd = 2)

lines(medias.1[,2], type="l", col = rgb(1, 0,0, alpha = 0.5))
lines(smooth.basis(argvals = medias.1[,1],y= medias.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(medias.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(medias.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(medias.0),y= medias.0, 
                   create.bspline.basis(rangeval = c(0,length(medias.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### standard deviation


sddata <- apply(M, 1, function(x) sd(x, na.rm=T))
sddata.1 <- apply(M[,which(y==1)], 1, function(x) sd(x, na.rm=T))
sddata.1 <- na.omit(cbind(1:1289, sddata.1))
sddata.0 <- apply(M[,which(y!=1)], 1, function(x) sd(x, na.rm=T))



plot(sddata, type="l", col = rgb(0, 0,0, alpha = 0.5), xlab =expression("t"["j"]), ylab = bquote(hat(sigma)(t)))
lines(smooth.basis(argvals = 1:length(sddata),y= sddata, 
                   create.bspline.basis(rangeval = c(0,length(sddata)), nbasis = 10, norder = 3)
),
lwd = 2)

lines(sddata.1[,2], type="l", col = rgb(1, 0,0, alpha = 0.5))
lines(smooth.basis(argvals = sddata.1[,1],y= sddata.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(sddata.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(sddata.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(sddata.0),y= sddata.0, 
                   create.bspline.basis(rangeval = c(0,length(sddata.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### Covariance function

eval.grid <- as.integer(seq(1,1289,length.out = 100))
grid <- expand.grid(eval.grid, eval.grid)
grid <-cbind(grid,Var3 =NA)

M.2 <- M
M.2 <- M.2 - rowMeans(M, na.rm = T)
rowMeans(M.2, na.rm = T)
grid

for(i in 1:nrow(grid)){
  grid[i,3] <- sum(M.2[grid[i,1],]*M.2[grid[i,2],], na.rm = T)/112
  
}
grid
grid <- xtabs(Var3 ~ Var1+Var2, data = grid)


persp(grid, theta = 20, phi = 20, col="deepskyblue", box = F)
image(grid)
contour(grid, add=T)

############################################## Temperature ################################################
########################################################################################################

Temp <- myfiles[[11]]

Temp <- Temp %>% mutate(valor = as.numeric(as.character(valor)))
Temp <- na.omit(Temp)

sum(is.na(Temp$valor))

Temp$fecha <- as.Date(Temp$fecha, format = "%d/%m/%Y")


Temp <- Temp[Temp$fecha > "2009-01-01",] ## 61


Temp <- Temp %>% 
  dplyr::group_by(tel_p, fecha) %>% 
  dplyr::summarise(valor = mean(valor)) %>% 
  dplyr::mutate(fecha = as.numeric(fecha))


Temp <- as.data.frame(Temp)

plot(Temp$valor)


splited_Temp <- split(Temp, Temp$tel_p)

a <- as.numeric(min(Temp$fecha));
b <- as.numeric(max(Temp$fecha));
A <- Temp %>% dplyr::select(fecha, valor, tel_p)
splited <- split(A, Temp$tel_p)
splited <- lapply(splited, "[", -c(3))
for(i in 1:112){ splited[[i]]$fecha <- splited[[i]]$fecha - min(splited[[i]]$fecha) +a}


M.Temp <- fda::sparse.mat(splited)
M <- M.Temp
M
medias <- apply(M, 1, function(x) mean(x, na.rm=T))
medias.1 <- apply(M[,which(y==1)], 1, function(x) mean(x, na.rm=T))
medias.1 <- na.omit(cbind(1:1289, medias.1))
medias.0 <- apply(M[,which(y!=1)], 1, function(x) mean(x, na.rm=T))



plot(medias, type="l", col = rgb(0, 0,0, alpha = 0.5),
     ylim=c(35,36.7),xlab =expression("t"["j"]), ylab = bquote(bar(x)(t)))
lines(smooth.basis(argvals = 1:length(medias),y= medias, 
                   create.bspline.basis(rangeval = c(0,length(medias)), nbasis = 10, norder = 3)
),
lwd = 2)

lines(medias.1[,2], type="l", col = rgb(1, 0,0, alpha = 0.5))
lines(smooth.basis(argvals = medias.1[,1],y= medias.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(medias.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(medias.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(medias.0),y= medias.0, 
                   create.bspline.basis(rangeval = c(0,length(medias.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### standard deviation


sddata <- apply(M, 1, function(x) sd(x, na.rm=T))
sddata.1 <- apply(M[,which(y==1)], 1, function(x) sd(x, na.rm=T))
sddata.1 <- na.omit(cbind(1:1289, sddata.1))
sddata.0 <- apply(M[,which(y!=1)], 1, function(x) sd(x, na.rm=T))



plot(sddata, type="l", col = rgb(0, 0,0, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(sddata),y= sddata, 
                   create.bspline.basis(rangeval = c(0,length(sddata)), nbasis = 10, norder = 3)
),
lwd = 2)

lines(sddata.1[,2], type="l", col = rgb(1, 0,0, alpha = 0.5))
lines(smooth.basis(argvals = sddata.1[,1],y= sddata.1[,2], 
                   create.bspline.basis(rangeval = c(0,max(sddata.1[,1])), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0.85, 0,0, alpha = 1))

lines(sddata.0, type="l", col = rgb(0, 0,1, alpha = 0.5))
lines(smooth.basis(argvals = 1:length(sddata.0),y= sddata.0, 
                   create.bspline.basis(rangeval = c(0,length(sddata.0)), nbasis = 10, norder = 3)
),
lwd = 2, col = rgb(0, 0,0.85, alpha = 1))


### Covariance function

eval.grid <- as.integer(seq(1,1289,length.out = 100))
grid <- expand.grid(eval.grid, eval.grid)
grid <-cbind(grid,Var3 =NA)

M.2 <- M
M.2 <- M.2 - rowMeans(M, na.rm = T)
rowMeans(M.2, na.rm = T)
grid

for(i in 1:nrow(grid)){
  grid[i,3] <- sum(M.2[grid[i,1],]*M.2[grid[i,2],], na.rm = T)/112
  
}
grid
grid <- xtabs(Var3 ~ Var1+Var2, data = grid)


persp(grid, theta = 20, phi = 20, col="deepskyblue", box = F)
image(grid)
contour(grid, add=T)


###################### Figure 3.11 ###################################################


max.fecha <- sapply(splited, function(x) max(x$fecha))

par(mfrow=c(1,1))
plot(rep(-6,length(splited$`1`$fecha))~splited$`1`$fecha, type="l", ylim= c(-6,5),
     ylab = "", xlab= expression("t"["j"]), yaxt="n", xaxt="n")
axis(side = 1, at = c(0, 72, 213, 365, 547, 800, 1000, 1200), labels = c(0, 72, 213, 365, 547, 800, 1000, 1200))
pos.y = -6.1
for(i in order(max.fecha, decreasing = T)){
  pos.y <- pos.y +0.1
  lines(rep(pos.y,length(splited[[i]]$fecha))~splited[[i]]$fecha, col = (y+1)[i])
}

y.pos <-seq(5.1, by = -0.1, length.out = 112)
abline(v = 72)
abline(v = 213)
abline(h=y.pos[which(max.fecha[order(max.fecha, decreasing = F)] == 213)], col = rgb(0.5,0.5,0.5, alpha = 0.5))
abline(v = 365)
abline(h=y.pos[which(max.fecha[order(max.fecha, decreasing = F)] == 382)],col = rgb(0.5,0.5,0.5, alpha = 0.5))
abline(v=547)
abline(h = y.pos[which(max.fecha[order(max.fecha, decreasing = F)] == 560)], col = rgb(0.5,0.5,0.5, alpha = 0.5))

#################################### Figure 3.12. ###################################################

tpower <- function(x, t, p){
  # Truncated p-th power function
  (x-t)^p*(x>t)
}
bbase <- function(x, xl, xr, ndx, deg){
  # Construct a B-spline basis of degree ’deg’
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B
  return(B)
}

P.bspline.smoothing <- function(M, xl=0, xr=1290, argvals = 1:1289 , nbasis=180, norder=3, lambda = 0.1,
                                do.plot = F, ...) {
  
  basis <- create.bspline.basis(rangeval = c(xl,xr), nbasis = nbasis+norder, norder = norder)
  coefs <- matrix(0, ncol = ncol(M), nrow = nbasis+norder)
  B <- bbase(argvals, xl, xr, nbasis, norder)
  D <- ncol(B); D <- (diff(diff(diag(D))))
  P <- t(D)%*%D
  
  for(i in 1:(ncol(M))){
    data=data.frame(tiempo = argvals, curve = rep(0, nrow(M)))
    data$curve <- M[,i]
    iota <- diag(as.numeric(!is.na(data$curve)))
    data$curve[is.na(data$curve)] <- 0
    
    c<-solve(t(B)%*%iota%*%B + lambda*P)%*%t(B)%*%iota%*%data$curve
    
    coefs[,i] <- c
  }
  
  if(do.plot == T){
    plot(NA, xlim = c(1, nrow(M)), ylim = c(min(coefs), max(coefs)), xlab = expression(t[j]), ylab = "")
    for(i in 1:ncol(M)) lines(B%*%coefs[,i])
  }
  
  return(fd(coefs, basisobj = basis))
}


gcv.pen <-function(curve, lambda, argvals){
  data=data.frame(tiempo = argvals, curve = rep(0, length(curve)))
  data$curve <- curve
  data.2=na.omit(data)
  B <- bbase(data.2$tiempo, 0,1, ndx = 180, deg=3)
  D <- ncol(B); D <- ((diff(diag(D))))
  lambda = lambda
  P <- t(D)%*%D
  a <- solve(t(B)%*%B+lambda*P)%*%t(B)%*%data.2$curve
  H = B%*%solve(t(B)%*%B+lambda*P)%*%t(B)
  B.2 <- bbase(seq(0,1,length.out = 1289), 0,1, ndx = 180, deg=3)
  y1 <- B.2%*%a
  trh <- sum(diag(H))
  
  ((1/nrow(data.2))*sum((data$curve-y1)^2, na.rm = T))/((1/nrow(data.2))*trh)
  #2*sum((data$curve-y1)^2, na.rm = T)-2*log(nrow(data.2))+2*log(trh) ## AIC uncomment
}

lambdas <- 0.001 * 10^(seq(0,6,1))
lambdas
gcv <- c()
for(i in lambdas) gcv <- c(gcv,gcv.pen(M[,6], lambda = i, argvals = seq(0,1, length.out = 1289)))
plot(gcv, xaxt="na", xlab = expression(lambda))
axis(1, at=1:7, labels = lambdas)


#################################### Figure 3.13. ###################################################

library(data.table);library(dplyr);library(plyr);library(stringr);library(fda)
archivos <- list.files(pattern = ".csv")
myfiles <- lapply(archivos, fread) # fread siempre es más rápido
pacientes <- readxl::read_xls(list.files(pattern = c(".xls")))
colnames(pacientes)[1] <- "tel_p";
podometro <- myfiles[[1]];podometro <- podometro %>% mutate(valor = as.numeric(as.character(valor)))
podometro <- na.omit(podometro)
podometro$fecha <- as.Date(podometro$fecha, format = "%d/%m/%Y")
podometro <- podometro[podometro$fecha > "2009-01-01",] ## 61
fechas <- podometro$fecha
podometro <- podometro[podometro$valor<=40000,]
podometro <- podometro %>% 
  dplyr::group_by(tel_p, fecha) %>% 
  dplyr::summarise(valor = sum(valor)) %>% 
  dplyr::mutate(fecha = as.numeric(fecha))
y <- pacientes$mortalidad_2years[pacientes$tel_p %in% unique(podometro$tel_p)]
podometro <- as.data.frame(podometro)

splited_podometro <- split(podometro, podometro$tel_p)

a <- as.numeric(min(podometro$fecha));
b <- as.numeric(max(podometro$fecha));
A <- podometro %>% dplyr::select(fecha, valor, tel_p)
splited <- split(A, podometro$tel_p)

splited$`41` <- splited$`41`[-1,] # 2nd value three months later!
splited$`46` <- splited$`46`[-1,] # 2nd value one months later!
splited$`109` <- splited$`109`[-c(1:7),]
splited$`110` <- splited$`110`[-c(1:5),]
splited$`70` <- splited$`70`[-c(1:2),]

splited <- lapply(splited, "[", -c(3))
for(i in 1:112){ splited[[i]]$fecha <- splited[[i]]$fecha - min(splited[[i]]$fecha) }
max.fecha <- sapply(splited, function(x) max(x$fecha))

M <- fda::sparse.mat(splited)
rm(splited_podometro,myfiles,A)

M[which(M[,4]>9000),4] <- 317
M[which(M[,6]>5000),6] <- c(596, 710, 978) ## Curve 6 change
M[which(M[,7]>30000),7] <- 3010 ## Curve 7 change
M[which(M[,12] >15000),12] <- c(2191, 2183, 2150) ## Curve 12
M[which(M[,15] > 25000),15] <- 8935 ## Curve 15
M[which(M[,16]>25000),16] <- 9512 ## Curve 16
M[which(M[,19]>15000),19] <- mean(M[,19], na.rm=T) ## Curve 19
M[which(M[,20] < 1000),20] <- mean(M[,20], na.rm=T)
M[which(M[,23] <1000)[-c(1,2)],23] <- mean(M[,23], na.rm = T)
M[which(M[,27]>15000),27] <- c(5208,3545,2692)
M[which(M[,33]>15000),33] <- 5303
M[which(M[,34]>30000),34] <- 3521
M[which(M[,36]>60000),36] <- NA
M[which(M[,65]>10000),65] <- mean(M[,65], na.rm = T)
M[which(M[,67]>9000),67] <- mean(M[,67], na.rm = T)
M[which(M[,70]>6000),70] <- mean(M[,70], na.rm = T)
M[which(M[,85]==540),85] <- 5400
M[which(M[,90] > 25000),90] <- 2754
M[which(M[,94] > 30000),94] <- 7008
M[which(M[,96]>40000),96] <- 8793
M[which(M[,110]>40000),110] <- mean(M[,110], na.rm=T)
M[1,112] <- 3480



n <- 72
M.72 <- M[1:n,]

n <- 213
M.213 <- M[1:n,]
M.213 <- M.213[, max.fecha>=n]

n <- 365
M.365 <- M[1:n,]
M.365 <- M.365[, max.fecha>=n]

library(imputeTS)

n <- 547
M.2 <- M
M.2[,51] <- na_mean(M[,51])

M.547 <- M.2[1:n,]
M.547 <- M.547[, max.fecha>=n]

n <- 730
M.730 <- M.2[1:n,]
M.730 <- M.730[, max.fecha>=n]



smooth.72 <- P.bspline.smoothing(M.72, xl = 0, xr = 1, argvals = seq(0,1, length.out = 72), nbasis = 40, lambda = 10,norder=3, do.plot=F)
smooth.213 <- P.bspline.smoothing(M.213, xl = 0, xr = 1, argvals = seq(0,1, length.out = 213), nbasis = 50, lambda = 10,norder=3, do.plot=F)
smooth.365 <- P.bspline.smoothing(M.365, xl = 0, xr = 1, argvals = seq(0,1, length.out = 365), nbasis = 50, lambda = 10, norder=3,do.plot=F)
smooth.547 <- P.bspline.smoothing(M.547, xl = 0, xr = 1, argvals = seq(0,1, length.out = 547), nbasis = 60, lambda = 10,norder=3, do.plot=F)

par(mfrow=c(2,2))
B = bbase(seq(0,1, length.out = 72), 0, 1, 40, 3)
plot(NA, xlim = c(1, nrow(M.72)), ylim = c(min(smooth.72$coefs), max(smooth.72$coefs)), xlab = expression(t[j]), ylab = expression(hat(X)(t)[72]))
for(i in 1:ncol(M)) lines(B%*%smooth.72$coefs[,i], col = (y+1)[i])

B = bbase(seq(0,1, length.out = 213), 0, 1, 50, 3)
plot(NA, xlim = c(1, nrow(M.213)), ylim = c(min(smooth.213$coefs), max(smooth.213$coefs)), xlab = expression(t[j]), ylab = expression(hat(X)(t)[213]))
for(i in 1:ncol(M)) lines(B%*%smooth.213$coefs[,i], col = (y[max.fecha>=213]+1)[i])

B = bbase(seq(0,1, length.out = 365), 0, 1, 50, 3)
plot(NA, xlim = c(1, nrow(M.365)), ylim = c(min(smooth.365$coefs), max(smooth.365$coefs)), xlab = expression(t[j]), ylab = expression(hat(X)(t)[365]))
for(i in 1:ncol(M)) lines(B%*%smooth.365$coefs[,i], col = (y[max.fecha>=365]+1)[i])

B = bbase(seq(0,1, length.out = 547), 0, 1, 60, 3)
plot(NA, xlim = c(1, nrow(M.547)), ylim = c(min(smooth.365$coefs), max(smooth.547$coefs)), xlab = expression(t[j]), ylab = expression(hat(X)(t)[547]))
for(i in 1:ncol(M)) lines(B%*%smooth.547$coefs[,i], col = (y[max.fecha>=547]+1)[i])

################################# FUNCTIONAL FPCA ###################################


y.72 <- y[max.fecha >= 72]
y.213 <- y[max.fecha >= 213]
y.365 <- y[max.fecha >= 365]
y.547 <- y[max.fecha >= 547]

pca.72 <- pca.fd(smooth.72, nharm = 3)
pca.213 <- pca.fd(smooth.213, nharm = 3)
pca.365 <- pca.fd(smooth.365, nharm = 3)
pca.547 <- pca.fd(smooth.547, nharm = 3)

print(xtable::xtable(rbind(
  pca.72$varprop,
  pca.213$varprop,
  pca.365$varprop,
  pca.547$varprop
)))



par(mfrow=c(1,4))
plot(pca.72$harmonics[,1], ylab = "")
plot(pca.213$harmonics[,1], ylab = "")
plot(pca.365$harmonics[,1], ylab = "")
plot(pca.547$harmonics[,1], ylab = "")


plot(pca.72$scores[,1:2], xlab = "PC1", ylab = "PC2", pch = 16, col = y+1)
plot(pca.213$scores[,1:2], xlab = "PC1", ylab = "PC2", pch=16, col = y.213+1)
plot(pca.365$scores[,1:2], xlab = "PC1", ylab = "PC2", pch = 16, col = y.365+1)
plot(pca.547$scores[,1:2], xlab = "PC1", ylab = "PC2", pch=16, col = y.547+1)


################################# Functional Logistic Regression Analysis #########################################

fpca.log.reg <- function(fd, argvals = seq(0,1,length.out = 1289), y, method = "FPC_log", do.plot =F, n.pcs = 2){
  require(expm); require(pracma); require(caret)
  A <- t(fd$coefs)
  basis <- fd$basis
  psi <- inprod(basis, basis)
  G <- eigen(cov(A%*%expm::sqrtm(psi)))$vectors
  pca.apsi <- prcomp(A%*%expm::sqrtm(psi), center = F, scale. = F)
  phimat <- bbase(argvals, xl=fd$basis$rangeval[1],xr=fd$basis$rangeval[2], fd$basis$nbasis-3,3)
  if(method=="basis_expansion"){
    X = A%*%psi
    model<- glm(y~X, family = "binomial")
    
    hat.beta <- coefficients(model)[-1]
    alpha <- coefficients(model)[1]
    hat.b.t <- phimat%*%hat.beta
    if(do.plot==T) plot(hat.b.t, type = "l")
    
    basis.coef.2 <- smooth.basis(argvals = argvals, y = hat.b.t[,1], basis)
    beta.basis.coef <- basis.coef.2$fd$coefs
    
    l_0 <- alpha + A%*%psi%*%beta.basis.coef
    pi_0 <- exp(l_0)/(1+exp(l_0))
    library(caret)
    conf.mat <- confusionMatrix(table(as.numeric(pi_0>0.5), y))
  } else if(method =="FPC_log"){
    F.1 <- solve(expm::sqrtm(psi))%*%G
    X.2 = pca.apsi$x[,1:n.pcs]
    model <- glm(y~X.2, family = "poisson")
    gamma <- coefficients(model)[-1]
    f.j <- phimat%*%F.1
    hat.b.t <- f.j[,1:n.pcs]%*%as.matrix(gamma)
    
    F.1 <- solve(expm::sqrtm(psi))%*%G[,1:n.pcs]
    if(do.plot==T) plot(hat.b.t, type = "l")
    
    if(fd$fdnames$funs == "values") {
      name <- "non-centered FPC_log"
      alpha <- coefficients(model)[1] - trapz(argvals,hat.b.t*phimat%*%mean(fd)$coefs)
    } else alpha <- coefficients(model)[1]
    
    beta.basis.coef <- solve(t(phimat)%*%phimat)%*%t(phimat)%*%hat.b.t
    
    l_0 <- alpha + A%*%psi%*%beta.basis.coef
    pi_0 <- exp(l_0)/(1+exp(l_0))
    conf.mat <- confusionMatrix(data = as.factor(as.numeric(pi_0>0.5)), reference = as.factor(y), positive= "1")
  }
  
  
  res <- list(model, hat.b.t , pi_0,conf.mat)  
  class(res) <- name
  return(res)
}


predict.flog <- function(model, train, test, argvals) {
  A <- t(test$coefs)
  basis <- train$basis
  psi <- inprod(basis, basis)
  xl=train$basis$rangeval[1];xr=train$basis$rangeval[2]
  argvals = argvals
  phimat <- bbase(argvals, xl=xl,xr=xr, train$basis$nbasis-3,3)
  
  glm.model <- model[[1]]
  hat.b.t <- model[[2]]
  
  if(class(model) == "non-centered FPC_log"){
    alpha <- coefficients(glm.model)[1] - trapz(argvals,hat.b.t*phimat%*%mean(train)$coefs)
    
    beta.basis.coef <- solve(t(phimat)%*%phimat)%*%t(phimat)%*%hat.b.t
    
    l_0 <- alpha + A%*%psi%*%beta.basis.coef
    pi_0 <- exp(l_0)/(1+exp(l_0))
  }
  return(pi_0) 
}


smooth.72 <- P.bspline.smoothing(M.72, xl = 0, xr = 1, argvals = seq(0,1, length.out = 72), nbasis = 30, lambda = 10,norder=3, do.plot=T)

fanova.onefactor(fdata(smooth.72), group = factor(y))
par(mfrow=c(1,1))
log.res.72 <- fpca.log.reg(smooth.72, argvals = seq(0,1, length.out = 72), y=y, do.plot=T, n.pcs = 1)
log.res.72[[2]]
ROC(log.res.72[[3]], y)
MKmisc::HLgof.test(log.res.72[[3]], y)


fdata.72 <- fdata(smooth.72)
y.df.72 <- data.frame(y=factor(y))

train.pairs <- list("df"=y.df.72,"x"=fdata.72)

basis.b = create.bspline.basis(rangeval=c(0,1),nbasis=50)
basis.x = create.bspline.basis(rangeval = c(0,1), nbasis = 10)

fda.class.glm <- fregre.glm(y~x,family=binomial(),train.pairs, basis.x = basis.x,
                            basis.b=basis.b)
ROC(fda.class.glm$fitted.values, factor(y))
MKmisc::HLgof.test(fda.class.glm$fitted.values, y)

plot(fda.class.glm$beta.l$x)

#####################

##################################################################################################################################
set.seed(12)
n <- ncol(M.72)
b=150

# d <- sample(which(y ==1), size = length(which(y ==1))*0.6)
# d <- c(d, sample(which(y ==0), size = length(which(y ==0))*0.6))
# length(d)
# 
# predictions.2 <- matrix(NA, nrow = n-length(d), ncol = b)
# observations.2 <- predictions.2

set.seed(12)
n <- ncol(M.72)
d <- sample(1:n, size = n*0.6)
b=150
predictions.2 <- matrix(NA, nrow = n-length(d), ncol = b)
observations.2 <- predictions.2

for(i in 1:b){
  d <- sample(1:n, size = n*0.6)
  
  train <- fd(smooth.72$coefs[,d], basisobj = smooth.72$basis); y.train <- y[d]
  test <- fd(smooth.72$coefs[,-d], basisobj = smooth.72$basis); y.test <- y[-d]
  
  res.log.reg <- fpca.log.reg(train, argvals = seq(0,1,length.out = 72), y=y.train, "FPC_log", do.plot=F)
  
  predictions.2[,i] <- predict.flog(res.log.reg, train, test, argvals = seq(0,1, length.out = 72))
  observations.2[,i] <- y.test
}

library(cvAUC)
out.2 <- cvAUC(predictions.2, observations.2)
plot(out.2$perf, col="grey82", lty=3)
plot(out.2$perf, col="red", avg="vertical", add=TRUE)
##################################################################################################################################

set.seed(12)
n <- ncol(M.72)
d <- sample(1:n, size = n*0.6)
b=150
predictions.2 <- matrix(NA, nrow = n-length(d), ncol = b)
observations.2 <- predictions.2
basis.b = create.bspline.basis(rangeval=c(0,1),nbasis=50)
basis.x = create.bspline.basis(rangeval = c(0,1), nbasis = 10)

for(i in 1:b){
  d <- sample(1:n, size = n*0.6)
  
  train <- fdata(fd(smooth.72$coefs[,d], basisobj = smooth.72$basis)); y.train <- data.frame(y=factor(y[d]))
  test <- fdata(fd(smooth.72$coefs[,-d], basisobj = smooth.72$basis)); y.test <- data.frame(y=y[-d])
  
  train.pairs <- list("df" = y.train, "x" = train)
  test.pairs <- list("df" = y.test, "x" = test)
  
  fda.class.glm <- fregre.glm(y~x,family=binomial(),train.pairs, basis.x = basis.x,
                              basis.b=basis.b)
  predictions.2[,i] <- predict(fda.class.glm, test.pairs)
  observations.2[,i] <- y.test$y
}

library(cvAUC)
out.2 <- cvAUC(predictions.2, observations.2)
plot(out.2$perf, col="grey82", lty=3)
plot(out.2$perf, col="red", avg="vertical", add=TRUE)
out.2$cvAUC

##################################################################################################################################
##################################################################################################################################

smooth.213 <- P.bspline.smoothing(M.213, xl = 0, xr = 1, argvals = seq(0,1, length.out = 213), nbasis = 50, lambda = 10, do.plot=T)
smooth.213$coefs
y.213 <- y[max.fecha >= 213]
fanova.onefactor(fdata(smooth.213), group = factor(y.213), plot=T)

par(mfrow=c(1,1))
log.res.213 <- fpca.log.reg(smooth.213, argvals = seq(0,1, length.out = 213), y=y.213, do.plot=T, n.pcs = 2)
plot(log.res.213[[2]], type = "l")
ROC(log.res.213[[3]], y.213)

set.seed(12)
n <- ncol(M.213)
d <- sample(1:n, size = n*0.6)
b=150

predictions.2 <- matrix(NA, nrow = n-length(d), ncol = b)
observations.2 <- predictions.2

for(i in 1:b){
  d <- sample(1:n, size = n*0.6)
  
  train <- fd(smooth.213$coefs[,d], basisobj = smooth.213$basis); y.train <- y.213[d]
  test <- fd(smooth.213$coefs[,-d], basisobj = smooth.213$basis); y.test <- y.213[-d]
  
  res.log.reg <- fpca.log.reg(train, argvals = seq(0,1,length.out = 213), y=y.train, "FPC_log", do.plot=F)
  
  predictions.2[,i] <- predict.flog(res.log.reg, train, test, argvals = seq(0,1, length.out = 213))
  observations.2[,i] <- y.test
}

library(cvAUC)
out.2 <- cvAUC(predictions.2, observations.2)
plot(out.2$perf, col="grey82", lty=3)
plot(out.2$perf, col="red", avg="vertical", add=TRUE)
quantile(out.2$fold.AUC, c(0.1,0.9))
out.2$cvAUC

fdata.213 <- fdata(smooth.213)
y.df.213 <- data.frame(y=factor(y.213))

train.pairs <- list("df"=y.df.213,"x"=fdata.213)

basis.b = create.bspline.basis(rangeval=c(0,1),nbasis=50)
basis.x = create.bspline.basis(rangeval = c(0,1), nbasis = 10)

fda.class.glm <- fregre.glm(y~x,family=binomial(),train.pairs, basis.x = basis.x,
                            basis.b=basis.b)
ROC(fda.class.glm$fitted.values, factor(y.df.213$y))
plot(fda.class.glm$beta.l$x)

##################################################################################################################################
##################################################################################################################################

smooth.365 <- P.bspline.smoothing(M.365, xl = 0, xr = 1, argvals = seq(0,1, length.out = 365), nbasis = 50, lambda = 10, do.plot=T)
smooth.365$coefs
y.365 <- y[max.fecha >= 365]

fanova.onefactor(fdata(smooth.365), group = factor(y.365))

log.res.365 <- fpca.log.reg(smooth.365, argvals = seq(0,1, length.out = 365), y=y.365, do.plot=T, n.pcs = 1)
ROC(log.res.365[[3]], y.365 , plot ="ROC")

fdata.365 <- fdata(smooth.365)
y.df.365 <- data.frame(y=factor(y.365))

train.pairs <- list("df"=y.df.365,"x"=fdata.365)

basis.b = create.bspline.basis(rangeval=c(0,1),nbasis=10)
basis.x = create.bspline.basis(rangeval = c(0,1), nbasis = 10)

fda.class.glm <- fregre.glm(y~x,family=binomial(),train.pairs, 
                            basis.b=basis.b)
ROC(fda.class.glm$fitted.values, factor(y.df.365$y))
plot(fda.class.glm$beta.l$x)
fda.class.glm$fitted.values


##################################################################################################################################
##################################################################################################################################

smooth.547 <- P.bspline.smoothing(M.547, xl = 0, xr = 1, argvals = seq(0,1, length.out = 547), nbasis = 60, lambda = 10, do.plot=T)
y.547 <- y[max.fecha >= 547]

fanova.onefactor(fdata(smooth.547), group = factor(y.547))

log.res.547 <- fpca.log.reg(smooth.547, argvals = seq(0,1, length.out = 547), y=y.547, do.plot=T, n.pcs = 1)
ROC(log.res.547[[3]], y.547 , plot ="ROC")

fdata.547 <- fdata(smooth.547)
y.df.547 <- data.frame(y=factor(y.547))

train.pairs <- list("df"=y.df.547,"x"=fdata.547)

basis.b = create.bspline.basis(rangeval=c(0,1),nbasis=10)
basis.x = create.bspline.basis(rangeval = c(0,1), nbasis = 10)

fda.class.glm <- fregre.glm(y~x,family=binomial(),train.pairs, 
                            basis.b=basis.b)
ROC(fda.class.glm$fitted.values, factor(y.df.547$y))
plot(fda.class.glm$beta.l$x)
fda.class.glm$fitted.values

par(mfrow=c(1,4))
plot(log.res.72[[2]], type="l", xlab = expression(t[j]), ylab = expression(hat(beta)(t)[72]))
plot(log.res.213[[2]], type="l", xlab = expression(t[j]), ylab = expression(hat(beta)(t)[213]))
plot(log.res.365[[2]], type="l", xlab = expression(t[j]), ylab = expression(hat(beta)(t)[365]))
plot(log.res.547[[2]], type="l", xlab = expression(t[j]), ylab = expression(hat(beta)(t)[547]))

par(mfrow = c(2,2))
pi=log.res.72[[3]]
ROC(pi, factor(y), plot="ROC")

pi=log.res.213[[3]]
ROC(pi, factor(y.213), plot="ROC")
pi=log.res.365[[3]]
ROC(pi, factor(y.365), plot="ROC")
pi=log.res.547[[3]]
ROC(pi, factor(y.547), plot="ROC")



set.seed(12)
n <- ncol(M.72)
d <- sample(1:n, size = n*0.6)

b=150
predictions.2 <- matrix(NA, nrow = n-length(d), ncol = b)
observations.2 <- predictions.2

for(i in 1:b){
  d <- sample(1:n, size = n*0.6)
  
  
  train <- fd(smooth.72$coefs[,d], basisobj = smooth.72$basis); y.train <- y[d]
  test <- fd(smooth.72$coefs[,-d], basisobj = smooth.72$basis); y.test <- y[-d]
  
  res.log.reg <- fpca.log.reg(train, argvals = seq(0,1,length.out = 72), y=y.train, "FPC_log", do.plot=F)
  
  predictions.2[,i] <- predict.flog(res.log.reg, train, test, argvals = seq(0,1, length.out = 72))
  observations.2[,i] <- y.test
}

library(cvAUC)
out.2 <- cvAUC(predictions.2, observations.2)
plot(out.2$perf, col="grey82", lty=3)
plot(out.2$perf, col="red", avg="vertical", add=TRUE)
out.2$fold.AUC; 
out.2$cvAUC

library(MKmisc)
HLgof.test(log.res.72[[3]], y)
##################################################################################################################
set.seed(12)
n <- ncol(M.213)
d <- sample(1:n, size = n*0.6)
b=150
predictions.2 <- matrix(NA, nrow = n-length(d), ncol = b)
observations.2 <- predictions.2

for(i in 1:b){
  d <- sample(1:n, size = n*0.6)
  
  train <- fd(smooth.213$coefs[,d], basisobj = smooth.213$basis); y.train <- y.213[d]
  test <- fd(smooth.213$coefs[,-d], basisobj = smooth.213$basis); y.test <- y.213[-d]
  
  res.log.reg <- fpca.log.reg(train, argvals = seq(0,1,length.out = 213), y=y.train, "FPC_log", do.plot=F)
  
  predictions.2[,i] <- predict.flog(res.log.reg, train, test, argvals = seq(0,1, length.out = 213))
  observations.2[,i] <- y.test
}

library(cvAUC)
out.2 <- cvAUC(predictions.2, observations.2)
plot(out.2$perf, col="grey82", lty=3)
plot(out.2$perf, col="red", avg="vertical", add=TRUE)
out.2$fold.AUC; out.2$cvAUC

HLgof.test(log.res.213[[3]],(y.213))

##################################################################################################################
set.seed(12)
n <- ncol(M.365)
d <- sample(1:n, size = n*0.6)
b=150
predictions.2 <- matrix(NA, nrow = n-length(d), ncol = b)
observations.2 <- predictions.2

for(i in 1:b){
  d <- sample(1:n, size = n*0.6)
  
  train <- fd(smooth.365$coefs[,d], basisobj = smooth.365$basis); y.train <- y.365[d]
  test <- fd(smooth.365$coefs[,-d], basisobj = smooth.365$basis); y.test <- y.365[-d]
  
  res.log.reg <- fpca.log.reg(train, argvals = seq(0,1,length.out = 365), y=y.train, "FPC_log", do.plot=F)
  
  predictions.2[,i] <- predict.flog(res.log.reg, train, test, argvals = seq(0,1, length.out = 365))
  observations.2[,i] <- y.test
}

library(cvAUC)
out.2 <- cvAUC(predictions.2, observations.2)
plot(out.2$perf, col="grey82", lty=3)
plot(out.2$perf, col="red", avg="vertical", add=TRUE)
out.2$fold.AUC; out.2$cvAUC
HLgof.test(log.res.365[[3]], y.365)

##################################################################################################################
set.seed(12)
n <- ncol(M.547)
d <- data.frame(x=1:n, y=y.547) %>% group_by(y) %>% sample_frac(0.8)
d <- d$x

b=150
predictions.2 <- matrix(NA, nrow = n-length(d), ncol = b)
observations.2 <- predictions.2

for(i in 1:b){
  d <- data.frame(x=1:n, y=y.547) %>% group_by(y) %>% sample_frac(0.8)
  d <- d$x
  
  train <- fd(smooth.547$coefs[,d], basisobj = smooth.547$basis); y.train <- y.547[d]
  test <- fd(smooth.547$coefs[,-d], basisobj = smooth.547$basis); y.test <- y.547[-d]
  
  res.log.reg <- fpca.log.reg(train, argvals = seq(0,1,length.out = 547), y=y.train, "FPC_log", do.plot=F)
  
  predictions.2[,i] <- predict.flog(res.log.reg, train, test, argvals = seq(0,1, length.out = 547))
  observations.2[,i] <- y.test
}
predictions.2
library(cvAUC)
out.2 <- cvAUC(predictions.2, observations.2)
plot(out.2$perf, col="grey82", lty=3)
plot(out.2$perf, col="red", avg="vertical", add=TRUE)
out.2$fold.AUC; out.2$cvAUC
HLgof.test(log.res.547[[3]], y.547)

fanova.onefactor(fdata(smooth.72), group = factor(y))
fanova.onefactor(fdata(smooth.213), group = factor(y.213))
fanova.onefactor(fdata(smooth.365), group = factor(y.365),nboot = 300, plot=T)
fanova.onefactor(fdata(smooth.547), group = factor(y.547),nboot=300, plot=T)


########################### FANOVA of SPO2 and Respiratory Rythm #################################################

library(data.table);library(dplyr);library(plyr);library(stringr);library(fda)
archivos <- list.files(pattern = ".csv")
myfiles <- lapply(archivos, fread) # fread siempre es más rápido

SPO2 <- myfiles[[10]]
SPO2 <- SPO2 %>% mutate(valor = as.numeric(as.character(valor)))
SPO2 <- na.omit(SPO2)
SPO2$fecha <- as.Date(SPO2$fecha, format = "%d/%m/%Y")


SPO2 <- SPO2[SPO2$fecha > "2009-01-01",] ## 61


SPO2 <- SPO2 %>% 
  dplyr::group_by(tel_p, fecha) %>% 
  dplyr::summarise(valor = mean(valor)) %>% 
  dplyr::mutate(fecha = as.numeric(fecha))

SPO2 <- as.data.frame(SPO2)
SPO2$valor[10689] <- 97

splited_SPO2 <- split(SPO2, SPO2$tel_p)

a <- as.numeric(min(SPO2$fecha));
b <- as.numeric(max(SPO2$fecha));
A <- SPO2 %>% dplyr::select(fecha, valor, tel_p)
splited <- split(A, SPO2$tel_p)

splited$`41` <- splited$`41`[-1,] # 2nd value three months later!
splited$`46` <- splited$`46`[-1,] # 2nd value one months later!
splited$`109` <- splited$`109`[-c(1:7),]
splited$`110` <- splited$`110`[-c(1:5),]
splited$`70` <- splited$`70`[-c(1:2),]

splited <- lapply(splited, "[", -c(3))
for(i in 1:112){ splited[[i]]$fecha <- splited[[i]]$fecha - min(splited[[i]]$fecha) +a}
M.spo2 <- fda::sparse.mat(splited)
M.spo2

n <- 72
M.spo2.72 <- M.spo2[1:n,]

n <- 213
M.spo2.213 <- M.spo2[1:n,]
M.spo2.213 <- M.spo2.213[, max.fecha>=n]


n <- 365
M.spo2.365 <- M.spo2[1:n,]
M.spo2.365 <- M.spo2.365[, max.fecha>=n]

n <- 547
M.spo2.2 <- M.spo2
M.spo2.2[,51] <- na_mean(M.spo2[,51])

M.spo2.547 <- M.spo2.2[1:n,]
M.spo2.547 <- M.spo2.547[, max.fecha>=n]

smooth.spo2.72 <- P.bspline.smoothing(M.spo2.72, xl=0, xr = 1, argvals = seq(0,1, length.out = 72),
                                      nbasis = 30, lambda = 10, do.plot=T)
smooth.spo2.213 <- P.bspline.smoothing(M.spo2.213, xl=0, xr = 1, argvals = seq(0,1, length.out = 213),
                                       nbasis = 50, lambda = 10, do.plot=T)
smooth.spo2.365 <- P.bspline.smoothing(M.spo2.365, xl=0, xr = 1, argvals = seq(0,1, length.out = 365),
                                       nbasis = 50, lambda = 10, do.plot=T)
smooth.spo2.547 <- P.bspline.smoothing(M.spo2.547, xl=0, xr = 1, argvals = seq(0,1, length.out = 547),
                                       nbasis = 60, lambda = 10, do.plot=T)

fanova.onefactor(fdata(smooth.spo2.72), group = factor(y), plot=T)
fanova.onefactor(fdata(smooth.spo2.213), group = factor(y.213), plot=T)
fanova.onefactor(fdata(smooth.spo2.365), group = factor(y.365), plot=T)
fanova.onefactor(fdata(smooth.spo2.547), group = factor(y.547), plot=T)
###################################################################################################################
###################################################################################################################

Ritmo.Resp <- myfiles[[9]]
Ritmo.Resp <- Ritmo.Resp %>% mutate(valor = as.numeric(as.character(valor)))
Ritmo.Resp <- na.omit(Ritmo.Resp)
Ritmo.Resp$fecha <- as.Date(Ritmo.Resp$fecha, format = "%d/%m/%Y")


Ritmo.Resp <- Ritmo.Resp[Ritmo.Resp$fecha > "2009-01-01",] ## 61


Ritmo.Resp <- Ritmo.Resp %>% 
  dplyr::group_by(tel_p, fecha) %>% 
  dplyr::summarise(valor = mean(valor)) %>% 
  dplyr::mutate(fecha = as.numeric(fecha))

Ritmo.Resp <- as.data.frame(Ritmo.Resp)


a <- as.numeric(min(Ritmo.Resp$fecha));
b <- as.numeric(max(Ritmo.Resp$fecha));
A <- Ritmo.Resp %>% dplyr::select(fecha, valor, tel_p)
splited <- split(A, Ritmo.Resp$tel_p)

splited$`41` <- splited$`41`[-1,] # 2nd value three months later!
splited$`46` <- splited$`46`[-1,] # 2nd value one months later!
splited$`109` <- splited$`109`[-c(1:7),]
splited$`110` <- splited$`110`[-c(1:5),]
splited$`70` <- splited$`70`[-c(1:2),]

splited <- lapply(splited, "[", -c(3))
for(i in 1:112){ splited[[i]]$fecha <- splited[[i]]$fecha - min(splited[[i]]$fecha) +a}
M.RitmoResp <- fda::sparse.mat(splited)
M.RitmoResp

n <- 72
M.RitmoResp.72 <- M.RitmoResp[1:n,]

n <- 213
M.RitmoResp.213 <- M.RitmoResp[1:n,]
M.RitmoResp.213 <- M.RitmoResp.213[, max.fecha>=n]


n <- 365
M.RitmoResp.365 <- M.RitmoResp[1:n,]
M.RitmoResp.365 <- M.RitmoResp.365[, max.fecha>=n]

n <- 547
M.RitmoResp.2 <- M.RitmoResp
M.RitmoResp.2[,51] <- na_mean(M.RitmoResp[,51])

M.RitmoResp.547 <- M.RitmoResp.2[1:n,]
M.RitmoResp.547 <- M.RitmoResp.547[, max.fecha>=n]

par(mfrow=c(1,1))
smooth.RitmoResp.72 <- P.bspline.smoothing(M.RitmoResp.72, xl=0, xr = 1, argvals = seq(0,1, length.out = 72),
                                           nbasis = 30, lambda = 10, do.plot=T)
smooth.RitmoResp.213 <- P.bspline.smoothing(M.RitmoResp.213, xl=0, xr = 1, argvals = seq(0,1, length.out = 213),
                                            nbasis = 50, lambda = 10, do.plot=T)
smooth.RitmoResp.365 <- P.bspline.smoothing(M.RitmoResp.365, xl=0, xr = 1, argvals = seq(0,1, length.out = 365),
                                            nbasis = 50, lambda = 10, do.plot=T)
smooth.RitmoResp.547 <- P.bspline.smoothing(M.RitmoResp.547, xl=0, xr = 1, argvals = seq(0,1, length.out = 547),
                                            nbasis = 60, lambda = 10, do.plot=T)

fanova.onefactor(fdata(smooth.RitmoResp.72), group = factor(y), plot=T)
fanova.onefactor(fdata(smooth.RitmoResp.213), group = factor(y.213), plot=T)
fanova.onefactor(fdata(smooth.RitmoResp.365), group = factor(y.365), plot=T)
fanova.onefactor(fdata(smooth.RitmoResp.547), group = factor(y.547), plot=T)


################################## FUNCTIONAL POISSON REGRESSION ########################################

y.prog <- pacientes$prog_adm_COPD[pacientes$tel_p %in% unique(podometro$tel_p)]
t_seguim <- pacientes$t_seguim[pacientes$tel_p %in% unique(podometro$tel_p)]
n <- 730
M.730 <- M.2[1:n,]
M.730 <- M.730[, max.fecha>=n]

y.730 <- y.prog[max.fecha>=n]
y.730 <- y.730
t_seguim.730 <- t_seguim[max.fecha>=n]

y.365 <- y.prog[max.fecha>=365]/2

par(mfrow=c(1,1))
smooth.730 <- P.bspline.smoothing(M.730, xl = 0, xr= 1, argvals = seq(0,1, length.out = 730), nbasis = 50 ,lambda = 1, do.plot=T)
smooth.365 <- P.bspline.smoothing(M.365, xl = 0, xr= 1, argvals = seq(0,1, length.out = 365), nbasis = 25 ,lambda = 10, do.plot=T)


fdata.730 <- fdata(smooth.730)
y.df.730 <- data.frame(y=y.730)



train.pairs <- list("df"=y.df.730,"x"=fdata.730)

basis.b = create.bspline.basis(rangeval=c(0,1),nbasis=50)
basis.x = create.bspline.basis(rangeval = c(0,1), nbasis = 50)


fda.poisson.glm <- fregre.glm(y~x,family=poisson(),train.pairs, 
                              basis.b=basis.b, CV=T)

cbind(fda.poisson.glm$fitted.values, y.730)
plot(fda.poisson.glm$beta.l$x)
par(mfrow=c(2,2))
plot(fda.poisson.glm)

summary(fda.poisson.glm)