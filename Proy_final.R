rm(list=ls())
gc()
require(forecast)
require(ggplot2)
require(mvtnorm)

##Cargar datos

Emigx<-read.csv("Emigx.csv",header=T)
fx5<-read.csv("fecu.csv",header=T)
Inmx<-read.csv("Inmx.csv",header=T)
mx<-read.csv("mx.csv",header=T)
Nx<-read.csv("Nx.csv",header=T)
Px<-read.csv("Px.csv",header=T)

#Renombrar variables

names(Emigx)[3:48]<-c(1970:2015)
names(Inmx)[3:48]<-c(1970:2015)
names(mx)[3:48]<-c(1970:2015)
names(Nx)[3:48]<-c(1970:2015)
names(Px)[3:49]<-c(1970:2016)

#quitar renglon columna de mx
mx <- mx [,-c(1,2)]

#Pasar edades quinquenales a edad simples por medio del método Heather Booth
asfr<-function (fx5,year){ 
  fx0<-fx5[fx5$Year==year,"fx"]
  
  Fx<-5*cumsum(fx0)
  #Es oigual al ultiimo elemento de la suma acumulada
  TGF<-Fx[7]
  #proporcion de las fecundidades acumuluadas con respecto a la TGF (hasta los 30 años se ha acumulado el 70% 
  FxF<-Fx/TGF
  
  #Definir las edades. Secuencia que comienza en 17 al 47, marca de clase de edades quinquenales el cual es el punto medio
  x5<-seq(17.5,47.5,5)
  
  #establecer funcion oara desagregar edades simples. Se multiplica por menos para que sea positivo
  Yx<-log(-log(FxF))
  #sale un infinito pero se omite de la regresion
  
  #Regresion de la varX menos el infinito con respecto a edad quinquenal menos la ultima (seran los primeros 6 valores)
  Yx.lm<-lm(Yx[-7]~x5[-7])
  Yx.lm
  
  #ordenada al origen 
  a<-Yx.lm$coefficients[1]
  #pendiente
  b<-Yx.lm$coefficients[2]
  
  #alfas y betas estimadas
  A<-exp(-exp(a))
  B<-exp(b)
  
  #F(x) estimadas para edades simples
  x1<-c(15:49)
  (Fx.estim<-TGF*A^(B^(x1)))
  
  #para desacumular las tasas especificas de fecundidad (nfx)
  fx<-Fx.estim[2:35]-Fx.estim[1:34]
  #para el calculo del primer valor
  (fx<-c(Fx.estim[1],fx))
}

#Matriz para edades simples
fx<-data.frame(matrix(0,35,46))
row.names(fx)<-c(15:49)
names(fx)<-c(1970:2015)

#Para rellenar matriz de edades simples
for(i in 1:46){
  fx[,i]<-asfr(fx5,1969+i)
}

#Tasas globales de fecundidad
colSums(fx)



#Inmigrantes a la edad x al tiempo t femeninas:
ixt.F <- as.matrix(Inmx[Inmx$Sexo=="Mujeres",-c(1:2)])/  
  as.matrix(Nx[Nx$Sexo=="Mujeres",-c(1:2)])
#Inmigrantes a la edad x al tiempo t masculinos: 
ixt.M <- as.matrix(Inmx[Inmx$Sexo=="Hombres",-c(1:2)])/
  as.matrix(Nx[Nx$Sexo=="Hombres",-c(1:2)])
#Emigrantes a la edad x al tiempo t femeninas:
ext.F<- as.matrix(Emigx[Emigx$Sexo=="Mujeres",-c(1:2)])/
  as.matrix(Nx[Nx$Sexo=="Mujeres",-c(1:2)])
#Emigrantes a la edad x al tiempo t masculinas: 
ext.M<- as.matrix(Emigx[Emigx$Sexo=="Hombres",-c(1:2)])/
  as.matrix(Nx[Nx$Sexo=="Hombres",-c(1:2)])

ixt.F<-ixt.F[1:(dim(ixt.F)[1]-35),]
ixt.M<-ixt.M[1:(dim(ixt.M)[1]-35),]
ext.F<-ext.F[1:(dim(ext.F)[1]-35),]
ext.M<-ext.M[1:(dim(ext.M)[1]-35),]

#Parametros para la función Lee-Carter
edades <- dim(mx) [1] #### igual a 220 grupos por edad
edades.fec <- dim(fx)[1]
tiempo.mort <- dim(mx)[2]
añoini.mort <- 1970
añoini.fec <- 1990
tiempo.fec <- dim(fx)[2] ## años para utilizar la proyeccion de la fecundidad
añobase <- 2015
horizonte <- 35
añofin <- añobase+horizonte
tiempo.tot <- tiempo.mort+horizonte
edades.mig <- dim(ixt.F)[1]
tiempo.mig <- dim(ixt.F)[2]
añoini.mig<-1995

#Función para proyectar por medio de Lee-Carter
lc.svd <- function(m, edades, tiempo1, tiempo2, ln){
  if (ln == TRUE){
    lm<-log(m)
  } else{
    lm <- m 
  }
  
  ax <- rowMeans(lm [, tiempo1:tiempo2])
  
  lm_a <- lm-ax 
  
  d <- matrix(0, nr = min(edades, tiempo2), 
              nc = min(edades, tiempo2))
  
  diag(d) <- svd(lm_a)$d
  
  
  kt <- (d%*%t(-svd(lm_a)$v))
  bx <- -svd(lm_a)$u
  
  lc.svd <- list(ax = ax, bx = bx, kt = kt, D = d)
  
}

#Tabla de mortalidad
tabmort <- function(m,edades,sex){
  
  mx <- m
  
  nax <- matrix(0.5,dim(mx)[1],dim(mx)[2])
  ## 1 MUJERES 2 HOMBRES
  if(sex==1){
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.01724){
        nax[1,i] <- 0.14903-2.05527*mx[1,i]
      }else if(mx[1,i]>=0.01724 & mx[1,i]<0.06891){
        nax[1,i] <- 0.04667+3.88089*mx[1,i]
      }else{nax[1,i] <- 0.31411}
    }
  }else{
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.023){
        nax[1,i] <- 0.14929-1.99545*mx[1,i]
      }else if(mx[1,i]>=0.023 & mx[1,i]<0.08307){
        nax[1,i] <- 0.02832+3.26021*mx[1,i]
      }else{nax[1,i] <- 0.29915}
    }
  }
  
  
  nax[edades,] <- 1/mx[edades,]
  
  qx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 1:(dim(mx)[1])){
    qx[i,]<-mx[i,]/(1+(1-nax[i,])*mx[i,])
  }
  
  px <- 1-qx
  
  lx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 2:dim(mx)[1]){
    lx[i,] <- lx[i-1,]*px[i-1,]
  }
  
  dx <- matrix(0,dim(mx)[1],dim(mx)[2])
  dx[dim(mx)[1],] <- lx[dim(mx)[1],]
  for(i in 1:(dim(mx)[1]-1)){
    dx[i,]<-lx[i,]-lx[i+1,]
  }
  
  
  Lx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Lx[1,] <- dx[1,]/mx[1,]
  Lx[edades,] <- dx[edades,]/mx[edades,]
  for(i in 2:(edades-1)){
    Lx[i,]<-(lx[i,]+lx[i+1,])/2
  }
  
  Tx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Tx[edades,]<-Lx[edades,]
  for(i in (edades-1):1){
    Tx[i,]<-Lx[i,]+Tx[i+1,]
  }
  
  ex <- Tx/lx
  
  Sx<-matrix(0,(dim(mx)[1]+1),dim(mx)[2])
  Sx[1,]<-Lx[1,]/lx[1,]
  Sx[(edades+1),] <- Tx[edades,]/Tx[(edades-1),]
  for(i in 2:edades){
    Sx[i,]<-Lx[i,]/Lx[i-1,]
  }
  
  tabmort <- list(Edad=c(0:(edades-1)),mx=mx, nax=nax, qx=qx, 
                  px=px, lx=lx, dx=dx, Lx=Lx, Tx=Tx, ex=ex, Sx=Sx)
}

#Función para cada componente
lc.mort <- lc.svd(mx, edades,
                  tiempo1 =1 ,
                  tiempo2 =tiempo.mort,
                  ln = TRUE)
lc.fec <- lc.svd(fx, edades = edades.fec, 
                 tiempo1 = 30, 
                 tiempo2 = 46, 
                 ln = TRUE)
lc.inmF <- lc.svd(ixt.F, edades = edades.mig,
                  tiempo1 = 1,
                  tiempo2 = tiempo.mig,
                  ln = TRUE)
lc.inmM <- lc.svd(ixt.M, edades = edades.mig,
                  tiempo1 = 1,
                  tiempo2 = tiempo.mig,
                  ln = TRUE)
lc.enmF <- lc.svd(ext.F, edades = edades.mig,
                  tiempo1 = 1,
                  tiempo2 = tiempo.mig,
                  ln = TRUE)
lc.enmM <- lc.svd(ext.M, edades = edades.mig,
                  tiempo1 = 1,
                  tiempo2 = tiempo.mig,
                  ln = TRUE)

#ARIMA
#mortalidad:
kt1.fit<-auto.arima(lc.mort$kt[1,], trace = T, d = 1)

#Fecundidad: 
ft1.fit <- auto.arima(lc.fec$kt[1,], trace = T,
                      allowdrift = F)
#Inmigracion mujeres:
it1F.fit <- auto.arima(lc.inmF$kt[1,], trace = T, 
                       allowdrift = F, d=1)
#Inmigracion hombres: 
it1M.fit <- auto.arima(lc.inmM$kt[1,], trace = T, 
                       allowdrift = F, d=1)
#Emigracion mujer: 
et1F.fit <- auto.arima(lc.enmF$kt[1,], trace = T, 
                       allowdrift = F, d=1)
#Emigracion hombres: 
et1M.fit <- auto.arima(lc.enmM$kt[1,], trace = T, 
                       allowdrift = F, d=1)


#Forecast por componete
h <- 35 
kt.for <- forecast(kt1.fit, h = h, c(95)) 
ft.for <- forecast(ft1.fit, h = h, c(95)) 
itF.for <- forecast(it1F.fit, h = h, c(95)) 
itM.for <- forecast(it1M.fit, h = h, c(95)) 
etF.for <- forecast(et1F.fit, h = h, c(95)) 
etM.for <- forecast(et1M.fit, h = h, c(95)) 


#Tasa central de mortalidad
mx.for <- exp(lc.mort$ax + lc.mort$bx[,1]%*%t(kt.for$mean))

#Probabilidad de suoervivencia
SxF.for <- tabmort(mx.for[111:220,], edades = 110, sex = 1)$Sx
SxM.for <- tabmort(mx.for[1:110,], edades = 110, sex = 2)$Sx

#Tasas especificas de fecundidad: 
fx.for <- exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$mean))

#Tasas globales de fecundidad proyectada
TGF.for<-colSums(fx.for)


#Calcular las tasas especificas de inmigrantes y emigrantes por sexo: 
ixF.for<-rbind(exp(lc.inmF$ax + lc.inmF$bx[,1]%*%t(itF.for$mean)), matrix(0,35,35))
ixM.for<-rbind(exp(lc.inmM$ax + lc.inmM$bx[,1]%*%t(itM.for$mean)), matrix(0,35,35))
exF.for<-rbind(exp(lc.enmF$ax + lc.enmF$bx[,1]%*%t(etF.for$mean)), matrix(0,35,35))
exM.for<-rbind(exp(lc.enmM$ax + lc.enmM$bx[,1]%*%t(etM.for$mean)), matrix(0,35,35))


#Poblaciones de mujeres conciliadas a inicio de año
PxF<-Px[Px$Sexo=="Mujeres",-c(1,2)]
PxM<-Px[Px$Sexo=="Hombres",-c(1,2)]

NxF<-Nx[Nx$Sexo=="Mujeres",-c(1,2)]
NxM<-Nx[Nx$Sexo=="Hombres",-c(1,2)]


#Matrices para proyecciones
PxF.for<-matrix(0,110,36)
PxM.for<-matrix(0,110,36)

NxF.for<-matrix(0,110,36)
NxM.for<-matrix(0,110,36)


#Para llenar 2015
PxF.for[,1]<-PxF[,"2016"]  
PxM.for[,1]<-PxM[,"2016"]
NxF.for[,1]<-NxF[,"2015"]
NxM.for[,1]<-NxM[,"2015"]
Bx<-matrix(0,35,35)
BF<-vector(length=35)
BM<-vector(length=35)

###Proyecciones para mujeres 

for(i in 2:36){
  
  ####Edades intermedias 1 a 108 años
  PxF.for[2:109,i]<-(PxF.for[1:108,i-1]+
                       0.5*NxF.for[1:108,i-1]*ixF.for[1:108,i-1])*SxF.for[1:108,i-1]+
    NxF.for[2:109,i-1]*0.5*ixF.for[2:109,i-1]- NxF.for[1:108,i-1]*exF.for[1:108,i-1]
  
  View(PxF.for)
  
  #ultimo grupo abierto 109+ (todos los que estaban en 108)
  PxF.for[110,i]<-(PxF.for[109,i-1]+
                     0.5*NxF.for[109,i-1]*ixF.for[109,i-1])*SxF.for[109,i-1]+
    NxF.for[110,i-1]*exF.for[109,i-1]+
    (PxF.for[110,i-1]+
       NxF.for[110,i-1]*0.5*ixF.for[110,i-1])*SxF.for[110,i-1]+
    NxF.for[110,i-1]*0.5*ixF.for[110,i-1]
  
  #Nacimientos
  Bx[,i-1]<-fx.for[,i-1]*(PxF.for[16:50,i-1]+
                            0.5*NxF.for[16:50,i-1]*ixF.for[16:50,i-1]+
                            PxF.for[16:50,i])*0.5
  
  #Nacimientos por sexo (mujeres)
  BF[i-1]<-(1/2.05)*sum(Bx[,i-1])
  
  #Primer grupo de edad
  PxF.for[1,i]<-BF[1]+SxF.for[1,i-1] +
    NxF.for[1,i-1]*0.5*ixF.for[1,i-1]-
    NxF.for[1,i-1]*exF.for[1,i-1]
  
  #Poblacion a mitad año
  NxF.for[,i]<-0.5*(PxF.for[,i-1]+PxF.for[,i])
  
}


###Proyecciones hombres

for(i in 2:36){
  
  ####Edades intermedias 1 a 108 años
  PxM.for[2:109,i]<-(PxM.for[1:108,i-1]+
                       0.5*NxM.for[1:108,i-1]*ixM.for[1:108,i-1])*SxM.for[1:108,i-1]+
    NxM.for[2:109,i-1]*0.5*ixM.for[2:109,i-1]- NxM.for[1:108,i-1]*exM.for[1:108,i-1]
  
  View(PxM.for)
  
  #ultimo grupo abierto 109+ (todos los que estaban en 108)
  PxM.for[110,i]<-(PxM.for[109,i-1]+
                     0.5*NxM.for[109,i-1]*ixM.for[109,i-1])*SxM.for[109,i-1]+
    NxM.for[110,i-1]*exM.for[109,i-1]+
    (PxM.for[110,i-1]+
       NxM.for[110,i-1]*0.5*ixM.for[110,i-1])*SxM.for[110,i-1]+
    NxM.for[110,i-1]*0.5*ixM.for[110,i-1]
  
  #Nacimientos
  Bx[,i-1]<-fx.for[,i-1]*(PxM.for[16:50,i-1]+
                            0.5*NxM.for[16:50,i-1]*ixM.for[16:50,i-1]+
                            PxM.for[16:50,i])*0.5
  
  #Nacimientos por sexo (mujeres)
  BM[i-1]<-(1/2.05)*sum(Bx[,i-1])
  
  #Primer grupo de edad
  PxM.for[1,i]<-BF[1]+SxM.for[1,i-1] +
    NxM.for[1,i-1]*0.5*ixM.for[1,i-1]-
    NxM.for[1,i-1]*exM.for[1,i-1]
  
  #Poblacion a mitad año
  NxM.for[,i]<-0.5*(PxM.for[,i-1]+PxM.for[,i])
  
}


#Redondear proyecciones
PxF.for<-round(PxF.for,0)
PxM.for<-round(PxM.for,0)
NxF.for<-round(NxF.for,0)
NxM.for<-round(NxM.for,0)

#Total Poblacion a mitad de año hombres y mujeres
colSums(NxM.for)+colSums(NxF.for)
colSums(PxM.for)+colSums(PxF.for)

#Pasar a excel
write.csv(PxM.for, file="PxM.for.csv")
write.csv(PxF.for, file="PxF.for.csv")


#Poblacion total
PxM.T<-colSums(PxM)
PxF.T<-colSums(PxF)

PxM.for.T<-colSums(PxM.for)
PxF.for.T<-colSums(PxF.for)

PxM.T<-(PxM.T[-c(46:47)])
PxF.T<-(PxF.T[-c(46:47)])

Px.dat<-data.frame(año=c(1970:2050),
                   hombres=c(PxM.T, PxM.for.T),
                   mujeres=c(PxF.T, PxF.for.T))

#poblacion proyectada
matplot(PxF.for,type="l")
matplot(PxM.for,type="l")


#Mortalidad
matplot(mx.for, type="l")
  #esperanza de vida
exF.for <- tabmort(mx.for[111:220,], edades = 110, sex = 1)$ex
exM.for<-tabmort(mx.for[1:110,], edades = 110, sex = 2)$ex
matplot(exF.for[1,], type="l")
ex.dat<-data.frame(year=c(2015:2050), hombres=c(exF.for), mujeres=c(exM.for))


#Tasas glabales de fecundidad
TGF<-colSums(fx)
TGF<-colSums(fx.for)
TGF.dat<-data.frame(year= c(1970:2050),
                    mean= c(TGF, TGF.for))
ggplot(TGF.dat, aes(x=year, y=mean))+
  geom_line(aes(y=mean), col= "blue")


#Inmigrantes
ixF.for.g<-ixF.for[1:(dim(ixF.for)[1]-35),]
ixM.for.g<-ixM.for[1:(dim(ixM.for)[1]-35),]
matplot(ixF.for.g,type="l")
matplot(ixM.for.g,type="l")




