#Algorithme de rejet pour f1 (REECRIRE LALGO POUR QUIL PRENNE F EN ARGUMENT POUR EVITER DE FAIRE DEUX ALGO POUR F1 F2)
#Question3
#on prendra n=1000
n<-10000
#Constante de normalisation
c<-1/((4*pi)*(pnorm(1)-pnorm(-1)))
#Temps d'execution trop lent (23s) probleme qq part 

rejet_f1<-function(n){
  ans<-NULL
  m<-n
  while (m>0) {
    x1<-rnorm(m%/%(4*pi*c)+1,0,2)
    x2<-rnorm(m%/%(4*pi*c)+1,0,1)
    u<-runif(m%/%(4*pi*c)+1)
    y<-(u <= (1/(4*pi))*(abs(x2)<=1))*cbind(x1,x2)
    ans<-rbind(ans,y[which(y[,1]!=0 | y[,2]!=0),])
    m<- n-length(ans[,1])
  }
  return(ans[1:n,])
}

#Densité de f1 
f1<-function(x1,x2){
  return(c*exp(-(x1^2)/8 - (x2^2)/2 )*(abs(x2)<1))
}
#Densité marginale f1(x1)
densite_marginale1 <-function(x){
  return(exp(-x^2/8)/(2*sqrt(2*pi)))
}
#Densité marginale f2(x2)
densite_marginale2 <-function(x){
  return(2*sqrt(2*pi)*c*exp(-x^2/2)*(abs(x<=1)))
}

#Echantillon pour n realisations 
x=rejet_f1(n)
dev.off()
par(mfrow=c(1,2)) #Pour afficher les deux graphiques ensembles

#Question4

#On valide maintenant notre algorithme du rejet en comparant la distribution de l’échantillon pour la 1ere variable avec la densité de la loi marginale de x1.
hist(x[,1], freq = FALSE, main = "Histogramme de x1", ylab = "Fréquences",xlab="x1")
t1 <- seq(min(x[,1]), max(x[,1]), length.out=nrow(x))
lines(t1, densite_marginale1(t1), col = "mediumseagreen", lwd = 3)
legend("topright", "f1(x1)", box.lty = 0, bg = "gray90", col = "mediumseagreen", lty = 1, lwd = 3, inset = 0.0001)

#On valide maintenant notre algorithme du rejet en comparant la distribution de l’échantillon pour la 2ere variable avec la densité de la loi marginale de x2.
hist(x[,2], freq = FALSE, main = "Histogramme de x2", ylab = "Fréquences",xlab="x2")
t2 <-seq(min(x[,2]), max(x[,2]), length.out=nrow(x))
lines(t2, densite_marginale2(t2), col = "mediumseagreen", lwd = 3)
legend("topright", "f1(x2)", box.lty = 0, bg = "gray90", col = "mediumseagreen", lty = 1, lwd = 3, inset = 0.01)
