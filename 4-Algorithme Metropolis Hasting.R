f1_tilde<-function(x1,x2){
  return(exp(-(x1^2)/8 - (x2^2)/2 )*(abs(x2)<1))
}
f2_tilde<-function(x1,x2){
  return(( cos(x1)^2 + (1/2)*(sin(3*x2)^2)*cos(x1)^4 )*exp(-(x1^2)/8 - (x2^2)/2 ))
}
densite_g<-function(x1,x2){
  return((1/4*pi)*exp((-1/2)*((x1^2)/4+x2^2)))
}
#Question1
algo_metroplis<-function(n,f){
  #Intialiser xo
  X<-matrix(0,n,ncol=2)
  #Simuler les x(t) revient a simuler des variables aleatoires suivant une loi binomiale de paramètre alpha(la deuxieme partie du min) par la methode de la fonction inverse
  for (i in 2:n) {
    eps1<-rnorm(1,0,2)
    eps2<-rnorm(1,0,1)
    u<-runif(1)
    p<-(f(eps1,eps2)*densite_g(X[i-1,1],X[i-1,2]))/(densite_g(eps1,eps2)*f(X[i-1,1],X[i-1,2]))
    X[i,]<- (u<=p)*c(eps1,eps2)+(u>p)*X[i-1,]
  }
  return(X)
}
chaine_markov_Y<-algo_metroplis(1000,f2_tilde)
chaine_markov_X<-algo_metroplis(1000,f1_tilde)

#Question2
convergence_MH_f1<-function(n,f){
  x<-algo_metroplis(n,f)
  v<-(exp(x[,1])+exp(x[,2]))>=5
  s<-seq(1,n,1)
  vect_estim<-cumsum(v)/s
  return(vect_estim)
}
convergence_MH_f2<-function(n,f){
  x<-algo_metroplis(n,f)
  v<-cos(x[,1]*x[,2])*sin(x[,1])*exp(sin(x[,1]+x[,2]))
  s<-seq(1,n,1)
  vect_estim<-cumsum(v)/s
  return(vect_estim)
}
dev.off()
par(mfrow=c(1,2))
plot(1:10000,convergence_MH_f1(10000,f1_tilde),type='l',col="green",main="Approximation (1) pour f1",xlab="T",ylab="")
plot(1:10000,convergence_MH_f2(10000,f2_tilde),type='l',col="blue",main="Approximation (1) pour f2",xlab="T",ylab="")

#Question3 : Commentaires sur le rapport 