#Question 13 

f_densite<-function(x,y,cst)
{
  c=1/(4*pi*((2*pnorm(1)-1))+cst*4*pi)
  res<-c*exp((-1/2)*((x^2)/4+y^2))*((abs(y)<1)+cst)
  return (res)
}


#Expliquer pourquoi on prend sqrt(2)
g_densite<-function(x,y)
{
  return(dnorm(x,0,2)*dnorm(y,0,sqrt(2)))
}


#Algortihme de rejet qui retourne les va acceptées mais aussi les va rejetées 
rejet_bis<-function(n,cst){
  M<-8*pi*sqrt(2)/(4*pi*((2*pnorm(1)-1))+cst*4*pi)
  x1<-rnorm(n,0,2)
  x2<-rnorm(n,0,sqrt(2))
  u<-runif(n)
  #Elements acceptés
  y<-(u<=f_densite(x1,x2,cst)/(M*g_densite(x1,x2)))*cbind(x1,x2)
  #Elements rejetés 
  z<-(u> f_densite(x1,x2,cst)/(M*g_densite(x1,x2)))*cbind(x1,x2)
  return(data.frame(Y=y,Z=z))
}

#Fonction qui retourne les estimateurs delta1,delta2,combinaison delta1,delta2
estim_rejet_recycle<-function(n,cst)
{
  M<-8*pi*sqrt(2)/(4*pi*((2*pnorm(1)-1))+cst*4*pi)
  
  x<-rejet_bis(n,cst)
  #delta1
  y1_accepted<-data.frame(x[which(x$Y.x1!=0|x$Y.x2!=0),])$Y.x1
  y2_accepted<-data.frame(x[which(x$Y.x1!=0|x$Y.x2!=0),])$Y.x2
  
  h_y_accepted<-(exp(y1_accepted)+exp(y2_accepted))>5

  delta1<-cumsum(h_y_accepted)/(1:length(y1_accepted))
 
  #Variance de delta1
  var_delta1<-var(h_y_accepted)/length(y1_accepted)
 

  #delta2
  z1_rejected<-data.frame(x[which(x$Z.x1!=0|x$Z.x2!=0),])$Z.x1
  z2_rejected<-data.frame(x[which(x$Z.x1!=0|x$Z.x2!=0),])$Z.x2
  
 
  f_z<-f_densite(z1_rejected,z2_rejected,cst)
  g_z<-g_densite(z1_rejected,z2_rejected)
  h_z_rejected<-(exp(z1_rejected)+exp(z2_rejected))>5
  
  delta2<-((M-1)*cumsum((h_z_rejected*f_z)/(M*g_z - f_z)))/(1:length(z2_rejected))
  
  #Variance de delta2
  var_delta2<-var(h_z_rejected)/length(z1_rejected)
  

  #alpha:estimée au cours d'une periode de chauffe pour eviter des problemes de dependances
  var_d2_echant<-var(h_z_rejected[1:min(length(z1_rejected),length(y1_accepted)/2)])
  var_d1_echant<-var(h_y_accepted[1:min(length(z1_rejected),length(y1_accepted)/2)])
  alpha<-var_d2_echant/(var_d1_echant+var_d2_echant)
  
  
  #delta3
  D1=delta1[1:min(length(delta1),length(delta2))]
  D2=delta2[1:min(length(delta1),length(delta2))]
  delta3=alpha*D1+(1-alpha)*D2
  
  
  return (list(delta1,delta2,delta3))
 
}


dev.off()
par(mfrow=c(1,2))
x=estim_rejet_recycle(100000,0)
plot(x[[1]],type='l',col="black",main="Convergence des estimateurs pour f1",xlab="n",ylab="")
lines(x[[2]],type='l',col="royalblue")
lines(x[[3]],type='l',col="firebrick")
legend("bottomright", c(expression(delta[1]),expression(delta[2]),expression(alpha*delta[1]+(1-alpha)*delta[2])),col=c("black","royalblue","firebrick"),lty = c(1, 1,1), lwd = c(3, 3,3), box.lty = 0, inset = 0.01, bg = "gray95",cex=0.45)
abline(0.25,0)
x=estim_rejet_recycle(100000,0.5)
plot(x[[1]],type='l',col="black",main="Convergence des estimateurs pour f3",xlab="n",ylab="")
lines(x[[2]],type='l',col="royalblue")
lines(x[[3]],type='l',col="firebrick")
legend("bottomright", c(expression(delta[1]),expression(delta[2]),expression(alpha*delta[1]+(1-alpha)*delta[2])),col=c("black","royalblue","firebrick"),lty = c(1, 1,1), lwd = c(3, 3,3), box.lty = 0, inset = 0.01, bg = "gray95",cex=0.45)