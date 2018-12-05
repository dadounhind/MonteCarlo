#Premier cas: p:=P(exp(x1)+exp(x2)>=5)

#Pour les différentes méthodes d’estimation ci-dessous, on tracera sur un même graphique l’évolution de l’estimateur 
#et l’évolution de l’intervalle de confiance au niveau 0.95 associé en fonction du nombre de simulation. 
#On prendra n = 10000.
n<-10000
level<-0.95
#Echantillon utilisé pour la methode de MC et va antithetique
x<-rejet_f1(n)

#Question 5a: Estimateur de MonteCarlo
estim_mc<-function(n,x){
  v<-exp(x[,1])+exp(x[,2])>=5
  s<-seq(1,n,1)
  vect_estim<-cumsum(v)/s
  vari<-var(v)
  print(vari)
  IC_min=vect_estim-qnorm((1+level)/2)*sqrt(vari/(1:length(x[,1])))
  IC_max=vect_estim+qnorm((1+level)/2)*sqrt(vari/(1:length(x[,1])))
  return(data.frame(estim=vect_estim,IC_min=IC_min,IC_max=IC_max))
}

#Question 5b: Estimateur par variable antithetique
#On prend A(X)=-X, A(X) suit la meme loi que X et (A et G) de monotonie differente 
estim_antith <-function(n,x){
  v1=((exp(x[,1])+exp(x[,2]))>=5) 
  v2=((exp(-x[,1])+exp(-x[,2]))>=5)
  vect_estim<-cumsum((v1+v2)/2)/seq(1,n,1)
  v<-(v1+v2)/2
  vari<-var(v)
  IC_min=vect_estim-qnorm((1+level)/2)*sqrt(vari/(1:length(x[,1])))
  IC_max=vect_estim+qnorm((1+level)/2)*sqrt(vari/(1:length(x[,1])))
  return(data.frame(estim=vect_estim,IC_min=IC_min,IC_max=IC_max))
}

#Question 5c:Estimateur par variable de controle
#l'integrale de la densité de x1+x2 avec x2 une normale tronquée depend de la fonction de repartition de x pour une normale 0,1
#Chose difficile à calculer, on pourrait eventuellement utiliser la fonction integrate pour trouver sa valeur numerique mais on perd bcp de precision(check livre)
#On pourrait aussi simuler un echantillon de cette densité par la méthode de rejet 
#On a choisit de considerer la loi de x1+x2 des normales car les informations sont concentrés sur -1,1 pour une normale 0,1 donc on perd pas bcp d'informations
#On choisit m = 200
m<-200
estim_va_controle<-function(n,k,m)
{
  #On commence par simuler un echantillon de taille n de X qui a pour densité f1
  x1<-rnorm(n,0,2)
  x2<-rnorm(n,0,1)
  #Estimation de la variable de controle
  #On prend un echantillon de taille m pour estimer la variable de controle b 
  
  g_burning=as.numeric(((exp(x1[1:m]) + exp(x2[1:m])) >5)*(abs(x2[1:m])<=1))/(pnorm(1)-pnorm(-1))
  g_burning=g_burning-mean(g_burning)
  
  #On sait calculer la loi de X1+X2 avec X1 qui suit une Normale(0,4) et X2 qui suit une Normale(0,1)
  #On va donc utiliser h(x)=(X1+X2>2*log(k)) pour trouver notre variable de controle
  
  esp_h<-1-pnorm(2*log(k),0,sqrt(5)) #Esperance de h(x)
  h_burning=as.numeric(sqrt(exp(x1[1:m]+x2[1:m]))>k)
  h_burning<-h_burning-esp_h
  
  #On trouve b :
  b<-mean(g_burning*h_burning)/mean(h_burning^2)
  
  #On prend ensuite un echantillon de taille n-m pour estimer p 
  g=as.numeric(((exp(x1[(m+1):n]) + exp(x2[(m+1):n]) )>5)*(abs(x2[(m+1):n])<=1))/(pnorm(1)-pnorm(-1))
  h=as.numeric(sqrt(exp(x1[(m+1):n]+x2[(m+1):n]))>k)
  #Coeffiction de correlation
  print(cov(g,h)/sqrt(var(h)*var(g)))
  
  #Finalement on a par le cours que l'estimateur à la forme suivante
  #p<-mean(g-b*(h-esp_h))
  vect_estim<-cumsum(g-b*(h-esp_h))/seq(1,n-m,1)
  v<-g-b*(h-esp_h)
  vari<-var(v)
  IC_min=vect_estim-qnorm((1+level)/2)*sqrt(vari/(1:length(v)))
  IC_max=vect_estim+qnorm((1+level)/2)*sqrt(vari/(1:length(v)))
  return(data.frame(estim=vect_estim,IC_min=IC_min,IC_max=IC_max))
}
#Question 5d:Estimateur par methode de stratification
estim_strat<-function (n,L){
  Y=matrix(0,nrow=L,ncol=n/L)
  V=rep(0,n)
  IC_min=rep(0,n)
  IC_max=rep(0,n)
  #la proba que X1 soit dans un des intervalles verifie 
  P_dl<-1/L
  nL <-floor(n/L)
  for (l in 2:(L-1)) {
    #x1 definit sur un ensemble D qu'on va partionné de la maniere suivante (comme tp3)
    dl<-c(2*qnorm((l-1)/L),2*qnorm(l/L))
    #On simule la loi conditionnelle de X1
    u=runif(nL)
    X1_l=qnorm( pnorm(dl[2],0,2) + u*(pnorm(dl[2],0,2) - pnorm(dl[1],0,2) ),0,2)
    X2_l=X2<-rejet_f1(nL)[,2]
    Y[l,]=cumsum((exp(X1_l)+exp(X2_l))>=5)/(1:nL)
    V[l]=var((exp(X1_l)+exp(X2_l))>=5)
  }
  estim=colSums(Y)*P_dl
  V<-P_dl*sum(V)
  print(V)
  IC_min=estim-qnorm((1+level)/2)*sqrt((V)/(1:nL))
  IC_max=estim+qnorm((1+level)/2)*sqrt((V)/(1:nL))
  return(data.frame(estim=estim,IC_min=IC_min,IC_max=IC_max))
}

#Annexe 


#Fonction qui donne un graphique de l'evolution de l'estimateur et son intervalle de confiance en fonction de n 
graphique<-function(x,methode){
  n<-length(x$estim)
  plot(1:n,x$estim,type="l",col = "royalblue", lwd = 2,
       main ="",
       xlab = "n",ylab="")
  title(main = list(methode, cex = 0.8,
                    col = "black", font = 2))
  lines(1:n,x$IC_min,col="mediumseagreen")
  lines(1:n,x$IC_max,col="firebrick")
  legend("bottomright", c( "Etimateur","Borne sup IC","Borne inf IC"),col=c("royalblue","firebrick","mediumseagreen"),lty = c(1, 1,1), lwd = c(3, 3,3), box.lty = 0, inset = 0.05, bg = "gray95",cex=0.45)
}
#Graphiques 
dev.off()
par(mar=c(0.7,0.7,0.7,0.7))
par(mfrow=c(2,2))
graphique(estim_mc(n,x),'Monte Carlo classique')
graphique(estim_antith(n,x),'Variable antithétique')
graphique(estim_va_controle(n,2,200),'Variable de controle')
graphique(estim_strat(n,3),'Stratification par allocation proportionnelle')

#Temps d'execution des fonctions 
sys_time_mc<-system.time(estim_mc(n,x))
sys_time_antith<-system.time(estim_antith(n,x))
sys_time_va_controle<-system.time(estim_va_controle(n,2,200)) 
sys_time_strat<-system.time(estim_strat(n,10))


#Seconde partie 
estim_mc_f2<-function(n){
  x1<-rnorm(n,0,2)
  x2<-rnorm(n,0,1)
  v<-(cos(x1)^2+0.5*sin(3*x2)^2*cos(x1)^4)
  return(mean(v))
}
rejet_f2<-function(n){
  c<-estim_mc_f2(n)
  ans<-NULL
  m<-n
  while (m>0) {
    x1<-rnorm(m%/%(6*pi*c)+1,0,2)
    x2<-rnorm(m%/%(6*pi*c)+1,0,1)
    u<-runif(m%/%(6*pi*c)+1)
    y<-(u <= (1/(6*pi))*(cos(x1)^2+0.5*sin(3*x2)^2*cos(x1)^4))*cbind(x1,x2)
    ans<-rbind(ans,y[which(y[,1]!=0 | y[,2]!=0),])
    m<- n-length(ans[,1])
  }
  return(ans[1:n,])
}
x=rejet_f2(n)
fonction_int<-function(x)
{
  return (cos(x[,1]*x[,2])*sin(x[,1])*exp(sin(x[,1]+x[,2])))
}


estim_int_mc<-function(n,x)
{
  y=cumsum(fonction_int(x))/(1:n)
  return (y)
}

estim_int_anth<-function(n,x)
{
  z=cumsum((fonction_int(x)+fonction_int(-x))/2)/(1:n)
  return(z)
  
}


