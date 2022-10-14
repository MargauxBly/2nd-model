library(shiny)
library(MASS) 
library(rockchalk)
library(ggplot2)
library(SuppDists)
library(gridExtra)
library(png) 
library(doParallel)
library(dplyr)



Nb_maint<-function(n,tau_periode){
  return(floor((n - 1) / (tau_periode + 1)))
}

ARD1 <-  function(al,ale,mu1,mu2,sig21,sig22,rho,r, n,p,tau_periode,perio) { #r coef de correlation
    Yt <- NULL 
    nb_maint = floor((n - 1) / (tau_periode + 1))
    #mu=rep(mu,2) ; sigma2=rep(sigma2,2) #on simule 2 processsus
    n = n - nb_maint
    if (perio) {
      #t <- seq(0, (n * p - 1), by = p)
       t<-seq(0, by = p,length.out = n)
    } else{
      set.seed(al)
      t <- c(0, sort(sample(1:(n * 5), (n - 1))))
      set.seed(NULL)
    }
      if(perio){ #si periodique
      covar=r*sqrt(sig21*sig22)
      COV<-matrix(c(sig21*p,rep(covar,2)*p,sig22*p),2,2) #COV (X_a,X_b)
      set.seed(ale)
      inc<-rockchalk::mvrnorm(n-1,mu=c(mu1*p,mu2*p),Sigma=COV)
      
      }
    Xt<-rbind(0,apply(inc,2,cumsum))
    
    if (nb_maint > 0) {
      taus0<-t[seq(tau_periode + 1, by = tau_periode, length.out = nb_maint)]
      taus <- c(0, taus0, t[length(t)])
      Xt0 <- Xt[, 1][t <= taus[2]]
      ty <- sort(c(t, taus[-c(1, length(taus))]))
      for (i in 2:(length(taus) - 1)) {
        Xt2 <- Xt[, 2][t <= taus[i + 1] & t >= taus[i]]
        Xt1 <- Xt[, 1][t <= taus[i + 1] & t >= taus[i]]
        Xtf = Xt1 - rho * Xt2[1]
        Yt <- c(Yt, Xtf)
        
      }
      Yt <- c(Xt0, Yt)

    } else{
      taus <- 0
      Yt = Xt[, 2]
    }
    df_Yt <- data.frame(Temps = ty, Degradation = Yt)
    list(dfY = df_Yt, t = t, taus = taus,Xt=Xt)
}

Xtab.bis<-function(Xt){
  Xa=Xt[,1]-Xt[,2]
  Xb=Xt[,2]
  return(cbind(Xa,Xb))
}

Xt12.bis<-function(Xt){
  X1=Xt[,1]+Xt[,2]
  X2=Xt[,2]
  return(cbind(X1,X2))
}

Xtab<-function(ale,mua,mub,sig2a,sig2b,corrab, n,p,tau_periode){
  nb_maint = floor((n - 1) / (tau_periode + 1))
  
  n = n - nb_maint
 
  covar=corrab*sqrt(sig2a*sig2b)
  COV<-matrix(c(sig2a*p,rep(covar,2)*p,sig2b*p),2,2) #COV (X_a,X_b)
  set.seed(ale)
  inc<-rockchalk::mvrnorm(n-1,mu=c(mua*p,mub*p),Sigma=COV)
  
  Xt<-rbind(0,apply(inc,2,cumsum))
  
  return(Xt)
}

ARD1_12bis <-  function(Xt,p,rho,tau_periode,n) { #r coef de correlation
  Yt <- NULL 
  nb_maint = floor((n - 1) / (tau_periode + 1))
  t<-seq(0, by = p,length.out = n-nb_maint)
  if (nb_maint > 0) {
    taus <-
      c(0, t[seq(tau_periode + 1, by = tau_periode, length.out = nb_maint)], t[length(t)])
    Xt0 <- Xt[, 1][t <= taus[2]]
    ty <- sort(c(t, taus[-c(1, length(taus))]))
    for (i in 2:(length(taus) - 1)) {
      Xt2 <- Xt[, 2][t <= taus[i + 1] & t >= taus[i]]
      Xt1 <- Xt[, 1][t <= taus[i + 1] & t >= taus[i]]
      Xtf = Xt1- rho * Xt2[1]
      Yt <- c(Yt, Xtf)
      
    }
    Yt <- c(Xt0, Yt)
    
  } else{
    taus <- 0
    Yt = Xt[, 2]
  }
  df_Yt <- data.frame(Temps = ty, Degradation = Yt)
  list(dfY = df_Yt, t = t, taus = taus)
}



ARD1_AB <-  function(Xt,p,rho,tau_periode,n) { #r coef de correlation
  Yt <- NULL 
  nb_maint = floor((n - 1) / (tau_periode + 1))
  t<-seq(0, by = p,length.out = n-nb_maint)
  if (nb_maint > 0) {
    taus <-
      c(0, t[seq(tau_periode + 1, by = tau_periode, length.out = nb_maint)], t[length(t)])
    Xt0 <- Xt[, 1][t <= taus[2]]+Xt[, 2][t <= taus[2]]
    ty <- sort(c(t, taus[-c(1, length(taus))]))
    for (i in 2:(length(taus) - 1)) {
      Xtb <- Xt[, 2][t <= taus[i + 1] & t >= taus[i]]
      Xta <- Xt[, 1][t <= taus[i + 1] & t >= taus[i]]
      Xtf = Xta+Xtb - rho * Xtb[1]
      Yt <- c(Yt, Xtf)
      
    }
    Yt <- c(Xt0, Yt)
    
  } else{
    taus <- 0
    Yt = Xt[, 2]
  }
  df_Yt <- data.frame(Temps = ty, Degradation = Yt)
  list(dfY = df_Yt, t = t, taus = taus)
}



manip_df<-function(df){
  Dt=diff(df[,"Temps"]) #Delta t_{j,i} et temps des sauts
  Dy=diff(df[,"Degradation"]) #Delta y_{j,i} et sauts
  df2=data.frame(DTemps=Dt,DDegradation=Dy)
  nb_maint=sum(Dt==0)
  tau_periode=(which(Dt==0)[1])-1
  df2_NOj=df2[df2[,1]!=0,] #sans les sauts : Delta t_{j,i} et  Delta y_{j,i}
  djumps=df2[df2[,1]==0,2]  #juste les sauts
  Mrt=matrix(df2_NOj[1:(tau_periode*nb_maint),1],tau_periode,nb_maint)
  Dtauj=apply(Mrt,2,sum) # Delta tau_j + on convertit en vecteur
  Mrt2=matrix(df2_NOj[1:(tau_periode*nb_maint),2],tau_periode,nb_maint)
  Dyj=apply(Mrt2,2,sum) # Delta y_j
  dtji=df2_NOj[,1] ; dyji=df2_NOj[,2]
  list(Dt=Dt, Dy=Dy, df2=df2, nb_maint=nb_maint, tau_periode=tau_periode, 
       df2_NOj=df2_NOj, djumps=djumps, Dtauj=Dtauj, Dyj=Dyj, dtji=dtji, dyji=dyji)
}

LV1<-function(x,df){
  mdf=manip_df(df)
  mu1=x[[1]] ; mu2=mu1 ; sigma21=x[[2]] ; sigma22=x[[3]] ; rho=x[[4]] ; cov=x[[5]]*sqrt(sigma21*sigma22)
  djumps=mdf[[7]]; Dtauj=mdf[[8]]; Dyj=mdf[[9]]; dtji=mdf[[10]]; dyji=mdf[[11]]
  logvrais=0.5*(sum(log(2*pi *sigma21 *dtji )+ ((dyji-mu1*dtji)^2)/(sigma21*dtji))
             +sum(log(2*pi *rho^2 *Dtauj*(sigma22-(cov^2/sigma21)))
                  +((djumps+rho*(mu2*Dtauj+(cov/sigma21)*(Dyj-mu1*Dtauj)))^2)/(rho^2*Dtauj*(sigma22-(cov^2/sigma21)))))
  return(logvrais)
}  

# LV1.corr<-function(x,df){
#   mdf=manip_df(df)
#   mu1=x[[1]] ; mu2=x[[2]] ; sigma21=x[[3]] ; sigma22=x[[4]] ; rho=x[[5]] ; r=x[[6]]
#   djumps=mdf[[7]]; Dtauj=mdf[[8]]; Dyj=mdf[[9]]; dtji=mdf[[10]]; dyji=mdf[[11]]
#   logvrais=0.5*(sum(log(2*pi *sigma21 *dtji )+ ((dyji-mu1*dtji)^2)/(sigma21*dtji))
#                 +sum(log(2*pi *rho^2 *Dtauj*sigma22*(1-r^2))
#                      +((djumps+rho*(mu2*Dtauj+(r*sqrt(sigma22/sigma21))*(Dyj-mu1*Dtauj)))^2)/(rho^2*Dtauj*sigma22*(1-r^2))))
#   return(logvrais)
# }  



estim_init<-function(df){
  
  Dt=manip_df(df)[[1]];   djumps=manip_df(df)[[7]]; Dtauj=manip_df(df)[[8]]; Dyj=manip_df(df)[[9]] ;
  dtji=manip_df(df)[[10]]; dyji=manip_df(df)[[11]]


  mu10<-sum(dyji)/sum(dtji)
  
  sigma210<-sum(((dyji-mu10*dtji)^2)/dtji)/length(dtji)
  #sigma220<-(-((mean(djumps)+rho0*mu10*Dtauj[1])*sqrt(sigma210))/(rho0*r0*(Dyj[length(Dyj)]-mu10*Dtauj[1])))^2
  sigma220<-sigma210-(sigma210/10)
  
  # if(length(djumps)>1){
  #   sigma220<-var(djumps)/((rho0^2)*((1-r0)^2)*Dtauj[1])
  # }else{
  #   sigma220<-0.1
  # }
  # 
  list(mu10=mu10,mu20=mu10,sigma210=sigma210,sigma220=sigma220)
  
}

  
  


est_para<-function(al,ale, mu1,mu2,sig21,sig22,rho,r, n,p,tau_periode,perio,nb_simu){ # flv=fonction log vraisemblance
  bounds <- matrix(c(
    -Inf,Inf,
    0,Inf,
    0,Inf,
    -1,2,
    -0.99,0.99
  ), nc=2, byrow=TRUE)
  colnames(bounds) <- c("lower", "upper")
  
  # Convert the constraints to the ui and ci matrices
  nr <- nrow(bounds)
  ui <- rbind( diag(nr), -diag(nr) )
  ci <- c( bounds[,1], - bounds[,2] )
  
  # Remove the infinite values
  i <- as.vector(is.finite(bounds))
  ui <- ui[i,]
  ci <- ci[i]
  
  
  alea=sort(sample((floor(ale)+1):(ale+nb_simu),size=nb_simu,replace = F))
  Mopt=matrix(0,length(alea),5)
  k=1

  for (i in alea){
    df=ARD1(al,i, mu1,mu2,sig21,sig22,rho,r, n,p,tau_periode,perio)$dfY
    mu10=estim_init(df)$mu10 
    sig210=estim_init(df)$sigma210 ; sig220=estim_init(df)$sigma220
    rho0<-runif(1,0.001,0.999) ; r0<-runif(1,-0.98,0.98)
    opt_ard<-constrOptim(c(mu10,sig210,sig220,rho0,r0),f=LV1,df=df,grad=NULL,ui=ui,ci=ci,method = "Nelder-Mead") #avec contraintes
    #opt_ard<-optim(c(2,2,4,5,0.5,0.5),LV1,df=df) #sans contraintes
    Mopt[k,]<-c(opt_ard$par)
    k=k+1
  }
  return(Mopt)
}


est_bxp<-function(Mopt,par,para_vec){
  para_v=para_vec
  para_c=c("mu1","mu2","sigma21","sigma22","rho","correlation")
  
  gg=ggplot(as.data.frame(Mopt),aes(y=Mopt[,par]))+geom_boxplot(fill="light blue")+
    theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())+theme_bw()+
    theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
    labs(y=para_c[par])+
    geom_hline(yintercept = para_v[par],color="red")
  
      
     clim_mu<-c(para_vec[par]-3,para_vec[par]+3)
     clim_sig<-c(para_vec[par]-5,para_vec[par]+5)
     clim_rho<-c(para_vec[par]-0.5,para_vec[par]+0.5)
     
     if(par==1 |par==2){
     gg<-gg+coord_cartesian(ylim = clim_mu)
     }
     if(par==3 | par==4){
     gg<-gg+coord_cartesian(ylim = clim_sig)
       
     }
     if(par==5){
       gg<-gg+coord_cartesian(ylim = clim_rho)
       
     }
     if(par==6){
       gg<-gg+scale_y_continuous(limits=c(-1,1))
       
     }
   
    #coord_cartesian(ylim = c(0,1))
    #scale_y_continuous(limits=c(-1,1))
   

  return(gg)
}

#Mopt=est_para(al=24,ale=24,mu1=4,mu2=3,sig21=5,sig22=3,rho=0.5,r=0.7, n=20,p=1,tau_periode=3,perio=T,nb_simu=5000)
#est_bxp(Mopt,par=5,para_vec=c(4,3,5,3,0.5,0.7))




temps_seuil<-function(s,nb_simu,al,ale,mu1,mu2,sig21,sig22,rho,r, n,p,tau_periode,perio){
   st<-NULL ; t<-NULL ;alea<- sort(sample((floor(ale)+1):(ale+nb_simu),size=nb_simu,replace = F))
  for (i in 1:nb_simu){
    df=ARD1(al,ale=alea[i],mu1,mu2,sig21,sig22,rho,r, n,p,tau_periode,perio)$dfY
    if(sum(df$Degradation >= s)!=0){
      t<-df$Temps[df$Degradation>=s][1]
      st<-c(st,t)
    }
  }
  return(st)
}
 

RUL1 <-function(mu1,mu2,sig21,sig22,rho,r,p,tau_periode,s) {
  df <- data.frame(Temps = 0, Degradation = 0)
  time <- df[, 1]
  deg <- df[, 2]
  incX <- c(0, 0)
  incY <- 0
  Xt <- incX
  t<-0
  j <- 2
  tau_p<-tau_periode
  while ((deg[length(deg)] < s) & (j<10000)) {
    covar <- r * sqrt(sig21 * sig22)
    COV <- matrix(c(sig21 * p, rep(covar, 2) * p, sig22 * p), 2, 2)
    t <- t + p
    inc <- rockchalk::mvrnorm(1, mu = c(mu1 * p, mu2 * p), Sigma = COV)
    incY <- incY + inc[1]
    incX <- incX + inc
    Xt <- rbind(Xt, incX)
    df[j, ] <- c(t, incY)
    if (as.character(t)==as.character(p * tau_p) & (incY<s)) { #temps de maintenance
      j <- j + 1
      incY <- Xt[dim(Xt)[1], 1] - rho * Xt[dim(Xt)[1], 2]
      df[j, ] <- c(t, incY)
      tau_p<-tau_p+tau_periode
    }
    time <- df[, 1]
    deg <- df[, 2]
    j <- j + 1
  }
  return (t-p)
}




temps_seuil2<-function(mu1,mu2,sig21,sig22,rho,r,p,tau_periode,s,nb_simu){
  ts<-c()
  if(nb_simu>0){
    withProgress(message = "Making Plot" , {
      n<-nb_simu
    for(i in 1:nb_simu){
     ts<-c(ts,RUL1(mu1,mu2,sig21,sig22,rho,r,p,tau_periode,s))
     incProgress((nb_simu*10^(-2))/nb_simu)
    }
      
    })
  }
  return(ts)
}


plot_rul<-function(TS,nb_simu,liss){
  lts<-length(TS)
ggplot(data=data.frame(rul=TS),aes(x=rul))+geom_histogram(aes(y=..density..),bins=1+ceiling(log(nb_simu)/log(2)),color="black",fill="light grey")+
  geom_density(alpha=0.3,fill="#FF6666",adjust=liss)+
  theme_bw()+ggtitle(paste0("Nombre de trajectoires : ",as.character(lts)))

}


plot_rul2<-function(TS,nb_simu,liss,bins2){
  lts<-length(TS)
  ggplot(data=data.frame(rul=TS),aes(x=rul))+geom_histogram(aes(y=..density..),bins=bins2,color="black",fill="light grey")+
    geom_density(alpha=0.3,fill="#FF6666",adjust=liss)+
    theme_bw()+ggtitle(paste0("Nombre de trajectoires : ",as.character(lts)))
}


plotARD1 <- function(df,taus,Xt,t,mod) {
  
  gg = ggplot(df, aes(x = Temps, y = Degradation)) + geom_point() + geom_line() +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = "Time") + scale_x_continuous(breaks = taus)
  
  if(mod){
    vec_traj<-rep(1:2,each=length(t))
    couleur<-c( "#e06cf0","#56B4E9")
  }else{vec_traj<-rep(c("a","b"),each=length(t))
    couleur<-c( "#E69F00", "#56B4E9")
  }
  
  df_Xt<-data.frame(Temps=rep(t,2),Degradation=c(Xt),
                    traj=vec_traj)
  
    df_Xt$traj<-as.factor(df_Xt$traj)
    gg=gg+geom_line(data=df_Xt,aes(x=Temps,y=Degradation,color=traj),alpha=0.2)+
      geom_point(data=df_Xt,aes(x=Temps,y=Degradation,color=traj),alpha=0.2)+
      scale_color_manual(name = "X(t)",values=couleur)
    
      
  
  #gg=gg+scale_color_discrete(name = "X(t)")
  return(gg)
}



###############
#### CAS 4 ####
###############

Nb_maint4<-function(n,tau_periode){
  return(floor((n-2)/(tau_periode-1)))
}

traj<-function(df){
  t<-df[,1]
  nb_maint=length(which(diff(t)==0))
  n=length(t)
  tau_periode=floor((n/(nb_maint+1))-1)
  if(n==((tau_periode+1)*(nb_maint+1))){ #on enleve aussi la derniere valeur si besoin
    new_df<-df[-c(which(diff(t)==0),which(diff(t)==0)+1,length(t)),]
  }else{
    new_df<-df[-c(which(diff(t)==0),which(diff(t)==0)+1),]
  }
  return(new_df)
}




plotARD1_4 <- function(df1,taus,Xt,t,mod) { #df1 doit correspondre au df du cas 1
  dfNA<-df1
  taus.b<-taus[-c(1,length(taus))]
  ind_t<-pmatch(rep(taus.b,each=2),df1[,1])
  ind_t<-ind_t[!is.na(ind_t)]
  dfNA[ind_t,2]<-NA
  if( (t[length(t)])==(taus.b[length(taus.b)]+unique(diff(taus[-length(taus)])))){
    dfNA[dim(df1)[1],2]<-NA} # derniere valeur en NA si tfinal instant de maint
  gg = ggplot(df1, aes(x = Temps, y = Degradation))+ geom_line(linetype="dotted") +
    geom_line(data=dfNA)+geom_point(data=dfNA)+
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = "Time") + scale_x_continuous(breaks = taus)
  
  if(mod){
    vec_traj<-rep(1:2,each=length(t))
  }else{vec_traj<-rep(c("a","b"),each=length(t))}
  
  df_Xt<-data.frame(Temps=rep(t,2),Degradation=c(Xt),
                    traj=vec_traj)
  
  df_Xt$traj<-as.factor(df_Xt$traj)
  gg=gg+geom_line(data=df_Xt,aes(x=Temps,y=Degradation,color=traj),alpha=0.2)+
    geom_point(data=df_Xt,aes(x=Temps,y=Degradation,color=traj),alpha=0.2)
  
  
  gg=gg+scale_color_discrete(name = "X(t)")
  
  return(gg)
}


manip_df4<-function(df){
  Dt=diff(df[,"Temps"]) #Delta t_{j,i} et temps des sauts
  Dy=diff(df[,"Degradation"]) #Delta y_{j,i} et sauts
  p<-Dt[1]
  df2=data.frame(DTemps=Dt,DDegradation=Dy)
  nb_maint=sum((Dt)==(max(Dt))) #si perio, ne marche qu avec un operateur comme + comprends pas pq...
  tau_periode=which((Dt)==(max(Dt)))[1]
  df2_NOj=df2[(df2[,1])!=(max(Dt)),] #sans les sauts : Delta t_{j,i} et  Delta y_{j,i}
  djumps=df2[(df2[,1])==(max(Dt)),2]  #juste les sauts
  Dtauj=rep(p*tau_periode,nb_maint) 
  Mrt1=matrix(df2_NOj[1:(tau_periode-1),2],tau_periode-1,1)# \Delta y_{0,i}
  Mrt2=matrix(df2_NOj[seq(tau_periode,length.out=(nb_maint-1)*(tau_periode-2)),2],tau_periode-2,nb_maint-1)
  Dyj=c(sum(Mrt1),apply(Mrt2,2,sum)) # Delta y_j
  dtji=df2_NOj[,1] ; dyji=df2_NOj[,2]
  list(Dt=Dt, Dy=Dy, df2=df2, nb_maint=nb_maint, tau_periode=tau_periode, 
       df2_NOj=df2_NOj, djumps=djumps, Dtauj=Dtauj, Dyj=Dyj, dtji=dtji, dyji=dyji, p=p)
}

LV4<-function(x,df){
  mdf=manip_df4(df)
  mu1=x[[1]] ; mu2=x[[2]] ; sigma21=x[[3]] ; sigma22=x[[4]] ; rho=x[[5]] ; cov=x[[6]]*sqrt(sigma21*sigma22)
  djumps=mdf[[7]]; Dtauj=mdf[[8]]; Dyj=mdf[[9]]; dtji=mdf[[10]]; dyji=mdf[[11]] ; p=mdf[[12]]
  
  nb_maint<-length(djumps)
  pvec<-rep(p,nb_maint)
  Dt_obs<-c(Dtauj[1]-p,Dtauj[-1]-2*p)
  
  mu.zj<-mu1*2*pvec-rho*mu2*Dtauj #length=nb_maint
  sig2.zj<-sigma21*2*pvec+rho^2*sigma22*Dtauj-2*rho*cov*pvec
  

  logvrais=(1/2)*(sum(log(2*pi *sigma21 *dtji )+((dyji-mu1*dtji)^2)/(sigma21*dtji))
             +sum(log(2*pi *(sig2.zj-rho^2 *cov^2*((Dt_obs/sigma21)+(c(0,(pvec[-1])^2)/c(1,sig2.zj[-(length(sig2.zj))]))))) #1 au pif
                  +(djumps-mu.zj+((rho*cov)/sigma21)*(Dyj-mu1*Dt_obs)
                 +(rho*cov*c(0,pvec[-1])/c(1,sig2.zj[-(length(sig2.zj))]))*(c(0,djumps[-length(djumps)])-c(0,mu.zj[-length(mu.zj)])))^2/(sig2.zj-rho^2 *cov^2*((Dt_obs/sigma21)+((c(0,pvec[-1])^2)/c(1,sig2.zj[-length(sig2.zj)]))))))
  return(logvrais)
}  




est_para4<-function(ale, mu1,mu2,sig21,sig22,rho,r, n,p,tau_periode,perio,nb_simu){ # flv=fonction log vraisemblance

  bounds <- matrix(c(
    -Inf,Inf,
    -Inf,Inf,
    0,Inf,
    0,Inf,
    -1,2,
    -0.99,0.99
  ), nc=2, byrow=TRUE)
  colnames(bounds) <- c("lower", "upper")
  
  # Convert the constraints to the ui and ci matrices
  nr <- nrow(bounds)
  ui <- rbind( diag(nr), -diag(nr) )
  ci <- c( bounds[,1], - bounds[,2] )
  
  # Remove the infinite values
  i <- as.vector(is.finite(bounds))
  ui <- ui[i,]
  ci <- ci[i]
  
  
  alea=sort(sample((floor(ale)+1):(ale+nb_simu),size=nb_simu,replace = F))
  Mopt=matrix(0,length(alea),6)
  k=1
  
  for (i in alea){
    df0=ARD1(al,i, mu1,mu2,sig21,sig22,rho,r, n,p,tau_periode,perio)$dfY
    df4<-traj(df0)
    opt_ard<-constrOptim(c(mu1,mu2,sig21,sig22,rho,r),f=LV4,df=df4,grad=NULL,ui=ui,ci=ci,method = "Nelder-Mead") #avec contraintes
    #opt_ard<-optim(c(2,2,4,5,0.5,0.5),LV1,df=df) #sans contraintes
    Mopt[k,]<-c(opt_ard$par)
    k=k+1
  }
  return(Mopt)
}


est_bxp4<-function(Mopt,par,para_vec){
  para_c=c("mu1","mu2","sigma21","sigma22","rho","correlation")
  
  gg=ggplot(as.data.frame(Mopt),aes(y=Mopt[,par]))+geom_boxplot(fill="light blue")+
    theme(panel.grid.major.x =element_blank(),panel.grid.minor.x=element_blank())+theme_bw()+
    theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
    labs(y=para_c[par])+
    geom_hline(yintercept = para_vec[par],color="red")
    
  clim_mu<-c(para_vec[par]-3,para_vec[par]+3)
  clim_sig<-c(para_vec[par]-5,para_vec[par]+5)
  clim_rho<-c(para_vec[par]-0.5,para_vec[par]+0.5)
  
  if(par==1 |par==2){
    gg<-gg+coord_cartesian(ylim = clim_mu)
  }
  if(par==3 | par==4){
    gg<-gg+coord_cartesian(ylim = clim_sig)
    
  }
  if(par==5){
    gg<-gg+coord_cartesian(ylim = clim_rho)
    
  }
  if(par==6){
    gg<-gg+scale_y_continuous(limits=c(-1,1))
    
  }
  
  return(gg)

  }






shinyServer(function(input, output,session) {
  val<-reactiveValues(al=24,ale=24)
  
    
    output$p1<-renderUI({
      if(input$model=="X^{(1)}, X^{(2)}"){
        numericInput(
          inputId = "mu1",
          label = "mu_1",
          value = 4
        )
      }else{
        numericInput(
          inputId = "mua",
          label = "mu_a",
          value = 4)
      }
    })

    output$p2<-renderUI({
      if(input$model=="X^{(1)}, X^{(2)}"){
        numericInput(
          inputId = "mu2",
          label = "mu_2",
          value = 2
        )
      }else{
        numericInput(
          inputId = "mub",
          label = "mu_b",
          value = 2)
      }
    })

    output$p3<-renderUI({
      if(input$model=="X^{(1)}, X^{(2)}"){
        numericInput(
          inputId = "var1",
          label = "sigma2_1",
          value = 7
        )
      }else{
        numericInput(
          inputId = "sig2a",
          label = "sigma2_a",
          value = 7)
      }
    })

    output$p4<-renderUI({
      if(input$model=="X^{(1)}, X^{(2)}"){
        numericInput(
          inputId = "var2",
          label = "sigma2_2",
          value = 3
        )
      }else{
        numericInput(
          inputId = "sig2b",
          label = "sigma2_b",
          value = 3)
      }
    })

    output$p5<-renderUI({
      if(input$model=="X^{(1)}, X^{(2)}"){

        numericInput(
          inputId = "corr",
          label = "r",
          value = 0.7,min=-1,max=1,step=0.1)

      }else{
        numericInput(
          inputId = "corrab",
          label = "corr",
          #value = (sqrt(input$var1*input$var2))/(input$var1+2*input$var2-2*input$corr*sqrt(input$var1*input$var2)))
          value = 0.7,min=-1,max=1,step=0.1)

          }
    })


    observeEvent(input$button, {
        b<-val$ale+input$nbsimu
        val$ale<-b
        
    })
      
    output$myplot <- renderPlotly({
      validate()
      
      if (input$model == "X^{(1)}, X^{(2)}") {
        
        if (input$but == "CAS 1") {
          ard1 <-
            ARD1(
              val$al,
              val$ale,
              input$mu1,
              input$mu2,
              input$var1,
              input$var2,
              input$rho,
              input$corr,
              input$n,
              input$p,
              input$tau_periode,
              perio = T
            )
          df <- ard1$dfY
          taus <- ard1$taus
          Xt <- ard1$Xt
          t <- ard1$t
          Xtab = Xtab.bis(Xt)
          
          plotARD1(df, taus, Xt, t, mod = T)
        } else{
          ard1 <-
            ARD1(
              val$al,
              val$ale,
              input$mu1,
              input$mu2,
              input$var1,
              input$var2,
              input$rho,
              input$corr,
              input$n+2*Nb_maint4(input$n,input$tau_periode),
              input$p,
              input$tau_periode,
              perio = T
            )
          df <- ard1$dfY
          taus <- ard1$taus
          Xt <- ard1$Xt
          t <- ard1$t
          Xtab = Xtab.bis(Xt)
          
          plotARD1_4(df, taus, Xt, t, mod = T)
        }
      } else{
       
        if (input$but == "CAS 1") {
          Xt <-
            Xtab(
              val$ale,
              input$mua,
              input$mub,
              input$sig2a,
              input$sig2b,
              input$corrab,
              input$n,
              input$p,
              input$tau_periode
            )
          ard1ab <-
            ARD1_AB(Xt, input$p, input$rho, input$tau_periode, input$n)
          df <- ard1ab$dfY
          taus <- ard1ab$taus
          t <- ard1ab$t
          
          
          plotARD1(df, taus, Xt, t, mod = F)
        } else{
          
          Xt <-
            Xtab(
              val$ale,
              input$mua,
              input$mub,
              input$sig2a,
              input$sig2b,
              input$corrab,
              input$n+2*Nb_maint4(input$n,input$tau_periode),
              input$p,
              input$tau_periode
            )
          ard1ab <-
            ARD1_AB(Xt, input$p, input$rho, input$tau_periode, input$n)
          df <- ard1ab$dfY
          taus <- ard1ab$taus
          t <- ard1ab$t
          
          
          plotARD1_4(df, taus, Xt, t, mod = F)
          
        }
      }
    })
    
    output$tablepara<-renderTable({
      if(input$but == "CAS 1"){
        if(input$model=="X^{(1)}, X^{(2)}"){
        Mp<-data.frame(mu_a=input$mu1-input$mu2,mu_b=input$mu2,sigma2_a=input$var1+input$var2-2*(input$corr*sqrt(input$var1*input$var2)),
                     sigma2_b=input$var2,r_ab=(input$corr*sqrt(input$var1)-sqrt(input$var2))/sqrt(input$var1+input$var2-2*input$corr*sqrt(input$var1*input$var2)))
        withMathJax() 
        #colnames(Mp)<-c("$$\\mu_a$$","$$\\mu_b$$","$$\\sigma_a^2$$","$$\\sigma_b^2$$","$$r_{ab}$$")
        Mp
        }else{
          data.frame(mu1=input$mua+input$mub,mu2=input$mub,sigma21=input$sig2a+input$sig2b+2*(input$corrab*sqrt(input$sig2a*input$sig2b)),
                         sigma22=input$sig2b,r_12=(input$corrab*sqrt(input$sig2a)+sqrt(input$sig2b))/sqrt(input$sig2a+input$sig2b+2*input$corrab*sqrt(input$sig2a*input$sig2b)))
          
        }
      }
    })
    
    output$myplotAB<-renderPlotly({
      validate()
      if(input$but == "CAS 1"){
        if(input$model=="X^{(1)}, X^{(2)}"){
          ard1<-ARD1(val$al,val$ale,input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr, input$n,input$p,input$tau_periode,perio=T)
          Xt<-ard1$Xt ; Xtab<-Xtab.bis(Xt)
          ard1.ab<-ARD1_AB(Xtab,input$p, input$rho,input$tau_periode,input$n)
          plotARD1(ard1.ab$dfY,ard1.ab$taus,Xtab,ard1.ab$t,mod=F)
        }else{
          Xt0<-Xtab(val$ale,input$mua,input$mub,input$sig2a,input$sig2b,input$corrab, input$n,input$p,input$tau_periode)
          Xt<-Xt12.bis(Xt0)
          ard1ab<-ARD1_12bis(Xt,input$p, input$rho,input$tau_periode,input$n)
          df<-ard1ab$dfY ; taus<-ard1ab$taus ; t<-ard1ab$t
          plotARD1(ard1ab$dfY,ard1ab$taus,Xt,ard1ab$t,mod=T)
          
        }
        
      }
    })
    
    output$table<-renderTable({
      if(input$but == "CAS 1"){
        M1<-data.frame(Maintenances=round(Nb_maint(input$n,input$tau_periode),0),n_j=input$tau_periode-1)
        colnames(M1)<-c("Nombre de maintenances","n_j")
        M1
      }else{
        
        M1<-data.frame(Maintenances=floor((input$n-2)/(input$tau_periode-1)),n_j=input$tau_periode-1)
        colnames(M1)<-c("Nombre de maintenances","n_j")
        M1
      }
    })
    
    val2<-reactiveValues(al=24,ale=24)
    
    observeEvent(input$button2, {
      b<-val2$ale+input$Nbtr
      val2$ale<-b
      
    })  
    
    
        
    output$cov<-renderTable({
      if(input$model=="X^{(1)}, X^{(2)}"){
      data.frame(cov=input$corr*sqrt(input$var1*input$var2))}
      else{
          data.frame(cov=input$corrab*sqrt(input$sig2a*input$sig2b))
      }
    })
      

    output$density<-renderPlotly({
      validate()
      TS<-temps_seuil(input$s,input$Nbtr,val2$ale,val2$ale,
                      input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr, 
                      input$n,input$p,input$tau_periode,perio=T)
      plot_rul(TS=TS,input$Nbtr,input$lissage)
      
    })
    
    val3<-reactiveValues(ale=1)
    observeEvent(input$button3, {
      e<-val3$ale+input$Nbtr2
      val3$ale<-e
      
    }) 
    
    # output$TS2<-renderTable({
    #   TS<-temps_seuil2(
    #     input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr, 
    #     input$p,input$tau_periode,input$s2,input$Nbtr2)
    #   data.frame(TS=TS)
    #   
    # })
    
    output$density2<-renderPlotly({
      set.seed(val3$ale)
      
      TS<-temps_seuil2(
                      input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr,
                      input$p,input$tau_periode,input$s2,input$Nbtr2)

      
      plot_rul2(TS=TS,input$Nbtr2,input$lissage2,input$bins2)
      
      
      
    })
    

    

    
    
    
    
    # observeEvent(input$button3, {
    #   last.seed=temps_seuil2(val3$ale,
    #                        input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr, 
    #                        input$p2,tau_periode=1,input$s2,input$Nbtr2)$alea
    #   c3<-val3$ale+last.seed+1
    #   val3$ale<-c3
    #   
    # })  
    
    
    # output$density2<-renderPlotly({
    #   validate()
    #   TS<-temps_seuil2(val3$ale,
    #                   input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr, 
    #                 input$p2,tau_periode=1,input$s2,input$Nbtr2)$trul
    #   plot_rul2(TS=TS,input$Nbtr2,input$lissage)
    # })
    #   
    
    
    val4<-reactiveValues(al=24,ale=24)
    
    observeEvent(input$button4, {
      d<-val4$ale+input$nbsimu
      val4$ale<-d
      
    })  
    
    
    
    output$bxARD1<-renderPlot({
      validate()
      if(input$but=="CAS 1"){
  
      Mopt<-est_para(val4$al,val4$ale,input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr, input$n,input$p,input$tau_periode,perio=T,input$nbsimu)
      
      p1<-est_bxp(Mopt,1,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
      p2<-est_bxp(Mopt,2,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
      p3<-est_bxp(Mopt,3,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
      p4<-est_bxp(Mopt,4,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
      p5<-est_bxp(Mopt,5,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
      p6<-est_bxp(Mopt,6,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
      grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2, nrow = 3)
      }else{
        Mopt4<-est_para4(val4$ale,input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr, input$n,input$p,input$tau_periode,perio=T,input$nbsimu)
        
        p1<-est_bxp4(Mopt4,1,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
        p2<-est_bxp4(Mopt4,2,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
        p3<-est_bxp4(Mopt4,3,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
        p4<-est_bxp4(Mopt4,4,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
        p5<-est_bxp4(Mopt4,5,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
        p6<-est_bxp4(Mopt4,6,c(input$mu1,input$mu2,input$var1,input$var2,input$rho,input$corr))
        grid.arrange(p1,p2,p3,p4,p5,p6,ncol=2, nrow = 3)
      }
    })
    

}
)


