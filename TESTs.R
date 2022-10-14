###############
#### TESTS ####
###############

 setwd("C:/Users/Margaux/OneDrive/Documents/DEGRADATION MULTIVARIEE/Inference/Applis Shiny/New_ARD1_3")
 source("server.R")

##### REPRESENTATION GRAPHIQUE ####


al=24;ale=25;mu1=5;mu2=3;sig21=9;sig22=4;rho=0.5;r=0.7;n=20;p=1;tau_periode=3;perio=T
# #df=myDF[,2:3]
# setwd("C:/Users/Margaux/OneDrive/Documents/DEGRADATION MULTIVARIEE/Inference/Applis Shiny/New_ARD1_3")
# source("server.R")
# mua=mu1-mu2 ; mub=mu2
# sig2a=sig21+sig22-2*(r*sqrt(sig21*sig22))
# sig2b=sig22
# corrab=(r*sqrt(sig21)-sqrt(sig22))/(sqrt(sig21+sig22-2*r*sqrt(sig21)*sqrt(sig22)))
#  
# ARD1_AB(al=24,ale=24,mua,mub,sig2a,sig2b,rho,corrab , n=20,p=1,tau_periode=3,perio=T)$Xt
 Xt=Xtab.bis(ard1$Xt)
#ard1ab=ARD1_AB(Xt,p=1,rho=0.5,tau_periode=3,n=200)
ard1=ARD1(al,ale=17,mu1=5,mu2=5,sig21=10,sig22=7,rho=0.5,r=0.7, n=20,p=1,tau_periode=3,perio=T)

df=ard1$dfY  
taus=ard1$taus
t=ard1$t
manip_df(df)


plotARD1(df,taus,Xt12.bis(ard1$Xt),t,mod=F)

TS=temps_seuil2(40,30,10,7,0.8,0.7,0.001,50,2.6,1)

plotARD1.rul <- function(df,taus,Xt,t,mod,seuil,ts) {
  
  gg = ggplot(df, aes(x = Temps, y = Degradation))  + geom_line() +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = "Time") + scale_x_continuous(breaks = as.numeric(sort(c(ts,taus[-length(taus)]))))+
    scale_y_continuous(breaks=sort(c(seuil,seq(0,12,4))))+
    theme(axis.text.x = element_text(face=NULL, size=10, angle=90))+
  
    geom_segment(x=0,y=seuil,xend=ts,yend=seuil,color="red",linetype="dashed")+ 
    geom_segment(x=ts,y=0,xend=ts,yend=seuil,color="red",linetype="dashed")

  return(gg)
}

plotARD1.rul(ard1$df,ard1$taus,ard1$Xt,ard1$t,mod=T,seuil=2.6,ts=0.092)


LV1(c(2,2,10,9,0,0.7),df=ard1$dfY)
LV1(c(4,4,10,36,0,0.7),df=ard1$dfY)


LV1(c(mu1=5,mu2=3,sig21=9,sig22=4,rho=0.5,r=0.7),df1)



corrab=0.8 ;sig2a=15 ; sig2b=10 ; mua=4 ; mub=3
Xt=Xtab.bis(ard1$Xt)
 
#Xt=Xt12.bis(Xtab(ale=4,mua,mub,sig2a,sig2b,corrab, n,p,tau_periode))
#ard=ARD1_12bis(Xt,p,rho,tau_periode,n)

Xt=Xtab(ale=30,mua,mub,sig2a,sig2b,corrab, n,p,tau_periode)


plotARD1 <- function(df,taus,Xt,t,mod) {
  
  gg = ggplot(df, aes(x = Temps, y = Degradation))  + geom_line() +
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
    scale_color_manual(name = "X(t)",values=couleur)
    #scale_y_continuous(limits=c(0,58.55))
  
  #gg=gg+scale_color_discrete(name = "X(t)")
  return(gg)
}



#ard=ARD1_AB(Xt,p,rho=0.6,tau_periode,n)
plotARD1(ard1$dfY,ard1$taus,ard1$Xt,ard1$t,mod=T)
plotARD1_4(ard1$dfY,ard1$taus,ard1$Xt,ard1$t,mod=T)

# 
# 
# traj(ard$dfY)


############################
#### INFERENCE CAS 1 #######
############################
mv<-NULL
for (i in sample(100:100000,100)){
ard1=ARD1(al,ale=i,mu1=5,mu2=5,sig21=10,sig22=1,rho=0.8,r=-0.9, n=500,p=1,tau_periode=6,perio=T)
df=ard1$dfY
mv<-c(mv,(manip_df(df)[[7]]))
}
print(var(mv)/((0.8^2)*6*(1-(-0.9^2))))
print((var(mv)/((0.8^2)*6)))




Mopt6=est_para(al=24,ale=24,mu1=5,mu2=5,sig21=10,sig22=7,rho=0.8,r=0.7, n=20,p=1,tau_periode=3,perio=T,nb_simu=5000)
save(Mopt4, file="Mopt4bis.RData")

Mopt=est_para(al=24,ale=24,mu1=5,mu2=5,sig21=10,sig22=7,rho=0.5,r=0.7, n=20,p=1,tau_periode=3,perio=T,nb_simu=50,rho0=0.7,r0=0.5)

boxplot(Mopt6[,5])
abline(h=0.7,col="red")



# M=rbind(Mopt1,Mopt2)
# estim_bxp2<-function(M,par){
#   L=list()
#   for(i in 1:5){ # car 5 parametres
#     estim_para<-c(Mopt1[,i],Mopt2[,i])
#     L[[i]]<-data.frame(para=estim_para,situation=rep(1:2,each=nb_simu)) #nb de situations
#   }
#   ggplot(dfmu,aes(x=as.factor(situation),y=mu))+geom_boxplot(fill="light blue")+
#     labs(x="Situation")+theme_bw()+geom_hline(yintercept = para_v[par],color="red")
#   
# }

bxp_par(par,nb_simu,valtheo)

bxp_par<-function(par,nb_simu,valtheo){ #par numero du parametre entre 1 et 5 valtheo valeur theorique du parametre
  estim<-c(Mopt1[,par],Mopt2[,par],Mopt3[,par],Mopt4[,par],Mopt5[,par],Mopt6[,par]) #1 situation par Mopt
  df=data.frame(para=estim,situation=rep(1:6,each=nb_simu)) #rep nb de situations
  char_para<-c("mu","sigma_1^2","sigma_2^2","rho","r_{12}")
  
  ggplot(df,aes(x=as.factor(situation),y=para))+geom_boxplot(fill="light blue")+
  labs(x="Situation",y=char_para[par])+theme_bw()+
  geom_hline(yintercept=valtheo,color="red") 
  #scale_y_continuous(breaks=c(0,0.2,0.5,0.8,1))+
  #coord_cartesian(ylim=c(0,30))
   #geom_segment(x=0,y=valtheo[1],xend=4.5,yend=valtheo[1],color="red")+
  #geom_segment(x=4.5,y=valtheo[2],xend=5.5,yend=valtheo[2],color="red")+
    #geom_segment(x=5.5,y=valtheo[3],xend=6.6,yend=valtheo[3],color="red")
}

bxp_par(5,nb_simu=5000,valtheo=0.7)


bxp_par(5,5000,0.7)

est_bxp(Mopt1,par=5,para_vec=c(5,10,7,0.5,0.7))

for(i in c(20,200)){
  if(i==20){
    Mopt=est_para(al=24,ale=24,mu1=4,mu2=3,sig21=5,sig22=3,rho=0.5,r=0, n=i,p=1,tau_periode=3,perio=T,nb_simu=5000)
    save(Mopt,file="Mopt_n20_corr0_tp2.RData")
  }else{
    Mopt=est_para(al=24,ale=24,mu1=4,mu2=3,sig21=5,sig22=3,rho=0.5,r=0, n=i,p=1,tau_periode=3,perio=T,nb_simu=5000)
    save(Mopt,file="Mopt_n200_corr0_tp2.RData")
  }
}


###############
#### RUL ######
###############

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
      if (as.character(t)==as.character(p * tau_p) & (incY<s)) {
        j <- j + 1
        incY <- Xt[dim(Xt)[1], 1] - rho * Xt[dim(Xt)[1], 2]
        df[j, ] <- c(t, incY)
        tau_p<-tau_p+tau_periode
      }
      time <- df[, 1]
      deg <- df[, 2]
      j <- j + 1
    }
    #return (df[-(dim(df)[1]),])
    return (t-p)
}

plot(RUL1(40,30,7,3,1,0.9,0.001,50,2.6),type="l")


ts<-RUL1(40,30,7,3,1,0.9,0.001,50,2.6)
set.seed(10)

temps_seuil2<-function(mu1,mu2,sig21,sig22,rho,r,p,tau_periode,s,nb_simu){
  ts<-c()
  if(nb_simu>0){
  for(i in 1:nb_simu){
    ts<-c(ts,RUL1(mu1,mu2,sig21,sig22,rho,r,p,tau_periode,s))
    }
  }
  return(ts)
}

TS=temps_seuil2(40,30,10,7,0.8,0.7,0.001,50,2.6,5000)

hist(ts,breaks=1+ceiling(log(5000)/log(2)))

hist(TS,breaks=30)

plot_rul3<-function(TS,liss,bins2){
  lts<-length(TS)
  ggplot(data=data.frame(rul=TS),aes(x=rul))+geom_histogram(aes(y=..density..),bins=bins2,color="black",fill="light grey")+
    geom_density(alpha=0.3,fill="#FF6666",adjust=liss)+
    theme_bw()+
    scale_x_continuous(breaks=seq(0,max(TS),by=tau_periode*p))
    
  }

TS=RUL1(40,30,15,3,1,0.9,0.001,50,2.6)

plot_rul3(TS,1.6,15)



sort(TS)+((1:nb_simu)*10^{-7})
nbrec=10
sort(TS)[c(seq(1,by=floor(length(TS)/nbrec),length(TS)),length(TS))]

plot_rul2<-function(TS,liss,nbrec, tau_periode,p){
  nb.rec=nbrec-2
  TS<-sort(TS)+((1:nb_simu)*10^{-7})
  lts<-length(TS)
  ggplot(data=data.frame(rul=TS),aes(x=rul))+geom_histogram(aes(y=..density..),breaks=sort(TS)[c(seq(1,by=floor(length(TS)/nbrec),length(TS)),length(TS))],color="black",fill="light grey")+
    geom_density(alpha=0.3,fill="#FF6666",adjust=liss)+
    theme_bw()+
    scale_x_continuous(breaks=seq(0,max(TS),by=tau_periode*p))
  #+ggtitle(paste0("Nombre de trajectoires : ",as.character(lts)))
}


#breaks=seq(0,max(TS),by=tau_periode*p)
plot_rul2(TS2,1.7,10,50)
TS=temps_seuil2(40,30,7,3,1,0.9,0.001,50,2.6,5000)

plot_rul2(TS,liss=1.7,12)

sort(TS)
TS2<-TS+p*0.1*runif(length(TS),0,1)
length(c(0,TS[seq(1,by=round(500/8),length(TS))],max(TS)))

new.seq[1]<-1
hist(TS2,breaks=sort(TS2)[c(seq(1,by=floor(length(TS2)/10),length(TS2)),length(TS2))])
boxplot(TS)

hist(TS,nclass=50)
hist(TS*3,nclass=33)


coup<-sort(TS2)[c(seq(1,by=floor(length(TS2)/15),length(TS2)),length(TS2))]

len<-NULL
for(i in 1:(length(coup)-1)){
len<-c(len,length(sort(TS2)[coup[i]<=sort(TS2) & sort(TS2)<coup[(i+1)]]))
}
len

length(c(seq(1,by=floor(length(TS)/10),length(TS)),length(TS)))

############################
#### INFERENCE CAS 4 #######
############################

df=ARD1(al=24,ale=24,mu1=5,mu2=3,sig21=9,sig22=4,rho=0.5,r=0.7, n=20,p=1,tau_periode=3,perio=T)$dfY 
df4<-traj(df)

ARD1(al=24,ale=24,mu1=4,mu2=3,sig21=5,sig22=3,rho=0.5,r=0.7, n=3000,p=0.0001,tau_periode=999,perio=T)

LV4(c(5,3,9,4,0.5,0.7),df4)
LV1(c(4,3,5,3,0,0.7),df)
nb_maint=floor((n-2)/(tau_periode-1))

manip_df(df)$djumps
manip_df4(df4)

for(i in c(20,200)){
  if(i==20){
    Mopt4=est_para4(ale=24,mu1=4,mu2=3,sig21=5,sig22=3,rho=0.5,r=0, n=i,p=1,tau_periode=3,perio=T,nb_simu=5000)
    save(Mopt4,file="Mopt4_equiv_n20_cas1_corr0.RData")
  }else{
    Mopt4=est_para4(ale=24,mu1=4,mu2=3,sig21=5,sig22=3,rho=0.5,r=0, n=i,p=1,tau_periode=3,perio=T,nb_simu=5000)
    save(Mopt4,file="Mopt4_equiv_n200_cas1_corr0.RData")
  }
}

Mopt4=est_para4(ale=24,mu1=4,mu2=3,sig21=5,sig22=3,rho=0.5,r=0, n=20,p=1,tau_periode=3,perio=T,nb_simu=200)

for(i in c(20,200)){
  if(i==0){
    Mopt4=est_para4(ale=24,mu1=4,mu2=3,sig21=5,sig22=3,rho=0.5,r=0, n=1000,p=1,tau_periode=19,perio=T,nb_simu=5000)
    save(Mopt4,file="Mopt4_n1000_corr0_tp19.RData")
  }else{
    Mopt4=est_para4(ale=24,mu1=4,mu2=3,sig21=5,sig22=3,rho=0.5,r=i, n=1000,p=1,tau_periode=19,perio=T,nb_simu=5000)
    save(Mopt4,file="Mopt4_n1000_corr07_tp19.RData")
  }
}

est_bxp4(Mopt4,6,c(4,3,5,3,0.5,0))

print(RUL1(4,3,7,6,rho=0.5,r=0.7,p=0.1,tau_periode=3,s=10))


