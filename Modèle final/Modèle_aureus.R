
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

#Code modèle page

Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+IRa+ISa+S+CR+CS
    
    
    new_teta<-teta-log(1-vec_virus(t)*ATB)
    
    
    dSa <- -Sa*((betaR*ct*(CRa+IRa+CR)/N)+betaS*(CSa+ISa+CS)/N)-omega*Sa+new_teta*S+(gammaR+alpha*(1-sigmaR))*CRa+(gammaS+alpha*(1-sigmaS))*CSa
    dCRa <- Sa*(betaR*ct*(CRa+IRa+CR)/N)-(gammaR+alpha*(1-sigmaR))*CRa-rhoR*CRa-omega*CRa+new_teta*CR
    dCSa <- Sa*(betaS*(CSa+ISa+CS)/N)-(gammaS+alpha*(1-sigmaS))*CSa-rhoS*CSa-omega*CSa+new_teta*CS
    dIRa <- rhoR*CRa-deltaRa*IRa+rhoR*CR
    dISa <- rhoS*CSa-deltaSa*ISa+rhoS*CS
    
    dS <- -S*((betaR*ct*(CRa+IRa+CR)/N)+betaS*(CSa+ISa+CS)/N)+omega*Sa-new_teta*S+gammaR*CR+gammaS*CS+deltaRa*IRa+deltaSa*ISa
    dCR <- S*(betaR*ct*(CRa+IRa+CR)/N)-gammaR*CR-rhoR*CR+omega*CRa-new_teta*CR
    dCS <- S*(betaS*(CSa+ISa+CS)/N)-gammaS*CS-rhoS*CS+omega*CSa-new_teta*CS
    
    
    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS)
    
    
    exp<-(Sa+CRa+CSa+IRa+ISa)*100/N
    non_exp<-(S+CR+CS)*100/N
    
    list(res,N=N)
    
  })
  
}


create_params<-function(betaR=0.0147,betaS=0.01428,ct=0.95,deltaRa=0,deltaSa=0,gammaR=9.80*10^-3,gammaS=1.02*10^-2,rhoR=8.22*10^-6,rhoS=1.20*10^-4,teta=0.0014,omega=0.08, alpha=0.33, sigmaR=1,sigmaS=0, ATB=0.3)
{
  list(betaR=betaR,betaS=betaS,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gammaR=gammaR,gammaS=gammaS,rhoR=rhoR,rhoS=rhoS,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB)
}

create_initial_cond<-function(Sa0=1000,CRa0=1,CSa0=1,IRa0=0,ISa0=0,S0=1000,CR0=1,CS0=1){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0)
}

run<-function(Init.cond,param,Tmax=400,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model, param,vec_virus=vec_virus))
  return(result)
  
}

graph<- function(data,filter_values,title){
  #data_name<-as.character(substitute(data))
  
  
  if(!is.null(filter_values))
  {
    p<-data %>%
      melt(id = "time") %>%
      filter(variable %in% filter_values) %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 8, face = "bold"),
            legend.text = element_text(size = 6)) +
      labs(title=title,x = "Time", y = "Value", colour = "Population:")
    
    
  }
  else{
    p<-data %>%
      melt(id = "time") %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 8, face = "bold"),
            legend.text = element_text(size = 6)) +
      labs(title=title,x = "Time", y = "Value", colour = "Population:")
    
    
  }
  
  return(p)
}
# pas d'épidémie
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond()
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL,title="S.Aureus colonization without a virus epidemic")
graph(run1,c("CRa","CSa","CR","CS"),title=NULL)
graph(run1,c("IRa","ISa"),title=NULL)


# Epidémie de grippe mais pas de vaccination
param<-create_params(rhoR=0,rhoS=0)
Init.cond<-create_initial_cond()
run0<-run(Init.cond,param)

vec_virus=I_vac_0
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL,title="S.Aureus colonization with influenza epidemic, no vaccination")
graph(run2,c("IRa","ISa"), title=NULL)
graph(run2,c("CRa","CSa","CR","CS"),title=NULL)

# Epidémie de grippe vaccination 50%
vec_virus=I_vac_50
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run3<-run(Init.cond,param)
run3_g<-graph(run3,NULL,title="S.Aureus colonization with influenza epidemic, vaccination 50%")
graph(run3,c("IRa","ISa"),title=NULL)
graph(run3,c("CRa","CSa","CR","CS"),title=NULL)


# Epidémie de grippe vaccination 80%
vec_virus=I_vac_80
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run4<-run(Init.cond,param)
run4_g<-graph(run4,NULL,title="S.Pneumoniae colonization with influenza epidemic, vaccination 80%")
graph(run4,c("IRa","ISa"),title=NULL)
graph(run4,c("CRa","CSa","CR","CS"),title=NULL)

all_res <- data.frame(
  time = run2$time,
  IR_no_vaccination = run2$IRa,
  IS_no_vaccination = run2$ISa,
  IR_50_vaccination = run3$IRa,
  IS_50_vaccination = run3$ISa,
  IR_80_vaccination = run4$IRa,
  IS_80_vaccination = run4$ISa,
  I_no_vaccination= run2$IRa + run2$ISa,
  I_50_vaccination= run3$IRa + run3$ISa,
  I_80_vaccination= run4$IRa + run4$ISa
)
IR_g<-graph(all_res,c("IR_no_vaccination","IR_50_vaccination","IR_80_vaccination"),title=NULL)
IS_IR_g<-graph(all_res,c("IR_no_vaccination","IR_50_vaccination","IR_80_vaccination","IS_no_vaccination","IS_50_vaccination","IS_80_vaccination"),title=NULL)
I_tot_g<-graph(all_res,c("I_no_vaccination","I_50_vaccination","I_80_vaccination"),title=NULL)

grid.arrange(run1_g,run2_g,run3_g,IR_g,ncol=2)

res <- data.frame(time = r1$time)
for (i in seq(1,19,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                 IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                                 CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
  runt<-run(Init.cond,param)
  
  col<-paste("vaccination",results_df[i,1])
  
  res[[col]]<- runt$IRa
  
}

graph(res,NULL,title=NULL)
graph(res,c("vaccination 0.1","vaccination 0.95"),title=NULL)
