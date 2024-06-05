library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

#Code modèle page 12

Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=CRa+CSa+IRa+ISa+CR+CS
    
    
    new_teta<-teta-log(1-vec_virus(t)*ATB)
    
    
    dCRa <- -CRa*(betaS*(CSa+ISa+CS)/N)+CSa*(betaR*ct*(CR+IRa+CRa)/N)-CRa*(gamma+alpha*(1-sigmaR))+CSa*(gamma+alpha*(1-sigmaS))-omega*CRa+new_teta*CR-rho*CRa
    dCSa <- -CSa*(betaR*ct*(CR+IRa+CRa)/N)+CRa*(betaS*(CSa+ISa+CS)/N)-CSa*(gamma+alpha*(1-sigmaS))+CRa*(gamma+alpha*(1-sigmaR))-omega*CSa+new_teta*CS-rho*CSa
    dIRa <- rho*CRa+rho*CR
    dISa <- rho*CSa+rho*CS
    
    dCR <- -CR*(betaS*(CSa+ISa+CS)/N)+CS*(betaR*ct*(CR+IRa+CRa)/N)-new_teta*CR+omega*CRa-rho*CR-gamma*CR+gamma*CS
    dCS <- -CS*(betaR*ct*(CR+IRa+CRa)/N)+CR*(betaS*(CSa+ISa+CS)/N)-new_teta*CS+omega*CSa-rho*CS-gamma*CS+gamma*CR
    
    
    res<-c(dCRa,dCSa,dIRa,dISa,dCR,dCS)
    
    
    exp<-(CRa+CSa+IRa+ISa)*100/N
    non_exp<-(CR+CS)*100/N
    
    list(res,new_teta=new_teta,N=N)
    
  })
  
}


create_params<-function(betaS=0.012,betaR=0.015,ct=1,deltaRa=0,deltaSa=0,gamma=1/120,rho=1/5,teta=0.0014,omega=0.08, alpha=0.33, sigmaR=1,sigmaS=0, ATB=0.3)
{
  list(betaS=betaS,betaR=betaR,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rho=rho,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB)
}

create_initial_cond<-function(CRa0=1000,CSa0=100,IRa0=0,ISa0=0,CR0=100,CS0=1000){
  c(CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,CR=CR0,CS=CS0)
}

run<-function(Init.cond,param,Tmax=300,dt=1){
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

#Pas d'épidémie 
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond()
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL,title="E.Coli colonization without a virus epidemic")
graph(run1,c("CRa","CSa","CR","CS"),title=NULL)
graph(run1,c("IRa","ISa"),title=NULL)

# Epidémie de rotavirus mais pas de vaccination
param<-create_params(rho=0)
Init.cond<-create_initial_cond()
run0<-run(Init.cond,param)

vec_virus=I_vac_0
param<-create_params()
Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL,title="E.Coli colonization with rotavirus epidemic, no vaccination")
graph(run2,c("IRa","ISa"), title=NULL)
graph(run2,c("CRa","CSa","CR","CS"),title=NULL)

# Epidémie de rotavirus vaccination 50%
vec_virus=I_vac_50
param<-create_params()
Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run3<-run(Init.cond,param)
run3_g<-graph(run3,NULL,title="E.Coli colonization with rotavirus epidemic, vaccination 50%")
graph(run3,c("IRa","ISa"),title=NULL)
graph(run3,c("CRa","CSa","CR","CS"),title=NULL)


# # Epidémie de rotavirus vaccination 80%
vec_virus=I_vac_80
param<-create_params()
Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run4<-run(Init.cond,param)
run4_g<-graph(run4,NULL,title="E.Coli colonization with rotavirus epidemic, vaccination 80%")
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
  Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                 IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                 CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
  runt<-run(Init.cond,param)
  
  col<-paste("vaccination",results_df[i,1])
  
  res[[col]]<- runt$IRa
  
}

graph(res,NULL,title=NULL)
graph(res,c("vaccination 0.1","vaccination 0.95"),title=NULL)


