library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

#Code modèle page8 n°2

Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+IRa+ISa+S+CR+CS
    
    
    new_teta<-teta+vec_virus(t)

    
    dSa <- -Sa*((beta*ct*(CRa+IRa+CR)/N)+beta*(CSa+ISa+CS)/N)-omega*Sa+new_teta*S+(gamma+alpha*(1-sigmaR))*CRa+(gamma+alpha*(1-sigmaS))*CSa
    dCRa <- Sa*(beta*ct*(CRa+IRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-rhoRa*CRa-omega*CRa+new_teta*CR
    dCSa <- Sa*(beta*(CSa+ISa+CS)/N)-(gamma+alpha*(1-sigmaS))*CSa-rhoSa*CSa-omega*CSa+new_teta*CS
    dIRa <- rhoRa*CRa-deltaRa*IRa+rho*CR
    dISa <- rhoSa*CSa-deltaSa*ISa+rho*CS
    
    dS <- -S*((beta*ct*(CRa+IRa+CR)/N)+beta*(CSa+ISa+CS)/N)+omega*Sa-new_teta*S+gamma*CR+gamma*CS+deltaRa*IRa+deltaSa*ISa
    dCR <- S*(beta*ct*(CRa+IRa+CR)/N)-gamma*CR-rho*CR+omega*CRa-teta*CR
    dCS <- S*(beta*(CSa+ISa+CS)/N)-gamma*CS-rho*CS+omega*CSa-teta*CS
    
    
    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS)
    
    
    exp<-(Sa+CRa+CSa+IRa+ISa)*100/N
    non_exp<-(S+CR+CS)*100/N
    
    list(res)
    
  })
  
}


create_params<-function(beta=0.056,ct=0.9652,deltaRa=0,deltaSa=0,gamma=0.05,rho=3*10^-6,rhoRa=3*10^-6,rhoSa=3*10^-6,teta=0.0014,omega=0.08, alpha=0.33, sigmaR=1,sigmaS=0)
{
  list(beta=beta,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS)
}

create_initial_cond<-function(Sa0=100,CRa0=1,CSa0=1,IRa0=0,ISa0=0,S0=100,CR0=1,CS0=1){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0)
}

run<-function(Init.cond,param,Tmax=365,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model, param,vec_virus=vec_virus))
  return(result)
  
}

graph<- function(data,filter_values){
  data_name<-as.character(substitute(data))
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
      labs(title=data_name,x = "Time", y = "Value", colour = "Population:")
    
    
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
      labs(title=data_name,x = "Time", y = "Value", colour = "Population:")
    
    
  }
  
  return(p)
}

vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond()
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL)
graph(run1,c("CRa","CSa","CR","CS"))
graph(run1,c("IRa","ISa"))

param<-create_params(rho=0,rhoRa=0,rhoSa=0)
Init.cond<-create_initial_cond()
run0<-run(Init.cond,param)

vec_virus=vec_virus_v
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL)
graph(run2,c("IRa","ISa"))
graph(run2,c("CRa","CSa","CR","CS"))


vec_virus=vec_virus_nv
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run3<-run(Init.cond,param)
run3_g<-graph(run3,NULL)
graph(run3,c("IRa","ISa"))
graph(run3,c("CRa","CSa","CR","CS"))

grid.arrange(graph(run2,c("IRa","ISa")),graph(run3,c("IRa","ISa")),ncol=1)
grid.arrange(graph(run2,c("CRa","CR")),graph(run3,c("CRa","CR")),ncol=1)
grid.arrange(run1_g,run2_g,run3_g,ncol=2)



merge_run<-function(data1,data2){
  new_run <- merge(data1, data2,by = "time", suffixes = c(".vaccine", ".non_vaccine"))
  return(new_run)

}
graph(merge_run(run2,run3),c("new_teta.vaccine","new_teta.non_vaccine"))
graph_Iv_Inv<-graph(merge_run(run2,run3),c("IRa.vaccine","IRa.non_vaccine","ISa.vaccine","ISa.non_vaccine"))

grid.arrange(run1_g,run2_g,run3_g,graph_Iv_Inv,ncol=2)


runI <- run2
for (col_name in colnames(run2)) {
  if (col_name %in% colnames(run3)) {
    runI[[col_name]] <- run2[[col_name]] + run3[[col_name]]
  }
}
runI_g<-graph(runI,c("ISa","IRa"))

runtest <- run3
for (col_name in colnames(run3)) {
  if (col_name %in% colnames(run2)) {
    runtest[[col_name]] <- run3[[col_name]] - run2[[col_name]]
  }
}
