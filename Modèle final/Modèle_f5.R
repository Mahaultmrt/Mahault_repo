library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

#Code modèle page8 n°1

Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+IRa+ISa+S+CR+CS+IR+IS
    
    
    new_teta<-teta+vec_virus(t)
    new_tetaI<-tetaI+vec_virus(t)
    
    
    dSa <- -Sa*((beta*ct*(CRa+IRa+CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+delta*IRa+deltaSa*ISa-omega*Sa+new_teta*S+gamma*CRa+(gamma+alpha)*CSa
    dCRa <- Sa*(beta*ct*(CRa+IRa+CR+IR)/N)-gamma*CRa-rhoRa*CRa-omega*CRa+new_teta*CR
    dCSa <- Sa*(beta*(CSa+ISa+CS+IS)/N)-(gamma+alpha)*CSa-rhoSa*CSa-omega*CSa+new_teta*CS
    dIRa <- rhoRa*CRa-delta*IRa+new_tetaI*IR
    dISa <- rhoSa*CSa-deltaSa*ISa+new_tetaI*IS
    
    dS <- -S*((beta*ct*(CRa+IRa+CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+delta*IR+delta*IS+omega*Sa-new_teta*S+gamma*CR+gamma*CS
    dCR <- S*(beta*ct*(CRa+IRa+CR+IR)/N)-gamma*CR-rho*CR+omega*CRa-new_teta*CR
    dCS <- S*(beta*(CSa+ISa+CS+IS)/N)-gamma*CS-rho*CS+omega*CSa-new_teta*CS
    dIR <- rho*CR-delta*IR-new_tetaI*IR
    dIS <- rho*CS-delta*IS-new_tetaI*IS
    
    
    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS,dIR,dIS)
    
    
    exp<-(Sa+CRa+CSa+IRa+ISa)*100/N
    non_exp<-(S+CR+CS+IR+IS)*100/N
    
    list(res,new_teta=new_teta)
    
  })
  
}


create_params<-function(beta=0.065,ct=0.95,delta=0.14,deltaSa=0.33,gamma=0.03,rho=0.03,rhoRa=0.03,rhoSa=0.03,teta=0.0014,omega=0.14,tetaI=1, alpha=0.5)
{
  list(beta=beta,ct=ct,delta=delta,deltaSa=deltaSa,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,tetaI=tetaI,alpha=alpha)
}

create_initial_cond<-function(Sa0=100,CRa0=1,CSa0=1,IRa0=0,ISa0=0,S0=100,CR0=1,CS0=1,IR0=0,IS0=0){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0,IR=IR0,IS=IS0)
}

run<-function(Init.cond,param,Tmax=400,dt=1){
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
            legend.text = element_text(size = 8)) +
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
            legend.text = element_text(size = 8)) +
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
graph(run1,c("IRa","ISa","IR","IS"))

# vec_virus=vec_virus_v
# param<-create_params()
# Init.cond<-create_initial_cond(Sa0=tail(run1$Sa, n = 1),CRa0=tail(run1$CRa, n = 1),CSa0=tail(run1$CSa, n = 1),
#                                IRa0=tail(run1$IRa, n = 1),ISa0=tail(run1$ISa, n = 1),S0=tail(run1$S, n = 1),
#                                CR0=tail(run1$CR, n = 1),CS0=tail(run1$CS, n = 1),IR0=tail(run1$IR, n = 1),IS0=tail(run1$IS, n = 1))
# run2<-run(Init.cond,param)
# run2_g<-graph(run2,NULL)
# graph(run2,c("IRa","ISa","IR","IS"))
# graph(run2,c("CRa","CSa","CR","CS"))
# 
# 
# vec_virus=vec_virus_nv
# param<-create_params()
# Init.cond<-create_initial_cond(Sa0=tail(run1$Sa, n = 1),CRa0=tail(run1$CRa, n = 1),CSa0=tail(run1$CSa, n = 1),
#                                IRa0=tail(run1$IRa, n = 1),ISa0=tail(run1$ISa, n = 1),S0=tail(run1$S, n = 1),
#                                CR0=tail(run1$CR, n = 1),CS0=tail(run1$CS, n = 1),IR0=tail(run1$IR, n = 1),IS0=tail(run1$IS, n = 1))
# run3<-run(Init.cond,param)
# run3_g<-graph(run3,NULL)
# graph(run3,c("IRa","ISa","IR","IS"))
# graph(run3,c("CRa","CSa","CR","CS"))
# 
# merge_run<-function(data1,data2){
#   new_run <- merge(data1, data2,by = "time", suffixes = c(".vaccine", ".non_vaccine"))
#   return(new_run)
#   
# }
# graph(merge_run(run2,run3),c("new_teta.vaccine","new_teta.non_vaccine"))
# 
# grid.arrange(run1_g,run2_g,ncol=2)
# grid.arrange(graph(run1,c("IRa","ISa","IR","IS")),graph(run2,c("IRa","ISa","IR","IS")),ncol=2,graph(run3,c("IRa","ISa","IR","IS")))
# 
