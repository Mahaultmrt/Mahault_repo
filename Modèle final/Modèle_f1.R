library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

#code du modèle page 2

Res_model <- function(t, pop, param,vec_virus) {

  with(as.list(c(pop, param)), {


    N=Sa+CRa+CSa+IRa+ISa+S+CR+CS+IR+IS

    infection<- vec_virus(t)
    
    dSa <- -Sa*((beta*ct*(CRa+IRa+CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+delta*IRa+delta*ISa-omega*Sa+teta*S+gamma*CRa+(gamma+alpha)*CSa+infection*atb
    dCRa <- Sa*(beta*ct*(CRa+IRa+CR+IR)/N)-gamma*CRa-rhoRa*CRa-omega*CRa+teta*CR
    dCSa <- Sa*(beta*(CSa+ISa+CS+IS)/N)-(gamma+alpha)*CSa-rhoSa*CSa-omega*CSa+teta*CS
    dIRa <- rhoRa*CRa-delta*IRa-omega*IRa+teta*IR
    dISa <- rhoSa*CSa-delta*ISa-omega*ISa+teta*IS

    dS <- -S*((beta*ct*(CRa+IRa+CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+delta*IR+delta*IS+omega*Sa-teta*S+gamma*CR+gamma*CS
    dCR <- S*(beta*ct*(CRa+IRa+CR+IR)/N)-gamma*CR-rho*CR+omega*CRa-teta*CR
    dCS <- S*(beta*(CSa+ISa+CS+IS)/N)-gamma*CS-rho*CS+omega*CSa-teta*CS
    dIR <- rho*CR-delta*IR+omega*IRa-teta*IR
    dIS <- rho*CS-delta*IS+omega*ISa-teta*IS


    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS,dIR,dIS)



    list(res)

  })

}

# On part du principe qu'il n'y a pas de suscpetible non exposé au début (seulement ceux qui ne seront plus exposés aux ATB)

create_params<-function(beta=1,ct=0.8,delta=0.14,gamma=0.03,rho=0.1,rhoRa=0.08,rhoSa=0,teta=0.022,omega=0.07,alpha=0.5,atb=0.7)
{
  list(beta=beta,ct=ct,delta=delta,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,atb=atb)
}

create_initial_cond<-function(Sa0=0,CRa0=1,CSa0=1,IRa0=0,ISa0=0,S0=0,CR0=1,CS0=1,IR0=0,IS0=0){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0,IR=IR0,IS=IS0)
}

run<-function(Init.cond,param,Tmax=200,dt=1){
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

vec_virus<- approxfun(r1$time, r1$Iv)
param<-create_params()
Init.cond<-create_initial_cond()
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL)


vec_virus<- approxfun(r1$time, r1$Inv)
param<-create_params()
Init.cond<-create_initial_cond()
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL)

merge_run<-function(data1,data2){
  new_run <- merge(data1, data2, by = "time", suffixes = c(".vaccine", ".non_vaccine"))
  return(new_run)
  
}
graph(merge_run(run1,run2),c("CR.vaccine","CR.non_vaccine","CRa.vaccine","CRa.non_vaccine"))

grid.arrange(run1_g,run2_g,ncol=1)
