library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

#code du modèle page 5

Res_model <- function(t, pop, param,vec_virus_v,vec_virus_nv) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Snva+CRnva+CSnva+IRnva+ISnva+Snv+CRnv+CSnv+IRnv+ISnv+
      Sva+CRva+CSva+IRva+ISva+Sv+CRv+CSv+IRv+ISv
    
    infection_v<- vec_virus_v(t)
    infection_nv<- vec_virus_nv(t)
    
    
    dSnva <- -Snva*((beta*ct*(CRnva+CRnv+CRva+CRv+IRnva+IRnv+IRva+IRv)/N)+
                      beta*(CSnva+CSnv+CSva+CSv+ISnva+ISnv+ISva+ISv)/N)+
                      delta*IRnva+delta*ISnva-omega*Snva+teta*Snv+gamma*CRnva+(gamma+alpha)*CSnva+infection_nv*atb #*100
    dCRnva <- Snva*(beta*ct*(CRnva+CRnv+CRva+CRv+IRnva+IRnv+IRva+IRv)/N)-gamma*CRnva-rhoRa*CRnva-omega*CRnva+teta*CRnv
    dCSnva <- Snva*(beta*(CSnva+CSnv+CSva+CSv+ISnva+ISnv+ISva+ISv)/N)-(gamma+alpha)*CSnva-rhoSa*CSnva-omega*CSnva+teta*CSnv
    dIRnva <- rhoRa*CRnva-delta*IRnva-omega*IRnva+teta*IRnv
    dISnva <- rhoSa*CSnva-delta*ISnva-omega*ISnva+teta*ISnv
    
    dSnv <- -Snv*((beta*ct*(CRnva+CRnv+CRva+CRv+IRnva+IRnv+IRva+IRv)/N)+
                      beta*(CSnva+CSnv+CSva+CSv+ISnva+ISnv+ISva+ISv)/N)+
      delta*IRnv+delta*ISnv+omega*Snva-teta*Snv+gamma*CRnv+gamma*CSnv
    dCRnv <- Snv*(beta*ct*(CRnva+CRnv+CRva+CRv+IRnva+IRnv+IRva+IRv)/N)-gamma*CRnv-rho*CRnv+omega*CRnva-teta*CRnv
    dCSnv <- Snv*(beta*(CSnva+CSnv+CSva+CSv+ISnva+ISnv+ISva+ISv)/N)-gamma*CSnv-rho*CSnv+omega*CSnva-teta*CSnv
    dIRnv <- rho*CRnv-delta*IRnv+omega*IRnva-teta*IRnv
    dISnv <- rho*CSnv-delta*ISnv+omega*ISnva-teta*ISnv
    
    dSva <- -Sva*((beta*ct*(CRnva+CRnv+CRva+CRv+IRnva+IRnv+IRva+IRv)/N)+
                      beta*(CSnva+CSnv+CSva+CSv+ISnva+ISnv+ISva+ISv)/N)+
      delta*IRva+delta*ISva-omega*Sva+teta*Sv+gamma*CRva+(gamma+alpha)*CSva+infection_v*atb #*100
    dCRva <- Sva*(beta*ct*(CRnva+CRnv+CRva+CRv+IRnva+IRnv+IRva+IRv)/N)-gamma*CRva-rhoRa*CRva-omega*CRva+teta*CRv
    dCSva <- Snva*(beta*(CSnva+CSnv+CSva+CSv+ISnva+ISnv+ISva+ISv)/N)-(gamma+alpha)*CSnva-rhoSa*CSnva-omega*CSnva+teta*CSnv
    dIRva <- rhoRa*CRva-delta*IRva-omega*IRva+teta*IRv
    dISva <- rhoSa*CSva-delta*ISva-omega*ISva+teta*ISv
    
    dSv <- -Sv*((beta*ct*(CRnva+CRnv+CRva+CRv+IRnva+IRnv+IRva+IRv)/N)+
                    beta*(CSnva+CSnv+CSva+CSv+ISnva+ISnv+ISva+ISv)/N)+
      delta*IRv+delta*ISv+omega*Sva-teta*Sv+gamma*CRv+gamma*CSv
    dCRv <- Sv*(beta*ct*(CRnva+CRnv+CRva+CRv+IRnva+IRnv+IRva+IRv)/N)-gamma*CRv-rho*CRv+omega*CRva-teta*CRv
    dCSv <- Snv*(beta*(CSnva+CSnv+CSva+CSv+ISnva+ISnv+ISva+ISv)/N)-gamma*CSnv-rho*CSnv+omega*CSnva-teta*CSnv
    dIRv <- rho*CRv-delta*IRv+omega*IRva-teta*IRv
    dISv <- rho*CSv-delta*ISv+omega*ISva-teta*ISv
    
 
    
    res<-c(dSnva,dCRnva,dCSnva,dIRnva,dISnva,dSnv,dCRnv,dCSnv,dIRnv,dISnv,
             dSva,dCRva,dCSva,dIRva,dISva,dSv,dCRv,dCSv,dIRv,dISv)
    
    
    
    list(res,c(infection_v,infection_nv))
    
  })
  
}

# On part du principe qu'il n'y a pas de suscpetible non exposé au début (seulement ceux qui ne seront plus exposés aux ATB)

create_params<-function(beta=1,ct=0.8,delta=0.14,gamma=0.03,rho=0.1,rhoRa=0.08,rhoSa=0,teta=0.022,omega=0.07,alpha=0.5,atb=0.7)
{
  list(beta=beta,ct=ct,delta=delta,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,atb=atb)
}

create_initial_cond<-function(Snva0=0,CRnva0=1,CSnva0=1,IRnva0=0,ISnva0=0,Snv0=0,CRnv0=1,
                              CSnv0=1,IRnv0=0,ISnv0=0,Sva0=0,CRva0=1,CSva0=1,IRva0=0,ISva0=0,Sv0=0,
                              CRv0=1,CSv0=1,IRv0=0,ISv0=0){
  c(Snva=Snva0,CRnva=CRnva0,CSnva=CSnva0,IRnva=IRnva0,ISnva=ISnva0,Snv=Snv0,CRnv=CRnv0,
    CSnv=CSnv0,IRnv=IRnv0,ISnv=ISnv0,Sva=Sva0,CRva=CRva0,CSva=CSva0,IRva=IRva0,ISva=ISva0,Sv=Sv0,
    CRv=CRv0,CSv=CSv0,IRv=IRv0,ISv=ISv0)
}

run<-function(Init.cond,param,Tmax=200,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model, param,vec_virus_v=vec_virus_v,vec_virus_nv=vec_virus_nv))
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


param<-create_params()
Init.cond<-create_initial_cond()
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL)

graph(run1,c("IRnva","ISnva","IRnv","ISnv","IRva","ISva","IRv","ISv"))
graph(run1,c("CRnva","CRnv","CRva","CRv"))
graph(run1,c("Snva","Sva","Snv","Sv"))
graph(run1,c("IRnva","IRnv","IRva","IRv"))


