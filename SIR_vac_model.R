library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)


SIR_model_vacc <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    
    N=S+Iv+Inv+Rv+Rnv
    
    
    dS <- -S*((beta*Pv*Iv/N)+beta*(1-Pv)*Inv/N)
    dIv <- (S*beta*Pv*Iv/N)-gamma*Iv
    dInv <- (S*beta*(1-Pv)*Inv/N)-gamma*Inv
    dRv <- gamma*Iv
    dRnv <- gamma*Inv
    
    res <-c(dS,dIv,dInv,dRv,dRnv)
    
    list(res)
    
  })
  
}

create_params<-function(beta=4,gamma=0.2,Pv=(1-vc*vf),vc=0.9,vf=1)
{
  list(beta=beta,gamma=gamma,Pv=Pv,vc=vc,vf=vf)
}

create_initial_cond<-function(S0=100,Iv0=1,Inv0=1,Rv0=0,Rnv0=0){
  c(S=S0,Iv=Iv0,Inv=Inv0,Rv=Rv0,Rnv=Rnv0)
}

run<-function(Init.cond,param,Tmax=50,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, SIR_model_vacc, param))
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
r1<-run(Init.cond,param)
r1_g<- graph(r1,NULL)

param<-create_params(vf=0.7)
Init.cond<-create_initial_cond()
r2<-run(Init.cond,param)
r2_g<- graph(r2,NULL)

param<-create_params(vc=0.7,vf=0.5)
Init.cond<-create_initial_cond()
r3<-run(Init.cond,param)
r3_g<- graph(r3,NULL)

param<-create_params(vc=1,vf=0.1)
Init.cond<-create_initial_cond()
r4<-run(Init.cond,param)
r4_g<- graph(r4,NULL)

grid.arrange(r1_g,r2_g,r3_g,r4_g,ncol=2)


