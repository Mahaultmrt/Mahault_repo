library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)


SIR_model_vacc_2 <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    N=S+Sv+Snv+Iv+Inv+Rv+Rnv
    
    dS<- -S*(Pv+(1-Pv))
    dSv<-S*Pv-Sv*((beta*PI*Iv/N)+beta*Inv/N)
    dSnv<-S*(1-Pv)-Snv*((beta*PI*Iv/N)+beta*Inv/N)
    dIv<- Sv*((beta*PI*Iv/N)+beta*Inv/N) -gamma*Iv
    dInv<-Snv*((beta*PI*Iv/N)+beta*Inv/N)-gamma*Inv
    dRv<-gamma*Iv
    dRnv<-gamma*Inv
    res <-c(dS,dSv,dSnv,dIv,dInv,dRv,dRnv)
    
    list(res)
    
  })
  
}

create_params<-function(beta=4,gamma=0.2,Pv=0.8,PI=(1-vf),vf=0.95)
{
  list(beta=beta,gamma=gamma,Pv=Pv,PI=PI,vf=vf)
}

create_initial_cond<-function(S0=100,Sv0=0,Snv0=0,Iv0=1,Inv0=1,Rv0=0,Rnv0=0){
  c(S=S0,Sv=Sv0,Snv=Snv0,Iv=Iv0,Inv=Inv0,Rv=Rv0,Rnv=Rnv0)
}

run<-function(Init.cond,param,Tmax=50,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, SIR_model_vacc_2, param))
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

param<-create_params()
Init.cond<-create_initial_cond(Iv=50,Inv=50)
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


 