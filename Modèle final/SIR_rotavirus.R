library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)


SIR_model_vacc_2 <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    N=Sv+Snv+Iv+Inv+Rv+Rnv
    
    
    dSv<- -Sv*beta*(1-vf)*((Iv+Inv)/N)+gamma*Ps*Iv
    dSnv<- -Snv*beta*((Iv+Inv)/N)+gamma*Ps*Inv
    dIv<- Sv*beta*(1-vf)*((Iv+Inv)/N) -gamma*(1-Ps)*Iv-gamma*Ps*Iv
    dInv<-Snv*beta*((Iv+Inv)/N)-gamma*(1-Ps)*Inv-gamma*Ps*Inv
    dRv<-gamma*Iv*(1-Ps)
    dRnv<-gamma*Inv*(1-Ps)
    res <-c(dSv,dSnv,dIv,dInv,dRv,dRnv)
    list(res)
    
    
  })
  
}


create_params<-function(beta=2.5,gamma=0.2,vf=0.643,Ps=0.75)
{
  list(beta=beta,gamma=gamma,vf=vf,Ps=Ps)
}

create_initial_cond<-function(Sv0=100,Snv0=100,Iv0=0,Inv0=1,Rv0=0,Rnv0=0){
  c(Sv=Sv0,Snv=Snv0,Iv=Iv0,Inv=Inv0,Rv=Rv0,Rnv=Rnv0)
}

run<-function(Init.cond,param,Tmax=300,dt=1){
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
r1$Iv_Inv<-r1$Iv+r1$Inv
r1_g<- graph(r1,NULL)
graph(r1,"Iv_Inv")
graph(r1,c("Iv","Inv"))
grid.arrange(graph(r1,c("Iv","Inv")),graph(r1,"Iv_Inv"),ncol=1)


# vec_virus_v <- approxfun(r1$time, r1$Iv)
# vec_virus_nv<- approxfun(r1$time,r1$Inv)

vec_virus_v<-approxfun(r1$time,r1%>%
                         mutate(propIv=Iv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                         select(propIv)%>%
                         pull)

vec_virus_nv<-approxfun(r1$time,r1%>%
                          mutate(propInv=Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                          select(propInv)%>%
                          pull)

vec_virus_0<-function(time){
  return(0)
}
