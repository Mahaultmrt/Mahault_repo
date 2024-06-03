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

create_initial_cond<-function(Sv0=0,Snv0=1000,Iv0=0,Inv0=1,Rv0=0,Rnv0=0){
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



I_vac_0<-approxfun(r1$time,r1%>%
                     mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                     select(propI)%>%
                     pull)
param<-create_params()
Init.cond<-create_initial_cond(Sv0=500,Snv0=500)
r2<-run(Init.cond,param)
r2$Iv_Inv<-r2$Iv+r2$Inv
r2_g<- graph(r2,NULL)
graph(r2,"Iv_Inv")
graph(r2,c("Iv","Inv"))
grid.arrange(graph(r2,c("Iv","Inv")),graph(r2,"Iv_Inv"),ncol=1)


I_vac_50<-approxfun(r2$time,r2%>%
                      mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                      select(propI)%>%
                      pull)

param<-create_params()
Init.cond<-create_initial_cond(Sv0=800,Snv0=200)
r3<-run(Init.cond,param)
r3$Iv_Inv<-r3$Iv+r3$Inv
r3_g<- graph(r3,NULL)
graph(r3,"Iv_Inv")
graph(r3,c("Iv","Inv"))
grid.arrange(graph(r3,c("Iv","Inv")),graph(r3,"Iv_Inv"),ncol=1)

I_vac_80<-approxfun(r3$time,r3%>%
                      mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                      select(propI)%>%
                      pull)

vec_virus_0<-function(time){
  return(0)
}

grid.arrange(r1_g,r2_g,r3_g,ncol=2)


results_df <- data.frame(vacc = numeric(), Max_I = numeric())
I_vac<-list()
for (i in seq(0.1,1,by=0.05)){
  param<-create_params()
  Init.cond<-create_initial_cond(Sv0=1000*i,Snv0=1000*(1-i))
  r<-run(Init.cond,param)
  r$Iv_Inv<-r$Iv+r$Inv
  Max_I=max(r$Iv_Inv,na.rm=TRUE)
  new_row=data.frame(vacc=i, Max_I, n = 1)
  results_df <- bind_rows(results_df, new_row)
  virus<-approxfun(r$time,r%>%
                     mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                     select(propI)%>%
                     pull)
  I_vac<-append(I_vac,virus)
  
  
}

ggplot(data   = results_df,              
       mapping = aes(x = vacc,    
                     y = Max_I)) +   
  geom_point(color="blue")+
  geom_line(color="purple")+
  theme_bw()+
  labs(title="Epidemic spike according to vaccination  ")