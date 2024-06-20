library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)


SIR_model_vacc_2 <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    N=Sv+Snv+Iv+Inv+Rv+Rnv
    
   
    dSv<- -Sv*beta*(1-vf)*((Iv+Inv)/N)
    dSnv<- -Snv*beta*((Iv+Inv)/N)
    dIv<- Sv*beta*(1-vf)*((Iv+Inv)/N) -gamma*Iv
    dInv<-Snv*beta*((Iv+Inv)/N)-gamma*Inv
    dRv<-gamma*Iv
    dRnv<-gamma*Inv
    res <-c(dSv,dSnv,dIv,dInv,dRv,dRnv)
    list(res,N=N)

    
  })
  
}

create_params<-function(beta=0.28,gamma=0.14,vf=0.6)
{
  list(beta=beta,gamma=gamma,vf=vf)
}

create_initial_cond<-function(Sv0=0,Snv0=100000,Iv0=0,Inv0=1,Rv0=0,Rnv0=0){
  c(Sv=Sv0,Snv=Snv0,Iv=Iv0,Inv=Inv0,Rv=Rv0,Rnv=Rnv0)
}

run<-function(Init.cond,param,Tmax=365,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, SIR_model_vacc_2, param))
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
#pas de vaccination
param<-create_params()
Init.cond<-create_initial_cond()
r1<-run(Init.cond,param)
r1$Iv_Inv<-r1$Iv+r1$Inv
r1_g<- graph(r1,NULL,title="Epidemic Dynamics of Influenza per 100,000 Population without vaccination")
I_g1<-graph(r1,"Iv_Inv", title="Cumulative Incidence of Infected Individuals per 100,000 without vaccination")
Iv_Inv_g1<-graph(r1,c("Iv","Inv"),title="Infected Individuals per 100,000 without vaccination")
grid.arrange(I_g1,Iv_Inv_g1,ncol=1)


I_vac_0<-approxfun(r1$time,r1%>%
                        mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                        select(propI)%>%
                        pull)

# vaccination 50%
param<-create_params()
Init.cond<-create_initial_cond(Sv0=100000*0.5,Snv0=100000*0.5)
r2<-run(Init.cond,param)
r2$Iv_Inv<-r2$Iv+r2$Inv
r2_g<- graph(r2,NULL,title="Epidemic Dynamics of Influenza per 100,000 Population with 50% vaccine coverage for Influenza")
I_g2<-graph(r2,"Iv_Inv", title="Cumulative Incidence of Infected Individuals per 100,000 with 50% vaccine covergage")
Iv_Inv_g2<-graph(r2,c("Iv","Inv"),title="Infected Individuals per 100,000 with 50% vaccine coverage")
prop_I2=r2%>%
  mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
  select(propI)%>%
  pull
r2$propI<-prop_I2
graph(r2,"propI","Cumulative Incidence of Infected Individuals per 100,000 with 50% Vaccination Coverage for Influenza")
grid.arrange(I_g2,Iv_Inv_g2,ncol=1)


I_vac_50<-approxfun(r2$time,r2%>%
                     mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                     select(propI)%>%
                     pull)

# vaccination 80%
param<-create_params()
Init.cond<-create_initial_cond(Sv0=100000*0.8,Snv0=100000*0.2)
r3<-run(Init.cond,param)
r3$Iv_Inv<-r3$Iv+r3$Inv
r3_g<- graph(r3,NULL,title="Epidemic Dynamics of Influenza per 100,000 Population with 80% vaccine coverage for Influenza")
I_g3<-graph(r3,"Iv_Inv", title="Cumulative Incidence of Infected Individuals per 100,000 with 80% vaccine covergage")
Iv_Inv_g3<-graph(r3,c("Iv","Inv"),title="Infected Individuals per 100,000 with 80% vaccine coverage")
grid.arrange(I_g3,Iv_Inv_g3,ncol=1)

I_vac_80<-approxfun(r3$time,r3%>%
                      mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                      select(propI)%>%
                      pull)

vec_virus_0<-function(time){
  return(0)
}

grid.arrange(Iv_Inv_g1,Iv_Inv_g2,Iv_Inv_g3,ncol=2)

results_df <- data.frame(vacc = numeric(), max_propI = numeric(), last_propR=numeric())
I_vac<-list()
for (i in seq(0.1,1,by=0.05)){
  param<-create_params()
  Init.cond<-create_initial_cond(Sv0=1000*i,Snv0=1000*(1-i))
  r<-run(Init.cond,param)
  r$Iv_Inv<-r$Iv+r$Inv
  r$Rv_Rnv<-r$Rv+r$Rnv
  prop_I=r%>%
    mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
    select(propI)%>%
    pull
  max_propI=max(prop_I,na.rm=TRUE)
  prop_R<-r%>%
    mutate(propR=Rv_Rnv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
    select(propR)%>%
    pull
  last_propR=tail(prop_R, n = 1)
  new_row=data.frame(vacc=i, max_propI,last_propR, n = 1)
  results_df <- bind_rows(results_df, new_row)
  virus<-approxfun(r$time,prop_I)
  I_vac<-append(I_vac,virus)
 

}

ggplot() +   
  geom_point(data=results_df, aes(x=vacc,y=max_propI,colour="Infected people at Epidemic peak"))+
  geom_line(data=results_df, aes(x=vacc,y=max_propI, colour="Infected people at Epidemic peak"))+
  geom_point(data=results_df, aes(x=vacc,y=last_propR, colour="Annual recovery"))+
  geom_line(data=results_df, aes(x=vacc,y=last_propR, colour="Annual recovery"))+
  labs(title = "Infected people at Epidemic peak and annual recovery according to vaccination", y = "Population (per 100,000)",
       x = "Vaccine coverage") +
  scale_colour_manual(name = "Legend", values = c("Infected people at Epidemic peak" = "purple", "Annual recovery" = "orange")) +
  theme_bw()




