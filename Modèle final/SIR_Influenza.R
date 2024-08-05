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
    
    Incidence<-Sv*beta*(1-vf)*((Iv+Inv)/N)+Snv*beta*((Iv+Inv)/N)
    
    res <-c(dSv,dSnv,dIv,dInv,dRv,dRnv)
    list(res,Incidence=Incidence)

    
  })
  
}

create_params<-function(beta=0.28,gamma=0.14,vf=0.6)
{
  list(beta=beta,gamma=gamma,vf=vf)
}

create_initial_cond<-function(Sv0=0,Snv0=99500,Iv0=0,Inv0=500,Rv0=0,Rnv0=0){
  c(Sv=Sv0,Snv=Snv0,Iv=Iv0,Inv=Inv0,Rv=Rv0,Rnv=Rnv0)
}

run<-function(Init.cond,param,Tmax=365,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, SIR_model_vacc_2, param))
  tot <- sum(result[1, !(colnames(result) %in% c("time", "Incidence"))])
  result$Incidence<-result$Incidence/tot
  return(result)
  
}


#pas de vaccination
param<-create_params()
Init.cond<-create_initial_cond()
r1<-run(Init.cond,param)
r1$Iv_Inv<-r1$Iv+r1$Inv
r1_g<- graph2(r1,NULL,title="Epidemic Dynamics of Influenza per 100,000 Population without vaccination")
I_g1<-graph2(r1,"Iv_Inv", title="Infected Individuals without vaccination")
Iv_Inv_g1<-graph2(r1,c("Iv","Inv"),title="Infected Individuals per 100,000 without vaccination")
prop_I1=r1%>%
  mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
  select(propI)%>%
  pull
r1$propI<-prop_I1
propI1_g<-graph2(r1,"propI","Porpotion of Infected Individuals \nwithout Vaccination for Influenza")
Incidence1<-graph2(r1,"Incidence","Incidence of infected people without vaccination")
grid.arrange(I_g1,Iv_Inv_g1,ncol=1)


I_vac_0<-approxfun(r1$time,r1%>%
                        mutate(propInc=Incidence)%>%
                        select(propInc)%>%
                        pull)

# vaccination 50%
param<-create_params()
Init.cond<-create_initial_cond(Sv0=99500*0.5,Snv0=99500*0.5)
r2<-run(Init.cond,param)
r2$Iv_Inv<-r2$Iv+r2$Inv
r2_g<- graph2(r2,NULL,title="Epidemic Dynamics of Influenza per 100,000 Population with 50% vaccine coverage")
I_g2<-graph2(r2,"Iv_Inv", title="Infected Individuals with 50% vaccine covergage")
Iv_Inv_g2<-graph2(r2,c("Iv","Inv"),title="Infected Individuals with 50% vaccine coverage")
prop_I2=r2%>%
  mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
  select(propI)%>%
  pull
r2$propI<-prop_I2
propI2_g<-graph2(r2,"propI","Proportion of People infected by Influenza with 50% vaccine Coverage")
Incidence2<-graph2(r2,"Incidence","Incidence of infected people \nwith 50% vaccine coverage")
grid.arrange(I_g2,Iv_Inv_g2,ncol=1)


I_vac_50<-approxfun(r2$time,r2%>%
                     mutate(propInc=Incidence)%>%
                     select(propInc)%>%
                     pull)



# vaccination 80%
param<-create_params()
Init.cond<-create_initial_cond(Sv0=99500*0.8,Snv0=99500*0.2)
r3<-run(Init.cond,param)
r3$Iv_Inv<-r3$Iv+r3$Inv
r3_g<- graph2(r3,NULL,title="Epidemic Dynamics of Influenza per 100,000 Population with 80% vaccine coverage for Influenza")
I_g3<-graph2(r3,"Iv_Inv", title="Infected Individuals with 80% vaccine covergage")
Iv_Inv_g3<-graph2(r3,c("Iv","Inv"),title="Infected Individuals with 80% vaccine coverage")
prop_I3=r3%>%
  mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
  select(propI)%>%
  pull
r3$propI<-prop_I3
propI3_g<-graph2(r3,"propI","Porpotion of Infected Individuals \nwith 80% vaccine Coverage for Influenza")
Incidence3<-graph2(r3,"Incidence","Incidence of infected people \nwith 80% vaccine coverage")
grid.arrange(I_g3,Iv_Inv_g3,ncol=1)

I_vac_80<-approxfun(r3$time,r3%>%
                      mutate(propInc=Incidence)%>%
                      select(propInc)%>%
                      pull)

grid.arrange(I_g3,Iv_Inv_g3,ncol=1)


# vaccination 70% efficacité vaccinale 30%
param<-create_params(vf=0.3)
Init.cond<-create_initial_cond(Sv0=99500*0.7,Snv0=99500*0.3)
r70<-run(Init.cond,param)
I_vac_70<-approxfun(r70$time,r70%>%
                      mutate(propInc=Incidence)%>%
                      select(propInc)%>%
                      pull)

# vaccination 30% efficacité vaccinale 70%
param<-create_params(vf=0.7)
Init.cond<-create_initial_cond(Sv0=99500*0.3,Snv0=99500*0.7)
r30<-run(Init.cond,param)
I_vac_30<-approxfun(r30$time,r30%>%
                      mutate(propInc=Incidence)%>%
                      select(propInc)%>%
                      pull)


# vaccination 46% efficacité vaccinale 46%
param<-create_params(vf=0.46)
Init.cond<-create_initial_cond(Sv0=99500*0.46,Snv0=99500*0.54)
r46<-run(Init.cond,param)
I_vac_46<-approxfun(r46$time,r46%>%
                      mutate(propInc=Incidence)%>%
                      select(propInc)%>%
                      pull)


vec_virus_0<-function(time){
  return(0)
}


results_df <- data.frame(vacc = numeric(), max_propI = numeric(), last_propR=numeric())
I_vac<-list()
for (i in seq(0,1,by=0.05)){
  param<-create_params()
  Init.cond<-create_initial_cond(Sv0=99500*i,Snv0=99500*(1-i))
  r<-run(Init.cond,param)
  r$Iv_Inv<-r$Iv+r$Inv
  r$Rv_Rnv<-r$Rv+r$Rnv
  prop_I=r%>%
    mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
    select(propI)%>%
    pull
  prop_Inc=r%>%
    mutate(prop_Inc=Incidence)%>%
    select(prop_Inc)%>%
    pull
  max_propI=max(prop_I,na.rm=TRUE)
  prop_R<-r%>%
    mutate(propR=Rv_Rnv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
    select(propR)%>%
    pull
  last_propR=tail(prop_R, n = 1)
  new_row=data.frame(vacc=i, max_propI,last_propR, n = 1)
  results_df <- bind_rows(results_df, new_row)
  virus<-approxfun(r$time,prop_Inc)
  I_vac<-append(I_vac,virus)
 

}


I_R<-graph_I_R(results_df)
  

combined_Incidence<-data.frame(time=seq(from=0,to=365,by=1),no_vaccination=r1$Incidence,vaccination_50=r2$Incidence,vaccination_80=r3$Incidence)
I_flu<-graph2(combined_Incidence,NULL,NULL)


