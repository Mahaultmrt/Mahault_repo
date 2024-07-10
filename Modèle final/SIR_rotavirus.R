library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)


SIR_model_vacc_2 <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    N=Sv+Snv+Iv+Inv+Rv+Rnv
    
    
    dSv<- -Sv*beta*(1-vf)*((Iv+Inv)/N)+gamma*Ps*Iv*(1-vf)
    dSnv<- -Snv*beta*((Iv+Inv)/N)+gamma*Ps*Inv
    dIv<- Sv*beta*(1-vf)*((Iv+Inv)/N) -gamma*(1-Ps*(1-vf))*Iv-gamma*Ps*Iv*(1-vf)
    dInv<-Snv*beta*((Iv+Inv)/N)-gamma*(1-Ps)*Inv-gamma*Ps*Inv
    dRv<-gamma*Iv*(1-(Ps*(1-vf)))
    dRnv<-gamma*Inv*(1-Ps)
    
    Incidence<-Sv*beta*(1-vf)*((Iv+Inv)/N)+Snv*beta*((Iv+Inv)/N)
    
    res <-c(dSv,dSnv,dIv,dInv,dRv,dRnv)
    list(res,Incidence=Incidence,N=N)
    
    
  })
  
}


create_params<-function(beta=2.5,gamma=0.2,vf=0.643,Ps=0.75)
{
  list(beta=beta,gamma=gamma,vf=vf,Ps=Ps)
}

create_initial_cond<-function(Sv0=0,Snv0=99500,Iv0=0,Inv0=500,Rv0=0,Rnv0=0){
  c(Sv=Sv0,Snv=Snv0,Iv=Iv0,Inv=Inv0,Rv=Rv0,Rnv=Rnv0)
}

run<-function(Init.cond,param,Tmax=365,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, SIR_model_vacc_2, param))
  return(result)
  
}



#no vaccincation
param<-create_params()
Init.cond<-create_initial_cond()
r1<-run(Init.cond,param)
r1$Iv_Inv<-r1$Iv+r1$Inv
r1_g<- graph2(r1,NULL,title="Epidemic Dynamics of Rotavirus per 100,000 Population without vaccination")
I_g1<-graph2(r1,"Iv_Inv", title="Infected Individuals without vaccination")
Iv_Inv_g1<-graph2(r1,c("Iv","Inv"),title="Infected Individuals without vaccination")
prop_I1=r1%>%
  mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
  select(propI)%>%
  pull
r1$propI<-prop_I1
propI1_g<-graph2(r1,"propI","Proportion od people infected by Rotavirus \nwithout Vaccination for Rotavirus")
grid.arrange(I_g1,Iv_Inv_g1,ncol=1)

I_vac_0<-approxfun(r1$time,r1%>%
                     mutate(propInc=Incidence/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                     select(propI)%>%
                     pull)

# vaccination 50%
param<-create_params()
Init.cond<-create_initial_cond(Sv0=99500*0.5,Snv0=99500*0.5)
r2<-run(Init.cond,param)
r2$Iv_Inv<-r2$Iv+r2$Inv
r2_g<- graph2(r2,NULL,title="Epidemic Dynamics of Rotavirus per 100,000 Population with 50% vaccine coverage")
I_g2<-graph2(r2,"Iv_Inv", title="Infected Individuals with 50% vaccine covergage")
Iv_Inv_g2<-graph2(r2,c("Iv","Inv"),title="Infected Individuals \nwith 50% vaccine coverage")
prop_I2=r2%>%
  mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
  select(propI)%>%
  pull
r2$propI<-prop_I2
propI2_g<-graph2(r2,"propI","Proportion od people infected by Rotavirus \nwith 50% Vaccination Coverage for Rotavirus")
grid.arrange(I_g2,Iv_Inv_g2,ncol=1)


I_vac_50<-approxfun(r2$time,r2%>%
                      mutate(propInc=Incidence/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                      select(propI)%>%
                      pull)



# vaccination 80%
param<-create_params()
Init.cond<-create_initial_cond(Sv0=99500*0.8,Snv0=99500*0.2)
r3<-run(Init.cond,param)
r3$Iv_Inv<-r3$Iv+r3$Inv
r3_g<- graph2(r3,NULL,title="Epidemic Dynamics of Rotavirus per 100,000 Population with 80% vaccine coverage")
I_g3<-graph2(r3,"Iv_Inv", title="Infected Individuals with 80% vaccine covergage")
Iv_Inv_g3<-graph2(r3,c("Iv","Inv"),title="Infected Individuals with 80% vaccine coverage")
prop_I3=r3%>%
  mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
  select(propI)%>%
  pull
r3$propI<-prop_I3
propI3_g<-graph2(r3,"propI","Proportion od people infected by Rotavirus \nwith 80% Vaccination Coverage for Rotavirus")
grid.arrange(I_g3,Iv_Inv_g3,ncol=1)

I_vac_80<-approxfun(r3$time,r3%>%
                      mutate(propInc=Incidence/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                      select(propI)%>%
                      pull)


# vaccination 100%
param<-create_params()
Init.cond<-create_initial_cond(Sv0=99500,Snv0=0)
r4<-run(Init.cond,param)
r4$Iv_Inv<-r4$Iv+r3$Inv
r4_g<- graph2(r4,NULL,title="Rotavirus epidemic 100% vaccination")
I_g4<-graph2(r4,"Iv_Inv", title="Rotavirus epidemic, 100% vaccination, infected people (Iv+InV)")
Iv_Inv_g4<-graph2(r4,c("Iv","Inv"),title="Rotavirus epidemic, 100% vaccination, infected people")
prop_I4=r4%>%
  mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
  select(propI)%>%
  pull
r4$propI<-prop_I4
propI4_g<-graph2(r4,"propI","Proportion od people infected by Rotavirus \nwith 100% Vaccination Coverage for Rotavirus")
grid.arrange(I_g4,Iv_Inv_g4,ncol=1)

I_vac_100<-approxfun(r4$time,r4%>%
                      mutate(propI=Iv_Inv/(Sv+Iv+Rv+Snv+Inv+Rnv))%>%
                      select(propI)%>%
                      pull)

vec_virus_0<-function(time){
  return(0)
}

grid.arrange(r1_g,r2_g,r3_g,ncol=2)
grid.arrange(I_g1,I_g2,I_g3,ncol=2)



results_df <- data.frame(vacc = numeric(), Max_I = numeric())
I_vac<-list()
for (i in seq(0,1,by=0.05)){
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


I_R<-ggplot() +   
  geom_point(data=results_df, aes(x=vacc,y=max_propI,colour="Infected people at Epidemic peak"))+
  geom_line(data=results_df, aes(x=vacc,y=max_propI, colour="Infected people at Epidemic peak"))+
  geom_point(data=results_df, aes(x=vacc,y=last_propR, colour="Cumulative incidence"))+
  geom_line(data=results_df, aes(x=vacc,y=last_propR, colour="Cumulative incidence"))+
  labs(title = "Infected people at Epidemic peak \nand cumulative incidence according to vaccination", y = "Proportion of Individuals",
       x = "Vaccine coverage",size=6) +
  scale_colour_manual(name = "Legend", values = c("Infected people at Epidemic peak" = "purple", "Cumulative incidence" = "orange")) +
  theme_bw()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold",hjust = 0.5))


grid.arrange(propI1_g,propI2_g,propI3_g,I_R,ncol=2)

combined_propI<-data.frame(time=seq(from=0,to=365,by=1),no_vaccination=r1$propI,vaccination_50=r2$propI,vaccination_80=r3$propI)
graph2(combined_propI,NULL,"Proportion of People infected by Rotavirus")
