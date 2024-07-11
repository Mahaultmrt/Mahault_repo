library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gt)
library(tidyr)

#Code modèle E.Coli

create_params<-function(beta=0.012,ct=0.96,deltaRa=0,deltaSa=0,gamma=0.01,rho=1.8*10^-6,rhoRa=1.8*10^-6,rhoSa=1.8*10^-6,teta=0.0014,omega=0.14, alpha=0.33, sigmaR=1, ATB=0.1,phi=9.83)
{
  list(beta=beta,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,ATB=ATB,phi=phi)
}

create_initial_cond<-function(CSa0=100000*0.01*0.8,CRa0=100000*0.01*0.2,CS0=100000*0.99*0.8,CR0=100000*0.99*0.2,IRa0=0,ISa0=0){
  c(CSa=CSa0,CRa=CRa0,CS=CS0,CR=CR0,IRa=IRa0,ISa=ISa0)
}

run<-function(Init.cond,param,Tmax=365,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model_Coli, param,vec_virus=vec_virus))
  tot <- sum(result[1, !(colnames(result) %in% c("time", "IRa", "ISa","N"))])
  proportion <- result
  proportion[ , -1] <- proportion[ , -1] / tot
  # proportion$CR_tot<-proportion$CRa+proportion$CR
  # proportion$CS_tot<-proportion$CSa+proportion$CS
  # proportion$C_tot<-proportion$CR_tot+proportion$CS_tot
  return(proportion)
  
  
}




#Pas d'épidémie pas de vaccination et pas d'infection
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond()
run0<-run(Init.cond,param)
run0_g<-graph(run0,c("CSa","CRa","CS","CR"),"E.Coli Colonization dynamics \nwithout virus epidemics")
CR_CS0<-graph2(run0,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people without a virus epidemic and without infection")


#Pas d'épidémie 
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL,title="E.Coli Colonization dynamics \nwithout virus épidemics")
propC1<-graph(run1,c("CRa","pCSa","CR","CS"),"E.Coli colonized people without a virus epidemic")
CR_CS1<-graph2(run1,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people without a virus epidemic")


# Epidémie de rotavirus mais pas de vaccination
vec_virus=I_vac_0
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL,title="E.Coli Colonization dynamics \nwith virus epidemic, no vaccination")
propC2<-graph(run2,c("CRa","CSa","CR","CS"),"E.Coli colonized people with rotavirus epidemic, no vaccination")
CR_CS2<-graph2(run2,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people with rotavirus epidemic, no vaccination")



# Epidémie de rotavirus vaccination 50%
vec_virus=I_vac_50
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run3<-run(Init.cond,param)
run3_g<-graph(run3,NULL,title="E.Coli Colonization dynamics \nwith rotavirus epidemic and vaccine coverage at 50%")
propC3<-graph(run3,c("CRa","CSa","CR","CS"),"E.Coli colonized people with rotavirus epidemic and vaccine coverage at 50%")
CR_CS3<-graph2(run3,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people with rotavirus epidemic and vaccine coverage at 50%")
tetas<-graph3(run3,c("teta","new_teta"),"Parameters teta for E.Coli colonization with 50% of vaccination for rotavirus")




# # Epidémie de rotavirus vaccination 80%
vec_virus=I_vac_80
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run4<-run(Init.cond,param)
run4_g<-graph(run4,NULL,title="E.Coli Colonization dynamics\nwith rotavirus epidemic and vaccine coverage at 80%")
propC4<-graph(run4,c("CRa","CSa","CR","CS"),"E.Coli colonized people with rotavirus epidemic and vaccine coverage at 80%")
CR_CS4<-graph2(run4,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people with rotavirus epidemic and vaccine coverage at 80%")



grid.arrange(run0_g,run2_g,run3_g,run4_g,ncol=2)


IR_final <- data.frame(vacc = numeric(), LastIR = numeric())
IS_final<- data.frame(vacc = numeric(), LastIS = numeric())
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                 IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                 CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
  runt<-run(Init.cond,param)
  LastIR=runt[["IRa"]][nrow(runt) - 1]
  new_row=data.frame(vacc=results_df[i,1], LastIR)
  IR_final <- bind_rows(IR_final, new_row)
  LastIS=runt[["ISa"]][nrow(runt) - 1]
  new_row2=data.frame(vacc=results_df[i,1], LastIS)
  IS_final <- bind_rows(IS_final, new_row2)
  col<-paste("vaccination",results_df[i,1])
  
}



I_final<-merge(IR_final,IS_final,by="vacc")
I_final <- pivot_longer(I_final, cols = c(LastIR,LastIS), names_to = "Strain", values_to = "Value")
I_final$Strain <- factor(I_final$Strain, levels = c("LastIS","LastIR"))

Cumulative_incidence<-graph_barplot(I_final)


data0<-percentage_final(run0)

run2$vaccination<-"vacc 0%"
run3$vaccination<-"vacc 50%"
run4$vaccination<-"vacc 80%"

all_run<-bind_rows(run2, run3, run4)
all_run <- melt(all_run, id.vars = c("time", "vaccination"))


res_graphs<-all_graph(all_run,NULL)+
  geom_hline(yintercept=tail(run0$CRa, n = 1), linetype="dashed",color="#DE6F00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CSa, n = 1), linetype="dashed",color="#2072BA",alpha=0.5)+
  geom_hline(yintercept=tail(run0$IRa, n = 1), linetype="dashed",color="#BD5E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$ISa, n = 1), linetype="dashed",color="#163F9E",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CR, n = 1), linetype="dashed",color="#FC7E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CS, n = 1), linetype="dashed",color="#2B9CFF",alpha=0.5)


diff<- data.frame(vacc = numeric(), diffIR=numeric(), diffIS = numeric())
for (i in seq(1,21,by=1)){
  
  diffIR=(IR_final$LastIR[i+1]-IR_final$LastIR[1])
  diffIS=(IS_final$LastIS[i+1]-IS_final$LastIS[1])
  new_row=data.frame(vacc=results_df[i,1], diffIR, diffIS)
  diff <- bind_rows(diff, new_row)
  
  
}

diff_graph(diff)
