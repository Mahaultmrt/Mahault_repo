library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(gt)
library(tidyr)
library(forcats)

#Code modèle S.Aureus

create_params<-function(beta=1.65*10^-2,ct=0.90,deltaRa=0,deltaSa=0,gamma=1.02*10^-2,rhoR=8.22*10^-6,rhoS=1.20*10^-4,rhoRa=8.22*10^-6,rhoSa=1.20*10^-4,teta=0.0014,omega=0.08, alpha=0.33, sigmaR=1,sigmaS=0, ATB=0.1)
{
  list(beta=beta,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rhoR=rhoR,rhoS=rhoS,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB)
}

create_initial_cond<-function(Sa0=100000*0.02*0.7,CRa0=100000*0.02*0.03,CSa0=100000*0.02*0.27,IRa0=0,ISa0=0,S0=100000*0.98*0.7,CR0=100000*0.98*0.03,CS0=100000*0.98*0.27){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0)
}

run<-function(Init.cond,param,Tmax=365,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model, param,vec_virus=vec_virus))
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
run0_g<-graph(run0,c("Sa","CRa","CSa","S","CR","CS"),"S.Aureus Colonization dynamics \n(without virus epidemic)")
CR_CS0<-graph2(run0,c("CR_tot","CS_tot","C_tot"),"S.Aureus colonized people without a virus epidemic and without IPD")

# pas d'épidémie
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond()
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL,title="S.Aureus Colonization dynamics \nwithout virus épidemics")
propC1<-graph(run1,c("CRa","CSa","CR","CS"),"S.Aureus colonized people without virus epidemic")
CR_CS1<-graph2(run1,c("CR_tot","CS_tot"),"S.Aureus colonized people without a virus epidemic")


# Epidémie de grippe mais pas de vaccination
param<-create_params(rhoR=0,rhoS=0)
Init.cond<-create_initial_cond()
run0<-run(Init.cond,param)

vec_virus=I_vac_0
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL,title="S.Aureus Colonization dynamics \nwith influenza epidemic, no vaccination")
propC2<-graph(run2,c("CRa","CSa","CR","CS"),"S.Aureus colonized people with influenza epidemic, no vaccination")
CR_CS2<-graph2(run2,c("CR_tot","CS_tot"),"S.Aureus colonized people with influenza epidemic, no vaccination")


# Epidémie de grippe vaccination 50%
vec_virus=I_vac_50
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run3<-run(Init.cond,param)
run3_g<-graph(run3,NULL,title="S.Aureus Colonization dynamics \nwith influenza epidemic and vaccine coverage at 50%")
propC3<-graph(run3,c("CRa","CSa","CR","CS"),"S.Aureus colonized people with influenza epidemic and vaccine coverage at 50%")
CR_CS3<-graph2(run3,c("CR_tot","CS_tot"),"S.Aureus colonized people with influenza epidemic and vaccine coverage at 50%")
tetas<-graph2(run3,c("teta","new_teta"),"Parameters teta for S.Aureus colonization with 50% of vaccination for influenza")




# Epidémie de grippe vaccination 80%
vec_virus=I_vac_80
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run4<-run(Init.cond,param)
run4_g<-graph(run4,NULL,title="S.Aureus Colonization dynamics \nwith influenza epidemic and vaccine coverage at 80%")
propC4<-graph(run4,c("CRa","CSa","CR","CS"),"S.Aureus colonized people with influenza epidemic and vaccine coverage at 80%")
CR_CS4<-graph2(run4,c("CR_tot","CS_tot"),"S.Aureus colonized people with influenza epidemic and vaccine coverage at 80%")

grid.arrange(run0_g,run2_g,run3_g,run4_g,ncol=2)



IR_final <- data.frame(vacc = numeric(), LastIR = numeric())
IS_final<- data.frame(vacc = numeric(), LastIS = numeric())
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                 IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
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


corr_vacc_ATB_ISIR<- data.frame(vacc = numeric(), ATB=numeric(), LastISIR = numeric(), LastpropIR=numeric())
for (i in seq(1,21,by=1)){
  for(j in seq(0,0.5,by=0.1)){
    
    vec_virus=I_vac[[i]]
    param<-create_params(ATB=j)
    Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runt<-run(Init.cond,param)
    runt$ISIR=runt$ISa+runt$IRa
    LastIRa=runt[["IRa"]][nrow(runt) - 1]
    LastISIR=runt[["ISIR"]][nrow(runt) - 1]
    LastpropIR=round(LastIRa*100/LastISIR,digits=0)
    new_row=data.frame(vacc=results_df[i,1], ATB=j, LastISIR,LastpropIR)
    corr_vacc_ATB_ISIR <- bind_rows(corr_vacc_ATB_ISIR, new_row)
    
    
    
  }
  
}


h1<-heatmap(corr_vacc_ATB_ISIR,"vacc","ATB","LastISIR","vaccine coverage","Antibiotics",
        "Cumulative incidence of infection (per 100,000)", "Cumulative incidence of infection \ndepending on the vaccine coverage and the proportion of ATB",values=TRUE,var_text="LastpropIR")



data0<-percentage_final(run0)

run2$vaccination<-"vacc 0%"
run3$vaccination<-"vacc 50%"
run4$vaccination<-"vacc 80%"

all_run<-bind_rows(run2, run3, run4)
all_run <- melt(all_run, id.vars = c("time", "vaccination"))


res_graphs<-all_graph(all_run,NULL)+
  geom_hline(yintercept=tail(run0$Sa, n = 1), linetype="dashed",color="#499124",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CRa, n = 1), linetype="dashed",color="#DE6F00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CSa, n = 1), linetype="dashed",color="#2072BA",alpha=0.5)+
  geom_hline(yintercept=tail(run0$IRa, n = 1), linetype="dashed",color="#BD5E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$ISa, n = 1), linetype="dashed",color="#163F9E",alpha=0.5)+
  geom_hline(yintercept=tail(run0$S, n = 1), linetype="dashed",color="#68CF33",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CR, n = 1), linetype="dashed",color="#FC7E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CS, n = 1), linetype="dashed",color="#2B9CFF",alpha=0.5)
