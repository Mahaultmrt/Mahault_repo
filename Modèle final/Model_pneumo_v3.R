library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(tidyr)
library(forcats)
library(ppcor)
library(tibble)

library(epiR)

#Code modèle S.Pneumoniae

create_params<-function(beta=0.065,ct=0.96,deltaRa=0,deltaSa=0,gamma=0.05,rhoR=3*10^-6,rhoS=3*10^-6,rhoRa=3.0*10^-6,rhoSa=3*10^-6,teta=0.0014,omega=0.08, alpha=0.33, sigmaR=1,sigmaS=0,ATB=0.1){
  list(beta=beta,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rhoR=rhoR,rhoS=rhoS,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB)
}

create_initial_cond<-function(Sa0=100000*0.02*0.8,CRa0=100000*0.02*0.04,CSa0=100000*0.02*0.16,IRa0=0,ISa0=0,S0=100000*0.98*0.8,CR0=100000*0.98*0.04,CS0=100000*0.98*0.16){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0)
}

run<-function(Init.cond,param,Tmax=365,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model, param,vec_virus=vec_virus))
  tot <- sum(result[1, !(colnames(result) %in% c("time", "IRa", "ISa","N","new_teta"))])
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
run0_g<-graph(run0,c("Sa","CRa","CSa","S","CR","CS"),"S.Pneumoniae Colonization dynamics \n(without virus epidemic)")
CR_CS0<-graph2(run0,c("CR_tot","CS_tot","C_tot"),"S.Pneumoniae colonized people without a virus epidemic and without IPD")

#Pas d'épidémie pas de vaccination et infection
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL,title="S.Pneumoniae Colonization dynamics \nwithout virus épidemic")
propC1<-graph(run1,c("CRa","CSa","CR","CS"),"S.Pneumoniae colonized people without virus epidemic")
CR_CS1<-graph2(run1,c("CR_tot","CS_tot","C_tot"),"S.Pneumoniae colonized people without a virus epidemic")


# Epidémie de grippe mais pas de vaccination
vec_virus=I_vac_0
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL,title="S.Pneumoniae Colonization dynamics \nwith influenza epidemic, no vaccination")
propC2<-graph(run2,c("CRa","CSa","CR","CS"),"S.Pneumoniae colonized people with influenza epidemic, no vaccination")
CR_CS2<-graph2(run2,c("CR_tot","CS_tot","C_tot"),"S.Pneumoniae colonized people with influenza epidemic, no vaccination")


# Epidémie de grippe mais vaccination 50%
vec_virus=I_vac_50
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run3<-run(Init.cond,param)
run3_g<-graph(run3,NULL,title="S.Pneumoniae Colonization dynamics \nwith influenza epidemic and vaccine coverage at 50%")
propC3<-graph(run3,c("CRa","CSa","CR","CS"),"S.Pneumoniae colonized people with influenza epidemic and vaccine coverage at 50%")
CR_CS3<-graph2(run3,c("CR_tot","CS_tot","C_tot"),"S.Pneumoniae colonized peoplewith influenza epidemic and vaccine coverage at 50%")
tetas<-graph3(run3,c("teta","new_teta"),"Parameter teta for S.pneumonia colonization with 50% of vaccination for influenza")


# Epidémie de grippe mais vaccination 80%
vec_virus=I_vac_80
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run4<-run(Init.cond,param)
run4_g<-graph(run4,NULL,title="S. pneumoniae Colonization dynamics \nwith influenza epidemic and vaccine coverage at 80%")
propC4<-graph(run4,c("CRa","CSa","CR","CS"),"S.Pneumoniae colonized people with influenza epidemicand vaccine coverage at 80%")
CR_CS4<-graph2(run4,c("CR_tot","CS_tot","C_tot"),"S.Pneumoniae colonized people with influenza epidemic and vaccine coverage at 80%")

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
    LastpropIR=round(LastIRa*100/LastISIR,digits =0)
    new_row=data.frame(vacc=results_df[i,1], ATB=j, LastISIR,LastpropIR)
    corr_vacc_ATB_ISIR <- bind_rows(corr_vacc_ATB_ISIR, new_row)
    
    
    
  }
  
}


h1<-heatmap(corr_vacc_ATB_ISIR,"vacc","ATB","LastISIR","Vaccine coverage","Antibiotics",
            "Cumulative incidence of IPDs (per 100,000)", "Cumulative incidence of IPDs depending on the vaccine coverage \nand the proportion of ATB",values=TRUE,var_text="LastpropIR")


data0<-percentage_final(run0)

run2bis<-run2
run3bis<-run3
run4bis<-run4

run2bis$vaccination<-"vacc 0%"
run3bis$vaccination<-"vacc 50%"
run4bis$vaccination<-"vacc 80%"

all_run<-bind_rows(run2bis, run3bis, run4bis)
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


diff<- data.frame(vacc = numeric(), diffIR=numeric(), diffIS = numeric(),diffISIR=numeric())
for (i in seq(1,21,by=1)){

  diffIR=(IR_final$LastIR[i+1]-IR_final$LastIR[1])/tail(run0$IRa, n=1)
  diffIS=(IS_final$LastIS[i+1]-IS_final$LastIS[1])/tail(run0$ISa, n=1)
  diffISIR=((IR_final$LastIR[i+1]+IS_final$LastIS[i+1])-(IR_final$LastIR[1]+IS_final$LastIS[1]))/(tail(run0$IRa, n=1)+tail(run0$ISa, n=1))
  new_row=data.frame(vacc=results_df[i,1], diffIR, diffIS,diffISIR)
  diff <- bind_rows(diff, new_row)
  
}



diff_graph(diff)


psa<-data.frame(beta = numeric(), ct=numeric(), gamma = numeric(),alpha = numeric(), teta=numeric(), omega = numeric(),ATB = numeric(), incidenceR=numeric(),incidenceS=numeric())

psa=data.frame(beta=runif(1000,0.9*0.065,1.1*0.065),
               ct=runif(1000,0.9*0.96,1.1*0.96),
               gamma=runif(1000,0.9*0.05,1.1*0.05),
               alpha=runif(1000,0.9*0.33,1.1*0.33),
               teta=runif(1000,0.9*0.0014,1.1*0.0014),
               omega=runif(1000,0.9*0.08,1.1*0.08),
               ATB=runif(1000,0.9*0.1,1.1*0.1),
               incidenceR=0,
               incidenceS=0)


for (i in seq(1,1000,by=1)){
  beta<-psa$beta[i]
  ct<-psa$ct[i]
  gamma<-psa$gamma[i]
  alpha<-psa$alpha[i]
  teta<-psa$teta[i]
  omega<-psa$omega[i]
  ATB<-psa$ATB[i]
  
  vec_virus=I_vac_50
  param<-create_params(beta=beta,ct=ct,gamma=gamma,alpha=alpha,teta=teta,omega=omega,ATB=ATB)
  Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                 IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                                 CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
  runi<-run(Init.cond,param)
  psa$incidenceR[i]=runi[["IRa"]][nrow(runi) - 1]
  psa$incidenceS[i]=runi[["ISa"]][nrow(runi) - 1]
  
}



psa_graph<-density_graph(psa)+
  geom_vline(xintercept=run3[["IRa"]][nrow(run3) - 1],linetype="dashed",color="#BD5E00")+
  geom_vline(xintercept=run3[["ISa"]][nrow(run3) - 1],linetype="dashed",color="#163F9E")+
  geom_vline(xintercept=run3[["IRa"]][nrow(run3) - 1]+run3[["ISa"]][nrow(run3) - 1],linetype="dashed",color="#4B0082")


psa_R<-psa[,-c(9)]
psa_R = psa_R %>%
  dplyr::select(beta,ct,gamma,alpha,teta,omega,ATB,incidenceR) %>%
  epi.prcc() %>%
  rename(param = var)

# and we plot the results
ggplot(psa_R) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with cumulative incidence") +
  scale_x_discrete(labels = c(bquote(beta),
                              bquote(ct),
                              bquote(gamma),
                              bquote(alpha),
                              bquote(teta),
                              bquote(omega),
                              bquote(ATB)))


psa_S<-psa[,-c(8)]
psa_S = psa_S %>%
  dplyr::select(beta,ct,gamma,alpha,teta,omega,ATB,incidenceS) %>%
  epi.prcc() %>%
  rename(param = var)

ggplot(psa_S) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with cumulative incidence") +
  scale_x_discrete(labels = c(bquote(beta),
                              bquote(ct),
                              bquote(gamma),
                              bquote(alpha),
                              bquote(teta),
                              bquote(omega),
                              bquote(ATB)))

p_cor<-pcor(psa)

estimate_pcor<-p_cor$estimate
names<- c("beta","ct","gamma","alpha","teta","omega","ATB","incidenceR","incidenceS")
dimnames(estimate_pcor)<-list(names,names)

estimate_pcor <- as.data.frame(estimate_pcor)

cor_incidenceR <- data.frame(
  Parameters =c("beta","ct","gamma","alpha","teta","omega","ATB","incidenceR","incidenceS") ,
  Value = estimate_pcor$incidenceR,
  Error = 0.1 
)

cor_incidenceR<-cor_incidenceR[-c(8,9),-c(8,9)]


p_cor_graph<-function(data){
  ggplot(data, aes(x=Parameters,y=Value))+
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Value-Error,ymax=Value+Error),width=0.3)+
    labs(title="Correlation between incidence of infection (resistant strain) and parameters", x="Parameters",y="Value")+
    theme_minimal()+
    theme_bw()

}

p_cor_graph(cor_incidenceR)



cor_incidenceS <- data.frame(
  Parameters =c("beta","ct","gamma","alpha","teta","omega","ATB","incidenceR","incidenceS") ,
  Value = estimate_pcor$incidenceS,
  Error = 0.1 
)

cor_incidenceS<-cor_incidenceS[-c(8,9),-c(8,9)]

p_cor_graph(cor_incidenceS)




