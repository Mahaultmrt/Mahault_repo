library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(tidyr)
library(forcats)
library(tibble)
library(epiR)
library(cowplot)

source(here::here("Scripts", "Equations.R"))
source(here::here("Scripts", "Graphs.R"))
source(here::here("Scripts", "SIR_Influenza.R"))

#Code modèle S.Pneumoniae
create_params<-function(beta=0.065,fitness=0.974,deltaRa=0,deltaSa=0,gamma=0.05,rhoR=3*10^-6,rhoS=3*10^-6,rhoRa=3.0*10^-6,rhoSa=3*10^-6,theta=0.0015,omega=0.08, alpha=1/7, sigmaR=1,sigmaS=0,ATB=0.1,extra_rho=50){
  list(beta=beta,fitness=fitness,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rhoR=rhoR,rhoS=rhoS,rhoRa=rhoRa,rhoSa=rhoSa,theta=theta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB,extra_rho=extra_rho)
}

create_initial_cond<-function(Sa0=100000*0.02*0.8,CRa0=100000*0.02*0.04,CSa0=100000*0.02*0.16,IRa0=0,ISa0=0,S0=100000*0.98*0.8,CR0=100000*0.98*0.04,CS0=100000*0.98*0.16,Inc0=0){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0,Inc=Inc0)
}

run<-function(Init.cond,param,Tmax=140,dt=1, vec_virus=vec_virus){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model_pneumo, param,vec_virus=vec_virus))
  tot <- sum(result[1, !(colnames(result) %in% c("time", "IRa", "ISa","N","new_theta"))])
  proportion <- result
  proportion[ , -1] <- proportion[ , -1] / tot
  
  return(proportion)
}


#Pas d'épidémie pas de vaccination et pas d'infection
# calibrate fitness to have an equilibrium
is_it_stable = function(f_to_test){
  param<-create_params(fitness=f_to_test)
  Init.cond<-create_initial_cond()
  run_to_test<-run(Init.cond,param, vec_virus = vec_virus_0, Tmax = 10000)
  
  return(abs(diff(c(run_to_test$CR[9999], run_to_test$CR[10000]))))
  # return(abs(diff(c(run_to_test$CR[4999], run_to_test$CR[5000])))/run_to_test$CR[4999])
}

f_to_use = optim(0.974, is_it_stable, method="L-BFGS-B", lower=0.95, upper=0.999)

create_params<-function(beta=0.065,fitness=f_to_use$par,deltaRa=0,deltaSa=0,gamma=0.05,rhoR=3*10^-6,rhoS=3*10^-6,rhoRa=3.0*10^-6,rhoSa=3*10^-6,theta=0.0015,omega=0.08, alpha=1/7, sigmaR=1,sigmaS=0,ATB=0.1,extra_rho=50){
  list(beta=beta,fitness=fitness,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rhoR=rhoR,rhoS=rhoS,rhoRa=rhoRa,rhoSa=rhoSa,theta=theta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB,extra_rho=extra_rho)
}

param<-create_params()
Init.cond<-create_initial_cond()
run0p<-run(Init.cond,param,vec_virus=vec_virus_0, Tmax=10000)
Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                               IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                               CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
run0p<-run(Init.cond,param,vec_virus=vec_virus_0, Tmax=140)
run0_g<-graph(run0p,c("Sa","CRa","CSa","S","CR","CS"),"S.Pneumoniae Colonization dynamics \n(without virus epidemic)")
run0_g

#Pas d'épidémie pas de vaccination et infection
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                               IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                               CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
run1p<-run(Init.cond,param,vec_virus=vec_virus_0)
run1_g<-graph(run1p,NULL,title="S.Pneumoniae Colonization dynamics \nwithout virus épidemic")


# Epidémie de grippe mais pas de vaccination
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                               IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                               CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
run2p<-run(Init.cond,param,vec_virus=I_vac_0_flu)
run2_g<-graph(run2p,NULL,title="S.Pneumoniae Colonization dynamics \nwith influenza epidemic, no vaccination")


# Epidémie de grippe mais vaccination 25%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                               IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                               CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
run3p<-run(Init.cond,param,vec_virus=I_vac_25_flu)
run3_g<-graph(run3p,NULL,title="S.Pneumoniae Colonization dynamics \nwith influenza epidemic and vaccine coverage at 25%")


# Epidémie de grippe mais vaccination 50%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                               IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                               CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
run4p<-run(Init.cond,param,vec_virus=I_vac_50_flu)
run4_g<-graph(run4p,NULL,title="S. pneumoniae Colonization dynamics \nwith influenza epidemic and vaccine coverage at 50%")

# Tableau avec les différentes incidences cumulées des infections en fonction de la couverture vaccinale
IR_final <- data.frame(vacc = numeric(), LastIR = numeric(), ratioIR=numeric()) #Souche résistante
IS_final<- data.frame(vacc = numeric(), LastIS = numeric()) # Souche sensible
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                                 IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                                 CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
  runt<-run(Init.cond,param,vec_virus=vec_virus)
  LastIR=runt[["IRa"]][nrow(runt) - 1]
  LastIS=runt[["ISa"]][nrow(runt) - 1]
  new_row2=data.frame(vacc=results_df[i,1], LastIS)
  IS_final <- bind_rows(IS_final, new_row2)
  col<-paste("vaccination",results_df[i,1])
  ratioIR=LastIR/(LastIR+LastIS)
  new_row=data.frame(vacc=results_df[i,1], LastIR,ratioIR)
  IR_final <- bind_rows(IR_final, new_row)
  
  
}

# Fusion des deux tableaux pour avoir les deux souches dans un même tableau
I_final<-merge(IR_final,IS_final,by="vacc")
I_final <- pivot_longer(I_final, cols = c(LastIR,LastIS), names_to = "Strain", values_to = "Value")
I_final$Strain <- factor(I_final$Strain, levels = c("LastIS","LastIR"))

# Tableau des incidences cumulées des infections (S+R) en fonction de la vaccinaion et ATB
corr_vacc_ATB_ISIR<- data.frame(vacc = numeric(), ATB=numeric(), LastISIR = numeric(), LastpropIR=numeric())
for (i in seq(1,21,by=1)){
  for(j in seq(0,0.5,by=0.1)){
    
    vec_virus=I_vac[[i]]
    param<-create_params(ATB=j)
    Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                                   IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                                   CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
    runt<-run(Init.cond,param,vec_virus=vec_virus)
    runt$ISIR=runt$ISa+runt$IRa
    LastIRa=runt[["IRa"]][nrow(runt) - 1]
    LastISIR=runt[["ISIR"]][nrow(runt) - 1]
    LastpropIR=round(LastIRa*100/LastISIR,digits =0)
    new_row=data.frame(vacc=results_df[i,1], ATB=j, LastISIR,LastpropIR)
    corr_vacc_ATB_ISIR <- bind_rows(corr_vacc_ATB_ISIR, new_row)
    
  }
  
}


# Tableau des pourcentages d'individus susceptibles, colonisés, exposés et non exposés aux antibiotiques
data0<-percentage_final(run0p)

# Création du tableau pour afficher les sorties du modèle sur un graphique
run2bis<-run2p
run3bis<-run3p
run4bis<-run4p

run2bis$vaccination<-"vacc 0%"
run3bis$vaccination<-"vacc 25%"
run4bis$vaccination<-"vacc 50%"

all_run<-bind_rows(run2bis, run3bis, run4bis)
all_run <- melt(all_run, id.vars = c("time", "vaccination"))

all_runp=all_run%>% filter(variable!="Inc")

# Graphique avec les trois sorties du modèle en fonction des vaccination
res_graphs_P<-all_graph(all_runp,NULL)+
  geom_hline(yintercept=tail(run0p$Sa, n = 1), linetype="dashed",color="#499124",alpha=0.5)+
  geom_hline(yintercept=tail(run0p$CRa, n = 1), linetype="dashed",color="#DE6F00",alpha=0.5)+
  geom_hline(yintercept=tail(run0p$CSa, n = 1), linetype="dashed",color="#2072BA",alpha=0.5)+
  geom_hline(yintercept=tail(run0p$IRa, n = 1), linetype="dashed",color="#BD5E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0p$ISa, n = 1), linetype="dashed",color="#163F9E",alpha=0.5)+
  geom_hline(yintercept=tail(run0p$S, n = 1), linetype="dashed",color="#68CF33",alpha=0.5)+
  geom_hline(yintercept=tail(run0p$CR, n = 1), linetype="dashed",color="#FC7E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0p$CS, n = 1), linetype="dashed",color="#2B9CFF",alpha=0.5)


# Tableau des différence relative selon la vaccination (S et R)
diff<- data.frame(vacc = numeric(), diffIR=numeric(), diffIS = numeric(),diffISIR=numeric(),diffratio=numeric())
for (i in seq(0,20,by=1)){
  
  diffIR=(IR_final$LastIR[i+1]-IR_final$LastIR[1])/IR_final$LastIR[1]
  diffIS=(IS_final$LastIS[i+1]-IS_final$LastIS[1])/IS_final$LastIS[1]
  diffISIR=((IR_final$LastIR[i+1]+IS_final$LastIS[i+1])-(IR_final$LastIR[1]+IS_final$LastIS[1]))/(IR_final$LastIR[1]+IS_final$LastIS[1])
  diffratio=(IR_final$ratioIR[i+1]-IR_final$ratioIR[1])/IR_final$ratioIR[1]
  new_row=data.frame(vacc=results_df[i+1,1], diffIR, diffIS,diffISIR,diffratio)
  diff <- bind_rows(diff, new_row)
  
  
}


# Calcule de l'exposition finale aux antibiotiques selon la vaccination
exp_finalp<- data.frame(vacc = numeric(), exp = numeric())
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                                 IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                                 CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
  runt<-run(Init.cond,param,vec_virus=vec_virus)
  exp=runt[["Inc"]][nrow(runt) - 1]
  new_row=data.frame(vacc=results_df[i,1], exp)
  exp_finalp <- bind_rows(exp_finalp, new_row)
  col<-paste("vaccination",results_df[i,1])
  
}

# Différences relatives des expostions finales aux antibiotiques
diff_expp<- data.frame(vacc = numeric(), diffexp=numeric())
for (i in seq(0,20,by=1)){
  
  diffexp=(exp_finalp$exp[i+1]-exp_finalp$exp[1])/exp_finalp$exp[1]
  new_row=data.frame(vacc=results_df[i+1,1], diffexp)
  diff_expp <- bind_rows(diff_expp, new_row)
  
  
}

# Courbes des différences relatives d'expositions aux antibiotiques en fonction de la vaccination
graph_diff_expP<-graph_exp(diff_expp,diff)

# Calcules des ratio infection/exposition aux antibiotiques
ratio_P=data.frame(vacc=seq(0,1,by=0.05),ratio_exp_R=NA,ratioR_RS_exp=NA)
ratio_P$ratio_exp_R=diff$diffIR/diff_expp$diffexp
ratio_P$ratioR_RS_exp=(diff$diffratio)/diff_expp$diffexp
# ratio_P$ratio_exp_R=diff_exp$diffexp/diff$diffIR



# Analyse de sensibilité probabiliste
psa_vacc_50p=data.frame(beta=runif(500,0.8*0.065,1.2*0.065),
                        fitness=runif(500,0.8*f_to_use$par,min(1,1.2*f_to_use$par)),
                        gamma=runif(500,0.8*0.05,1.2*0.05),
                        rho=runif(500,0.8*3*10^-6,1.2*3*10^-6),
                        extra_rho=runif(500,0.8*50,1.2*50),
                        alpha=runif(500,0.8*1/7,1.2*1/7),
                        theta=runif(500,0.8*0.0015,1.2*0.0015),
                        omega=runif(500,0.8*0.08,1.2*0.08),
                        ATB=runif(500,0.8*0.1,1.2*0.1),
                        incidenceR=0,
                        incidenceS=0,
                        incidenceT=0)

for (i in seq(1,500,by=1)){
  beta<-psa_vacc_50p$beta[i]
  fitness<-psa_vacc_50p$fitness[i]
  gamma<-psa_vacc_50p$gamma[i]
  rho<-psa_vacc_50p$rho[i]
  extra_rho<-psa_vacc_50p$extra_rho[i]
  alpha<-psa_vacc_50p$alpha[i]
  theta<-psa_vacc_50p$theta[i]
  omega<-psa_vacc_50p$omega[i]
  ATB<-psa_vacc_50p$ATB[i]
  
  param<-create_params(beta=beta,fitness=fitness,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB,rhoR=rho,rhoS=rho,rhoRa=rho,rhoSa=rho,extra_rho=extra_rho)
  Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                                 IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                                 CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
  runi0<-run(Init.cond,param, vec_virus = I_vac_0_flu)
  runi<-run(Init.cond,param, vec_virus = I_vac_50_flu)
  psa_vacc_50p$incidenceR[i]=(runi[["IRa"]][nrow(runi) - 1] - runi0[["IRa"]][nrow(runi0) - 1]) / runi0[["IRa"]][nrow(runi0) - 1]
  psa_vacc_50p$incidenceS[i]=(runi[["ISa"]][nrow(runi) - 1] - runi0[["ISa"]][nrow(runi0) - 1]) / runi0[["ISa"]][nrow(runi0) - 1]
  psa_vacc_50p$incidenceT[i]=( (runi[["IRa"]][nrow(runi) - 1]+runi[["ISa"]][nrow(runi) - 1]) - (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1]) ) / (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1])
  
}

univar_vacc_50p = expand.grid(alpha = seq(0.11, 0.17, 0.01),
                              # attr_exp = seq(0.12, 0.18, 0.01),
                              ATB = seq(0.07, 0.13, 0.01),
                              # ATB = 0,
                              incidenceR=0,
                              incidenceS=0,
                              incidenceT=0)

# find_ATB_for_exp = function(x, target_exp = 0.15){
#   
#   (target_exp - mean(-log(1-I_vac_0_flu(0:140)*x)/(-log(1-I_vac_0_flu(0:140)*x)+0.0015)))^2
#   
# }
# 
# univar_vacc_50p$ATB = sapply(univar_vacc_50p$attr_exp, function(x) optim(0.1, find_ATB_for_exp, method="BFGS", target_exp=x)$par)

for(i in c(1:nrow(univar_vacc_50p))){
  param<-create_params(alpha=univar_vacc_50p$alpha[i],ATB=univar_vacc_50p$ATB[i])
  Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                                 IRa0=0,ISa0=0,S0=tail(run0a$S, n = 1),
                                 CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
  runi0<-run(Init.cond,param, vec_virus = I_vac_0_flu)
  runi<-run(Init.cond,param, vec_virus = I_vac_50_flu)
  univar_vacc_50p$incidenceR[i]=(runi[["IRa"]][nrow(runi) - 1] - runi0[["IRa"]][nrow(runi0) - 1]) / runi0[["IRa"]][nrow(runi0) - 1]
  univar_vacc_50p$incidenceS[i]=(runi[["ISa"]][nrow(runi) - 1] - runi0[["ISa"]][nrow(runi0) - 1]) / runi0[["ISa"]][nrow(runi0) - 1]
  univar_vacc_50p$incidenceT[i]=( (runi[["IRa"]][nrow(runi) - 1]+runi[["ISa"]][nrow(runi) - 1]) - (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1]) ) / (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1])
}


# Tableau pour créer les courbes de différences relatives avec plusieurs simulations

# Souche résistante
IR_vaccp<-data.frame(matrix(ncol=20),nrow=0)
colnames(IR_vaccp)<-seq(0, 1, by = 0.05)
IR_vacc_bis<-IR_vaccp

IS_vaccp<-data.frame(matrix(ncol=20),nrow=0)
colnames(IS_vaccp)<-seq(0, 1, by = 0.05)
IS_vacc_bis<-IS_vaccp

ISR_vaccp<-data.frame(matrix(ncol=20),nrow=0)
colnames(ISR_vaccp)<-seq(0, 1, by = 0.05)
ISR_vacc_bis<-ISR_vaccp

abx_vaccp<-data.frame(matrix(ncol=20),nrow=0)
colnames(abx_vaccp)<-seq(0, 1, by = 0.05)

for (i in seq(1,21,by=1)){
  for (j in seq(1,500,by=1)){
    beta<-psa_vacc_50p$beta[j]
    fitness<-psa_vacc_50p$fitness[j]
    gamma<-psa_vacc_50p$gamma[j]
    rho<-psa_vacc_50p$rho[j]
    extra_rho<-psa_vacc_50p$extra_rho[j]
    alpha<-psa_vacc_50p$alpha[j]
    theta<-psa_vacc_50p$theta[j]
    omega<-psa_vacc_50p$omega[j]
    ATB<-psa_vacc_50p$ATB[j]
    
    vec_virus=I_vac[[i]]
    param<-create_params(beta=beta,fitness=fitness,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB,rhoR=rho,rhoS=rho,rhoRa=rho,rhoSa=rho,extra_rho=extra_rho)
    Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                                   IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                                   CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
    runi<-run(Init.cond,param,vec_virus = vec_virus)
    IR_vaccp[j,i]<-runi[["IRa"]][nrow(runi) - 1]
    IS_vaccp[j,i]<-runi[["ISa"]][nrow(runi) - 1]
    ISR_vaccp[j,i]<-runi[["IRa"]][nrow(runi) - 1]+runi[["ISa"]][nrow(runi) - 1]
    abx_vaccp[j,i]<-runi[["Inc"]][nrow(runi) - 1]
  }
  
}


for (j in seq(1,500,by=1)){
  for (i in seq(0,20,by=1)){
    IR_vacc_bis[j,i+1]<-(IR_vaccp[j,i+1]-IR_vaccp[j,1])/IR_vaccp[j,1]
    IS_vacc_bis[j,i+1]<-(IS_vaccp[j,i+1]-IS_vaccp[j,1])/IS_vaccp[j,1]
    ISR_vacc_bis[j,i+1]<-(ISR_vaccp[j,i+1]-ISR_vaccp[j,1])/ISR_vaccp[j,1]
  }
}

IR_vacc_bis <- gather(IR_vacc_bis, key = "vacc", value = "incidence")
IR_vacc_bis$vacc <- as.numeric(as.character(IR_vacc_bis$vacc))

IS_vacc_bis <- gather(IS_vacc_bis, key = "vacc", value = "incidence")
IS_vacc_bis$vacc <- as.numeric(as.character(IS_vacc_bis$vacc))

ISR_vacc_bis <- gather(ISR_vacc_bis, key = "vacc", value = "incidence")
ISR_vacc_bis$vacc <- as.numeric(as.character(ISR_vacc_bis$vacc))

diff_graph(diff,IR_vacc_bis,IS_vacc_bis,ISR_vacc_bis)

diff_simp<-data.frame(vacc=seq(0,1,by=0.05))

# Calculs des valeurs moyennes
mean_values_IR <- IR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_IR = mean(incidence, na.rm = TRUE))

diff_simp$mean_incidence_IR<-mean_values_IR$mean_incidence_IR

mean_values_IS <- IS_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_IS = mean(incidence, na.rm = TRUE))
diff_simp$mean_incidence_IS<-mean_values_IS$mean_incidence_IS

mean_values_ISR <- ISR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_ISR = mean(incidence, na.rm = TRUE))
diff_simp$mean_incidence_ISR<-mean_values_ISR$mean_incidence_ISR

# Calculs des écarts types
sd_values_IR <- IR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_IR = sd(incidence, na.rm = TRUE))
diff_simp$sd_incidence_IR<-sd_values_IR$sd_incidence_IR

sd_values_IS <- IS_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_IS = sd(incidence, na.rm = TRUE))
diff_simp$sd_incidence_IS<-sd_values_IS$sd_incidence_IS

sd_values_ISR <- ISR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_ISR = sd(incidence, na.rm = TRUE))
diff_simp$sd_incidence_ISR<-sd_values_ISR$sd_incidence_ISR

# Calculs des intervals de confiance
diff_simp$IR_ic_l=diff_simp$mean_incidence_IR-diff_simp$sd_incidence_IR
diff_simp$IR_ic_u=diff_simp$mean_incidence_IR+diff_simp$sd_incidence_IR
diff_simp$IS_ic_l=diff_simp$mean_incidence_IS-diff_simp$sd_incidence_IS
diff_simp$IS_ic_u=diff_simp$mean_incidence_IS+diff_simp$sd_incidence_IS
diff_simp$ISR_ic_l=diff_simp$mean_incidence_ISR-diff_simp$sd_incidence_ISR
diff_simp$ISR_ic_u=diff_simp$mean_incidence_ISR+diff_simp$sd_incidence_ISR

# Courbes des différences relatives avec les valeurs moyennes et intervalles de confiance
diff_sim_P<-diff_graph_sim(diff_simp)



## COMPARING EQUAL COVERAGE x EFFICACY #####

# Couverture vaccinale 70% efficacité 30%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                               IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                               CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
run70<-run(Init.cond,param,vec_virus=I_vac_70_flu)
IR70<-run70[["IRa"]][nrow(run70) - 1]
IS70<-run70[["ISa"]][nrow(run70) - 1]
ISR70<-run70[["ISa"]][nrow(run70) - 1]+run70[["IRa"]][nrow(run70) - 1]

# Couverture vaccinale 30% efficacité 70%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                               IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                               CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
run30<-run(Init.cond,param,vec_virus=I_vac_30_flu)
IR30<-run30[["IRa"]][nrow(run30) - 1]
IS30<-run30[["ISa"]][nrow(run30) - 1]
ISR30<-run30[["ISa"]][nrow(run30) - 1]+run30[["IRa"]][nrow(run30) - 1]

# Couverture vaccinale 46% efficacité 46%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0p$Sa, n = 1),CRa0=tail(run0p$CRa, n = 1),CSa0=tail(run0p$CSa, n = 1),
                               IRa0=0,ISa0=0,S0=tail(run0p$S, n = 1),
                               CR0=tail(run0p$CR, n = 1),CS0=tail(run0p$CS, n = 1))
run46<-run(Init.cond,param,vec_virus=I_vac_46_flu)
IR46<-run46[["IRa"]][nrow(run46) - 1]
IS46<-run46[["ISa"]][nrow(run46) - 1]
ISR46<-run46[["ISa"]][nrow(run46) - 1]+run46[["IRa"]][nrow(run46) - 1]
grid.arrange(graph(run70,NULL,"vaccine coverage 70%, vaccine efficacy 30%"),
             graph(run30,NULL,"vaccine coverage 30%, vaccine efficacy 70%"),graph(run46,NULL,"vaccine converage 46%, vaccine efficacy 46%"),ncol=2)

test_vacc<-data.frame(vacc=c(0.70,0.30,0.46),incidenceR=c(IR70,IR30,IR46),incidenceS=c(IS70,IS30,IS46))
test_vacc <- pivot_longer(test_vacc, cols = c(incidenceR,incidenceS), names_to = "Strain", values_to = "Value")


