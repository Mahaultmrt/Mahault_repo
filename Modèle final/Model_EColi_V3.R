library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gt)
library(tidyr)
library(epiR)


#Code modèle E.Coli

create_params<-function(beta=0.012,fcost=0.96,deltaRa=0,deltaSa=0,gamma=0.01,rho=1.8*10^-6,rhoRa=1.8*10^-6,rhoSa=1.8*10^-6,theta=0.0014,omega=0.14, alpha=0.33, sigmaR=1, ATB=0.1,phi=9.83)
{
  list(beta=beta,fcost=fcost,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,theta=theta,omega=omega,alpha=alpha,sigmaR=sigmaR,ATB=ATB,phi=phi)
}

create_initial_cond<-function(CSa0=100000*0.01*0.8,CRa0=100000*0.01*0.2,CS0=100000*0.99*0.8,CR0=100000*0.99*0.2,IRa0=0,ISa0=0,Inc0=0){
  c(CSa=CSa0,CRa=CRa0,CS=CS0,CR=CR0,IRa=IRa0,ISa=ISa0,Inc=Inc0)
}

run<-function(Init.cond,param,Tmax=365,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model_Coli, param,vec_virus=vec_virus))
  tot <- sum(result[1, !(colnames(result) %in% c("time", "IRa", "ISa","N","new_theta"))])
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
run0_E<-graph(run0,c("CSa","CRa","CS","CR"),NULL)
CR_CS0<-graph2(run0,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people without a virus epidemic and without infection")


#Pas d'épidémie 
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL,title="E.Coli Colonization dynamics \nwithout virus epidemics")
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
thetas<-graph3(run3,c("theta","new_theta"),"Parameters theta for E.Coli colonization with 50% of vaccination for rotavirus")


# Epidémie de rotavirus vaccination 80%
vec_virus=I_vac_80
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run4<-run(Init.cond,param)
run4_g<-graph(run4,NULL,title="E.Coli Colonization dynamics\nwith rotavirus epidemic and vaccine coverage at 80%")
propC4<-graph(run4,c("CRa","CSa","CR","CS"),"E.Coli colonized people with rotavirus epidemic and vaccine coverage at 80%")
CR_CS4<-graph2(run4,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people with rotavirus epidemic and vaccine coverage at 80%")



# Tableau avec les différentes incidences cumulées des infections en fonction de la couverture vaccinale
IR_final <- data.frame(vacc = numeric(), LastIR = numeric(), ratioIR=numeric())
IS_final<- data.frame(vacc = numeric(), LastIS = numeric())
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                 IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                 CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
  runt<-run(Init.cond,param)
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

# Barplot des différentes incidences cumulées des infections en fonction de la vaccination
Cumulative_incidence_E<-graph_barplot(I_final)

# Tableau des pourcentages d'individus susceptibles, colonisés, exposés et non exposés aux antibiotiques
data0<-percentage_final(run0)

# Création du tableau pour afficher les sorties du modèle sur un graphique
run2bis<-run2
run3bis<-run3
run4bis<-run4

run2bis$vaccination<-"vacc 0%"
run3bis$vaccination<-"vacc 50%"
run4bis$vaccination<-"vacc 80%"

all_run<-bind_rows(run2bis, run3bis, run4bis)
all_run <- melt(all_run, id.vars = c("time", "vaccination"))

all_run=all_run%>% filter(variable!="Inc")

# Graphique avec les trois sorties du modèle en fonction des vaccination
res_graphs_E<-all_graph(all_run,NULL)+
  geom_hline(yintercept=tail(run0$CRa, n = 1), linetype="dashed",color="#DE6F00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CSa, n = 1), linetype="dashed",color="#2072BA",alpha=0.5)+
  geom_hline(yintercept=tail(run0$IRa, n = 1), linetype="dashed",color="#BD5E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$ISa, n = 1), linetype="dashed",color="#163F9E",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CR, n = 1), linetype="dashed",color="#FC7E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CS, n = 1), linetype="dashed",color="#2B9CFF",alpha=0.5)

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


# Courbes des différences relatives selon la vaccination
diff_E<-diff_graph(diff,NULL,NULL)

# Calcule de l'exposition finale aux antibiotiques selon la vaccination
exp_final<- data.frame(vacc = numeric(), exp = numeric())
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                 IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                 CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
  runt<-run(Init.cond,param)
  exp=runt[["Inc"]][nrow(runt) - 1]
  new_row=data.frame(vacc=results_df[i,1], exp)
  exp_final <- bind_rows(exp_final, new_row)
  col<-paste("vaccination",results_df[i,1])
  
  
  
}

# Différences relatives des expostions finales aux antibiotiques
diff_exp<- data.frame(vacc = numeric(), diffexp=numeric())
for (i in seq(0,20,by=1)){
  
  diffexp=(exp_final$exp[i+1]-exp_final$exp[1])/exp_final$exp[1]
  new_row=data.frame(vacc=results_df[i+1,1], diffexp)
  diff_exp <- bind_rows(diff_exp, new_row)
  
  
}

# Courbes des différences relatives d'expositions aux antibiotiques en fonction de la vaccination
graph_diff_expE<-graph_exp(diff_exp,diff)

# Calcules des ratio infection/exposition aux antibiotiques
ratio=data.frame(vacc=seq(0,1,by=0.05),ratio_exp_R=NA,ratioR_RS_exp=NA)
ratio$ratio_exp_R=diff$diffIR/diff_exp$diffexp
ratio$ratioR_RS_exp=(diff$diffratio)/diff_exp$diffexp

# Analyse de sensibilité probabiliste
psa<-data.frame(beta = numeric(), fcost=numeric(), gamma = numeric(),alpha = numeric(), theta=numeric(), omega = numeric(),ATB = numeric(), incidenceR=numeric(),incidenceS=numeric())

psa=data.frame(beta=runif(1000,0.8*0.065,1.2*0.065),
               fcost=runif(1000,0.8*0.96,1.2*0.96),
               gamma=runif(1000,0.8*0.05,1.2*0.05),
               alpha=runif(1000,0.8*0.33,1.2*0.33),
               theta=runif(1000,0.8*0.0014,1.2*0.0014),
               omega=runif(1000,0.8*0.08,1.2*0.08),
               ATB=runif(1000,0.8*0.1,1.2*0.1),
               incidenceR=0,
               incidenceS=0)

for (i in seq(1,1000,by=1)){
  beta<-psa$beta[i]
  fcost<-psa$fcost[i]
  gamma<-psa$gamma[i]
  alpha<-psa$alpha[i]
  theta<-psa$theta[i]
  omega<-psa$omega[i]
  ATB<-psa$ATB[i]
  
  vec_virus=I_vac_50
  param<-create_params(beta=beta,fcost=fcost,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB)
  Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                 IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                 CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
  runi<-run(Init.cond,param)
  psa$incidenceR[i]=runi[["IRa"]][nrow(runi) - 1]
  psa$incidenceS[i]=runi[["ISa"]][nrow(runi) - 1]
  
}



psa_graph_E<-density_graph(psa)+
  geom_vline(xintercept=run3[["IRa"]][nrow(run3) - 1],linetype="dashed",color="#BD5E00")+
  geom_vline(xintercept=run3[["ISa"]][nrow(run3) - 1],linetype="dashed",color="#163F9E")+
  geom_vline(xintercept=run3[["IRa"]][nrow(run3) - 1]+run3[["ISa"]][nrow(run3) - 1],linetype="dashed",color="#4B0082")

# Analyse de correlation partielle
psa_R<-psa[,-c(9)]
psa_R = psa_R %>%
  dplyr::select(beta,fcost,gamma,alpha,theta,omega,ATB,incidenceR) %>%
  epi.prcc() %>%
  rename(param = var)

pcorR_E<-graph_pcor(psa_R)

psa_S<-psa[,-c(8)]
psa_S = psa_S %>%
  dplyr::select(beta,fcost,gamma,alpha,theta,omega,ATB,incidenceS) %>%
  epi.prcc() %>%
  rename(param = var)

pcorS_E<-graph_pcor(psa_S)


psa$incidenceSR<-psa$incidenceR+psa$incidenceS
psa_SR<-psa[,-c(8,9)]
psa_SR = psa_SR %>%
  dplyr::select(beta,fcost,gamma,alpha,theta,omega,ATB,incidenceSR) %>%
  epi.prcc() %>%
  rename(param = var)
pcorSR_E<-graph_pcor(psa_SR)



# Tableau pour créer les courbes de différences relatives avec plusieurs simulations

# Souche résistante
IR_vacc<-data.frame(matrix(ncol=20),nrow=0)
colnames(IR_vacc)<-seq(0, 1, by = 0.05)
IR_vacc_bis<-IR_vacc


for (i in seq(1,21,by=1)){
  for (j in seq(1,100,by=1)){
    beta<-psa$beta[j]
    fcost<-psa$fcost[j]
    gamma<-psa$gamma[j]
    alpha<-psa$alpha[j]
    theta<-psa$theta[j]
    omega<-psa$omega[j]
    ATB<-psa$ATB[j]
    
    vec_virus=I_vac[[i]]
    param<-create_params(beta=beta,fcost=fcost,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB)
    Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runi<-run(Init.cond,param)
    IR_vacc[j,i]<-runi[["IRa"]][nrow(runi) - 1]
  }
  
}


for (j in seq(1,100,by=1)){
  for (i in seq(0,20,by=1)){
    
    IR_vacc_bis[j,i+1]<-(IR_vacc[j,i+1]-IR_vacc[j,1])/IR_vacc[j,1]
  }
  
}

IR_vacc_bis <- gather(IR_vacc_bis, key = "vacc", value = "incidence")
IR_vacc_bis$vacc <- as.numeric(as.character(IR_vacc_bis$vacc))

# Souche sensible

IS_vacc<-data.frame(matrix(ncol=20),nrow=0)
colnames(IS_vacc)<-seq(0, 1, by = 0.05)
IS_vacc_bis<-IS_vacc

for (i in seq(1,21,by=1)){
  for (j in seq(1,100,by=1)){
    beta<-psa$beta[j]
    fcost<-psa$fcost[j]
    gamma<-psa$gamma[j]
    alpha<-psa$alpha[j]
    theta<-psa$theta[j]
    omega<-psa$omega[j]
    ATB<-psa$ATB[j]
    
    vec_virus=I_vac[[i]]
    param<-create_params(beta=beta,fcost=fcost,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB)
    Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runi<-run(Init.cond,param)
    IS_vacc[j,i]<-runi[["ISa"]][nrow(runi) - 1]
  }
  
}


for (j in seq(1,100,by=1)){
  for (i in seq(0,20,by=1)){
    
    IS_vacc_bis[j,i+1]<-(IS_vacc[j,i+1]-IS_vacc[j,1])/IS_vacc[j,1]
  }
  
}

IS_vacc_bis <- gather(IS_vacc_bis, key = "vacc", value = "incidence")
IS_vacc_bis$vacc <- as.numeric(as.character(IS_vacc_bis$vacc))

# Souche résistante + sensible

ISR_vacc<-data.frame(matrix(ncol=20),nrow=0)
colnames(ISR_vacc)<-seq(0, 1, by = 0.05)
ISR_vacc_bis<-ISR_vacc

for (i in seq(1,21,by=1)){
  for (j in seq(1,100,by=1)){
    beta<-psa$beta[j]
    fcost<-psa$fcost[j]
    gamma<-psa$gamma[j]
    alpha<-psa$alpha[j]
    theta<-psa$theta[j]
    omega<-psa$omega[j]
    ATB<-psa$ATB[j]
    
    vec_virus=I_vac[[i]]
    param<-create_params(beta=beta,fcost=fcost,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB)
    Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runi<-run(Init.cond,param)
    ISR_vacc[j,i]<-runi[["ISa"]][nrow(runi) - 1]+runi[["IRa"]][nrow(runi) - 1]
  }
  
}


for (j in seq(1,100,by=1)){
  for (i in seq(0,20,by=1)){
    
    ISR_vacc_bis[j,i+1]<-(ISR_vacc[j,i+1]-ISR_vacc[j,1])/ISR_vacc[j,1]
  }
  
}

ISR_vacc_bis <- gather(ISR_vacc_bis, key = "vacc", value = "incidence")
ISR_vacc_bis$vacc <- as.numeric(as.character(ISR_vacc_bis$vacc))

diff_graph(diff,IR_vacc_bis,IS_vacc_bis,ISR_vacc_bis)

diff_sim<-data.frame(vacc=seq(0,1,by=0.05))

# Calculs des valeurs moyennes
mean_values_IR <- IR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_IR = mean(incidence, na.rm = TRUE))

diff_sim$mean_incidence_IR<-mean_values_IR$mean_incidence_IR

mean_values_IS <- IS_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_IS = mean(incidence, na.rm = TRUE))
diff_sim$mean_incidence_IS<-mean_values_IS$mean_incidence_IS

mean_values_ISR <- ISR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_ISR = mean(incidence, na.rm = TRUE))
diff_sim$mean_incidence_ISR<-mean_values_ISR$mean_incidence_ISR

# Calculs des écarts types
sd_values_IR <- IR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_IR = sd(incidence, na.rm = TRUE))
diff_sim$sd_incidence_IR<-sd_values_IR$sd_incidence_IR

sd_values_IS <- IS_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_IS = sd(incidence, na.rm = TRUE))
diff_sim$sd_incidence_IS<-sd_values_IS$sd_incidence_IS

sd_values_ISR <- ISR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_ISR = sd(incidence, na.rm = TRUE))
diff_sim$sd_incidence_ISR<-sd_values_ISR$sd_incidence_ISR

# Calculs des intervals de confiance
diff_sim$IR_ic_l=diff_sim$mean_incidence_IR-diff_sim$sd_incidence_IR
diff_sim$IR_ic_u=diff_sim$mean_incidence_IR+diff_sim$sd_incidence_IR
diff_sim$IS_ic_l=diff_sim$mean_incidence_IS-diff_sim$sd_incidence_IS
diff_sim$IS_ic_u=diff_sim$mean_incidence_IS+diff_sim$sd_incidence_IS
diff_sim$ISR_ic_l=diff_sim$mean_incidence_ISR-diff_sim$sd_incidence_ISR
diff_sim$ISR_ic_u=diff_sim$mean_incidence_ISR+diff_sim$sd_incidence_ISR

# Courbes des différences relatives avec les valeurs moyennes et intervalles de confiance
diff_sim_E<-diff_graph_sim(diff_sim)




# Tests couverture vaccinale et efficacité vaccinale

# Couverture vaccinale 70% effecicité 30%
vec_virus=I_vac_70
param<-create_params()
Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run70<-run(Init.cond,param)
IR70<-run70[["IRa"]][nrow(run70) - 1]
IS70<-run70[["ISa"]][nrow(run70) - 1]
ISR70<-run70[["ISa"]][nrow(run70) - 1]+run70[["IRa"]][nrow(run70) - 1]

# Couverture vaccinale 30% effecicité 70%
vec_virus=I_vac_30
param<-create_params()
Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run30<-run(Init.cond,param)
IR30<-run30[["IRa"]][nrow(run30) - 1]
IS30<-run30[["ISa"]][nrow(run30) - 1]
ISR30<-run30[["ISa"]][nrow(run30) - 1]+run30[["IRa"]][nrow(run30) - 1]

# Couverture vaccinale 46% effeicité 46%
vec_virus=I_vac_46
param<-create_params()
Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run46<-run(Init.cond,param)
IR46<-run46[["IRa"]][nrow(run46) - 1]
IS46<-run46[["ISa"]][nrow(run46) - 1]
ISR46<-run46[["ISa"]][nrow(run46) - 1]+run46[["IRa"]][nrow(run46) - 1]
grid.arrange(graph(run70,NULL,"vaccine coverage 70%, vaccine efficacy 30%"),
             graph(run30,NULL,"vaccine coverage 30%, vaccine efficacy 70%"),graph(run46,NULL,"vaccine converage 46%, vaccine efficacy 46%"),ncol=2)

test_vacc<-data.frame(vacc=c(0.70,0.30,0.46),incidenceR=c(IR70,IR30,IR46),incidenceS=c(IS70,IS30,IS46))
test_vacc <- pivot_longer(test_vacc, cols = c(incidenceR,incidenceS), names_to = "Strain", values_to = "Value")


