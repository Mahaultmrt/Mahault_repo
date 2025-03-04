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

source(here::here("Impact of vaccination on ABR", "Equations.R"))
source(here::here("Impact of vaccination on ABR", "Graphs.R"))
source(here::here("Impact of vaccination on ABR", "SIR_rotavirus.R"))

#Code modèle E.Coli

create_params<-function(beta=0.012,fitness=0.96,deltaRa=0,deltaSa=0,gamma=0.01,rho=1.8*10^-6,rhoRa=1.8*10^-6,rhoSa=1.8*10^-6,theta=0.0014,omega=0.14, alpha=0.33, sigmaR=1, ATB=0.1,phi=9.83)
{
  list(beta=beta,fitness=fitness,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,theta=theta,omega=omega,alpha=alpha,sigmaR=sigmaR,ATB=ATB,phi=phi)
}

create_initial_cond<-function(CSa0=100000*0.01*0.8,CRa0=100000*0.01*0.2,CS0=100000*0.99*0.8,CR0=100000*0.99*0.2,IRa0=0,ISa0=0,Inc0=0){
  c(CSa=CSa0,CRa=CRa0,CS=CS0,CR=CR0,IRa=IRa0,ISa=ISa0,Inc=Inc0)
}

run<-function(Init.cond,param,Tmax=365,dt=1,vec_virus=vec_virus){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model_Coli, param,vec_virus=vec_virus))
  tot <- sum(result[1, !(colnames(result) %in% c("time", "IRa", "ISa","N","new_theta"))])
  proportion <- result
  proportion[ , -1] <- proportion[ , -1] / tot
  
  return(proportion)
}




#Pas d'épidémie pas de vaccination et pas d'infection
param<-create_params()
Init.cond<-create_initial_cond()
run0e<-run(Init.cond,param,vec_virus=vec_virus_0)
run0_g<-graph(run0e,c("CSa","CRa","CS","CR"),"E.Coli Colonization dynamics \nwithout virus epidemics")


#Pas d'épidémie 
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0e$CSa, n = 1),CRa0=tail(run0e$CRa, n = 1),CS0=tail(run0e$CS, n = 1),
                               CR0=tail(run0e$CR, n = 1),IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1))
run1e<-run(Init.cond,param,vec_virus=vec_virus_0)
run1_g<-graph(run1e,NULL,title="E.Coli Colonization dynamics \nwithout virus epidemics")


# Epidémie de rotavirus mais pas de vaccination
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0e$CSa, n = 1),CRa0=tail(run0e$CRa, n = 1),CS0=tail(run0e$CS, n = 1),
                               CR0=tail(run0e$CR, n = 1),IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1))
run2e<-run(Init.cond,param,vec_virus=I_vac_0)
run2_g<-graph(run2e,NULL,title="E.Coli Colonization dynamics \nwith virus epidemic, no vaccination")



# Epidémie de rotavirus vaccination 50%
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0e$CSa, n = 1),CRa0=tail(run0e$CRa, n = 1),CS0=tail(run0e$CS, n = 1),
                               CR0=tail(run0e$CR, n = 1),IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1))
run3e<-run(Init.cond,param,vec_virus=I_vac_50)
run3_g<-graph(run3e,NULL,title="E.Coli Colonization dynamics \nwith rotavirus epidemic and vaccine coverage at 50%")


# Epidémie de rotavirus vaccination 80%
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0e$CSa, n = 1),CRa0=tail(run0e$CRa, n = 1),CS0=tail(run0e$CS, n = 1),
                               CR0=tail(run0e$CR, n = 1),IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1))
run4e<-run(Init.cond,param,vec_virus=I_vac_80)
run4_g<-graph(run4e,NULL,title="E.Coli Colonization dynamics\nwith rotavirus epidemic and vaccine coverage at 80%")



# Tableau avec les différentes incidences cumulées des infections en fonction de la couverture vaccinale
IR_final <- data.frame(vacc = numeric(), LastIR = numeric(), ratioIR=numeric())
IS_final<- data.frame(vacc = numeric(), LastIS = numeric())
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(CRa0=tail(run0e$CRa, n = 1),CSa0=tail(run0e$CSa, n = 1),
                                 IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1),
                                 CR0=tail(run0e$CR, n = 1),CS0=tail(run0e$CS, n = 1))
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

# Tableau des pourcentages d'individus susceptibles, colonisés, exposés et non exposés aux antibiotiques
data0<-percentage_final(run0e)

# Création du tableau pour afficher les sorties du modèle sur un graphique
run2bis<-run2e
run3bis<-run3e
run4bis<-run4e

run2bis$vaccination<-"vacc 0%"
run3bis$vaccination<-"vacc 50%"
run4bis$vaccination<-"vacc 80%"

all_run<-bind_rows(run2bis, run3bis, run4bis)
all_run <- melt(all_run, id.vars = c("time", "vaccination"))

all_rune=all_run%>% filter(variable!="Inc")

# Graphique avec les trois sorties du modèle en fonction des vaccination
res_graphs_E<-all_graph(all_rune,NULL)+
  geom_hline(yintercept=tail(run0e$CRa, n = 1), linetype="dashed",color="#DE6F00",alpha=0.5)+
  geom_hline(yintercept=tail(run0e$CSa, n = 1), linetype="dashed",color="#2072BA",alpha=0.5)+
  geom_hline(yintercept=tail(run0e$IRa, n = 1), linetype="dashed",color="#BD5E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0e$ISa, n = 1), linetype="dashed",color="#163F9E",alpha=0.5)+
  geom_hline(yintercept=tail(run0e$CR, n = 1), linetype="dashed",color="#FC7E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0e$CS, n = 1), linetype="dashed",color="#2B9CFF",alpha=0.5)

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
exp_finale<- data.frame(vacc = numeric(), exp = numeric())
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(CRa0=tail(run0e$CRa, n = 1),CSa0=tail(run0e$CSa, n = 1),
                                 IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1),
                                 CR0=tail(run0e$CR, n = 1),CS0=tail(run0e$CS, n = 1))
  runt<-run(Init.cond,param,vec_virus=vec_virus)
  exp=runt[["Inc"]][nrow(runt) - 1]
  new_row=data.frame(vacc=results_df[i,1], exp)
  exp_finale <- bind_rows(exp_finale, new_row)
  col<-paste("vaccination",results_df[i,1])
  
  
  
}

# Différences relatives des expostions finales aux antibiotiques
diff_expe<- data.frame(vacc = numeric(), diffexp=numeric())
for (i in seq(0,20,by=1)){
  
  diffexp=(exp_finale$exp[i+1]-exp_finale$exp[1])/exp_finale$exp[1]
  new_row=data.frame(vacc=results_df[i+1,1], diffexp)
  diff_expe <- bind_rows(diff_expe, new_row)
  
  
}

# Courbes des différences relatives d'expositions aux antibiotiques en fonction de la vaccination
graph_diff_expE<-graph_exp(diff_expe,diff)

# Calcules des ratio infection/exposition aux antibiotiques
ratio_E=data.frame(vacc=seq(0,1,by=0.05),ratio_exp_R=NA,ratioR_RS_exp=NA)
ratio_E$ratio_exp_R=diff$diffIR/diff_expe$diffexp
ratio_E$ratioR_RS_exp=(diff$diffratio)/diff_expe$diffexp

# Analyse de sensibilité probabiliste
psa=data.frame(beta=runif(500,0.8*0.065,1.2*0.065),
               fitness=runif(500,0.8*0.96,1.2*0.96),
               gamma=runif(500,0.8*0.05,1.2*0.05),
               rho=runif(500,0.8*1.8*10^-6,1.2*1.8*10^-6),
               phi=runif(500,0.8*9.83,1.2*9.83),
               alpha=runif(500,0.8*0.33,1.2*0.33),
               theta=runif(500,0.8*0.0014,1.2*0.0014),
               omega=runif(500,0.8*0.08,1.2*0.08),
               ATB=runif(500,0.8*0.1,1.2*0.1),
               incidenceR=0,
               incidenceS=0,
               incidenceT=0)

for (i in seq(1,500,by=1)){
  beta<-psa$beta[i]
  fitness<-psa$fitness[i]
  gamma<-psa$gamma[i]
  rho<-psa$rho[i]
  phi<-psa$phi[i]
  alpha<-psa$alpha[i]
  theta<-psa$theta[i]
  omega<-psa$omega[i]
  ATB<-psa$ATB[i]
  
  param<-create_params(beta=beta,fitness=fitness,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB,rho=rho,rhoRa=rho,rhoSa=rho,phi=phi)
  Init.cond<-create_initial_cond(CRa0=tail(run0e$CRa, n = 1),CSa0=tail(run0e$CSa, n = 1),
                                 IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1),
                                 CR0=tail(run0e$CR, n = 1),CS0=tail(run0e$CS, n = 1))
  runi0<-run(Init.cond,param, vec_virus = I_vac_0)
  runi<-run(Init.cond,param, vec_virus = I_vac_50)
  psa$incidenceR[i]=(runi[["IRa"]][nrow(runi) - 1] - runi0[["IRa"]][nrow(runi0) - 1]) / runi0[["IRa"]][nrow(runi0) - 1]
  psa$incidenceS[i]=(runi[["ISa"]][nrow(runi) - 1] - runi0[["ISa"]][nrow(runi0) - 1]) / runi0[["ISa"]][nrow(runi0) - 1]
  psa$incidenceT[i]=( (runi[["IRa"]][nrow(runi) - 1]+runi[["ISa"]][nrow(runi) - 1]) - (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1]) ) / (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1])
  
}


# Analyse de correlation partielle
psa_Re<-psa[,-c(11:12)]
psa_Re = psa_Re %>%
  dplyr::select(beta,fitness,gamma,rho,phi,alpha,theta,omega,ATB,incidenceR) %>%
  mutate(incidenceR=abs(incidenceR)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorR_E<-graph_pcor(psa_Re)

psa_Se<-psa[,-c(10,12)]
psa_Se = psa_Se %>%
  dplyr::select(beta,fitness,gamma,rho,phi,alpha,theta,omega,ATB,incidenceS) %>%
  mutate(incidenceS=abs(incidenceS)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorS_E<-graph_pcor(psa_Se)


psa_SRe<-psa[,-c(10,11)]
psa_SRe = psa_SRe %>%
  dplyr::select(beta,fitness,gamma,rho,phi,alpha,theta,omega,ATB,incidenceT) %>%
  mutate(incidenceT=abs(incidenceT)) %>%
  epi.prcc() %>%
  rename(param = var)
pcorSR_E<-graph_pcor(psa_SRe)



# Tableau pour créer les courbes de différences relatives avec plusieurs simulations

# Souche résistante
IR_vacce<-data.frame(matrix(ncol=20),nrow=0)
colnames(IR_vacce)<-seq(0, 1, by = 0.05)
IR_vacc_bis<-IR_vacce

IS_vacce<-data.frame(matrix(ncol=20),nrow=0)
colnames(IS_vacce)<-seq(0, 1, by = 0.05)
IS_vacc_bis<-IS_vacce

ISR_vacce<-data.frame(matrix(ncol=20),nrow=0)
colnames(ISR_vacce)<-seq(0, 1, by = 0.05)
ISR_vacc_bis<-ISR_vacce

abx_vacce<-data.frame(matrix(ncol=20),nrow=0)
colnames(abx_vacce)<-seq(0, 1, by = 0.05)

for (i in seq(1,21,by=1)){
  for (j in seq(1,500,by=1)){
    beta<-psa$beta[j]
    fitness<-psa$fitness[j]
    gamma<-psa$gamma[j]
    rho<-psa$rho[j]
    phi<-psa$phi[j]
    alpha<-psa$alpha[j]
    theta<-psa$theta[j]
    omega<-psa$omega[j]
    ATB<-psa$ATB[j]
    
    vec_virus=I_vac[[i]]
    param<-create_params(beta=beta,fitness=fitness,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB,rho=rho,rhoRa=rho,rhoSa=rho,phi=phi)
    Init.cond<-create_initial_cond(CRa0=tail(run0e$CRa, n = 1),CSa0=tail(run0e$CSa, n = 1),
                                   IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1),
                                   CR0=tail(run0e$CR, n = 1),CS0=tail(run0e$CS, n = 1))
    runi<-run(Init.cond,param,vec_virus = vec_virus)
    IR_vacce[j,i]<-runi[["IRa"]][nrow(runi) - 1]
    IS_vacce[j,i]<-runi[["ISa"]][nrow(runi) - 1]
    ISR_vacce[j,i]<-runi[["IRa"]][nrow(runi) - 1]+runi[["ISa"]][nrow(runi) - 1]
    abx_vacce[j,i]<-runi[["Inc"]][nrow(runi) - 1]
  }
  
}


for (j in seq(1,500,by=1)){
  for (i in seq(0,20,by=1)){
    IR_vacc_bis[j,i+1]<-(IR_vacce[j,i+1]-IR_vacce[j,1])/IR_vacce[j,1]
    IS_vacc_bis[j,i+1]<-(IS_vacce[j,i+1]-IS_vacce[j,1])/IS_vacce[j,1]
    ISR_vacc_bis[j,i+1]<-(ISR_vacce[j,i+1]-ISR_vacce[j,1])/ISR_vacce[j,1]
  }
}

IR_vacc_bis <- gather(IR_vacc_bis, key = "vacc", value = "incidence")
IR_vacc_bis$vacc <- as.numeric(as.character(IR_vacc_bis$vacc))

IS_vacc_bis <- gather(IS_vacc_bis, key = "vacc", value = "incidence")
IS_vacc_bis$vacc <- as.numeric(as.character(IS_vacc_bis$vacc))

ISR_vacc_bis <- gather(ISR_vacc_bis, key = "vacc", value = "incidence")
ISR_vacc_bis$vacc <- as.numeric(as.character(ISR_vacc_bis$vacc))

diff_sime<-data.frame(vacc=seq(0,1,by=0.05))

# Calculs des valeurs moyennes
mean_values_IR <- IR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_IR = mean(incidence, na.rm = TRUE))

diff_sime$mean_incidence_IR<-mean_values_IR$mean_incidence_IR

mean_values_IS <- IS_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_IS = mean(incidence, na.rm = TRUE))
diff_sime$mean_incidence_IS<-mean_values_IS$mean_incidence_IS

mean_values_ISR <- ISR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_ISR = mean(incidence, na.rm = TRUE))
diff_sime$mean_incidence_ISR<-mean_values_ISR$mean_incidence_ISR

# Calculs des écarts types
sd_values_IR <- IR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_IR = sd(incidence, na.rm = TRUE))
diff_sime$sd_incidence_IR<-sd_values_IR$sd_incidence_IR

sd_values_IS <- IS_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_IS = sd(incidence, na.rm = TRUE))
diff_sime$sd_incidence_IS<-sd_values_IS$sd_incidence_IS

sd_values_ISR <- ISR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_ISR = sd(incidence, na.rm = TRUE))
diff_sime$sd_incidence_ISR<-sd_values_ISR$sd_incidence_ISR

# Calculs des intervals de confiance
diff_sime$IR_ic_l=diff_sime$mean_incidence_IR-diff_sime$sd_incidence_IR
diff_sime$IR_ic_u=diff_sime$mean_incidence_IR+diff_sime$sd_incidence_IR
diff_sime$IS_ic_l=diff_sime$mean_incidence_IS-diff_sime$sd_incidence_IS
diff_sime$IS_ic_u=diff_sime$mean_incidence_IS+diff_sime$sd_incidence_IS
diff_sime$ISR_ic_l=diff_sime$mean_incidence_ISR-diff_sime$sd_incidence_ISR
diff_sime$ISR_ic_u=diff_sime$mean_incidence_ISR+diff_sime$sd_incidence_ISR

# Courbes des différences relatives avec les valeurs moyennes et intervalles de confiance
diff_sim_E<-diff_graph_sim(diff_sime)




## COMPARING EQUAL COVERAGE x EFFICACY #####

# Couverture vaccinale 70% efficacité 30%
param<-create_params()
Init.cond<-create_initial_cond(CRa0=tail(run0e$CRa, n = 1),CSa0=tail(run0e$CSa, n = 1),
                               IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1),
                               CR0=tail(run0e$CR, n = 1),CS0=tail(run0e$CS, n = 1))
run70<-run(Init.cond,param,vec_virus=I_vac_70)
IR70<-run70[["IRa"]][nrow(run70) - 1]
IS70<-run70[["ISa"]][nrow(run70) - 1]
ISR70<-run70[["ISa"]][nrow(run70) - 1]+run70[["IRa"]][nrow(run70) - 1]

# Couverture vaccinale 30% efficacité 70%
param<-create_params()
Init.cond<-create_initial_cond(CRa0=tail(run0e$CRa, n = 1),CSa0=tail(run0e$CSa, n = 1),
                               IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1),
                               CR0=tail(run0e$CR, n = 1),CS0=tail(run0e$CS, n = 1))
run30<-run(Init.cond,param,vec_virus=I_vac_30)
IR30<-run30[["IRa"]][nrow(run30) - 1]
IS30<-run30[["ISa"]][nrow(run30) - 1]
ISR30<-run30[["ISa"]][nrow(run30) - 1]+run30[["IRa"]][nrow(run30) - 1]

# Couverture vaccinale 46% efficacité 46%
param<-create_params()
Init.cond<-create_initial_cond(CRa0=tail(run0e$CRa, n = 1),CSa0=tail(run0e$CSa, n = 1),
                               IRa0=tail(run0e$IRa, n = 1),ISa0=tail(run0e$ISa, n = 1),
                               CR0=tail(run0e$CR, n = 1),CS0=tail(run0e$CS, n = 1))
run46<-run(Init.cond,param,vec_virus=I_vac_46)
IR46<-run46[["IRa"]][nrow(run46) - 1]
IS46<-run46[["ISa"]][nrow(run46) - 1]
ISR46<-run46[["ISa"]][nrow(run46) - 1]+run46[["IRa"]][nrow(run46) - 1]
grid.arrange(graph(run70,NULL,"vaccine coverage 70%, vaccine efficacy 30%"),
             graph(run30,NULL,"vaccine coverage 30%, vaccine efficacy 70%"),graph(run46,NULL,"vaccine converage 46%, vaccine efficacy 46%"),ncol=2)

test_vacc<-data.frame(vacc=c(0.70,0.30,0.46),incidenceR=c(IR70,IR30,IR46),incidenceS=c(IS70,IS30,IS46))
test_vacc <- pivot_longer(test_vacc, cols = c(incidenceR,incidenceS), names_to = "Strain", values_to = "Value")


