library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(gt)
library(tidyr)
library(forcats)
library(epiR)
library(cowplot)

source(here::here("Impact of vaccination on ABR", "Equations.R"))
source(here::here("Impact of vaccination on ABR", "Graphs.R"))
source(here::here("Impact of vaccination on ABR", "SIR_Influenza.R"))

#Code modèle S.Aureus

create_params<-function(beta=1.65*10^-2,fitness=0.90,deltaRa=0,deltaSa=0,gamma=1.02*10^-2,rhoR=3*10^-6,rhoS=3*10^-6,rhoRa=3*10^-6,rhoSa=3*10^-6,theta=0.0014,omega=0.08, alpha=0.33, sigmaR=1,sigmaS=0, ATB=0.1)
{
  list(beta=beta,fitness=fitness,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rhoR=rhoR,rhoS=rhoS,rhoRa=rhoRa,rhoSa=rhoSa,theta=theta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB)
}

create_initial_cond<-function(Sa0=100000*0.02*0.7,CRa0=100000*0.02*0.03,CSa0=100000*0.02*0.27,IRa0=0,ISa0=0,S0=100000*0.98*0.7,CR0=100000*0.98*0.03,CS0=100000*0.98*0.27,Inc0=0){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0,Inc=Inc0)
}

run<-function(Init.cond,param,Tmax=365,dt=1, vec_virus=vec_virus){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model_aureus, param,vec_virus=vec_virus))
  tot <- sum(result[1, !(colnames(result) %in% c("time", "IRa", "ISa","N","new_theta"))])
  proportion <- result
  proportion[ , -1] <- proportion[ , -1] / tot
  
  return(proportion)
}



#Pas d'épidémie pas de vaccination et pas d'infection
param<-create_params()
Init.cond<-create_initial_cond()
run0a<-run(Init.cond,param, vec_virus = vec_virus_0)
run0_g<-graph(run0a,c("Sa","CRa","CSa","S","CR","CS"),"S.Aureus Colonization dynamics \n(without virus epidemic)")

# Epidémie de grippe vaccination 0%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                               IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                               CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
run2a<-run(Init.cond,param, vec_virus = I_vac_0)
run2_g<-graph(run2a,NULL,title="S.Aureus Colonization dynamics \nwith influenza epidemic, no vaccination")


# Epidémie de grippe vaccination 50%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                               IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                               CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
run3a<-run(Init.cond,param,vec_virus=I_vac_50)
run3_g<-graph(run3a,NULL,title="S.Aureus Colonization dynamics \nwith influenza epidemic and vaccine coverage at 50%")


# Epidémie de grippe vaccination 80%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                               IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                               CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
run4a<-run(Init.cond,param, vec_virus = I_vac_80)
run4_g<-graph(run4a,NULL,title="S.Aureus Colonization dynamics \nwith influenza epidemic and vaccine coverage at 80%")


# Tableau avec les différentes incidences cumulées des infections en fonction de la couverture vaccinale
IR_final <- data.frame(vacc = numeric(), LastIR = numeric(), ratioIR=numeric())
IS_final<- data.frame(vacc = numeric(), LastIS = numeric())
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                                 IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                                 CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
  runt<-run(Init.cond,param, vec_virus = vec_virus)
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
I_finala<-merge(IR_final,IS_final,by="vacc")
I_finala <- pivot_longer(I_finala, cols = c(LastIR,LastIS), names_to = "Strain", values_to = "Value")
I_finala$Strain <- factor(I_finala$Strain, levels = c("LastIS","LastIR"))

# Tableau des incidences cumulées des infections (S+R) en fonction de la vaccinaion et ATB
corr_vacc_ATB_ISIR<- data.frame(vacc = numeric(), ATB=numeric(), LastISIR = numeric(), LastpropIR=numeric())
for (i in seq(1,21,by=1)){
  for(j in seq(0,0.5,by=0.1)){
    
    vec_virus=I_vac[[i]]
    param<-create_params(ATB=j)
    Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                                   IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                                   CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
    runt<-run(Init.cond,param, vec_virus = vec_virus)
    runt$ISIR=runt$ISa+runt$IRa
    LastIRa=runt[["IRa"]][nrow(runt) - 1]
    LastISIR=runt[["ISIR"]][nrow(runt) - 1]
    LastpropIR=round(LastIRa*100/LastISIR,digits=0)
    new_row=data.frame(vacc=results_df[i,1], ATB=j, LastISIR,LastpropIR)
    corr_vacc_ATB_ISIR <- bind_rows(corr_vacc_ATB_ISIR, new_row)
    
  }
  
}

# Tableau des pourcentages d'individus susceptibles, colonisés, exposés et non exposés aux antibiotiques
data0<-percentage_final(run0a)

run2bis<-run2a
run3bis<-run3a
run4bis<-run4a

run2bis$vaccination<-"vacc 0%"
run3bis$vaccination<-"vacc 50%"
run4bis$vaccination<-"vacc 80%"

all_run<-bind_rows(run2bis, run3bis, run4bis)
all_run <- melt(all_run, id.vars = c("time", "vaccination"))

all_runa=all_run%>% filter(variable!="Inc")

# Graphique avec les trois sorties du modèle en fonction des vaccination
res_graphs_A<-all_graph(all_runa,NULL)+
  geom_hline(yintercept=tail(run0a$Sa, n = 1), linetype="dashed",color="#499124",alpha=0.5)+
  geom_hline(yintercept=tail(run0a$CRa, n = 1), linetype="dashed",color="#DE6F00",alpha=0.5)+
  geom_hline(yintercept=tail(run0a$CSa, n = 1), linetype="dashed",color="#2072BA",alpha=0.5)+
  geom_hline(yintercept=tail(run0a$IRa, n = 1), linetype="dashed",color="#BD5E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0a$ISa, n = 1), linetype="dashed",color="#163F9E",alpha=0.5)+
  geom_hline(yintercept=tail(run0a$S, n = 1), linetype="dashed",color="#68CF33",alpha=0.5)+
  geom_hline(yintercept=tail(run0a$CR, n = 1), linetype="dashed",color="#FC7E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0a$CS, n = 1), linetype="dashed",color="#2B9CFF",alpha=0.5)


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
exp_finala<- data.frame(vacc = numeric(), exp = numeric())
for (i in seq(1,21,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                                 IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                                 CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
  runt<-run(Init.cond,param, vec_virus = vec_virus)
  exp=runt[["Inc"]][nrow(runt) - 1]
  new_row=data.frame(vacc=results_df[i,1], exp)
  exp_finala <- bind_rows(exp_finala, new_row)
  col<-paste("vaccination",results_df[i,1])
  
}

# Différences relatives des expostions finales aux antibiotiques
diff_expa<- data.frame(vacc = numeric(), diffexp=numeric())
for (i in seq(0,20,by=1)){
  
  diffexp=(exp_finala$exp[i+1]-exp_finala$exp[1])/exp_finala$exp[1]
  new_row=data.frame(vacc=results_df[i+1,1], diffexp)
  diff_expa <- bind_rows(diff_expa, new_row)
  
  
}

# Courbes des différences relatives d'expositions aux antibiotiques en fonction de la vaccination
graph_diff_expA<-graph_exp(diff_expa,diff)


# Calcules des ratio infection/exposition aux antibiotiques
ratio_A=data.frame(vacc=seq(0,1,by=0.05),ratio_exp_R=NA,ratioR_RS_exp=NA)
ratio_A$ratio_exp_R=diff$diffIR/diff_expa$diffexp
ratio_A$ratioR_RS_exp=(diff$diffratio)/diff_expa$diffexp

# Analyse de sensibilité univariée
# Analyse de sensibilité probabiliste
univar=data.frame(param=c("beta", "beta",
                          "fitness", "fitness",
                          "gamma", "gamma",
                          "rho", "rho",
                          "alpha", "alpha",
                          "theta", "theta",
                          "omega", "omega",
                          "ATB", "ATB"),
                  val = c(0.8*0.065,1.2*0.065,
                          0.8*0.96,1.2*0.96,
                          0.8*0.05,1.2*0.05,
                          0.8*3*10^-6,1.2*3*10^-6,
                          0.8*0.33,1.2*0.33,
                          0.8*0.0014,1.2*0.0014,
                          0.8*0.08,1.2*0.08,
                          0.8*0.1,1.2*0.1),
                  incidenceR = 0, incidenceS = 0, incidenceT = 0)

for (i in 1:nrow(univar)){

  param<-create_params()
  
  param[grepl(univar$param[i], names(param))] = univar$val[i]
  
  Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                                 IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                                 CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
  runi0<-run(Init.cond,param, vec_virus = I_vac_0)
  runi<-run(Init.cond,param, vec_virus = I_vac_50)
  univar$incidenceR[i]=(runi[["IRa"]][nrow(runi) - 1] - runi0[["IRa"]][nrow(runi0) - 1]) / runi0[["IRa"]][nrow(runi0) - 1]
  univar$incidenceS[i]=(runi[["ISa"]][nrow(runi) - 1] - runi0[["ISa"]][nrow(runi0) - 1]) / runi0[["ISa"]][nrow(runi0) - 1]
  univar$incidenceT[i]=( (runi[["IRa"]][nrow(runi) - 1]+runi[["ISa"]][nrow(runi) - 1]) - (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1]) ) / (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1])
  
}

ggplot(univar) +
  geom_col(aes(x=incidenceR, y=param))

# Analyse de sensibilité probabiliste
psa=data.frame(beta=runif(500,0.8*0.065,1.2*0.065),
               fitness=runif(500,0.8*0.96,1.2*0.96),
               gamma=runif(500,0.8*0.05,1.2*0.05),
               rho=runif(500,0.8*3*10^-6,1.2*3*10^-6),
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
  alpha<-psa$alpha[i]
  theta<-psa$theta[i]
  omega<-psa$omega[i]
  ATB<-psa$ATB[i]
  
  param<-create_params(beta=beta,fitness=fitness,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB,rhoR=rho,rhoS=rho,rhoRa=rho,rhoSa=rho)
  Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                                 IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                                 CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
  runi0<-run(Init.cond,param, vec_virus = I_vac_0)
  runi<-run(Init.cond,param, vec_virus = I_vac_50)
  psa$incidenceR[i]=(runi[["IRa"]][nrow(runi) - 1] - runi0[["IRa"]][nrow(runi0) - 1]) / runi0[["IRa"]][nrow(runi0) - 1]
  psa$incidenceS[i]=(runi[["ISa"]][nrow(runi) - 1] - runi0[["ISa"]][nrow(runi0) - 1]) / runi0[["ISa"]][nrow(runi0) - 1]
  psa$incidenceT[i]=( (runi[["IRa"]][nrow(runi) - 1]+runi[["ISa"]][nrow(runi) - 1]) - (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1]) ) / (runi0[["IRa"]][nrow(runi0) - 1]+runi0[["ISa"]][nrow(runi0) - 1])
  
}

# I want to express the cor between param and impact of vaccine
# so + cor = impact increases if param increases
# but problem that vacc either increase or decrease inf
# if decrease inf, then greater impact is a smaller relative value...
# so flip around with 1- to have an impact that increases
# so it's the cor between relative increase or decrease and param

ggplot(psa) +
  geom_density(aes(incidenceR))
ggplot(psa) +
  geom_density(aes(incidenceS))
ggplot(psa) +
  geom_density(aes(incidenceT))


# Analyse de correlation partielle
psa_Ra<-psa[,-c(10:11)]
psa_Ra = psa_Ra %>%
  dplyr::select(beta,fitness,gamma,rho,alpha,theta,omega,ATB,incidenceR) %>%
  mutate(incidenceR=abs(incidenceR)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorR_A<-graph_pcor(psa_Ra)

psa_Sa<-psa[,-c(9,11)]
psa_Sa = psa_Sa %>%
  dplyr::select(beta,fitness,gamma,rho,alpha,theta,omega,ATB,incidenceS) %>%
  mutate(incidenceS=abs(incidenceS)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorS_A<-graph_pcor(psa_Sa)


psa_SRa<-psa[,-c(9:10)]
psa_SRa = psa_SRa %>%
  dplyr::select(beta,fitness,gamma,rho,alpha,theta,omega,ATB,incidenceT) %>%
  mutate(incidenceT=abs(incidenceT)) %>%
  epi.prcc() %>%
  rename(param = var)
pcorSR_A<-graph_pcor(psa_SRa)


# Tableau pour créer les courbes de différences relatives avec plusieurs simulations

# Souche résistante
IR_vacca<-data.frame(matrix(ncol=20),nrow=0)
colnames(IR_vacca)<-seq(0, 1, by = 0.05)
IR_vacc_bis<-IR_vacca

IS_vacca<-data.frame(matrix(ncol=20),nrow=0)
colnames(IS_vacca)<-seq(0, 1, by = 0.05)
IS_vacc_bis<-IS_vacca

ISR_vacca<-data.frame(matrix(ncol=20),nrow=0)
colnames(ISR_vacca)<-seq(0, 1, by = 0.05)
ISR_vacc_bis<-ISR_vacca

abx_vacca<-data.frame(matrix(ncol=20),nrow=0)
colnames(abx_vacca)<-seq(0, 1, by = 0.05)

for (i in seq(1,21,by=1)){
  for (j in seq(1,500,by=1)){
    beta<-psa$beta[j]
    fitness<-psa$fitness[j]
    gamma<-psa$gamma[j]
    rho<-psa$rho[j]
    alpha<-psa$alpha[j]
    theta<-psa$theta[j]
    omega<-psa$omega[j]
    ATB<-psa$ATB[j]
    
    vec_virus=I_vac[[i]]
    param<-create_params(beta=beta,fitness=fitness,gamma=gamma,alpha=alpha,theta=theta,omega=omega,ATB=ATB,rhoR=rho,rhoS=rho,rhoRa=rho,rhoSa=rho)
    Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                                   IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                                   CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
    runi<-run(Init.cond,param,vec_virus = vec_virus)
    IR_vacca[j,i]<-runi[["IRa"]][nrow(runi) - 1]
    IS_vacca[j,i]<-runi[["ISa"]][nrow(runi) - 1]
    ISR_vacca[j,i]<-runi[["IRa"]][nrow(runi) - 1]+runi[["ISa"]][nrow(runi) - 1]
    abx_vacca[j,i]<-runi[["Inc"]][nrow(runi) - 1]
  }
  
}


for (j in seq(1,500,by=1)){
  for (i in seq(0,20,by=1)){
    IR_vacc_bis[j,i+1]<-(IR_vacca[j,i+1]-IR_vacca[j,1])/IR_vacca[j,1]
    IS_vacc_bis[j,i+1]<-(IS_vacca[j,i+1]-IS_vacca[j,1])/IS_vacca[j,1]
    ISR_vacc_bis[j,i+1]<-(ISR_vacca[j,i+1]-ISR_vacca[j,1])/ISR_vacca[j,1]
  }
}

IR_vacc_bis <- gather(IR_vacc_bis, key = "vacc", value = "incidence")
IR_vacc_bis$vacc <- as.numeric(as.character(IR_vacc_bis$vacc))

IS_vacc_bis <- gather(IS_vacc_bis, key = "vacc", value = "incidence")
IS_vacc_bis$vacc <- as.numeric(as.character(IS_vacc_bis$vacc))

ISR_vacc_bis <- gather(ISR_vacc_bis, key = "vacc", value = "incidence")
ISR_vacc_bis$vacc <- as.numeric(as.character(ISR_vacc_bis$vacc))

diff_sima<-data.frame(vacc=seq(0,1,by=0.05))

# Calculs des valeurs moyennes
mean_values_IR <- IR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_IR = mean(incidence, na.rm = TRUE))

diff_sima$mean_incidence_IR<-mean_values_IR$mean_incidence_IR

mean_values_IS <- IS_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_IS = mean(incidence, na.rm = TRUE))
diff_sima$mean_incidence_IS<-mean_values_IS$mean_incidence_IS

mean_values_ISR <- ISR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(mean_incidence_ISR = mean(incidence, na.rm = TRUE))
diff_sima$mean_incidence_ISR<-mean_values_ISR$mean_incidence_ISR


sd_values_IR <- IR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_IR = sd(incidence, na.rm = TRUE))
diff_sima$sd_incidence_IR<-sd_values_IR$sd_incidence_IR

sd_values_IS <- IS_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_IS = sd(incidence, na.rm = TRUE))
diff_sima$sd_incidence_IS<-sd_values_IS$sd_incidence_IS

sd_values_ISR <- ISR_vacc_bis %>%
  group_by(vacc) %>%
  summarise(sd_incidence_ISR = sd(incidence, na.rm = TRUE))
diff_sima$sd_incidence_ISR<-sd_values_ISR$sd_incidence_ISR

# Calculs des intervals de confiance
diff_sima$IR_ic_l=diff_sima$mean_incidence_IR-diff_sima$sd_incidence_IR
diff_sima$IR_ic_u=diff_sima$mean_incidence_IR+diff_sima$sd_incidence_IR
diff_sima$IS_ic_l=diff_sima$mean_incidence_IS-diff_sima$sd_incidence_IS
diff_sima$IS_ic_u=diff_sima$mean_incidence_IS+diff_sima$sd_incidence_IS
diff_sima$ISR_ic_l=diff_sima$mean_incidence_ISR-diff_sima$sd_incidence_ISR
diff_sima$ISR_ic_u=diff_sima$mean_incidence_ISR+diff_sima$sd_incidence_ISR

# Courbes des différences relatives avec les valeurs moyennes et intervalles de confiance
diff_sim_A<-diff_graph_sim(diff_sima)



## COMPARING EQUAL COVERAGE x EFFICACY #####

# Couverture vaccinale 70% efficacité 30%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                               IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                               CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
run70<-run(Init.cond,param,vec_virus = I_vac_70)
IR70<-run70[["IRa"]][nrow(run70) - 1]
IS70<-run70[["ISa"]][nrow(run70) - 1]
ISR70<-run70[["ISa"]][nrow(run70) - 1]+run70[["IRa"]][nrow(run70) - 1]

# Couverture vaccinale 30% efficacité 70%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                               IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                               CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
run30<-run(Init.cond,param,vec_virus = I_vac_30)
IR30<-run30[["IRa"]][nrow(run30) - 1]
IS30<-run30[["ISa"]][nrow(run30) - 1]
ISR30<-run30[["ISa"]][nrow(run30) - 1]+run30[["IRa"]][nrow(run30) - 1]

# Couverture vaccinale 46% efficacité 46%
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0a$Sa, n = 1),CRa0=tail(run0a$CRa, n = 1),CSa0=tail(run0a$CSa, n = 1),
                               IRa0=tail(run0a$IRa, n = 1),ISa0=tail(run0a$ISa, n = 1),S0=tail(run0a$S, n = 1),
                               CR0=tail(run0a$CR, n = 1),CS0=tail(run0a$CS, n = 1))
run46<-run(Init.cond,param,vec_virus = I_vac_46)
IR46<-run46[["IRa"]][nrow(run46) - 1]
IS46<-run46[["ISa"]][nrow(run46) - 1]
ISR46<-run46[["ISa"]][nrow(run46) - 1]+run46[["IRa"]][nrow(run46) - 1]

test_vacc<-data.frame(vacc=c(0.70,0.30,0.46),incidenceR=c(IR70,IR30,IR46),incidenceS=c(IS70,IS30,IS46))
test_vacc <- pivot_longer(test_vacc, cols = c(incidenceR,incidenceS), names_to = "Strain", values_to = "Value")


