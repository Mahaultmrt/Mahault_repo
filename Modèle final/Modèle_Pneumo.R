library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(gt)

#Code modèle page8 n°2

Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+IRa+ISa+S+CR+CS
    
    
    new_teta<-teta-log(1-vec_virus(t)*ATB)

    
    dSa <- -Sa*((beta*ct*(CRa+IRa+CR)/N)+beta*(CSa+ISa+CS)/N)-omega*Sa+new_teta*S+(gamma+alpha*(1-sigmaR))*CRa+(gamma+alpha*(1-sigmaS))*CSa
    dCRa <- Sa*(beta*ct*(CRa+IRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-rhoRa*CRa-omega*CRa+new_teta*CR
    dCSa <- Sa*(beta*(CSa+ISa+CS)/N)-(gamma+alpha*(1-sigmaS))*CSa-rhoSa*CSa-omega*CSa+new_teta*CS
    dIRa <- rhoRa*CRa-deltaRa*IRa+rho*CR
    dISa <- rhoSa*CSa-deltaSa*ISa+rho*CS
    
    dS <- -S*((beta*ct*(CRa+IRa+CR)/N)+beta*(CSa+ISa+CS)/N)+omega*Sa-new_teta*S+gamma*CR+gamma*CS+deltaRa*IRa+deltaSa*ISa
    dCR <- S*(beta*ct*(CRa+IRa+CR)/N)-gamma*CR-rho*CR+omega*CRa-new_teta*CR
    dCS <- S*(beta*(CSa+ISa+CS)/N)-gamma*CS-rho*CS+omega*CSa-new_teta*CS
    
    
    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS)
    
    
    exp<-(Sa+CRa+CSa+IRa+ISa)*100/N
    non_exp<-(S+CR+CS)*100/N
    CR_tot<-(CRa+CR)/N
    CS_tot<-(CSa+CS)/N
    
    prop<-c(propSa=Sa/N,propCRa=CRa/N,propCSa=CSa/N,
            propIRa=IRa/N,propISa=ISa/N,propS=S/N,propCR=CR/N,propCS=CS/N)
    list(res,new_teta=new_teta,CR_tot=CR_tot,CS_tot=CS_tot,teta=teta,prop)
    
  })
  
}


create_params<-function(beta=0.056,ct=0.9652,deltaRa=0,deltaSa=0,gamma=0.05,rho=3*10^-6,rhoRa=3*10^-6,rhoSa=3*10^-6,teta=0.0014,omega=0.08, alpha=0.33, sigmaR=1,sigmaS=0,ATB=0.1)
{
  list(beta=beta,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB)
}

create_initial_cond<-function(Sa0=100000*0.5*0.8,CRa0=100000*0.5*0.04,CSa0=100000*0.5*0.16,IRa0=0,ISa0=0,S0=100000*0.5*0.8,CR0=100000*0.5*0.04,CS0=100000*0.5*0.16){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0)
}

run<-function(Init.cond,param,Tmax=400,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model, param,vec_virus=vec_virus))
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
            legend.text = element_text(size = 6),
            plot.title = element_text(size = 8, face = "bold",hjust = 0.5)) +
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
            legend.text = element_text(size = 6),
            plot.title = element_text(size = 8, face = "bold",hjust = 0.5)) +
      labs(title=title,x = "Time", y = "Value", colour = "Population:")
    
    
  }
  
  return(p)
}
heatmap<-function(data,x_var,y_var,fill_var,title=NULL,low_col="#377eb8",high_col="#e41a1c"){
  ggplot(data, aes_string(x = x_var, y = y_var, fill = fill_var)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = low_col, high = high_col) +
    labs(title = title,
         x = x_var,
         y = y_var,
         fill = fill_var)+
    theme_minimal()
}

#Pas d'épidémie pas de vaccination et pas d'infection
vec_virus=vec_virus_0
param<-create_params(rho=0,rhoRa=0,rhoSa=0)
Init.cond<-create_initial_cond()
run0<-run(Init.cond,param)

#Pas d'épidémie pas de vaccination et infection
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL,title="S.Pneumoniae colonization without a virus epidemic")
prop1<-graph(run1,c("propSa","propCRa","propCSa",
            "propIRa","propISa","propS","propCR","propCS"),
            "S.Pneumoniae colonization without a virus epidemic")
graph(run1,c("new_teta"),title=NULL)
propC1<-graph(run1,c("propCRa","propCSa","propCR","propCS"),"S.Pneumoniae colonized people without a virus epidemic")
CR_CS1<-graph(run1,c("CR_tot","CS_tot"),"S.Pneumoniae colonized people without a virus epidemic")


# Epidémie de grippe mais pas de vaccination

vec_virus=I_vac_0
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL,title="S.Pneumoniae colonization with influenza epidemic, no vaccination")
prop2<-graph(run2,c("propSa","propCRa","propCSa",
             "propIRa","propISa","propS","propCR","propCS"),
             "S.Pneumoniae colonization with influenza epidemic, no vaccination")
propC2<-graph(run2,c("propCRa","propCSa","propCR","propCS"),"S.Pneumoniae colonized people with influenza epidemic, no vaccination")
CR_CS2<-graph(run2,c("CR_tot","CS_tot"),"S.Pneumoniae colonized people with influenza epidemic, no vaccination")


# Epidémie de grippe mais vaccination 50%
vec_virus=I_vac_50
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run3<-run(Init.cond,param)
run3_g<-graph(run3,NULL,title="S.Pneumoniae colonization with influenza epidemic, vaccination 50%")
prop3<-graph(run3,c("propSa","propCRa","propCSa",
             "propIRa","propISa","propS","propCR","propCS"),
             "S.Pneumoniae colonization with influenza epidemic, vaccination 50%")
propC3<-graph(run3,c("propCRa","propCSa","propCR","propCS"),"S.Pneumoniae colonized peoplewith influenza epidemic, vaccination 50%")
CR_CS3<-graph(run3,c("CR_tot","CS_tot"),"S.Pneumoniae colonized peoplewith influenza epidemic, vaccination 50%")
tetas<-graph(run3,c("teta","new_teta"),"Parameters teta for S.pneumonia colonization with 50% of vaccination for influenza")


# Epidémie de grippe mais vaccination 80%
vec_virus=I_vac_80
param<-create_params()
Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                               IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                               CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
run4<-run(Init.cond,param)
run4_g<-graph(run4,NULL,title="S.Pneumoniae colonization with influenza epidemic, vaccination 80%")
prop4<-graph(run4,c("propSa","propCRa","propCSa",
             "propIRa","propISa","propS","propCR","propCS"),
             "S.Pneumoniae colonization with influenza epidemic, vaccination 80%")
propC4<-graph(run4,c("propCRa","propCSa","propCR","propCS"),"S.Pneumoniae colonized people with influenza epidemic, vaccination 80%")
CR_CS4<-graph(run4,c("CR_tot","CS_tot"),"S.Pneumoniae colonized people with influenza epidemic, vaccination 80%")


grid.arrange(propC1,propC2,propC3,propC4)
grid.arrange(CR_CS1,CR_CS2,CR_CS3,CR_CS4)


all_res <- data.frame(
  time = run2$time,
  IR_no_vaccination = run2$propIRa,
  IS_no_vaccination = run2$propISa,
  IR_50_vaccination = run3$propIRa,
  IS_50_vaccination = run3$propISa,
  IR_80_vaccination = run4$propIRa,
  IS_80_vaccination = run4$propISa,
  I_no_vaccination= run2$propIRa + run2$propISa,
  I_50_vaccination= run3$propIRa + run3$propISa,
  I_80_vaccination= run4$propIRa + run4$propISa
)
IR_g<-graph(all_res,c("IR_no_vaccination","IR_50_vaccination","IR_80_vaccination"),"Proportion of people infected by a resistant strain")
IS_IR_g<-graph(all_res,c("IR_no_vaccination","IR_50_vaccination","IR_80_vaccination","IS_no_vaccination","IS_50_vaccination","IS_80_vaccination"),title=NULL)
I_tot_g<-graph(all_res,c("I_no_vaccination","I_50_vaccination","I_80_vaccination"),title=NULL)

grid.arrange(prop1,prop2,prop3,IR_g,ncol=2)
grid.arrange(IR_g,IS_IR_g,I_tot_g,ncol=2)


res <- data.frame(time = r1$time)
IR_final <- data.frame(vacc = numeric(), LastIR = numeric())
I_relative<- data.frame(vacc = numeric(), value = numeric())
for (i in seq(1,19,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                 IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                                 CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
  runt<-run(Init.cond,param)
  LastIR=runt[["IRa"]][nrow(runt) - 1]
  new_row=data.frame(vacc=results_df[i,1], LastIR)
  IR_final <- bind_rows(IR_final, new_row)
  value<-LastIR/tail(run1$IRa, n = 1)
  new_row2=data.frame(vacc=results_df[i,1], value)
  I_relative <- bind_rows(I_relative, new_row2)
  col<-paste("vaccination",results_df[i,1])
 
  res[[col]]<- runt$IRa

}


graph(res,NULL,title=NULL)
graph(res,c("vaccination 0.1","vaccination 0.95"),title=NULL)

all_res <- all_res[-nrow(all_res), ]
tail(all_res$IR_no_vaccination, n = 1)-tail(all_res$IR_80_vaccination, n = 1)
all_res$s<-all_res$IR_no_vaccination - all_res$IR_80_vaccination

diff_IR <- all_res[seq(1, nrow(all_res), by = 50), ]
diff_IR[c("time","s")]

IR_final_table <- IR_final%>%
                    gt()

grid.newpage()
grid.draw(IR_final_table)

ggplot(IR_final, aes(x = vacc, y = LastIR)) +
  geom_bar(stat = "identity", fill = "#D8BFD8", color = "black") +
  labs(title = "Barplot of people infected depending on the vaccin coverage", x = "Vaccin coverage", y = "Infected people at the end of the season") +
  theme_minimal()


I_relative_table <- tableGrob(I_relative)


#heatmap pour voir le nombre de personnes infectées en fonction de la couverture vaccinale et proportion d'antibiotiques
corr_vacc_ATB<- data.frame(vacc = numeric(), ATB=numeric(), LastIR = numeric())
for (i in seq(1,19,by=1)){
  for(j in seq(0.1,0.5,by=0.1)){
  
    vec_virus=I_vac[[i]]
    param<-create_params(ATB=j)
    Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runt<-run(Init.cond,param)
    LastIR=runt[["IRa"]][nrow(runt) - 1]
    new_row=data.frame(vacc=results_df[i,1], ATB=j, LastIR)
    corr_vacc_ATB <- bind_rows(corr_vacc_ATB, new_row)
    
    
  }
  
}
heatmap(corr_vacc_ATB,"vacc","ATB","LastIR",NULL)

ggplot(all_res[seq(1, nrow(all_res), by = 25), ], aes(x = time)) +
  geom_bar(aes(y = IR_no_vaccination, fill = "0%"), stat = "identity", position = position_dodge(), width=14) +
  geom_bar(aes(y = IR_50_vaccination, fill = "50%"), stat = "identity", position = position_dodge(), width = 14) +
  geom_bar(aes(y = IR_80_vaccination, fill = "80%"), stat = "identity", position = position_dodge(), width = 14) +
  labs(title = "People infected by a resistant strain depending on vaccin coverage",
       x = "time",
       y = "infected people",
       fill = "vaccin coverage") +
  scale_fill_manual(values = c("0%" = "#77dd77", "50%" = "#ffcccb", "80%" = "#aec6cf")) + 
  theme_minimal()

#heatmap pour voir le nombre de personnes infectées en fonction de la couverture vaccinale 
# et du taux d'individus exposés aux antibiotiques
corr_vacc_teta<- data.frame(vacc = numeric(), teta=numeric(), LastIR = numeric())
for (i in seq(1,19,by=1)){
  for(j in seq(0.001,0.002,by=0.0001)){
    
    vec_virus=I_vac[[i]]
    param<-create_params(teta=j)
    Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runt<-run(Init.cond,param)
    LastIR=runt[["IRa"]][nrow(runt) - 1]
    new_row=data.frame(vacc=results_df[i,1], teta=j, LastIR)
    corr_vacc_teta <- bind_rows(corr_vacc_teta, new_row)
  }
  
}

heatmap(corr_vacc_teta,"vacc","teta","LastIR",NULL)



