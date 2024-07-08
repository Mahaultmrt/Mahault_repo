library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(gt)
library(tidyr)
library(forcats)
#library(ppcor)


#Code modèle page8 n°2

Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+S+CR+CS
    
    
    new_teta<-teta-log(1-vec_virus(t)*ATB)
    
    
    dSa <- -Sa*((beta*ct*(CRa+IRa+CR)/N)+beta*(CSa+ISa+CS)/N)-omega*Sa+new_teta*S+(gamma+alpha*(1-sigmaR))*CRa+(gamma+alpha*(1-sigmaS))*CSa
    dCRa <- Sa*(beta*ct*(CRa+IRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-omega*CRa+new_teta*CR
    dCSa <- Sa*(beta*(CSa+ISa+CS)/N)-(gamma+alpha*(1-sigmaS))*CSa-omega*CSa+new_teta*CS
    dIRa <- rhoRa*CRa-deltaRa*IRa+rho*CR
    dISa <- rhoSa*CSa-deltaSa*ISa+rho*CS
    
    dS <- -S*((beta*ct*(CRa+IRa+CR)/N)+beta*(CSa+ISa+CS)/N)+omega*Sa-new_teta*S+gamma*CR+gamma*CS+deltaRa*IRa+deltaSa*ISa
    dCR <- S*(beta*ct*(CRa+IRa+CR)/N)-gamma*CR+omega*CRa-new_teta*CR
    dCS <- S*(beta*(CSa+ISa+CS)/N)-gamma*CS+omega*CSa-new_teta*CS
    
    
    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS)

   
    list(res,new_teta=new_teta)
    
  })
  
}


create_params<-function(beta=0.065,ct=0.96,deltaRa=0,deltaSa=0,gamma=0.05,rho=3*10^-6,rhoRa=3.0*10^-6,rhoSa=3*10^-6,teta=0.0014,omega=0.08, alpha=0.33, sigmaR=1,sigmaS=0,ATB=0.1)
{
  list(beta=beta,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB)
}

create_initial_cond<-function(Sa0=100000*0.5*0.8,CRa0=100000*0.5*0.04,CSa0=100000*0.5*0.16,IRa0=0,ISa0=0,S0=100000*0.5*0.8,CR0=100000*0.5*0.04,CS0=100000*0.5*0.16){
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

graph<- function(data,filter_values,title){
  #data_name<-as.character(substitute(data))
  color_values <- c("Sa" = "#499124", "CRa" = "#DE6F00", "CSa" = "#2072BA", "IRa" = "#BD5E00", "ISa" = "#163F9E", 
                    "S" = "#68CF33", "CR"="#FC7E00", "CS"="#2B9CFF")
  
  
  
  if(!is.null(filter_values))
  {
    
    p<-data %>%
      melt(id = "time") %>%
      filter(variable %in% filter_values) %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")+
      scale_color_manual(values = color_values)
    
  }
  
  else{
    p<-data %>%
      melt(id = "time") %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")+
      scale_color_manual(values = color_values)
    
    
  }
  
  return(p)
}
graph2<- function(data,filter_values,title){
  #data_name<-as.character(substitute(data))
  
  
  if(!is.null(filter_values))
  {
    p<-data %>%
      melt(id = "time") %>%
      filter(variable %in% filter_values) %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")
    
    
  }
  else{
    p<-data %>%
      melt(id = "time") %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")
    
    
  }
  
  return(p)
}

heatmap <- function(data, x_var, y_var, fill_var, x_text = NULL, y_text = NULL, fill_text = NULL, title = NULL, low_col = "#377eb8", high_col = "#e41a1c", values = FALSE, var_text = NULL) {
  graph <- ggplot(data, aes_string(x = x_var, y = y_var, fill = fill_var)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = low_col, high = high_col) +
    labs(title = title,
         x = x_text,
         y = y_text,
         fill = fill_text) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
    theme_minimal()
  
  if (values & !is.null(var_text)) {
    graph <- graph + geom_text(aes_string(label = var_text), color = "white", size = 3)
  }
  
  return(graph)
}

#Pas d'épidémie pas de vaccination et pas d'infection
vec_virus=vec_virus_0
param<-create_params(rho=0,rhoRa=0,rhoSa=0)
Init.cond<-create_initial_cond()
run0<-run(Init.cond,param)
run0_g<-graph(run0,c("Sa","CRa","CSa","S","CR","CS"),"S.Pneumoniae Colonization dynamics \n(without virus épidemic)")
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
tetas<-graph2(run3,c("teta","new_teta"),"Parameters teta for S.pneumonia colonization with 50% of vaccination for influenza")


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


grid.arrange(propC1,propC2,propC3,propC4,CR_CS1,CR_CS2,CR_CS3,CR_CS4,ncol=2)
grid.arrange(CR_CS1,CR_CS2,CR_CS3,CR_CS4)


all_res <- data.frame(
  time =run0$time,
  IPD_resistant_no_vaccination = run2$IRa,
  IPD_sensitive_no_vaccination = run2$ISa,
  IPD_resistant_50_vaccination = run3$IRa,
  IPD_sensitive_50_vaccination = run3$ISa,
  IPD_resistant_80_vaccination = run4$IRa,
  IPD_sensitive_80_vaccination = run4$ISa,
  IPD_total_no_vaccination= run2$IRa + run2$ISa,
  IPD_total_50_vaccination= run3$IRa + run3$ISa,
  IPD_total_80_vaccination= run4$IRa + run4$ISa
)
IR_g<-graph2(all_res,c("IPD_resistant_no_vaccination","IPD_resistant_50_vaccination","IPD_resistant_80_vaccination"),"Annual incidence of IPD (resistant strain)")
IS_IR_g<-graph2(all_res,c("IPD_resistant_no_vaccination","IPD_resistant_50_vaccination","IPD_resistant_80_vaccination",
                         "IPD_sensitive_no_vaccination","IPD_sensitive_50_vaccination","IPD_sensitive_80_vaccination"),"Annual incidence of IPD (resistant and sensitive strain)")
I_tot_g<-graph2(all_res,c("IPD_total_no_vaccination","IPD_total_50_vaccination","IPD_total_80_vaccination"), "Total annual incidence of IPD")



grid.arrange(run0_g,run2_g,run3_g,run4_g,ncol=2)
grid.arrange(IR_g,IS_IR_g,I_tot_g,ncol=2)



res <- data.frame(time = r1$time)
IR_final <- data.frame(vacc = numeric(), LastIR = numeric())
IS_final<- data.frame(vacc = numeric(), LastIS = numeric())
IS_IR_final<- data.frame(vacc = numeric(), LastIS_IR = numeric())

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
  LastIS=runt[["ISa"]][nrow(runt) - 1]
  new_row2=data.frame(vacc=results_df[i,1], LastIS)
  IS_final <- bind_rows(IS_final, new_row2)
  runt$IR_IS=runt$IRa+runt$ISa
  LastIS_IR=runt[["IR_IS"]][nrow(runt) - 1]
  new_row3=data.frame(vacc=results_df[i,1], LastIS_IR)
  IS_IR_final<-bind_rows(IS_IR_final, new_row3)
  value<-LastIR/tail(run1$IRa, n = 1)
  new_row4=data.frame(vacc=results_df[i,1], value)
  I_relative <- bind_rows(I_relative, new_row4)
  col<-paste("vaccination",results_df[i,1])
  
  res[[col]]<- runt$IRa
  
}

graph2(res,NULL,title=NULL)
graph2(res,c("vaccination 0.1","vaccination 0.95"),title=NULL)


#heatmap pour voir le nombre de personnes infectées en fonction de la couverture vaccinale et proportion d'antibiotiques
corr_vacc_ATB<- data.frame(vacc = numeric(), ATB=numeric(), LastIR_relative = numeric())
for (i in seq(1,19,by=1)){
  for(j in seq(0.1,0.5,by=0.1)){
    
    vec_virus=I_vac[[i]]
    param<-create_params(ATB=j)
    Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runt<-run(Init.cond,param)
    LastIR=runt[["IRa"]][nrow(runt) - 1]
    LastIR_relative=LastIR/tail(run1$IRa, n = 1)
    new_row=data.frame(vacc=results_df[i,1], ATB=j, LastIR_relative)
    corr_vacc_ATB <- bind_rows(corr_vacc_ATB, new_row)
    
    
  }
  
}
h1<-heatmap(corr_vacc_ATB,"vacc","ATB","LastIR_relative","vaccine coverage","Antibiotics",
        "Annual relative incidence \nof IPD (resistant strain) per 100 000","Relative incidence of IPD \ndepending on the vaccine coverage and the proportion of ATB",values=FALSE)


I_final<-merge(IR_final,IS_final,by="vacc")
I_final <- pivot_longer(I_final, cols = c(LastIR,LastIS), names_to = "Strain", values_to = "Value")
I_final$Strain <- factor(I_final$Strain, levels = c("LastIS","LastIR"))
ggplot(I_final, aes(fill=Strain, y=Value, x=vacc)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(name=" ",labels = c("LastIS" = "Cumulative incidence of IPDs (senstive strain)", 
                               "LastIR" = "Cumulative incidence of IPDs (resistant strain)"),
                    values = c("LastIS" = "#163F9E", 
                               "LastIR" = "#BD5E00")) +
  labs(title = "Cumulative incidence of IPDs depending on the vaccine coverage", x = "Vaccine coverage", y = "Cumulative incidence of IPDs (per 100,000)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
  theme_minimal()


corr_vacc_ATB_ISIR<- data.frame(vacc = numeric(), ATB=numeric(), LastISIR = numeric(), LastpropIR=numeric())
for (i in seq(1,19,by=1)){
  for(j in seq(0.1,0.5,by=0.1)){
    
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


h2<-heatmap(corr_vacc_ATB_ISIR,"vacc","ATB","LastISIR","Vaccine coverage","Antibiotics",
        "Cumulative incidence of IPDs (per 100,000)", "Cumulative incidence of IPDs depending on the vaccine coverage \nand the proportion of ATB",values=TRUE,var_text="LastpropIR")


grid.arrange(h1,h2,ncol=2)


df_ISIR_barplot<-data.frame(vacc = numeric(), Incidence=numeric(), propI = numeric())
new_row=data.frame(vacc=0, Incidence=run2[["IRa"]][nrow(run2) - 1]+run2[["ISa"]][nrow(run2) - 1], 
                   propI=run2[["IRa"]][nrow(run2) - 1]/(run2[["IRa"]][nrow(run2) - 1]+run2[["ISa"]][nrow(run2) - 1]))
df_ISIR_barplot<-bind_rows(df_ISIR_barplot,new_row)
new_row2=data.frame(vacc=0.5, Incidence=run3[["IRa"]][nrow(run3) - 1]+run3[["ISa"]][nrow(run3) - 1], 
                    propI=run3[["IRa"]][nrow(run3) - 1]/(run3[["IRa"]][nrow(run3) - 1]+run3[["ISa"]][nrow(run3) - 1]))
df_ISIR_barplot<-bind_rows(df_ISIR_barplot,new_row2)
new_row3=data.frame(vacc=0.8, Incidence=run4[["IRa"]][nrow(run4) - 1]+run4[["ISa"]][nrow(run4) - 1], 
                    propI=run4[["IRa"]][nrow(run4) - 1]/(run2[["IRa"]][nrow(run4) - 1]+run4[["ISa"]][nrow(run4) - 1]))
df_ISIR_barplot<-bind_rows(df_ISIR_barplot,new_row3)


ggplot(df_ISIR_barplot,aes(x=factor(vacc),y=Incidence))+
  geom_bar(stat="identity",fill="#00BFC4",width=0.3)+
  geom_text(aes(label=round(propI*100,4)),vjust=-0.5,color="black")+
  labs(title="Annual IPD incidence et percentage of IPD (resistant strain) depending on vaccination",
       x="vaccination",
       y="Annual Incidence")+
  theme_minimal()

corr_vacc_beta<- data.frame(vacc = numeric(), betas=numeric(), LastIR_relative=numeric())
for (i in seq(1,19,by=1)){
  for(j in seq(0.06,0.07,by=0.001)){
    
    vec_virus=I_vac[[i]]
    param<-create_params(beta=j)
    Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runt<-run(Init.cond,param)
    LastIR=runt[["IRa"]][nrow(runt) - 1]
    LastIR_relative=LastIR/tail(run1$IRa, n = 1)
    new_row=data.frame(vacc=results_df[i,1], betas=j, LastIR_relative)
    corr_vacc_beta <- bind_rows(corr_vacc_beta, new_row)
    
    
  }
  
}

h3<-heatmap(corr_vacc_beta,"vacc","betas","LastIR_relative","vaccine coverage","Beta",
            "Total annual incidence \nof IPD per 100,000", "Total incidence of IPD depending on the vaccine coverage \nand the rate of transmission",values=FALSE)

corr_vacc_alpha<- data.frame(vacc = numeric(), alphas=numeric(), LastIR_relative=numeric())
for (i in seq(1,19,by=1)){
  for(j in seq(0.3,0.45,by=0.01)){
    
    vec_virus=I_vac[[i]]
    param<-create_params(alpha=j)
    Init.cond<-create_initial_cond(Sa0=tail(run0$Sa, n = 1),CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),S0=tail(run0$S, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runt<-run(Init.cond,param)
    LastIR=runt[["IRa"]][nrow(runt) - 1]
    LastIR_relative=LastIR/tail(run1$IRa, n = 1)
    new_row=data.frame(vacc=results_df[i,1], alphas=j, LastIR_relative)
    corr_vacc_alpha <- bind_rows(corr_vacc_alpha, new_row)
    
    
  }
  
}

h4<-heatmap(corr_vacc_alpha,"vacc","alphas","LastIR_relative","vaccine coverage","Alpha",
            "Total annual incidence \nof IPD per 100,000", "Total incidence of IPD depending on the vaccine coverage \nand antibiotics clearance",values=FALSE)


tail(run0$Sa, n = 1)*100+tail(run0$S, n = 1)*100
tail(run0$CRa, n = 1)*100+tail(run0$CR, n = 1)*100
tail(run0$CSa, n = 1)*100+tail(run0$CS, n = 1)*100
tail(run0$ISa, n = 1)*100
tail(run0$IRa, n = 1)*100

tail(run0$Sa, n = 1)*100+tail(run0$CRa, n = 1)*100+tail(run0$CSa, n = 1)*100+tail(run0$ISa, n = 1)*100+tail(run0$IRa, n = 1)*100

run2$vaccination<-"vacc 0%"
run3$vaccination<-"vacc 50%"
run4$vaccination<-"vacc 80%"

all_run<-bind_rows(run2, run3, run4)
all_run <- melt(all_run, id.vars = c("time", "vaccination"))

color_values <- c("Sa" = "#499124", "CRa" = "#DE6F00", "CSa" = "#2072BA", "IRa" = "#BD5E00", "ISa" = "#163F9E", 
                  "S" = "#68CF33", "CR" = "#FC7E00", "CS" = "#2B9CFF")

all_graph <- ggplot(all_run, aes(x = time, y = value, colour = variable)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept=tail(run0$Sa, n = 1), linetype="dashed",color="#499124",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CRa, n = 1), linetype="dashed",color="#DE6F00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CSa, n = 1), linetype="dashed",color="#2072BA",alpha=0.5)+
  geom_hline(yintercept=tail(run0$IRa, n = 1), linetype="dashed",color="#BD5E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$ISa, n = 1), linetype="dashed",color="#163F9E",alpha=0.5)+
  geom_hline(yintercept=tail(run0$S, n = 1), linetype="dashed",color="#68CF33",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CR, n = 1), linetype="dashed",color="#FC7E00",alpha=0.5)+
  geom_hline(yintercept=tail(run0$CS, n = 1), linetype="dashed",color="#2B9CFF",alpha=0.5)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
  labs(title = "S.Pneumoniae Colonization dynamics", x = "Time", y = "Proportion of Individuals", colour = "Population:") +
  facet_grid(. ~ vaccination) + 
  scale_color_manual(values = color_values)+
  theme_bw() 




