library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(RColorBrewer)

#Code modèle page 12

Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=CSa+CRa+CS+CR+IRa+ISa
    
    
    new_teta<-teta-log(1-vec_virus(t)*ATB)
    
    dCSa<- - CSa*(beta*ct*(CRa+CR)/N)+(gamma+alpha*(1-sigmaR))*CRa-CSa*omega+new_teta*CS-rhoSa*CSa
    dCRa<- CSa*(beta*ct*(CRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-CRa*omega+new_teta*CR-rhoRa*CRa
    dCS<- -CS*(beta*ct*(CRa+CR)/N)+gamma*CR+CSa*omega-new_teta*CS-rho*CS
    dCR<- CS*(beta*ct*(CRa+CR)/N)-gamma*CR+CRa*omega-new_teta*CR-rho*CR
    dIRa<- rhoRa*CRa+rho*CR
    dISa<-rhoSa*CSa+rho*CS
   
    res<-c(dCSa,dCRa,dCS,dCR,dIRa,dISa)
    
    exp<-(CRa+CSa+IRa+ISa)*100/N
    non_exp<-(CR+CS)*100/N
    
    CR_tot<-(CRa+CR)/N
    CS_tot<-(CSa+CS)/N
    
    prop<-c(propCSa=CSa/N,propCRa=CRa/N,propCS=CS/N,propCR=CR/N,
            propIRa=IRa/N,propISa=ISa/N)
    list(res,new_teta=new_teta,teta=teta,CR_tot=CR_tot,CS_tot=CS_tot,prop,N=N)
    
  })
  
}


create_params<-function(beta=0.02,ct=0.9,deltaRa=0,deltaSa=0,gamma=0.01,rho=1.8*10^-6,rhoRa=1.8*10^-6,rhoSa=1.8*10^-6,teta=0.0014,omega=0.14, alpha=0.33, sigmaR=1,sigmaS=0, ATB=0.1)
{
  list(beta=beta,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB)
}

create_initial_cond<-function(CSa0=800,CRa0=200,CS0=800,CR0=200,IRa0=0,ISa0=0){
  c(CSa=CSa0,CRa=CRa0,CS=CS0,CR=CR0,IRa=IRa0,ISa=ISa0)
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

#Pas d'épidémie 
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond()
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL,title="E.Coli colonization without a virus epidemic")
graph(run1,c("CRa","CSa","CR","CS"),title=NULL)
graph(run1,c("IRa","ISa"),title=NULL)
prop1<-graph(run1,c("propCSa","propCRa","propCS","propCR","propIRa","propISa"),
             "E.Coli colonization without a virus epidemic")
propC1<-graph(run1,c("propCRa","propCSa","propCR","propCS"),"E.Coli colonized people without a virus epidemic")
CR_CS1<-graph(run1,c("CR_tot","CS_tot"),"E.Coli colonized people without a virus epidemic")


# Epidémie de rotavirus mais pas de vaccination
param<-create_params(rho=0,rhoRa=0,rhoSa=0)
Init.cond<-create_initial_cond()
run0<-run(Init.cond,param)

vec_virus=I_vac_0
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL,title="E.Coli colonization with rotavirus epidemic, no vaccination")
graph(run2,c("IRa","ISa"), title=NULL)
graph(run2,c("CRa","CSa","CR","CS"),title=NULL)
prop2<-graph(run2,c("propCSa","propCRa","propCS","propCR","propIRa","propISa"),
             "E.Coli colonization with rotavirus epidemic, no vaccination")
propC2<-graph(run2,c("propCRa","propCSa","propCR","propCS"),"E.Coli colonized people with rotavirus epidemic, no vaccination")
CR_CS2<-graph(run2,c("CR_tot","CS_tot"),"E.Coli colonized people with rotavirus epidemic, no vaccination")



# Epidémie de rotavirus vaccination 50%
vec_virus=I_vac_50
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run3<-run(Init.cond,param)
run3_g<-graph(run3,NULL,title="E.Coli colonization with rotavirus epidemic, vaccination 50%")
graph(run3,c("IRa","ISa"),title=NULL)
graph(run3,c("CRa","CSa","CR","CS"),title=NULL)
prop3<-graph(run3,c("propCSa","propCRa","propCS","propCR","propIRa","propISa"),
             "E.Coli colonization with rotavirus epidemic, vaccination 50%")
propC3<-graph(run3,c("propCRa","propCSa","propCR","propCS"),"E.Coli colonized peoplewith rotavirus epidemic, vaccination 50%")
CR_CS3<-graph(run3,c("CR_tot","CS_tot"),"E.Coli colonized peoplewith rotavirus epidemic, vaccination 50%")
tetas<-graph(run3,c("teta","new_teta"),"Parameters teta for E.Coli colonization with 50% of vaccination for rotavirus")




# # Epidémie de rotavirus vaccination 80%
vec_virus=I_vac_80
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run4<-run(Init.cond,param)
run4_g<-graph(run4,NULL,title="E.Coli colonization with rotavirus epidemic, vaccination 80%")
graph(run4,c("IRa","ISa"),title=NULL)
graph(run4,c("CRa","CSa","CR","CS"),title=NULL)
prop4<-graph(run4,c("propCSa","propCRa","propCS","propCR","propIRa","propISa"),
             "E.Coli colonization with rotavirus epidemic, vaccination 80%")
propC4<-graph(run4,c("propCRa","propCSa","propCR","propCS"),"E.Coli colonized people with rotavirus epidemic, vaccination 80%")
CR_CS4<-graph(run4,c("CR_tot","CS_tot"),"E.Coli colonized people with rotavirus epidemic, vaccination 80%")



grid.arrange(propC1,propC2,propC3,propC4)
grid.arrange(CR_CS1,CR_CS2,CR_CS3,CR_CS4)

grid.arrange(graph(run1,c("ISa","IRa"),NULL),graph(run2,c("ISa","IRa"),NULL),
             graph(run3,c("ISa","IRa"),NULL),graph(run4,c("ISa","IRa"),NULL))


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
IS_IR_g<-graph(all_res,c("IR_no_vaccination","IR_50_vaccination","IR_80_vaccination","IS_no_vaccination","IS_50_vaccination","IS_80_vaccination"),
               "Proportion of people infected by a resistant strain and sensitive strain")
I_tot_g<-graph(all_res,c("I_no_vaccination","I_50_vaccination","I_80_vaccination"),title=NULL)

grid.arrange(graph(run2,c("IRa","ISa"),NULL),graph(run3,c("IRa","ISa"),NULL),graph(run4,c("IRa","ISa"),NULL),ncol=2)
grid.arrange(prop1,prop2,prop3,IR_g,ncol=2)
grid.arrange(IR_g,IS_IR_g,I_tot_g,ncol=2)
grid.arrange(graph(run1,"new_teta",NULL),graph(run2,"new_teta",NULL),graph(run3,"new_teta",NULL),
             graph(run4,"new_teta",NULL),ncol=2)


resultat <- data.frame(time = r1$time)
IR_final <- data.frame(vacc = numeric(), LastIR = numeric())
I_relative<- data.frame(vacc = numeric(), value = numeric())
for (i in seq(1,19,by=1)){
  
  vec_virus=I_vac[[i]]
  param<-create_params()
  Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                                 CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
  runt<-run(Init.cond,param)
  LastIR=runt[["IRa"]][nrow(runt) - 1]
  new_row=data.frame(vacc=results_df[i,1], LastIR)
  IR_final <- bind_rows(IR_final, new_row)
  value<-LastIR/tail(run1$IRa, n = 1)
  new_row2=data.frame(vacc=results_df[i,1], value)
  I_relative <- bind_rows(I_relative, new_row2)
  col<-paste("vaccination",results_df[i,1])
  
  resultat[[col]]<- runt$IRa
  
}

graph(resultat,NULL,title=NULL)
graph(resultat,c("vaccination 0.1","vaccination 0.95"),title=NULL)


all_res <- all_res[-nrow(all_res), ]
tail(all_res$IR_no_vaccination, n = 1)-tail(all_res$IR_80_vaccination, n = 1)
all_res$s<-all_res$IR_no_vaccination - all_res$IR_80_vaccination

diff_IR <- all_res[seq(1, nrow(all_res), by = 50), ]
diff_IR[c("time","s")]






IR_final_table <- tableGrob(IR_final)

grid.newpage()
grid.draw(IR_final_table)

ggplot(IR_final, aes(x = vacc, y = LastIR)) +
  geom_bar(stat = "identity", fill = "#D8BFD8", color = "black") +
  labs(title = "Barplot of people infected depending on the vaccin coverage", x = "Vaccin coverage", y = "Infected people at the end of the season") +
  theme_minimal()

I_relative_table <- tableGrob(I_relative)

grid.newpage()
grid.draw(I_relative_table)

