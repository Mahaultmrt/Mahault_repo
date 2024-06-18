library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(gt)
library(tidyr)

#Code modèle page 12

Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=CSa+CRa+CS+CR+IRa+ISa
    
    
    new_teta<-teta-log(1-vec_virus(t)*ATB)
    
    dCSa<- - CSa*(beta*ct*phi*(CRa+CR)/N)+(gamma+alpha*(1-sigmaR))*CRa-CSa*omega+new_teta*CS-rhoSa*CSa
    dCRa<- CSa*(beta*ct*phi*(CRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-CRa*omega+new_teta*CR-rhoRa*CRa
    dCS<- -CS*(beta*ct*(CRa+CR)/N)+gamma*CR+CSa*omega-new_teta*CS-rho*CS
    dCR<- CS*(beta*ct*(CRa+CR)/N)-gamma*CR+CRa*omega-new_teta*CR-rho*CR
    dIRa<- rhoRa*CRa+rho*CR
    dISa<-rhoSa*CSa+rho*CS
    
    res<-c(dCSa,dCRa,dCS,dCR,dIRa,dISa)
  
    list(res,new_teta=new_teta)
    
  })
  
}


create_params<-function(beta=0.014,ct=0.96,deltaRa=0,deltaSa=0,gamma=0.01,rho=1.8*10^-6,rhoRa=1.8*10^-6,rhoSa=1.8*10^-6,teta=0.0014,omega=0.14, alpha=0.33, sigmaR=1,sigmaS=0, ATB=0.1,phi=1.1)
{
  list(beta=beta,ct=ct,deltaRa=deltaRa,deltaSa=deltaSa,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha,sigmaR=sigmaR,sigmaS=sigmaS,ATB=ATB,phi=phi)
}

create_initial_cond<-function(CSa0=800,CRa0=200,CS0=800,CR0=200,IRa0=0,ISa0=0){
  c(CSa=CSa0,CRa=CRa0,CS=CS0,CR=CR0,IRa=IRa0,ISa=ISa0)
}

run<-function(Init.cond,param,Tmax=400,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model, param,vec_virus=vec_virus))
  proportion=as.data.frame(t(apply(result, 1, function(row) row / sum(row[-1]))))
  proportion$CR_tot<-proportion$CRa+proportion$CR
  proportion$CS_tot<-proportion$CSa+proportion$CS
  proportion$C_tot<-proportion$CR_tot+proportion$CS_tot
  return(proportion)
  
  
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

heatmap<-function(data,x_var,y_var,fill_var,title=NULL,low_col="#377eb8",high_col="#e41a1c",values=FALSE,var_text=NULL){
  graph<-ggplot(data, aes_string(x = x_var, y = y_var, fill = fill_var)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = low_col, high = high_col) +
    labs(title = title,
         x = x_var,
         y = y_var,
         fill = fill_var)+
    theme_minimal()
  
  if (values &! is.null(var_text)){
    graph<-graph+geom_text(aes_string(label=var_text),color="white",size=3)
  }
  
  return(graph)
}


#Pas d'épidémie pas de vaccination et pas d'infection
vec_virus=vec_virus_0
param<-create_params(rho=0,rhoRa=0,rhoSa=0)
Init.cond<-create_initial_cond()
run0<-run(Init.cond,param)
run0_g<-graph(run0,NULL,"E.Coli colonized people without a virus epidemic and without infection")
CR_CS0<-graph(run0,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people without a virus epidemic and without infection")


#Pas d'épidémie 
vec_virus=vec_virus_0
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run1<-run(Init.cond,param)
run1_g<-graph(run1,NULL,title="E.Coli colonization without a virus epidemic")
propC1<-graph(run1,c("CRa","pCSa","CR","CS"),"E.Coli colonized people without a virus epidemic")
CR_CS1<-graph(run1,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people without a virus epidemic")


# Epidémie de rotavirus mais pas de vaccination
vec_virus=I_vac_0
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run2<-run(Init.cond,param)
run2_g<-graph(run2,NULL,title="E.Coli colonization with rotavirus epidemic, no vaccination")
propC2<-graph(run2,c("CRa","CSa","CR","CS"),"E.Coli colonized people with rotavirus epidemic, no vaccination")
CR_CS2<-graph(run2,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people with rotavirus epidemic, no vaccination")



# Epidémie de rotavirus vaccination 50%
vec_virus=I_vac_50
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run3<-run(Init.cond,param)
run3_g<-graph(run3,NULL,title="E.Coli colonization with rotavirus epidemic, vaccination 50%")
propC3<-graph(run3,c("CRa","CSa","CR","CS"),"E.Coli colonized peoplewith rotavirus epidemic, vaccination 50%")
CR_CS3<-graph(run3,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized peoplewith rotavirus epidemic, vaccination 50%")
tetas<-graph(run3,c("teta","new_teta"),"Parameters teta for E.Coli colonization with 50% of vaccination for rotavirus")




# # Epidémie de rotavirus vaccination 80%
vec_virus=I_vac_80
param<-create_params()
Init.cond<-create_initial_cond(CSa0=tail(run0$CSa, n = 1),CRa0=tail(run0$CRa, n = 1),CS0=tail(run0$CS, n = 1),
                               CR0=tail(run0$CR, n = 1),IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1))
run4<-run(Init.cond,param)
run4_g<-graph(run4,NULL,title="E.Coli colonization with rotavirus epidemic, vaccination 80%")
propC4<-graph(run4,c("CRa","CSa","CR","CS"),"E.Coli colonized people with rotavirus epidemic, vaccination 80%")
CR_CS4<-graph(run4,c("CR_tot","CS_tot","C_tot"),"E.Coli colonized people with rotavirus epidemic, vaccination 80%")


grid.arrange(propC1,propC2,propC3,propC4)
grid.arrange(CR_CS1,CR_CS2,CR_CS3,CR_CS4)

grid.arrange(graph(run1,c("ISa","IRa"),NULL),graph(run2,c("ISa","IRa"),NULL),
             graph(run3,c("ISa","IRa"),NULL),graph(run4,c("ISa","IRa"),NULL))


all_res <- data.frame(
  time =run0$time,
  IR_no_vaccination = run2$IRa,
  IS_no_vaccination = run2$ISa,
  IR_50_vaccination = run3$IRa,
  IS_50_vaccination = run3$ISa,
  IR_80_vaccination = run4$IRa,
  IS_80_vaccination = run4$ISa,
  I_no_vaccination= run2$IRa + run2$ISa,
  I_50_vaccination= run3$IRa + run3$ISa,
  I_80_vaccination= run4$IRa + run4$ISa
)
IR_g<-graph(all_res,c("IR_no_vaccination","IR_50_vaccination","IR_80_vaccination"),"Proportion of people infected by a resistant strain")
IS_IR_g<-graph(all_res,c("IR_no_vaccination","IR_50_vaccination","IR_80_vaccination","IS_no_vaccination","IS_50_vaccination","IS_80_vaccination"),
               "Proportion of people infected by a resistant strain and sensitive strain")
I_tot_g<-graph(all_res,c("I_no_vaccination","I_50_vaccination","I_80_vaccination"),title=NULL)

grid.arrange(run0_g,run2_g,run3_g,IR_g,ncol=2)
grid.arrange(IR_g,IS_IR_g,I_tot_g,ncol=2)


res <- data.frame(time = r1$time)
IR_final <- data.frame(vacc = numeric(), LastIR = numeric())
IS_final<- data.frame(vacc = numeric(), LastIS = numeric())
IS_IR_final<- data.frame(vacc = numeric(), LastIS_IR = numeric())

I_relative<- data.frame(vacc = numeric(), value = numeric())
for (i in seq(1,19,by=1)){
  
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

graph(res,NULL,title=NULL)
graph(res,c("vaccination 0.1","vaccination 0.95"),title=NULL)

all_res <- all_res[-nrow(all_res), ]
tail(all_res$IR_no_vaccination, n = 1)-tail(all_res$IR_80_vaccination, n = 1)
all_res$s<-all_res$IR_no_vaccination - all_res$IR_80_vaccination

diff_IR <- all_res[seq(1, nrow(all_res), by = 50), ]
diff_IR[c("time","s")]


IR_final_table <- IR_final%>%
  gt()

ggplot(IR_final, aes(x = vacc, y = LastIR)) +
  geom_bar(stat = "identity", fill = "#D8BFD8", color = "black") +
  labs(title = "Barplot of people infected depending on the vaccin coverage", x = "Vaccin coverage", y = "Infected people at the end of the season") +
  theme_minimal()




#heatmap pour voir le nombre de personnes infectées en fonction de la couverture vaccinale et proportion d'antibiotiques
corr_vacc_ATB<- data.frame(vacc = numeric(), ATB=numeric(), LastIR_relative = numeric())
for (i in seq(1,19,by=1)){
  for(j in seq(0.1,0.5,by=0.1)){
    
    vec_virus=I_vac[[i]]
    param<-create_params(ATB=j)
    Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runt<-run(Init.cond,param)
    LastIR=runt[["IRa"]][nrow(runt) - 1]
    LastIR_relative=LastIR/tail(run1$IRa, n = 1)
    new_row=data.frame(vacc=results_df[i,1], ATB=j, LastIR_relative)
    corr_vacc_ATB <- bind_rows(corr_vacc_ATB, new_row)
    
    
  }
  
}
heatmap(corr_vacc_ATB,"vacc","ATB","LastIR_relative","People infected depending on the vaccination coverage and the proportion of ATB (relative incidence)",values=FALSE)



I_final<-merge(IR_final,IS_final,by="vacc")
I_final<-merge(I_final,IS_IR_final,by="vacc")
I_final <- pivot_longer(I_final, cols = c(LastIR,LastIS,LastIS_IR), names_to = "Strain", values_to = "Value")
I_final$Strain <- factor(I_final$Strain, levels = c("LastIS_IR","LastIS","LastIR"))
ggplot(I_final, aes(fill=Strain, y=Value, x=vacc)) + 
  geom_bar(position="stack", stat="identity")+
  labs(title = "Barplot of people infected by a resistant strain depending on the vaccin coverage", x = "Vaccin coverage", y = "Infected people at the end of the season") +
  theme_minimal()


corr_vacc_ATB_ISIR<- data.frame(vacc = numeric(), ATB=numeric(), LastISIR = numeric(), LastpropIR=numeric())
for (i in seq(1,19,by=1)){
  for(j in seq(0.1,0.5,by=0.1)){
    
    vec_virus=I_vac[[i]]
    param<-create_params(ATB=j)
    Init.cond<-create_initial_cond(CRa0=tail(run0$CRa, n = 1),CSa0=tail(run0$CSa, n = 1),
                                   IRa0=tail(run0$IRa, n = 1),ISa0=tail(run0$ISa, n = 1),
                                   CR0=tail(run0$CR, n = 1),CS0=tail(run0$CS, n = 1))
    runt<-run(Init.cond,param)
    runt$ISIR=runt$ISa+runt$IRa
    LastIRa=runt[["IRa"]][nrow(runt) - 1]
    LastISIR=runt[["ISIR"]][nrow(runt) - 1]
    LastpropIR=round(LastIRa*100/LastISIR,2)
    new_row=data.frame(vacc=results_df[i,1], ATB=j, LastISIR,LastpropIR)
    corr_vacc_ATB_ISIR <- bind_rows(corr_vacc_ATB_ISIR, new_row)
    
    
    
  }
  
}


heatmap(corr_vacc_ATB_ISIR,"vacc","ATB","LastISIR","People infected depending on the vaccination coverage and the proportion of ATB (percentage of infected resistant)",values=TRUE,var_text="LastpropIR")




