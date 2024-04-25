library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

Res_model <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+IRa+ISa+S+CR+CS+IR+IS
    
    
    dSa <- -Sa*((beta*ct*(CRa+IRa)/N)+(beta*ct*(CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+delta*IRa+delta*ISa-omega*Sa+teta*S+gamma*CRa+(gamma+alpha)*CSa
    dCRa <- Sa*((beta*ct*(CRa+IRa)/N)+beta*ct*(CR+IR)/N)-gamma*CRa-rhoRa*CRa-omega*CRa+teta*CR
    dCSa <- Sa*(beta*(CSa+ISa+CS+IS)/N)-(gamma+alpha)*CSa-rhoSa*CSa-omega*CSa+teta*CS
    dIRa <- rhoRa*CRa-delta*IRa-omega*IRa+teta*IR
    dISa <- rhoSa*CSa-delta*ISa-omega*ISa+teta*IS
    
    dS <- -S*((beta*ct*(CRa+IRa)/N)+(beta*ct*(CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+delta*IR+delta*IS+omega*Sa-teta*S+gamma*CR+gamma*CS
    dCR <- S*((beta*ct*(CRa+IRa)/N)+beta*ct*(CR+IR)/N)-gamma*CR-rho*CR+omega*CRa-teta*CR
    dCS <- S*(beta*(CSa+ISa+CS+IS)/N)-gamma*CS-rho*CS+omega*CSa-teta*CS
    dIR <- rho*CR-delta*IR+omega*IRa-teta*IR
    dIS <- rho*CS-delta*IS+omega*ISa-teta*IS
    
    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS,dIR,dIS)
    
    list(res)
    
  })
  
}

create_params<-function(beta=1,ct=0.8,delta=0.14,gamma=0.03,rho=0.1,rhoRa=0.08,rhoSa=0,teta=0.022,omega=0.07,alpha=0.5)
{
  list(beta=beta,ct=ct,delta=delta,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha)
}

create_initial_cond<-function(Sa0=100,CRa0=1,CSa0=1,IRa0=0,ISa0=0,S0=120,CR0=1,CS0=1,IR0=0,IS0=0){
  c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0,IR=IR0,IS=IS0)
}

run<-function(Init.cond,param,Tmax=500,dt=1){
  Time=seq(from=0,to=Tmax,by=dt)
  result = as.data.frame(lsoda(Init.cond, Time, Res_model, param))
  return(result)
  
}

results_df <- data.frame(Alpha = numeric(), Last_CSa = numeric())

param<-create_params()
Init.cond<-create_initial_cond(CRa0=0,CSa0=0,CR0=0,CS0=0)
run_test<-run(Init.cond,param)


Init.cond<-create_initial_cond()
param<-create_params(rhoSa=0.08)
run1<-run(Init.cond,param)

param<-create_params()
run2<-run(Init.cond,param)

Init.cond<-create_initial_cond(CRa0=50,CSa0=50,CR0=50,CS0=50)
param<-create_params()
run3<-run(Init.cond,param)

graph<- function(data,filter_values){
  data_name<-as.character(substitute(data))
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
            legend.text = element_text(size = 8)) +
      labs(title=data_name,x = "Time", y = "Value", colour = "Population:")
    
    
  }
  else{
    p<-data %>%
      melt(id = "time") %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 8, face = "bold"),
            legend.text = element_text(size = 8)) +
      labs(title=data_name,x = "Time", y = "Value", colour = "Population:")
    
    
  }
  
  return(p)
}

run1_g<-graph(run1,NULL)
run2_g<-graph(run2,NULL)
run3_g<-graph(run3, NULL)
runtest_g<-graph(run_test,NULL)
grid.arrange(run1_g,run2_g,run3_g,runtest_g,ncol=2)

# on s'interesse aux individus colonisés par une souche résistante avec présence d'antibiotiques et sans
gCR<-graph(run1,c("CRa","CR"))

# on s'interesse aux individus colonisés par la souche sensible en présence et absence d'antibiotiques
gCS<-graph(run1,c("CSa","CS"))

#On va s'interesser à au nombres d'individus colonisés par la soucghe suscepitble selon le taux d'antibiotiques alpha
{
  results_df <- data.frame(Alpha = numeric(), Last_CSa = numeric())
  
  Init.cond<-create_initial_cond()
  param<-create_params(alpha=0.15)
  alpha1<-run(Init.cond,param)
  new_row=data.frame(Alpha=param$alpha, Last_CSa=tail(alpha1$CSa, n = 1))
  results_df <- bind_rows(results_df, new_row)
  
  param<-create_params(alpha=0.25)
  alpha2<-run(Init.cond,param)
  new_row=data.frame(Alpha=param$alpha, Last_CSa=tail(alpha2$CSa, n = 1))
  results_df <- bind_rows(results_df, new_row)
  
  param<-create_params(alpha=0.37)
  alpha3<-run(Init.cond,param)
  new_row=data.frame(Alpha=param$alpha, Last_CSa=tail(alpha3$CSa, n = 1))
  results_df <- bind_rows(results_df, new_row)
  
  param<-create_params(alpha=0.01)
  alpha4<-run(Init.cond,param)
  new_row=data.frame(Alpha=param$alpha, Last_CSa=tail(alpha4$CSa, n = 1))
  results_df <- bind_rows(results_df, new_row)
  
  param<-create_params(alpha=0.86)
  alpha5<-run(Init.cond,param)
  new_row=data.frame(Alpha=param$alpha, Last_CSa=tail(alpha5$CSa, n = 1))
  results_df <- bind_rows(results_df, new_row)
  
  param<-create_params(alpha=0.5)
  alpha6<-run(Init.cond,param)
  new_row=data.frame(Alpha=param$alpha, Last_CSa=tail(alpha6$CSa, n = 1))
  results_df <- bind_rows(results_df, new_row)
  
  param<-create_params(alpha=0.68)
  alpha7<-run(Init.cond,param)
  new_row=data.frame(Alpha=param$alpha, Last_CSa=tail(alpha7$CSa, n = 1))
  results_df <- bind_rows(results_df, new_row)
  
  param<-create_params(alpha=0.98)
  alpha8<-run(Init.cond,param)
  new_row=data.frame(Alpha=param$alpha, Last_CSa=tail(alpha8$CSa, n = 1))
  results_df <- bind_rows(results_df, new_row)
  
  param<-create_params(alpha=0)
  alpha9<-run(Init.cond,param)
  new_row=data.frame(Alpha=param$alpha, Last_CSa=tail(alpha9$CSa, n = 1))
  results_df <- bind_rows(results_df, new_row)
  
  ggplot(data   = results_df,              
         mapping = aes(x = Alpha,    
                       y = Last_CSa)) +   
    geom_point(color="blue")+
    geom_line(color="purple")+
    theme_bw()+
    labs(title="Number of people in CSa according to the rate of antiobiotics (alpha)")
  
  
  graph_alpha<-ggplot() + 
    geom_line(data=alpha9, aes(x=time,y=CSa), stat = "identity", colour="orange")+
    geom_line(data=alpha4, aes(x=time,y=CSa), stat = "identity", colour="blue")+
    geom_line(data=alpha6, aes(x=time,y=CSa), stat = "identity", colour="green")+
    geom_line(data=alpha8, aes(x=time,y=CSa),stat = "identity", colour="red")+
    theme_bw()
  
}

