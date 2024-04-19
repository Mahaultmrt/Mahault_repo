install.packages("gridExtra")

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)


Res_1 <- function(results_df,beta,ct,sigma, gamma,rho,rhoRa,rhoSa,teta,omega,alpha,Sa0,CRa0,CSa0,IRa0,ISa0,S0,CR0, CS0,IR0,IS0,Time) {
  require(deSolve) 
  
  Resistance_model_func <- function(t, pop, parameters) {
    with(as.list(c(pop, parameters)), {
      
      N=Sa+CRa+CSa+IRa+ISa+S+CR+CS+IR+IS
      
      dSa <- -Sa*((beta*ct*(CRa+IRa)/N)+(beta*ct*(CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+sigma*IRa+sigma*ISa-omega*Sa+teta*S+gamma*CRa+(gamma+alpha)*CSa
      dCRa <- Sa*((beta*ct*(CRa+IRa)/N)+beta*ct*(CR+IR)/N)-gamma*CRa-rhoRa*CRa-omega*CRa+teta*CR
      dCSa <- Sa*(beta*(CSa+ISa+CS+IS)/N)-(gamma+alpha)*CSa-rhoSa*CSa-omega*CSa+teta*CS
      dIRa <- rhoRa*CRa-sigma*IRa-omega*IRa+teta*IR
      dISa <- rhoSa*CSa-sigma*ISa-omega*ISa+teta*IS
      
      dS <- -S*((beta*ct*(CRa+IRa)/N)+(beta*ct*(CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+sigma*IR+sigma*IS+omega*Sa-teta*S+gamma*CR+gamma*CS
      dCR <- S*((beta*ct*(CRa+IRa)/N)+beta*ct*(CR+IR)/N)-gamma*CR-rho*CR+omega*CRa-teta*CR
      dCS <- S*(beta*(CSa+ISa+CS+IS)/N)-gamma*CS-rho*CS+omega*CSa-teta*CS
      dIR <- rho*CR-sigma*IR+omega*IRa-teta*IR
      dIS <- rho*CS-sigma*IS+omega*ISa-teta*IS
      
      res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS,dIR,dIS)
      
      return(list(res))
    })
  }
  
  parameters_values <- c(beta=beta,ct=ct,sigma=sigma,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha)
  
  initial_values <- c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0,IR=IR0,IS=IS0)
  
  out <- lsoda(initial_values, Time, Resistance_model_func, parameters_values)
  #as.data.frame(out)
  
  df_out <- as.data.frame(out)
  
  last_CSa <- tail(df_out$CSa, n = 1)
  
  new_row <- data.frame(Alpha = alpha, Last_CSa = last_CSa)
  results_df <- rbind(results_df, new_row)
  
  return(list(results=df_out,Result_DF=results_df))

}

results_df <- data.frame(Alpha = numeric(), Last_CSa = numeric())

{

  test<-Res_1(results_df,beta=1,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.1,
              omega=0.1,
              alpha=0.72,
              Sa0=100,
              CRa0=0,
              CSa0=0,
              IRa0=0,
              ISa0=0,
              S0=100,
              CR0=0,
              CS0=0,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
  )
  
r<-Res_1(results_df,beta=1,
         ct=0.1,
         sigma=0.14,
         gamma=0.03,
         rho=0.35,
         rhoRa=0.3,
         rhoSa=0,
         teta=0.022,
         omega=0.1,
         alpha=0.5,
         Sa0=100,
         CRa0=50,
         CSa0=50,
         IRa0=0,
         ISa0=0,
         S0=100,
         CR0=50,
         CS0=50,
         IR0=0,
         IS0=0,
         Time=seq(from=0,to=150,by=1)
)


r2<-Res_1(results_df,beta=1,
          ct=0.1,
          sigma=0.14,
          gamma=0.03,
          rho=0.35,
          rhoRa=0.3,
          rhoSa=0,
          teta=0.022,
          omega=0.1,
          alpha=0.10,
          Sa0=70,
          CRa0=1,
          CSa0=1,
          IRa0=0,
          ISa0=0,
          S0=65,
          CR0=1,
          CS0=1,
          IR0=0,
          IS0=0,
          Time=seq(from=0,to=150,by=1)
)

r3<-Res_1(results_df,beta=1,
         ct=0.1,
         sigma=0.14,
         gamma=0.03,
         rho=0.35,
         rhoRa=0.3,
         rhoSa=0,
         teta=0.022,
         omega=0.1,
         alpha=0.5,
         Sa0=70,
         CRa0=1,
         CSa0=1,
         IRa0=0,
         ISa0=0,
         S0=65,
         CR0=1,
         CS0=1,
         IR0=0,
         IS0=0,
         Time=seq(from=0,to=150,by=1)
)


}


graph<- function(data,filter_values){
  
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
      labs(x = "Time", y = "Value", colour = "Population:")
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
      labs(x = "Time", y = "Value", colour = "Population:")
    
  }
 
  return(p)
}

g1<-graph(r$results,NULL)
g2<-graph(r2$results,c("ISa","IS"))
g3<-graph(r3$results, NULL)
gtest<-graph(test$results,NULL)

grid.arrange(g1,g2,g3,gtest,ncol=2)
{
#On teste avec différents alpha
alpha1<-Res_1(results_df,beta=1,
          ct=0.1,
          sigma=0.14,
          gamma=0.03,
          rho=0.35,
          rhoRa=0.3,
          rhoSa=0,
          teta=0.022,
          omega=0.1,
          alpha=0.15,
          Sa0=70,
          CRa0=1,
          CSa0=1,
          IRa0=0,
          ISa0=0,
          S0=65,
          CR0=1,
          CS0=1,
          IR0=0,
          IS0=0,
          Time=seq(from=0,to=150,by=1)
)
alpha2<-Res_1(results_df,beta=1,
          ct=0.1,
          sigma=0.14,
          gamma=0.03,
          rho=0.35,
          rhoRa=0.3,
          rhoSa=0,
          teta=0.022,
          omega=0.1,
          alpha=0.5,
          Sa0=70,
          CRa0=1,
          CSa0=1,
          IRa0=0,
          ISa0=0,
          S0=65,
          CR0=1,
          CS0=1,
          IR0=0,
          IS0=0,
          Time=seq(from=0,to=150,by=1)
)
alpha3<-Res_1(results_df,beta=1,
          ct=0.1,
          sigma=0.14,
          gamma=0.03,
          rho=0.35,
          rhoRa=0.3,
          rhoSa=0,
          teta=0.022,
          omega=0.1,
          alpha=0.85,
          Sa0=70,
          CRa0=1,
          CSa0=1,
          IRa0=0,
          ISa0=0,
          S0=65,
          CR0=1,
          CS0=1,
          IR0=0,
          IS0=0,
          Time=seq(from=0,to=150,by=1)
)
alpha4<-Res_1(results_df,beta=1,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.022,
              omega=0.1,
              alpha=0.2,
              Sa0=70,
              CRa0=1,
              CSa0=1,
              IRa0=0,
              ISa0=0,
              S0=65,
              CR0=1,
              CS0=1,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
)
alpha5<-Res_1(results_df,beta=1,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.022,
              omega=0.1,
              alpha=0.32,
              Sa0=70,
              CRa0=1,
              CSa0=1,
              IRa0=0,
              ISa0=0,
              S0=65,
              CR0=1,
              CS0=1,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
)

alpha6<-Res_1(results_df,beta=1,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.022,
              omega=0.1,
              alpha=0.65,
              Sa0=70,
              CRa0=1,
              CSa0=1,
              IRa0=0,
              ISa0=0,
              S0=65,
              CR0=1,
              CS0=1,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
)
alpha7<-Res_1(results_df,beta=1,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.022,
              omega=0.1,
              alpha=0.98,
              Sa0=70,
              CRa0=1,
              CSa0=1,
              IRa0=0,
              ISa0=0,
              S0=65,
              CR0=1,
              CS0=1,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
)
alpha8<-Res_1(results_df,beta=1,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.022,
              omega=0.1,
              alpha=0.71,
              Sa0=70,
              CRa0=1,
              CSa0=1,
              IRa0=0,
              ISa0=0,
              S0=65,
              CR0=1,
              CS0=1,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
)
alpha9<-Res_1(results_df,beta=1,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.022,
              omega=0.1,
              alpha=0.56,
              Sa0=70,
              CRa0=1,
              CSa0=1,
              IRa0=0,
              ISa0=0,
              S0=65,
              CR0=1,
              CS0=1,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
)
}


# on s'interesse aux individus colonisés par une souche résistante avec présence d'antibiotiques et sans
gCR<-graph(r$results,c("CRa","CR"))

# on s'interesse aux individus colonisés par la souche sensible en présence et absence d'antibiotiques
gCS<-graph(r$results,c("CSa","CS"))



#On teste avec différents alpha

graph_alpha<-ggplot() + 
  geom_line(data=alpha1$results, aes(x=time,y=CSa), stat = "identity", colour="blue")+
  geom_line(data=alpha2$results, aes(x=time,y=CSa), stat = "identity", colour="green")+
  geom_line(data=alpha3$results, aes(x=time,y=CSa),stat = "identity", colour="red")+
  theme_bw()


#On teste avec différents beta
{
beta1<-Res_1(results_df,beta=1,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.022,
              omega=0.1,
              alpha=0.15,
              Sa0=70,
              CRa0=1,
              CSa0=1,
              IRa0=0,
              ISa0=0,
              S0=65,
              CR0=1,
              CS0=1,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
)

beta2<-Res_1(results_df,beta=3,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.022,
              omega=0.1,
              alpha=0.5,
              Sa0=70,
              CRa0=1,
              CSa0=1,
              IRa0=0,
              ISa0=0,
              S0=65,
              CR0=1,
              CS0=1,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
)
beta3<-Res_1(results_df,beta=6,
              ct=0.1,
              sigma=0.14,
              gamma=0.03,
              rho=0.35,
              rhoRa=0.3,
              rhoSa=0,
              teta=0.022,
              omega=0.1,
              alpha=0.85,
              Sa0=70,
              CRa0=1,
              CSa0=1,
              IRa0=0,
              ISa0=0,
              S0=65,
              CR0=1,
              CS0=1,
              IR0=0,
              IS0=0,
              Time=seq(from=0,to=150,by=1)
)
}
graph_beta<-ggplot() + 
  geom_line(data=beta1$results, aes(x=time,y=CRa+CSa+CR+CS),color="blue", stat = "identity")+
  geom_line(data=beta2$results, aes(x=time,y=CRa+CSa+CR+CS), color="red", stat = "identity")+
  geom_line(data=beta3$results, aes(x=time,y=CRa+CSa+CR+CS),color="green",stat = "identity")+
  theme_bw()


grid.arrange(graph_alpha,graph_beta,ncol=2)

data_point<-rbind(alpha1$Result_DF, alpha2$Result_DF,alpha3$Result_DF,alpha4$Result_DF,
                  alpha5$Result_DF,alpha6$Result_DF,alpha7$Result_DF,alpha8$Result_DF,alpha9$Result_DF)

#Permet de voir comment va varier le nombre d'individus colonisés susceptibles sous ATB en fonction de alpha
ggplot(data   = data_point,              
       mapping = aes(x = a,    
                     y = val)) +   
  geom_point()                   

