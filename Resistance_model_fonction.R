install.packages("gridExtra")

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

Res_1 <- function(beta,c,sigma, gamma,rho,rhoRa,rhoSa,teta,omega,alpha,Sa0,CRa0,CSa0,IRa0,ISa0,S0,CR0, CS0,IR0,IS0,times) {
  require(deSolve) 
  
  Resistance_model_func <- function(t, pop, parameters) {
    with(as.list(c(pop, parameters)), {
      
      N=Sa+CRa+CSa+IRa+ISa+S+CR+CS+IR+IS
      
      dSa <- -Sa*((beta*(1-c)*(CRa+IRa)/N)+(beta*c*(CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+sigma*IRa+sigma*ISa-omega*Sa+teta*S+gamma*CRa+(gamma+alpha)*CSa
      dCRa <- Sa*((beta*(1-c)*(CRa+IRa)/N)+beta*c*(CR+IR)/N)-gamma*CRa-rhoRa*CRa-omega*CRa+teta*CR
      dCSa <- Sa*(beta*(CSa+ISa+CS+IS)/N)-(gamma+alpha)*CSa-rhoSa*CSa-omega*CSa+teta*CS
      dIRa <- rhoRa*CRa-sigma*IRa-omega*IRa+teta*IR
      dISa <- rhoSa*CSa-sigma*ISa-omega*ISa+teta*IS
      
      dS <- -S*((beta*(1-c)*(CRa+IRa)/N)+(beta*c*(CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+sigma*IR+sigma*IS+omega*Sa-teta*S+gamma*CR+gamma*CS
      dCR <- S*((beta*(1-c)*(CRa+IRa)/N)+beta*c*(CR+IR)/N)-gamma*CR-rho*CR+omega*CRa-teta*CR
      dCS <- S*(beta*(CSa+ISa+CS+IS)/N)-gamma*CS-rho*CS+omega*CSa-teta*CS
      dIR <- rho*CR-sigma*IR+omega*IRa-teta*IR
      dIS <- rho*CS-sigma*IS+omega*ISa-teta*IS
      
      res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS,dIR,dIS)
      
      return(list(res))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta=beta,c=c,sigma=sigma,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha)
  
  # the initial values of variables:
  initial_values <- c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0,IR=IR0,IS=IS0)
  
  # solving
  out <- lsoda(initial_values, Time, Resistance_model_func, parameters_values)
  
  # returning the output:
  as.data.frame(out)
}

Res_1(beta=1,
  c=0.1,
  sigma=0.14,
  gamma=0.03,
  rho=0.35,
  rhoRa=0.3,
  rhoSa=0.15,
  teta=0.022,
  omega=0.1,
  alpha=0.72,
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
  times=seq(from=0,to=50,by=1)
)
r2<-Res_1(beta=1,
      c=0.1,
      sigma=0.14,
      gamma=0.03,
      rho=0.35,
      rhoRa=0.3,
      rhoSa=0.15,
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
      times=seq(from=0,to=50,by=1)
)

r<-Res_1(beta=1,
         c=0.1,
         sigma=0.14,
         gamma=0.03,
         rho=0.35,
         rhoRa=0.3,
         rhoSa=0.15,
         teta=0.022,
         omega=0.1,
         alpha=0.72,
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
         times=seq(from=0,to=50,by=1)
)
# méthode ggplot + reshape + joli

r3<-Res_1(beta=1,
         c=0.1,
         sigma=0.14,
         gamma=0.03,
         rho=0.35,
         rhoRa=0.3,
         rhoSa=0.15,
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
         times=seq(from=0,to=50,by=1)
)
graph1<-  r%>%
  melt(id = "time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x = "Time", y = "Value", colour = "Population:")

graph2<-r2 %>%
  melt(id = "time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x = "Time", y = "Value", colour = "Population:")

graph3<-r3 %>%
  melt(id = "time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x = "Time", y = "Value", colour = "Population:")


grid.arrange(graph1,graph2,graph3,ncol=2)

rho1<-Res_1(beta=1,
          c=0.1,
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
          times=seq(from=0,to=50,by=1)
)
rho2<-Res_1(beta=1,
            c=0.1,
            sigma=0.14,
            gamma=0.03,
            rho=0.5,
            rhoRa=0.3,
            rhoSa=0.5,
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
            times=seq(from=0,to=50,by=1)
)
rho3<-Res_1(beta=1,
            c=0.1,
            sigma=0.14,
            gamma=0.03,
            rho=0.5,
            rhoRa=0.3,
            rhoSa=0.15,
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
            times=seq(from=0,to=50,by=1)
)
ggplot() + 
  geom_bar(data=rho2, aes(x=time,y=ISa), stat = "identity", fill="blue",position=position_dodge(width=0.8),width=0.35)+
  geom_bar(data=rho3, aes(x=time,y=ISa), stat = "identity", fill="green",position=position_dodge(width=0.8),width=0.35)+
  geom_bar(data=rho1, aes(x=time,y=ISa),stat = "identity", fill="red",position=position_dodge(width=0.8),width=0.35)

