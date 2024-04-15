
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
# si les packages ne chargent pas, c'est qu'ils ne sont pas installés! Donc il faut faire install.packages("nomdupackage")

Resistance_model <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+IRa+ISa+S+CR+CS+IR+IS
    
    dSa <- -Sa*((beta*(1-c)*(CRa+IRa)/N)+(beta*c*(CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+sigma*IRa+sigma*ISa-omega*Sa+teta*S+gamma*CRa+(gamma+alpha)*CSa
    dCRa <- Sa*((beta*(1-c)*(CRa+IRa)/N)+beta*c*(CR+IR)/N)-gamma*CRa-rhoRa*CRa
    dCSa <- Sa*(beta*(CSa+ISa+CS+IS)/N)-(gamma+alpha)*CSa-rhoSa*CSa
    dIRa <- rhoRa*CRa-sigma*IRa
    dISa <- rhoSa*CSa-sigma*ISa
    
    dS <- -S*((beta*(1-c)*(CRa+IRa)/N)+(beta*c*(CR+IR)/N)+beta*(CSa+ISa+CS+IS)/N)+sigma*IR+sigma*IS+omega*Sa-teta*S+gamma*CR+gamma*CS
    dCR <- S*((beta*(1-c)*(CRa+IRa)/N)+beta*c*(CR+IR)/N)-gamma*CR-rho*CR
    dCS <- S*(beta*(CSa+ISa+CS+IS)/N)-gamma*CS-rho*CS
    dIR <- rho*CR-sigma*IR
    dIS <- rho*CS-sigma*IS
    
    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS,dIR,dIS)
    
    list(res)
    
  })
  
}

# on défini les paramètres
beta=1
c=0.2
sigma=0.14
gamma=0.03
rho=0.5
rhoRa=0.35
rhoSa=0
teta=0.022
omega=0.1
alpha=0.5
dt=0.1
Tmax=25


# on défini les conditions initiales
Sa0=50
CRa0=1
CSa0=1
IRa0=0
ISa0=0

S0=25
CR0=1
CS0=1
IR0=0
IS0=0

# on créé les vecteurs à partir des différentes valeurs
Time=seq(from=0,to=Tmax,by=dt)
Init.cond=c(Sa=Sa0,CRa=CRa0,CSa=CSa0,IRa=IRa0,ISa=ISa0,S=S0,CR=CR0,CS=CS0,IR=IR0,IS=IS0) 
param=c(beta=beta,c=c,sigma=sigma,gamma=gamma,rho=rho,rhoRa=rhoRa,rhoSa=rhoSa,teta=teta,omega=omega,alpha=alpha)

# on utilise la fonction lsoda, avec as.data.frame() autour pour un output plus pratique
result <- as.data.frame(lsoda(Init.cond, Time, Resistance_model, param))

# méthode ggplot + reshape + joli
result %>%
  melt(id = "time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11)) +
  labs(x = "Time", y = "Value", colour = "Population:")

# fonction pour enregistrer le plot (enregistre par défaut le dernier plot affiché)
# ggsave("nom_du_plot.png")
