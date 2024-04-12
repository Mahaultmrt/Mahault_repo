library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
# si les packages ne chargent pas, c'est qu'ils ne sont pas installés! Donc il faut faire install.packages("nomdupackage")

Seasonal_flu_model <- function(t, pop, param) {
  
  with((c(param, pop)), {
    
    Nc=Sc+Ic+Hc+Rc
    Na=Sa+Ia+Ha+Ra
    
    dSc <- -Sc*((matrice["c","c"]*Ic/Nc)+(matrice["a","c"]*Ia/Na))
    dIc <- Sc*((matrice["c","c"]*Ic/Nc)+(matrice["a","c"]*Ia/Na))-gamma*Ic
    dHc <- gamma*PH*Ic-mu*Hc 
    dRc <- gamma*(1-PH)*Ic-mu*(1-PD)*Hc 
    
    dSa <- -Sa*((matrice["c","a"]*Ic/Nc)+(matrice["a","a"]*Ia/Na))
    dIa <- Sa*((matrice["c","a"]*Ic/Nc)+(matrice["a","a"]*Ia/Na))-gamma*Ia
    dHa <- gamma*PH*Ia-mu*Ha 
    dRa <- gamma*(1-PH)*Ia-mu*(1-PD)*Ha 
    
    res<-c(dSc,dIc,dHc,dRc,dSa,dIa,dHa,dRa)
    
    list(res)
    
  })
  
}

# on défini les paramètres

gamma=0.14
dt=0.1
Tmax=75
mu=0.308
PH=0.0002
PD=0.11
c=c(0.3,0.2)
a=c(0.2,0.07)
matrice=matrix(data=c(f,m), nrow=2,byrow=TRUE)

#on crée la matrice de contact
rownames(matrice)=c("c","a")
colnames(matrice)=c("c","a")


# on défini les conditions initiales
Sc0=100
Ic0=1
Hc0=0
Rc0=0

Sa0=100
Ia0=1
Ha0=0
Ra0=0

# on créé les vecteurs à partir des différentes valeurs
Time=seq(from=0,to=Tmax,by=dt)
Init.cond=c(Sc=Sc0,Ic=Ic0,Hc=Hc0,Rc=Rc0,Sa=Sa0,Ia=Ia0,Ha=Ha0,Ra=Ra0) 
param=list(matrice=matrice,gamma=gamma,mu=mu)

# on utilise la fonction lsoda, avec as.data.frame() autour pour un output plus pratique
result <- as.data.frame(lsoda(Init.cond, Time, Seasonal_flu_model, param))

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

