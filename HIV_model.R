
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
# si les packages ne chargent pas, c'est qu'ils ne sont pas installés! Donc il faut faire install.packages("nomdupackage")

SEI_h_model <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    NF=SF+EF+IF
    NM=SM+EM+IM
    
    dSF <- -SF*((betaFF*IF/NF)+(betaMF*IM/NM))
    dEF <- SF*((betaFF*IF/NF)+(betaMF*IM/NM))-sigma*EF
    dIF <- sigma*EF 
    
    dSM <- -SM*((betaMM*IM/NM)+(betaFM*IF/NF))
    dEM <- SM*((betaMM*IM/NM)+(betaFM*IF/NF))-sigma*EM
    dIM <- sigma*EM 
    
    res<-c(dSF,dEF,dIF,dSM,dEM,dIM)
    
    list(res)
    
  })
  
}

# on défini les paramètres
betaMM=0.05
betaMF=0.10
betaFM=0.01
betaFF=0.01
sigma=0.0002
dt=0.1
Tmax=300
#NF=90000
#NM=100000
u=0.0001

# on défini les conditions initiales
IF0=0
SF0=100
EF0=1

IM0=0
SM0=100
EM0=1

# on créé les vecteurs à partir des différentes valeurs
Time=seq(from=0,to=Tmax,by=dt)
Init.cond=c(SF=SF0,EF=EF0,IF=IF0,SM=SM0,EM=EM0,IM=IM0) 
param=c(betaMM=betaMM,betaMF=betaMF,betaFM=betaFM,betaFF=betaFF,sigma=sigma,u=u)

# on utilise la fonction lsoda, avec as.data.frame() autour pour un output plus pratique
result <- as.data.frame(lsoda(Init.cond, Time, SEI_h_model, param))

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
