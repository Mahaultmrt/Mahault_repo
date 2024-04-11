
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
SEI_h_model_matrice <- function(t, pop, param) {
  
  SF <- pop[1]
  EF <- pop[2]
  IF <- pop[3]
  
  
  SM <- pop[4]
  EM <- pop[5]
  IM <- pop[6]
  
 
  
  matrice=param
  sigma=param[3,1]
  
  
  NF=SF+EF+IF
  NM=SM+EM+IM
  
  dSF <- -SF*((matrice["F","F"]*IF/NF)+(matrice["M","F"]*IM/NM))
  dEF <- SF*((matrice["F","F"]*IF/NF)+(matrice["M","F"]*IM/NM))-sigma*EF
  dIF <- sigma*EF 
  
  dSM <- -SM*((matrice["M","M"]*IM/NM)+(matrice["F","M"]*IF/NF))
  dEM <- SM*((matrice["M","M"]*IM/NM)+(matrice["F","M"]*IF/NF))-sigma*EM
  dIM <- sigma*EM 
  
  res<-c(dSF,dEF,dIF,dSM,dEM,dIM)
  
  list(res)
  
}
# on défini les paramètres

sigma=0.0002
dt=0.1
Tmax=300
u=0.0001
f=c(0.01,0.01)
m=c(0.10,0.05)
s=c(0.0002,0)
matrice=matrix(data=c(f,m), nrow=2,byrow=TRUE)
rownames(matrice)=c("F","M")
colnames(matrice)=c("F","M")


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
#param=c(matrice=matrice,sigma=sigma)
param=rbind(matrice,s)

# on utilise la fonction lsoda, avec as.data.frame() autour pour un output plus pratique
result <- as.data.frame(lsoda(Init.cond, Time, SEI_h_model_matrice, param))

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

