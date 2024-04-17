
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
# si les packages ne chargent pas, c'est qu'ils ne sont pas installés! Donc il faut faire install.packages("nomdupackage")

Resistance_model <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    
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
    
    list(res)
    
  })
  
}

# on défini les paramètres
beta=1
c=0.1
sigma=0.14
gamma=0.03
rho=0.35
rhoRa=0.3
rhoSa=0.15
teta=0.022
omega=0.1
alpha=0.72
dt=1
Tmax=50


# on défini les conditions initiales
Sa0=70
CRa0=1
CSa0=1
IRa0=0
ISa0=0

S0=65
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


# creating a matrix with 0 rows 
# and columns 
# mat = matrix(ncol = 2, nrow = 4)
# colnames(mat)=c("rho","res")

# converting the matrix to data 
# frame
# df=data.frame(mat)
# df["rho"]<-c("rho=1","rho=2","rho=3","rho=4")
# 
# df[1,2]<-result[251,"ISa"]
# df[2,2]<-result[251,"ISa"]
# df[3,2]<-result[251,"ISa"]
# df[4,2]<-result[251,"ISa"]
# ggplot(data=df, aes(x=rho, y=res)) +
#   geom_bar(stat="identity", fill="steelblue")+
#   theme_minimal()

# ggplot(result, aes(x=time)) + 
#   geom_line(aes(y = ISa), color = "red", linetype=3) + 
#   geom_line(aes(y = IS), color="blue", linetype=3) 

# voir l'evolution du nombre de cas infectés souche suceptible avec antibiotique lorsque rhoSA varie
p<-ggplot(result, aes(x=time)) + geom_line(aes(y = ISa, color="rhoSa1"), linetype=3)
p<- p+geom_line(data=result, aes(y = ISa, color="rhoSa2"),color="red", linetype=3)  + labs(color="rhoSa")


rhoSA1=result #rhoSA=0
rhoSA2= result #rhoSA=0.5
rhoSA3=result # rhoSA=0.15

# Barplot variation de rhoSA
ggplot() + 
  geom_bar(data=rhoSA2, aes(x=time,y=ISa), stat = "identity", fill="blue",position=position_dodge(width=0.8),width=0.35)+
  geom_bar(data=rhoSA3, aes(x=time,y=ISa), stat = "identity", fill="green",position=position_dodge(width=0.8),width=0.35)+
  geom_bar(data=rhoSA1, aes(x=time,y=ISa),stat = "identity", fill="red",position=position_dodge(width=0.8),width=0.35)
  
  
# barplot variation du nombre de personnes colonisées par la souche resistante en fonction du cout de resistance

c1=result #c=0
c2=result #c=0.15
c3=result #c=0.45

ggplot() + 
  geom_bar(data=c3, aes(x=time,y=CRa+CR), stat = "identity", fill="blue",position=position_dodge(width=0.8),width=0.35)+
  geom_bar(data=c2, aes(x=time,y=CRa+CR), stat = "identity", fill="green",position=position_dodge(width=0.8),width=0.35)+
  geom_bar(data=c1, aes(x=time,y=CRa+CR),stat = "identity", fill="red",position=position_dodge(width=0.8),width=0.35)

# barplot de variation du nombre de personnes colonisées (CSa) en fonction de l'efficacité de l'antibiotique

alpha1= result #alpha=0.12
alpha2=result #alpha=0.5
alpha3=result #alpha=0.72

ggplot() + 
  geom_bar(data=alpha1, aes(x=time,y=CSa), stat = "identity", fill="blue",position=position_dodge(width=0.8),width=0.35)+
  geom_bar(data=alpha2, aes(x=time,y=CSa), stat = "identity", fill="green",position=position_dodge(width=0.8),width=0.35)+
  geom_bar(data=alpha3, aes(x=time,y=CSa),stat = "identity", fill="red",position=position_dodge(width=0.8),width=0.35)