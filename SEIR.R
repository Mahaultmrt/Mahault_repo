
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)

SEIR.model <- function(t, pop, param) {
  
  with(as.list(c(pop, param)), {
    
    N=S+E+I+R
    
    
    dS <- -beta*S*I/N + u0*N - u1*S
    dE <- beta*I*S/N - u1*E-sigma*E
    dI <- sigma*E-u1*I-gamma*I
    dR <-  gamma*I - u1*R
    
    res<-c(dS, dI, dR,dE)
    
    list(res)
    
  })
  
}

# on défini les paramètres
beta=3
sigma=0.5
gamma=0.5
dt=1
Tmax=300
N=100000
u0=0.01
u1=0.01

# on défini les conditions initiales
I0=1
S0=N-10
R0=0
E0=1

# on créé les vecteurs à partir des différentes valeurs
Time=seq(from=0,to=Tmax,by=dt)
Init.cond=c(S=S0,I=I0,R=R0,E=E0) 
param=c(beta=beta,gamma=gamma,sigma=sigma,u0=u0,u1=u1)

# on utilise la fonction lsoda, avec as.data.frame() autour pour un output plus pratique
result <- as.data.frame(lsoda(Init.cond, Time, SEIR.model, param))

# méthode "basique" pour plot
plot(Time,result$S,type="l",col="green",xlab="Time",ylab="",ylim=c(0,N),bty="n")
lines(Time,result$I,type="l",col="red")
lines(Time,result$R,type="l",col="black")
legend("topright",c("S","I","R"),col=c("green","red","black"),lty=1,bty="n")

# méthode ggplot
ggplot(result) +
  geom_line(aes(time, S, colour = "S")) +
  geom_line(aes(time, I, colour = "I")) +
  geom_line(aes(time, R, colour = "R")) +
  theme_bw()

# méthode ggplot + reshape (pour réorganiser le dataframe, plus pratique)
result %>%
  melt(id = "time") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable)) +
  theme_bw()

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
