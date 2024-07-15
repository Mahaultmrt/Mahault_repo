#Code avec les équations des modèles
#Res_model est le modèle pour S.Pneumoniae et S.Aureus

Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+S+CR+CS
    
    
    new_teta<-teta-log(1-vec_virus(t)*ATB)
    
    
    dSa <- -Sa*((beta*ct*(CRa+CR)/N)+beta*(CSa+CS)/N)-omega*Sa+new_teta*S+(gamma+alpha*(1-sigmaR))*CRa+(gamma+alpha*(1-sigmaS))*CSa
    dCRa <- Sa*(beta*ct*(CRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-omega*CRa+new_teta*CR
    dCSa <- Sa*(beta*(CSa+CS)/N)-(gamma+alpha*(1-sigmaS))*CSa-omega*CSa+new_teta*CS
    dIRa <- rhoRa*CRa-deltaRa*IRa+rhoR*CR
    dISa <- rhoSa*CSa-deltaSa*ISa+rhoS*CS
    
    dS <- -S*((beta*ct*(CRa+CR)/N)+beta*(CSa+CS)/N)+omega*Sa-new_teta*S+gamma*CR+gamma*CS+deltaRa*IRa+deltaSa*ISa
    dCR <- S*(beta*ct*(CRa+CR)/N)-gamma*CR+omega*CRa-new_teta*CR
    dCS <- S*(beta*(CSa+CS)/N)-gamma*CS+omega*CSa-new_teta*CS
    
    
    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS)
    
    
    list(res)
    
  })
  
}

Res_model_Coli <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=CSa+CRa+CS+CR
    
    
    new_teta<-teta-log(1-vec_virus(t)*ATB)
    
    dCSa<- - CSa*(beta*ct*phi*(CRa+CR)/N)+(gamma+alpha*(1-sigmaR))*CRa-CSa*omega+new_teta*CS
    dCRa<- CSa*(beta*ct*phi*(CRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-CRa*omega+new_teta*CR
    dCS<- -CS*(beta*ct*(CRa+CR)/N)+gamma*CR+CSa*omega-new_teta*CS
    dCR<- CS*(beta*ct*(CRa+CR)/N)-gamma*CR+CRa*omega-new_teta*CR
    dIRa<- rhoRa*CRa+rho*CR
    dISa<-rhoSa*CSa+rho*CS
    
    res<-c(dCSa,dCRa,dCS,dCR,dIRa,dISa)
    
    list(res)
    
  })
  
}
