#Code avec les équations des modèles
#Res_model est le modèle pour S.Pneumoniae et S.Aureus
Res_model <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+S+CR+CS
    
    
    new_theta<-theta-log(1-vec_virus(t)*ATB)
    
    
    dSa <- -Sa*((beta*ct*(CRa+CR)/N)+beta*(CSa+CS)/N)-omega*Sa+new_theta*S+(gamma+alpha*(1-sigmaR))*CRa+(gamma+alpha*(1-sigmaS))*CSa
    dCRa <- Sa*(beta*ct*(CRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-omega*CRa+new_theta*CR
    dCSa <- Sa*(beta*(CSa+CS)/N)-(gamma+alpha*(1-sigmaS))*CSa-omega*CSa+new_theta*CS
    dIRa <- rhoRa*CRa-deltaRa*IRa+rhoR*CR
    dISa <- rhoSa*CSa-deltaSa*ISa+rhoS*CS
    
    dS <- -S*((beta*ct*(CRa+CR)/N)+beta*(CSa+CS)/N)+omega*Sa-new_theta*S+gamma*CR+gamma*CS+deltaRa*IRa+deltaSa*ISa
    dCR <- S*(beta*ct*(CRa+CR)/N)-gamma*CR+omega*CRa-new_theta*CR
    dCS <- S*(beta*(CSa+CS)/N)-gamma*CS+omega*CSa-new_theta*CS
    
    dInc<-new_theta*S+new_theta*CR+new_theta*CS
    

    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS,dInc)
    

    
    
    list(res)
    
  })
  
}

Res_model_Coli <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=CSa+CRa+CS+CR
    
    
    new_theta<-theta-log(1-vec_virus(t)*ATB)
    
    dCSa<- - CSa*(beta*ct*phi*(CRa+CR)/N)+(gamma+alpha*(1-sigmaR))*CRa-CSa*omega+new_theta*CS
    dCRa<- CSa*(beta*ct*phi*(CRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-CRa*omega+new_theta*CR
    dCS<- -CS*(beta*ct*(CRa+CR)/N)+gamma*CR+CSa*omega-new_theta*CS
    dCR<- CS*(beta*ct*(CRa+CR)/N)-gamma*CR+CRa*omega-new_theta*CR
    dIRa<- rhoRa*CRa+rho*CR
    dISa<-rhoSa*CSa+rho*CS
    
    dInc<-new_theta*CR+new_theta*CS
    
    res<-c(dCSa,dCRa,dCS,dCR,dIRa,dISa,dInc)
    
    list(res)
    
  })
  
}
