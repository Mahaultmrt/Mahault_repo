#Code avec les équations des modèles

Res_model_aureus <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+S+CR+CS
    
    
    new_theta<-theta-log(1-vec_virus(t)*ATB)
    
    
    dSa <- -Sa*((beta*fitness*(CRa+CR)/N)+beta*(CSa+CS)/N)-omega*Sa+new_theta*S+(gamma+alpha*(1-sigmaR))*CRa+(gamma+alpha*(1-sigmaS))*CSa
    dCRa <- Sa*(beta*fitness*(CRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-omega*CRa+new_theta*CR
    dCSa <- Sa*(beta*(CSa+CS)/N)-(gamma+alpha*(1-sigmaS))*CSa-omega*CSa+new_theta*CS
    dIRa <- rhoRa*CRa-deltaRa*IRa+rhoR*CR
    dISa <- rhoSa*CSa-deltaSa*ISa+rhoS*CS
    
    dS <- -S*((beta*fitness*(CRa+CR)/N)+beta*(CSa+CS)/N)+omega*Sa-new_theta*S+gamma*CR+gamma*CS+deltaRa*IRa+deltaSa*ISa
    dCR <- S*(beta*fitness*(CRa+CR)/N)-gamma*CR+omega*CRa-new_theta*CR
    dCS <- S*(beta*(CSa+CS)/N)-gamma*CS+omega*CSa-new_theta*CS
    
    dInc<-new_theta*S+new_theta*CR+new_theta*CS
    

    res<-c(dSa,dCRa,dCSa,dIRa,dISa,dS,dCR,dCS,dInc)
    

    
    
    list(res)
    
  })
  
}

Res_model_pneumo <- function(t, pop, param,vec_virus) {
  
  with(as.list(c(pop, param)), {
    
    
    N=Sa+CRa+CSa+S+CR+CS
    
    
    new_theta<-theta-log(1-vec_virus(t)*ATB)
    # Assumin that there's a baseline incidence, and virus infected people have a higher risk extra_rho
    new_rhoRa = rhoRa*(1-vec_virus(t))+(rhoRa*extra_rho)*vec_virus(t)
    new_rhoSa = rhoSa*(1-vec_virus(t))+(rhoSa*extra_rho)*vec_virus(t)
    new_rhoR = rhoR*(1-vec_virus(t))+(rhoR*extra_rho)*vec_virus(t)
    new_rhoS = rhoS*(1-vec_virus(t))+(rhoS*extra_rho)*vec_virus(t)
    
    dSa <- -Sa*((beta*fitness*(CRa+CR)/N)+beta*(CSa+CS)/N)-omega*Sa+new_theta*S+(gamma+alpha*(1-sigmaR))*CRa+(gamma+alpha*(1-sigmaS))*CSa
    dCRa <- Sa*(beta*fitness*(CRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-omega*CRa+new_theta*CR
    dCSa <- Sa*(beta*(CSa+CS)/N)-(gamma+alpha*(1-sigmaS))*CSa-omega*CSa+new_theta*CS
    dIRa <- new_rhoRa*CRa-deltaRa*IRa+new_rhoR*CR
    dISa <- new_rhoSa*CSa-deltaSa*ISa+new_rhoS*CS
    
    dS <- -S*((beta*fitness*(CRa+CR)/N)+beta*(CSa+CS)/N)+omega*Sa-new_theta*S+gamma*CR+gamma*CS+deltaRa*IRa+deltaSa*ISa
    dCR <- S*(beta*fitness*(CRa+CR)/N)-gamma*CR+omega*CRa-new_theta*CR
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
    
    dCSa<- - CSa*(beta*fitness*phi*(CRa+CR)/N)+(gamma+alpha*(1-sigmaR))*CRa-CSa*omega+new_theta*CS
    dCRa<- CSa*(beta*fitness*phi*(CRa+CR)/N)-(gamma+alpha*(1-sigmaR))*CRa-CRa*omega+new_theta*CR
    dCS<- -CS*(beta*fitness*(CRa+CR)/N)+gamma*CR+CSa*omega-new_theta*CS
    dCR<- CS*(beta*fitness*(CRa+CR)/N)-gamma*CR+CRa*omega-new_theta*CR
    dIRa<- rhoRa*CRa+rho*CR
    dISa<-rhoSa*CSa+rho*CS
    
    dInc<-new_theta*CR+new_theta*CS
    
    res<-c(dCSa,dCRa,dCS,dCR,dIRa,dISa,dInc)
    
    list(res)
    
  })
  
}
