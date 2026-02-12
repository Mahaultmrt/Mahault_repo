# Analyse de correlation partielle
psa_Ra<-psa_vacc_50a[,-c(10:11)]
psa_Ra = psa_Ra %>%
  dplyr::select(beta,fitness,gamma,rho,alpha,theta,omega,ATB,incidenceR) %>%
  mutate(incidenceR=abs(incidenceR)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorR_A<-graph_pcor(psa_Ra)

psa_Sa<-psa_vacc_50a[,-c(9,11)]
psa_Sa = psa_Sa %>%
  dplyr::select(beta,fitness,gamma,rho,alpha,theta,omega,ATB,incidenceS) %>%
  mutate(incidenceS=abs(incidenceS)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorS_A<-graph_pcor(psa_Sa)


psa_SRa<-psa_vacc_50a[,-c(9:10)]
psa_SRa = psa_SRa %>%
  dplyr::select(beta,fitness,gamma,rho,alpha,theta,omega,ATB,incidenceT) %>%
  mutate(incidenceT=abs(incidenceT)) %>%
  epi.prcc() %>%
  rename(param = var)
pcorSR_A<-graph_pcor(psa_SRa)



# Analyse de correlation partielle
psa_Rp<-psa_vacc_50p[,-c(11:12)]
psa_Rp = psa_Rp %>%
  dplyr::select(beta,fitness,gamma,rho,extra_rho,alpha,theta,omega,ATB,incidenceR) %>%
  mutate(incidenceR=abs(incidenceR)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorR_P<-graph_pcor(psa_Rp)

psa_Sp<-psa_vacc_50p[,-c(10,12)]
psa_Sp = psa_Sp %>%
  dplyr::select(beta,fitness,gamma,rho,extra_rho,alpha,theta,omega,ATB,incidenceS) %>%
  mutate(incidenceS=abs(incidenceS)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorS_P<-graph_pcor(psa_Sp)

psa_SRp<-psa_vacc_50p[,-c(10,11)]
psa_SRp = psa_SRp %>%
  dplyr::select(beta,fitness,gamma,rho,extra_rho,alpha,theta,omega,ATB,incidenceT) %>%
  mutate(incidenceT=abs(incidenceT)) %>%
  epi.prcc() %>%
  rename(param = var)
pcorSR_P<-graph_pcor(psa_SRp)



# Analyse de correlation partielle
psa_Re<-psa_vacc_50e[,-c(10:11)]
psa_Re = psa_Re %>%
  dplyr::select(beta,fitness,gamma,rho,phi,theta,omega,ATB,incidenceR) %>%
  mutate(incidenceR=abs(incidenceR)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorR_E<-graph_pcor(psa_Re)

psa_Se<-psa_vacc_50e[,-c(9,11)]
psa_Se = psa_Se %>%
  dplyr::select(beta,fitness,gamma,rho,phi,theta,omega,ATB,incidenceS) %>%
  mutate(incidenceS=abs(incidenceS)) %>%
  epi.prcc() %>%
  rename(param = var)

pcorS_E<-graph_pcor(psa_Se)


psa_SRe<-psa_vacc_50e[,-c(9,10)]
psa_SRe = psa_SRe %>%
  dplyr::select(beta,fitness,gamma,rho,phi,theta,omega,ATB,incidenceT) %>%
  mutate(incidenceT=abs(incidenceT)) %>%
  epi.prcc() %>%
  rename(param = var)
pcorSR_E<-graph_pcor(psa_SRe)

