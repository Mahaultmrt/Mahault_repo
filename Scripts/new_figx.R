

# pour chaque run, noter force de lien entre antibio et selection
# c'est le paramètre alpha pour aureus et pneumo
# et le paramètre phi pour e coli

ggplot(univar_vacc_50a) +
  geom_point(aes(alpha, incidenceR))

ggplot(univar_vacc_50a) +
  geom_point(aes(ATB, incidenceR))

ggplot(univar_vacc_50a) +
  geom_point(aes(group = as.factor(alpha), x=ATB, colour=as.factor(alpha), incidenceR)) +
  theme_bw()

ggplot(univar_vacc_50a) +
  geom_tile(aes(x=as.factor(attr_exp), y=as.factor(alpha), fill=incidenceR)) +
  theme_bw()


ggplot(univar_vacc_50p) +
  geom_point(aes(alpha, incidenceR))

ggplot(univar_vacc_50p) +
  geom_point(aes(attr_exp, incidenceR))

ggplot(univar_vacc_50p) +
  geom_tile(aes(x=as.factor(attr_exp), y=as.factor(alpha), fill=incidenceR))


ggplot(univar_vacc_50e) +
  geom_point(aes(phi, incidenceR))

ggplot(univar_vacc_50e) +
  geom_point(aes(attr_exp, incidenceR))

ggplot(univar_vacc_50e) +
  geom_tile(aes(x=as.factor(attr_exp), y=as.factor(phi), fill=incidenceR))

