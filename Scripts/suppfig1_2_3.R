

p1a=ggplot(psa_Ra) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with\nabsolute change in R incidence")+
  scale_x_discrete(labels = c("alpha"=bquote(alpha),
                              "ATB"=bquote(ATB),
                              "beta"=bquote(beta),
                              "fitness"=bquote(f),
                              "gamma"=bquote(gamma),
                              "omega"=bquote(omega),
                              "rho"=bquote(rho),
                              "theta"=bquote(theta[0])))

p2a=ggplot(psa_Sa) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with\nabsolute change in S incidence")+
  scale_x_discrete(labels = c("alpha"=bquote(alpha),
                              "ATB"=bquote(ATB),
                              "beta"=bquote(beta),
                              "fitness"=bquote(f),
                              "gamma"=bquote(gamma),
                              "omega"=bquote(omega),
                              "rho"=bquote(rho),
                              "theta"=bquote(theta[0])))

p3a=ggplot(psa_SRa) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with\nabsolute change in total incidence")+
  scale_x_discrete(labels = c("alpha"=bquote(alpha),
                              "ATB"=bquote(ATB),
                              "beta"=bquote(beta),
                              "fitness"=bquote(f),
                              "gamma"=bquote(gamma),
                              "omega"=bquote(omega),
                              "rho"=bquote(rho),
                              "theta"=bquote(theta[0])))


plot_grid(p1a,p2a,p3a,ncol=2,labels=c("a)", "b)", "c)"), hjust=0)
ggsave(here::here("Figures", "suppfig1.png"), height = 6, width = 7, bg = "white")




p1p=ggplot(psa_Rp) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with\nabsolute change in R incidence")+
  scale_x_discrete(labels = c("alpha"=bquote(alpha),
                              "ATB"=bquote(ATB),
                              "beta"=bquote(beta),
                              "extra_rho"=bquote(RM),
                              "fitness"=bquote(f),
                              "gamma"=bquote(gamma),
                              "omega"=bquote(omega),
                              "rho"=bquote(rho),
                              "theta"=bquote(theta[0])))

p2p=ggplot(psa_Sp) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with\nabsolute change in S incidence")+
  scale_x_discrete(labels = c("alpha"=bquote(alpha),
                              "ATB"=bquote(ATB),
                              "beta"=bquote(beta),
                              "extra_rho"=bquote(RM),
                              "fitness"=bquote(f),
                              "gamma"=bquote(gamma),
                              "omega"=bquote(omega),
                              "rho"=bquote(rho),
                              "theta"=bquote(theta[0])))

p3p=ggplot(psa_SRp) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with\nabsolute change in total incidence")+
  scale_x_discrete(labels = c("alpha"=bquote(alpha),
                              "ATB"=bquote(ATB),
                              "beta"=bquote(beta),
                              "extra_rho"=bquote(RM),
                              "fitness"=bquote(f),
                              "gamma"=bquote(gamma),
                              "omega"=bquote(omega),
                              "rho"=bquote(rho),
                              "theta"=bquote(theta[0])))


plot_grid(p1p,p2p,p3p,ncol=2,labels=c("a)", "b)", "c)"), hjust=0)
ggsave(here::here("Figures", "suppfig2.png"), height = 6, width = 7, bg = "white")





p1e=ggplot(psa_Re) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with\nabsolute change in R incidence")+
  scale_x_discrete(labels = c("alpha"=bquote(alpha),
                              "ATB"=bquote(ATB),
                              "beta"=bquote(beta),
                              "fitness"=bquote(f),
                              "gamma"=bquote(gamma),
                              "omega"=bquote(omega),
                              "phi"=bquote(phi),
                              "rho"=bquote(rho),
                              "theta"=bquote(theta[0])))

p2e=ggplot(psa_Se) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with\nabsolute change in S incidence")+
  scale_x_discrete(labels = c("alpha"=bquote(alpha),
                              "ATB"=bquote(ATB),
                              "beta"=bquote(beta),
                              "fitness"=bquote(f),
                              "gamma"=bquote(gamma),
                              "omega"=bquote(omega),
                              "phi"=bquote(phi),
                              "rho"=bquote(rho),
                              "theta"=bquote(theta[0])))

p3e=ggplot(psa_SRe) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid") +
  geom_hline(yintercept = -1, linetype = "solid") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(yintercept = -0.5, linetype = "dashed") +
  geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.3) +
  theme_bw() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12)) +
  labs(colour = "", x = "Parameters", y = "Correlation coefficient with\nabsolute change in total incidence")+
  scale_x_discrete(labels = c("alpha"=bquote(alpha),
                              "ATB"=bquote(ATB),
                              "beta"=bquote(beta),
                              "fitness"=bquote(f),
                              "gamma"=bquote(gamma),
                              "omega"=bquote(omega),
                              "phi"=bquote(phi),
                              "rho"=bquote(rho),
                              "theta"=bquote(theta[0])))


plot_grid(p1e,p2e,p3e,ncol=2,labels=c("a)", "b)", "c)"), hjust=0)
ggsave(here::here("Figures", "suppfig3.png"), height = 6, width = 7, bg = "white")


