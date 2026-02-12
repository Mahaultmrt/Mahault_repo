

dat_a = apply(abx_vacca, 2, function(x) (x-abx_vacca$`0`)/abx_vacca$`0`) %>%
  melt() %>%
  mutate(type="abx") %>%
  group_by(Var2) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup() %>%
  mutate(bac="S. aureus - influenza")

pa_abx = ggplot(dat_a) +
  geom_line(aes(Var2, mean), linewidth=1, colour = "black") +
  geom_point(aes(Var2, mean), size=2, colour = "black") +
  geom_errorbar(aes(Var2, ymin = mean-sd, ymax = mean+sd), linewidth=0.8, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  facet_grid(cols=vars(bac)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        strip.text = element_text(size=11, face="bold")) +
  labs(x = "Vaccine coverage", y = "Relative difference in\ncumulative antibiotic exposure")

dat_p = apply(abx_vaccp, 2, function(x) (x-abx_vaccp$`0`)/abx_vaccp$`0`) %>%
  melt() %>%
  mutate(type="abx") %>%
  group_by(Var2) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup() %>%
  mutate(bac="S. pneumoniae - influenza")

pp_abx = ggplot(dat_p) +
  geom_line(aes(Var2, mean), linewidth=1, colour = "black") +
  geom_point(aes(Var2, mean), size=2, colour = "black") +
  geom_errorbar(aes(Var2, ymin = mean-sd, ymax = mean+sd), linewidth=0.8, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  facet_grid(cols=vars(bac)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        strip.text = element_text(size=11, face="bold")) +
  labs(x = "Vaccine coverage", y = "Relative difference in\ncumulative antibiotic exposure")

dat_e = apply(abx_vacce, 2, function(x) (x-abx_vacce$`0`)/abx_vacce$`0`) %>%
  melt() %>%
  mutate(type="abx") %>%
  group_by(Var2) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup() %>%
  mutate(bac="E. coli - rotavirus")


pe_abx = ggplot(dat_e) +
  geom_line(aes(Var2, mean), linewidth=1, colour = "black") +
  geom_point(aes(Var2, mean), size=2, colour = "black") +
  geom_errorbar(aes(Var2, ymin = mean-sd, ymax = mean+sd), linewidth=0.8, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  facet_grid(cols=vars(bac)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        strip.text = element_text(size=11, face="bold")) +
  labs(x = "Vaccine coverage", y = "Relative difference in\ncumulative antibiotic exposure")


diff_sima$bac = "S. aureus - influenza"
pa = ggplot(diff_sima) +
  geom_line(aes(vacc, mean_incidence_IR, colour = "Resistant"), linewidth=1) +
  geom_line(aes(vacc, mean_incidence_IS, colour = "Sensitive"), linewidth=1) +
  geom_line(aes(vacc, mean_incidence_ISR, colour = "Total"), linewidth=1) +
  geom_point(aes(vacc, mean_incidence_IR, colour = "Resistant"), size=2) +
  geom_point(aes(vacc, mean_incidence_IS, colour = "Sensitive"), size=2) +
  geom_point(aes(vacc, mean_incidence_ISR, colour = "Total"), size=2) +
  geom_errorbar(aes(vacc, ymin = IR_ic_l, ymax = IR_ic_u, colour = "Resistant"), linewidth=0.8) +
  geom_errorbar(aes(vacc, ymin = IS_ic_l, ymax = IS_ic_u, colour = "Sensitive"), linewidth=0.8) +
  geom_errorbar(aes(vacc, ymin = ISR_ic_l, ymax = ISR_ic_u, colour = "Total"), linewidth=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  scale_colour_discrete(type = c("#BD5E00", "#163F9E","#4B0082")) +
  facet_grid(cols = vars(bac)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        strip.text = element_text(size=11, face = "bold")) +
  labs(x = "Vaccine coverage", y = "Relative difference in\ncumulative incidence of infection", colour = "Strain:")

vacca = rbind(data.frame(cov0 = IS_vacca$`0`, cov1 = IS_vacca$`1`, bac = "S"),
              data.frame(cov0 = IR_vacca$`0`, cov1 = IR_vacca$`1`, bac = "R"),
              data.frame(cov0 = ISR_vacca$`0`, cov1 = ISR_vacca$`1`, bac = "tot")) %>%
  mutate(cov0 = cov0*100000, cov1 = cov1*100000) %>%
  melt() %>%
  group_by(bac, variable) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup() %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = replace(variable, variable == "cov0", "0"),
         variable = replace(variable, variable == "cov1", "1"))

ggplot(vacca) +
  geom_pointrange(aes(variable, mean, ymin=mean-sd, ymax=mean+sd, colour = bac), position = position_dodge(0.1)) +
  scale_colour_discrete(type = c("#BD5E00", "#163F9E","#4B0082")) +
  theme_bw() +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold")) +
  labs(x = "Vaccine coverage", y = "Cumulative incidence of infection (per 100k)") +
  guides(colour="none")

diff_simp$bac = "S. pneumoniae - influenza"
pp = ggplot(diff_simp) +
  geom_line(aes(vacc, mean_incidence_IR, colour = "Resistant"), linewidth=1) +
  geom_line(aes(vacc, mean_incidence_IS, colour = "Sensitive"), linewidth=1) +
  geom_line(aes(vacc, mean_incidence_ISR, colour = "Total"), linewidth=1) +
  geom_point(aes(vacc, mean_incidence_IR, colour = "Resistant"), size=2) +
  geom_point(aes(vacc, mean_incidence_IS, colour = "Sensitive"), size=2) +
  geom_point(aes(vacc, mean_incidence_ISR, colour = "Total"), size=2) +
  geom_errorbar(aes(vacc, ymin = IR_ic_l, ymax = IR_ic_u, colour = "Resistant"), linewidth=0.8) +
  geom_errorbar(aes(vacc, ymin = IS_ic_l, ymax = IS_ic_u, colour = "Sensitive"), linewidth=0.8) +
  geom_errorbar(aes(vacc, ymin = ISR_ic_l, ymax = ISR_ic_u, colour = "Total"), linewidth=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  scale_colour_discrete(type = c("#BD5E00", "#163F9E","#4B0082")) +
  facet_grid(cols = vars(bac)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        strip.text = element_text(size=11, face = "bold")) +
  labs(x = "Vaccine coverage", y = "Relative difference in\ncumulative incidence of infection", colour = "Strain:")

vaccp = rbind(data.frame(cov0 = IS_vaccp$`0`, cov1 = IS_vaccp$`1`, bac = "S"),
              data.frame(cov0 = IR_vaccp$`0`, cov1 = IR_vaccp$`1`, bac = "R"),
              data.frame(cov0 = ISR_vaccp$`0`, cov1 = ISR_vaccp$`1`, bac = "tot")) %>%
  mutate(cov0 = cov0*100000, cov1 = cov1*100000) %>%
  melt() %>%
  group_by(bac, variable) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup() %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = replace(variable, variable == "cov0", "0"),
         variable = replace(variable, variable == "cov1", "1"))

ggplot(vaccp) +
  geom_pointrange(aes(variable, mean, ymin=mean-sd, ymax=mean+sd, colour = bac), position = position_dodge(0.1)) +
  scale_colour_discrete(type = c("#BD5E00", "#163F9E","#4B0082")) +
  theme_bw() +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold")) +
  labs(x = "Vaccine coverage", y = "Cumulative incidence of infection (per 100k)") +
  guides(colour="none")

diff_sime$bac = "E. coli - rotavirus"
pe = ggplot(diff_sime) +
  geom_line(aes(vacc, mean_incidence_IR, colour = "Resistant"), linewidth=1) +
  geom_line(aes(vacc, mean_incidence_IS, colour = "Sensitive"), linewidth=1) +
  geom_line(aes(vacc, mean_incidence_ISR, colour = "Total"), linewidth=1) +
  geom_point(aes(vacc, mean_incidence_IR, colour = "Resistant"), size=2) +
  geom_point(aes(vacc, mean_incidence_IS, colour = "Sensitive"), size=2) +
  geom_point(aes(vacc, mean_incidence_ISR, colour = "Total"), size=2) +
  geom_errorbar(aes(vacc, ymin = IR_ic_l, ymax = IR_ic_u, colour = "Resistant"), linewidth=0.8) +
  geom_errorbar(aes(vacc, ymin = IS_ic_l, ymax = IS_ic_u, colour = "Sensitive"), linewidth=0.8) +
  geom_errorbar(aes(vacc, ymin = ISR_ic_l, ymax = ISR_ic_u, colour = "Total"), linewidth=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  scale_colour_discrete(type = c("#BD5E00", "#163F9E","#4B0082")) +
  facet_grid(cols = vars(bac)) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        strip.text = element_text(size=11, face = "bold")) +
  labs(x = "Vaccine coverage", y = "Relative difference in\ncumulative incidence of infection", colour = "Strain:")




plot_grid(plot_grid(pa_abx, pp_abx+labs(y=NULL), pe_abx+labs(y=NULL),
                    pa+guides(colour="none"), pp+guides(colour="none")+labs(y=NULL), pe+guides(colour="none")+labs(y=NULL),
                    labels = c("a)", "", "", "b)", "", ""), hjust=0,
                    byrow = T, ncol = 3, rel_widths = c(1,0.85,0.85)),
          get_legend(pa+theme(legend.position="bottom")), nrow=2, rel_heights = c(1,0.1))

ggsave(here::here("Figures", "fig4.png"), bg="white", height=6, width=8.5)

dat_a %>% filter(Var2 == 1)
diff_sima %>% filter(vacc == 1)
0.0031/0.16
dat_p %>% filter(Var2 == 1)
diff_simp %>% filter(vacc == 1)
0.0647/0.16
dat_e %>% filter(Var2 == 1)
diff_sime %>% filter(vacc == 1)
0.0409/0.423
