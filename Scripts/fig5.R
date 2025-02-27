
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

plot_grid(pa+guides(colour="none"), pp+guides(colour="none"), pe+guides(colour="none"), get_legend(pa), labels = c("a)", "b)", "c)", ""), hjust=0)

ggsave(here::here("Figures", "fig5.png"), height = 6, width = 8, bg="white")

tail(diff_sima, 1)
tail(diff_simp, 1)
tail(diff_sime, 1)


