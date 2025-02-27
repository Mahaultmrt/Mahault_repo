

dat_a = rbind(apply(IR_vacca, 2, function(x) (x-IR_vacca$`0`)/IR_vacca$`0`) %>%
        melt() %>%
        mutate(type="R"),
      
      apply(abx_vacca, 2, function(x) (x-abx_vacca$`0`)/abx_vacca$`0`) %>%
        melt() %>%
        mutate(type="abx")) %>%
  group_by(type, Var2) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup() %>%
  mutate(bac="S. aureus - influenza")

pa = ggplot(dat_a) +
  geom_line(aes(Var2, mean, colour = type), linewidth=1) +
  geom_point(aes(Var2, mean, colour = type), size=2) +
  geom_errorbar(aes(Var2, ymin = mean-sd, ymax = mean+sd, colour = type), linewidth=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  facet_grid(cols=vars(bac)) +
  annotate("text", 0.25,-0.1,label=paste0(dat_a %>%
                                           dcast(Var2~type, value.var = "mean") %>%
                                           mutate(ratio = R/abx) %>%
                                           select(ratio) %>%
                                           pull %>%
                                           mean(na.rm = T) %>%
                                           round(2), 
                                         "% decrease in R infections\nper 1% decrease in abx exposure"), size=3) +
  scale_colour_discrete(type = c("black","#BD5E00"), labels = c("Resistant infections","Antibiotic exposure"), breaks = c("R", "abx")) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        strip.text = element_text(size=11, face="bold")) +
  labs(x = "Vaccine coverage", y = "Relative difference in\ncumulative incidence/exposure", colour = "")

dat_p = rbind(apply(IR_vaccp, 2, function(x) (x-IR_vaccp$`0`)/IR_vaccp$`0`) %>%
        melt() %>%
        mutate(type="R"),
      
      apply(abx_vaccp, 2, function(x) (x-abx_vaccp$`0`)/abx_vaccp$`0`) %>%
        melt() %>%
        mutate(type="abx")) %>%
  group_by(type, Var2) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup() %>%
  mutate(bac="S. pneumoniae - influenza")

pp = ggplot(dat_p) +
  geom_line(aes(Var2, mean, colour = type), linewidth=1) +
  geom_point(aes(Var2, mean, colour = type), size=2) +
  geom_errorbar(aes(Var2, ymin = mean-sd, ymax = mean+sd, colour = type), linewidth=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  facet_grid(cols=vars(bac)) +
  annotate("text", 0.25,-0.1,label=paste0(dat_p %>%
                                           dcast(Var2~type, value.var = "mean") %>%
                                           mutate(ratio = R/abx) %>%
                                           select(ratio) %>%
                                           pull %>%
                                           mean(na.rm = T) %>%
                                           round(2), 
                                         "% decrease in R infections\nper 1% decrease in abx exposure"), size=3) +
  scale_colour_discrete(type = c("black","#BD5E00"), labels = c("Resistant infections","Antibiotic exposure"), breaks = c("R", "abx")) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        strip.text = element_text(size=11, face="bold")) +
  labs(x = "Vaccine coverage", y = "Relative difference in\ncumulative incidence/exposure", colour = "")

dat_e = rbind(apply(IR_vacce, 2, function(x) (x-IR_vacce$`0`)/IR_vacce$`0`) %>%
        melt() %>%
        mutate(type="R"),
      
      apply(abx_vacce, 2, function(x) (x-abx_vacce$`0`)/abx_vacce$`0`) %>%
        melt() %>%
        mutate(type="abx")) %>%
  group_by(type, Var2) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ungroup() %>%
  mutate(bac="E. coli - rotavirus")
  

pe = ggplot(dat_e) +
  geom_line(aes(Var2, mean, colour = type), linewidth=1) +
  geom_point(aes(Var2, mean, colour = type), size=2) +
  geom_errorbar(aes(Var2, ymin = mean-sd, ymax = mean+sd, colour = type), linewidth=0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth=1) +
  facet_grid(cols=vars(bac)) +
  annotate("text", 0.25,-0.2,label=paste0(dat_e %>%
                                         dcast(Var2~type, value.var = "mean") %>%
                                         mutate(ratio = R/abx) %>%
                                         select(ratio) %>%
                                         pull %>%
                                         mean(na.rm = T) %>%
                                         round(2), 
                                         "% decrease in R infections\nper 1% decrease in abx exposure"), size=3) +
  scale_colour_discrete(type = c("black","#BD5E00"), labels = c("Resistant infections","Antibiotic exposure"), breaks = c("R", "abx")) +
  theme_bw() +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.text = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        strip.text = element_text(size=11, face="bold")) +
  labs(x = "Vaccine coverage", y = "Relative difference in\ncumulative incidence/exposure", colour = "")

plot_grid(pa+guides(colour="none"),pp+guides(colour="none"),pe+guides(colour="none"),get_legend(pa),
          labels = c("a)", "b)", "c)", ""), hjust=0, align="v", axis="l")

ggsave(here::here("Figures", "fig4.png"), bg="white", height=7, width=8.5)

dat_a %>% filter(Var2 == 1)
dat_p %>% filter(Var2 == 1)
dat_e %>% filter(Var2 == 1)
