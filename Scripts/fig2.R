
d_a = run0a %>%
  filter(time == max(time)) %>%
  mutate(Stot = Sa+S, CRtot = CR+CRa, CStot = CS+CSa,
         # abx = Sa+CRa+CSa, noabx=S+CR+CS,
         abx = Inc, noabx = 1-Inc,
         IRprop = IRa/(IRa+ISa), ISprop = ISa/(ISa+IRa)) %>%
  select(Stot, CRtot, CStot,
         abx, noabx,
         IRprop, ISprop) %>%
  melt() %>%
  mutate(panel_id = "Colonisation") %>%
  mutate(panel_id = replace(panel_id, variable %in% c("abx", "noabx"), "Abx exposure")) %>%
  mutate(panel_id = replace(panel_id, variable %in% c("IRprop", "ISprop"), "Infections")) %>%
  mutate(panel_id = factor(panel_id, levels = c("Colonisation", "Abx exposure", "Infections"))) %>%
  mutate(variable = factor(variable, levels = c("CRtot", "CStot", "Stot",
                                                "abx", "noabx",
                                                "IRprop", "ISprop"))) %>%
  mutate(bac = "S. aureus")

p_a = ggplot(d_a) +
  geom_bar(aes(panel_id, value, fill=variable), position="stack", stat="identity") +
  annotate("text", x="Colonisation", y=0.5, label = round(d_a$value[d_a$variable == "Stot"],2), colour = "white", size = 4) +
  annotate("text", x="Colonisation", y=0.85, label = round(d_a$value[d_a$variable == "CStot"],2), colour = "white", size = 4) +
  annotate("text", x="Abx exposure", y=0.9, label = round(d_a$value[d_a$variable == "abx"],2), colour = "white", size = 4) +
  annotate("text", x="Abx exposure", y=0.5, label = round(d_a$value[d_a$variable == "noabx"],2), colour = "white", size = 4) +
  annotate("text", x="Infections", y=0.5, label = paste0(round((tail(run0a$ISa, 1)) *100000), " / 100k"), colour = "white", size = 4) +
  annotate("text", x="Infections", y=0.95, label = paste0(round((tail(run0a$IRa, 1)) *100000), " / 100k"), colour = "white", size = 4) +
  scale_fill_discrete(type = c("#DE6F00", "#2072BA","#499124", 
                               "grey40", "grey60",
                               "#BD5E00", "#163F9E")) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  facet_wrap(vars(bac)) +
  labs(x = "", y = "Proportion") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        strip.text = element_text(size=11, face="italic"))

d_p = run0p %>%
  filter(time == max(time)) %>%
  mutate(Stot = Sa+S, CRtot = CR+CRa, CStot = CS+CSa,
         # abx = Sa+CRa+CSa, noabx=S+CR+CS,
         abx = Inc, noabx = 1-Inc,
         IRprop = IRa/(IRa+ISa), ISprop = ISa/(ISa+IRa)) %>%
  select(Stot, CRtot, CStot,
         abx, noabx,
         IRprop, ISprop) %>%
  melt() %>%
  mutate(panel_id = "Colonisation") %>%
  mutate(panel_id = replace(panel_id, variable %in% c("abx", "noabx"), "Abx exposure")) %>%
  mutate(panel_id = replace(panel_id, variable %in% c("IRprop", "ISprop"), "Infections")) %>%
  mutate(panel_id = factor(panel_id, levels = c("Colonisation", "Abx exposure", "Infections"))) %>%
  mutate(variable = factor(variable, levels = c("CRtot", "CStot", "Stot",
                                                "abx", "noabx",
                                                "IRprop", "ISprop"))) %>%
  mutate(bac = "S. pneumoniae")

p_p = ggplot(d_p) +
  geom_bar(aes(panel_id, value, fill=variable), position="stack", stat="identity") +
  annotate("text", x="Colonisation", y=0.5, label = round(d_p$value[d_p$variable == "Stot"],2), colour = "white", size = 4) +
  annotate("text", x="Colonisation", y=0.87, label = round(d_p$value[d_p$variable == "CStot"],2), colour = "white", size = 4) +
  annotate("text", x="Abx exposure", y=0.9, label = round(d_p$value[d_p$variable == "abx"],2), colour = "white", size = 4) +
  annotate("text", x="Abx exposure", y=0.5, label = round(d_p$value[d_p$variable == "noabx"],2), colour = "white", size = 4) +
  annotate("text", x="Infections", y=0.5, label = paste0(round((tail(run0p$ISa, 1)) *100000), " / 100k"), colour = "white", size = 4) +
  annotate("text", x="Infections", y=0.9, label = paste0(round((tail(run0p$IRa, 1)) *100000), " / 100k"), colour = "white", size = 4) +
  scale_fill_discrete(type = c("#DE6F00", "#2072BA","#499124", 
                               "grey40", "grey60",
                               "#BD5E00", "#163F9E")) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  facet_wrap(vars(bac)) +
  labs(x = "", y = "Proportion") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        strip.text = element_text(size=11, face="italic"))

d_e = run0e %>%
  filter(time == max(time)) %>%
  mutate(CRtot = CR+CRa, CStot = CS+CSa,
         # abx = CRa+CSa, noabx=CR+CS,
         abx = Inc, noabx = 1-Inc,
         IRprop = IRa/(IRa+ISa), ISprop = ISa/(ISa+IRa)) %>%
  select(CRtot, CStot,
         abx, noabx,
         IRprop, ISprop) %>%
  melt() %>%
  mutate(panel_id = "Colonisation") %>%
  mutate(panel_id = replace(panel_id, variable %in% c("abx", "noabx"), "Abx exposure")) %>%
  mutate(panel_id = replace(panel_id, variable %in% c("IRprop", "ISprop"), "Infections")) %>%
  mutate(panel_id = factor(panel_id, levels = c("Colonisation", "Abx exposure", "Infections"))) %>%
  mutate(variable = factor(variable, levels = c("CRtot", "CStot",
                                                "abx", "noabx",
                                                "IRprop", "ISprop"))) %>%
  mutate(bac = "E. coli")

p_e = ggplot(d_e) +
  geom_bar(aes(panel_id, value, fill=variable), position="stack", stat="identity") +
  annotate("text", x="Colonisation", y=0.5, label = round(d_e$value[d_e$variable == "CStot"],2), colour = "white", size = 4) +
  annotate("text", x="Colonisation", y=0.9, label = round(d_e$value[d_e$variable == "CRtot"],2), colour = "white", size = 4) +
  annotate("text", x="Abx exposure", y=0.9, label = round(d_e$value[d_e$variable == "abx"],2), colour = "white", size = 4) +
  annotate("text", x="Abx exposure", y=0.5, label = round(d_e$value[d_e$variable == "noabx"],2), colour = "white", size = 4) +
  annotate("text", x="Infections", y=0.5, label = paste0(round((tail(run0e$ISa, 1)) *100000), " / 100k"), colour = "white", size = 4) +
  annotate("text", x="Infections", y=0.9, label = paste0(round((tail(run0e$IRa, 1)) *100000), " / 100k"), colour = "white", size = 4) +
  scale_fill_discrete(type = c("#DE6F00", "#2072BA", 
                               "grey40", "grey60",
                               "#BD5E00", "#163F9E")) +
  scale_y_continuous(breaks = seq(0,1,0.2)) +
  facet_wrap(vars(bac)) +
  labs(x = "", y = "Proportion") +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        strip.text = element_text(size=11, face="italic"))

legend_col = get_plot_component(ggplot(d_a %>% filter(panel_id == "Colonisation")) +
                                  geom_bar(aes(panel_id, value, fill=variable), position="stack", stat="identity") +
                                  scale_fill_discrete(type = c("#DE6F00", "#2072BA","#499124"), 
                                                      labels = c("Colonised - R", "Colonised - S", "Uncolonised")) +
                                  labs(fill = "Colonisation:") +
                                  theme(legend.title = element_text(size=11, face="bold"),
                                        legend.text = element_text(size=11),
                                        legend.position = "bottom") +
                                  guides(fill=guide_legend(title.position="top", title.hjust=0.5)),
                                "guide-box-bot")


legend_abx = get_plot_component(ggplot(d_a %>% filter(panel_id == "Abx exposure")) +
                                  geom_bar(aes(panel_id, value, fill=variable), position="stack", stat="identity") +
                                  scale_fill_discrete(type = c("grey40", "grey60"), 
                                                      labels = c("Exposed", "Unexposed")) +
                                  labs(fill = "Antibiotic exposure:") +
                                  theme(legend.title = element_text(size=11, face="bold"),
                                        legend.text = element_text(size=11),
                                        legend.position = "bottom") +
                                  guides(fill=guide_legend(title.position="top", title.hjust=0.5)),
                                "guide-box-bot")

legend_inf = get_plot_component(ggplot(d_a %>% filter(panel_id == "Infections")) +
                                  geom_bar(aes(panel_id, value, fill=variable), position="stack", stat="identity") +
                                  scale_fill_discrete(type = c("#BD5E00", "#163F9E"), 
                                                      labels = c("Resistant", "Sensitive")) +
                                  labs(fill = "Infections:") +
                                  theme(legend.title = element_text(size=11, face="bold"),
                                        legend.text = element_text(size=11),
                                        legend.position = "bottom") +
                                  guides(fill=guide_legend(title.position="top", title.hjust=0.5)),
                                "guide-box-bot")

plot_1 = plot_grid(p_a + guides(fill="none"), p_p + guides(fill="none"), p_e + guides(fill="none"),
                   plot_grid(legend_col, legend_abx, legend_inf, ncol = 1))

p_flu = ggplot(combined_Incidence_flu %>% mutate(virus="Influenza") %>% filter(time <= 140)) +
  geom_point(aes(time, no_vaccination, colour = "0%"), size = 0.8) +
  geom_point(aes(time, vaccination_50, colour = "25%"), size = 0.8) +
  geom_point(aes(time, vaccination_80, colour = "50%"), size = 0.8) +
  scale_colour_discrete(type = c("darkorchid4", "darkorchid3", "darkorchid1")) +
  theme_bw() +
  facet_wrap(vars(virus)) +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        legend.text = element_text(size=11),
        strip.text = element_text(size=11)) +
  labs(x = "Time (days)", y = "Incidence rate per 100k", colour = "Vaccine coverage:")

p_rota = ggplot(combined_Incidence_rota %>% mutate(virus="Rotavirus") %>% filter(time <= 140)) +
  geom_point(aes(time, no_vaccination, colour = "0%"), size = 0.8) +
  geom_point(aes(time, vaccination_50, colour = "25%"), size = 0.8) +
  geom_point(aes(time, vaccination_80, colour = "50%"), size = 0.8) +
  scale_colour_discrete(type = c("darkorchid4", "darkorchid3", "darkorchid1")) +
  theme_bw() +
  facet_wrap(vars(virus)) +
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=11),
        legend.title = element_text(size=11, face = "bold"),
        legend.text = element_text(size=11),
        strip.text = element_text(size=11)) +
  labs(x = "Time (days)", y = "Incidence rate per 100k", colour = "Vaccine coverage:")

plot_2 = plot_grid(plot_grid(p_flu+guides(colour="none"), p_rota+guides(colour="none"), nrow=1), 
                   get_plot_component(p_flu+theme(legend.position = "bottom"), "guide-box-bot"),
                   nrow=2, rel_heights = c(1,0.1))

plot_grid(plot_1, plot_2, nrow=2,
          labels = c("a)", "b)"), hjust=0, rel_heights = c(1,0.7))

ggsave(here::here("Figures", "fig2.png"), bg = "white", height = 8, width = 8)
