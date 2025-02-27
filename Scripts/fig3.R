
pa = all_runa %>%
  mutate(variable = replace(variable, variable == "Sa", "S"),
         variable = replace(variable, variable == "CSa", "CS"),
         variable = replace(variable, variable == "CRa", "CR")) %>%
  mutate(vaccination = replace(vaccination, vaccination == "vacc 0%", "0% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 50%", "50% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 80%", "80% coverage")) %>%
  group_by(time, vaccination, variable) %>%
  summarise(value = sum(value)) %>%
  ungroup %>%
  mutate(variable = factor(variable, levels = c("S", "CS", "CR", "ISa", "IRa"))) %>%
  filter(variable %in% c("S", "CS", "CR")) %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 1) +
  geom_hline(yintercept=tail(run0a$Sa, 1)+tail(run0a$S, 1), linetype="dashed",color="#499124",alpha=0.5, linewidth = 1)+
  geom_hline(yintercept=tail(run0a$CRa, 1)+tail(run0a$CR, 1), linetype="dashed",color="#DE6F00",alpha=0.5, linewidth = 1)+
  geom_hline(yintercept=tail(run0a$CSa, 1)+tail(run0a$CS, 1), linetype="dashed",color="#2072BA",alpha=0.5, linewidth = 1)+
  scale_colour_discrete(type = c("#499124", "#2072BA","#DE6F00"), labels = c("Uncolonised", "Colonised - S", "Colonised - R")) +
  labs(x = "Time (days)", y = "Prop. colonised", colour = "Colonisation status:") +
  facet_grid(cols = vars(vaccination)) +
  theme_bw() +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        strip.text = element_text(size=11))

pa_res = all_runa %>%
  mutate(variable = replace(variable, variable == "Sa", "S"),
         variable = replace(variable, variable == "CSa", "CS"),
         variable = replace(variable, variable == "CRa", "CR")) %>%
  filter(variable %in% c("CS", "CR")) %>%
  mutate(vaccination = replace(vaccination, vaccination == "vacc 0%", "0% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 50%", "50% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 80%", "80% coverage")) %>%
  group_by(time, vaccination, variable) %>%
  summarise(value = sum(value)) %>%
  ungroup %>%
  group_by(time, vaccination) %>%
  mutate(value=value/(sum(value))) %>%
  ungroup %>%
  filter(variable=="CR") %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 1) +
  geom_hline(yintercept=(tail(run0a$CRa, 1)+tail(run0a$CR, 1))/(tail(run0a$CRa, 1)+tail(run0a$CR, 1)+tail(run0a$CSa, 1)+tail(run0a$CS, 1)), linetype="dashed",color="#DE6F00",alpha=0.5, linewidth = 1)+
  labs(x = "Time (days)", y = "Prop. colonised", colour = "Colonisation status:") +
  facet_grid(cols = vars(vaccination)) +
  theme_bw() +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        strip.text = element_text(size=11))
  

pa_abx = all_runa %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = replace(variable, variable == "Sa", "abx"),
         variable = replace(variable, variable == "CSa", "abx"),
         variable = replace(variable, variable == "CRa", "abx"),
         variable = replace(variable, variable == "S", "noabx"),
         variable = replace(variable, variable == "CS", "noabx"),
         variable = replace(variable, variable == "CR", "noabx")) %>%
  mutate(vaccination = replace(vaccination, vaccination == "vacc 0%", "0% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 50%", "50% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 80%", "80% coverage")) %>%
  group_by(time, vaccination, variable) %>%
  summarise(value = sum(value)) %>%
  ungroup %>%
  mutate(variable = factor(variable, levels = c("abx", "noabx", "ISa", "IRa"))) %>%
  filter(variable %in% c("abx")) %>%
  ggplot() +
  geom_line(aes(time, value), colour = "black", linewidth = 1) +
  geom_hline(yintercept=tail(run0a$Sa, 1)+tail(run0a$CRa, 1)+tail(run0a$CSa, 1), linetype="dashed",color="black",alpha=0.7, linewidth = 1)+
  labs(x = "Time (days)", y = "Prop. abx exposed", colour = "Abx status:") +
  facet_grid(cols = vars(vaccination)) +
  theme_bw() +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        strip.text = element_blank())

pa_all = plot_grid(pa + guides(colour = "none") + theme(axis.text.x = element_blank(), axis.title.x = element_blank()), pa_abx,
                   ncol=1, align = "v", rel_heights = c(1,0.8))



pp = all_runp %>%
  mutate(variable = replace(variable, variable == "Sa", "S"),
         variable = replace(variable, variable == "CSa", "CS"),
         variable = replace(variable, variable == "CRa", "CR")) %>%
  mutate(vaccination = replace(vaccination, vaccination == "vacc 0%", "0% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 50%", "50% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 80%", "80% coverage")) %>%
  group_by(time, vaccination, variable) %>%
  summarise(value = sum(value)) %>%
  ungroup %>%
  mutate(variable = factor(variable, levels = c("S", "CS", "CR", "ISa", "IRa"))) %>%
  filter(variable %in% c("S", "CS", "CR")) %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 1) +
  geom_hline(yintercept=tail(run0p$Sa, 1)+tail(run0p$S, 1), linetype="dashed",color="#499124",alpha=0.5, linewidth = 1)+
  geom_hline(yintercept=tail(run0p$CRa, 1)+tail(run0p$CR, 1), linetype="dashed",color="#DE6F00",alpha=0.5, linewidth = 1)+
  geom_hline(yintercept=tail(run0p$CSa, 1)+tail(run0p$CS, 1), linetype="dashed",color="#2072BA",alpha=0.5, linewidth = 1)+
  scale_colour_discrete(type = c("#499124", "#2072BA","#DE6F00")) +
  labs(x = "Time (days)", y = "Prop. colonised", colour = "Colonisation status:") +
  facet_grid(cols = vars(vaccination)) +
  theme_bw() +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        strip.text = element_text(size=11))



pp_abx = all_runp %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = replace(variable, variable == "Sa", "abx"),
         variable = replace(variable, variable == "CSa", "abx"),
         variable = replace(variable, variable == "CRa", "abx"),
         variable = replace(variable, variable == "S", "noabx"),
         variable = replace(variable, variable == "CS", "noabx"),
         variable = replace(variable, variable == "CR", "noabx")) %>%
  mutate(vaccination = replace(vaccination, vaccination == "vacc 0%", "0% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 50%", "50% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 80%", "80% coverage")) %>%
  group_by(time, vaccination, variable) %>%
  summarise(value = sum(value)) %>%
  ungroup %>%
  mutate(variable = factor(variable, levels = c("abx", "noabx", "ISa", "IRa"))) %>%
  filter(variable %in% c("abx")) %>%
  ggplot() +
  geom_line(aes(time, value), colour = "black", linewidth = 1) +
  geom_hline(yintercept=tail(run0p$Sa, 1)+tail(run0p$CRa, 1)+tail(run0p$CSa, 1), linetype="dashed",color="black",alpha=0.7, linewidth = 1)+
  labs(x = "Time (days)", y = "Prop. abx exposed", colour = "Abx status:") +
  facet_grid(cols = vars(vaccination)) +
  theme_bw() +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        strip.text = element_blank())

pp_all = plot_grid(pp + guides(colour = "none") + theme(axis.text.x = element_blank(), axis.title.x = element_blank()), pp_abx,
                   ncol=1, align = "v", rel_heights = c(1,0.8))


pe = all_rune %>%
  mutate(variable = replace(variable, variable == "CSa", "CS"),
         variable = replace(variable, variable == "CRa", "CR")) %>%
  mutate(vaccination = replace(vaccination, vaccination == "vacc 0%", "0% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 50%", "50% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 80%", "80% coverage")) %>%
  group_by(time, vaccination, variable) %>%
  summarise(value = sum(value)) %>%
  ungroup %>%
  mutate(variable = factor(variable, levels = c("CS", "CR", "ISa", "IRa"))) %>%
  filter(variable %in% c("CS", "CR")) %>%
  ggplot() +
  geom_line(aes(time, value, colour = variable), linewidth = 1) +
  geom_hline(yintercept=tail(run0e$CRa, 1)+tail(run0e$CR, 1), linetype="dashed",color="#DE6F00",alpha=0.5, linewidth = 1)+
  geom_hline(yintercept=tail(run0e$CSa, 1)+tail(run0e$CS, 1), linetype="dashed",color="#2072BA",alpha=0.5, linewidth = 1)+
  scale_colour_discrete(type = c("#2072BA","#DE6F00")) +
  labs(x = "Time (days)", y = "Prop. colonised", colour = "Colonisation status:") +
  facet_grid(cols = vars(vaccination)) +
  theme_bw() +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        strip.text = element_text(size=11))



pe_abx = all_rune %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = replace(variable, variable == "CSa", "abx"),
         variable = replace(variable, variable == "CRa", "abx"),
         variable = replace(variable, variable == "CS", "noabx"),
         variable = replace(variable, variable == "CR", "noabx")) %>%
  mutate(vaccination = replace(vaccination, vaccination == "vacc 0%", "0% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 50%", "50% coverage"),
         vaccination = replace(vaccination, vaccination == "vacc 80%", "80% coverage")) %>%
  group_by(time, vaccination, variable) %>%
  summarise(value = sum(value)) %>%
  ungroup %>%
  mutate(variable = factor(variable, levels = c("abx", "noabx", "ISa", "IRa"))) %>%
  filter(variable %in% c("abx")) %>%
  ggplot() +
  geom_line(aes(time, value), colour = "black", linewidth = 1) +
  geom_hline(yintercept=tail(run0e$CRa, 1)+tail(run0e$CSa, 1), linetype="dashed",color="black",alpha=0.7, linewidth = 1)+
  labs(x = "Time (days)", y = "Prop. abx exposed", colour = "Abx status:") +
  facet_grid(cols = vars(vaccination)) +
  theme_bw() +
  theme(axis.title = element_text(size=11),
        axis.text = element_text(size=11),
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        strip.text = element_blank())

pe_all = plot_grid(pe + guides(colour = "none") + theme(axis.text.x = element_blank(), axis.title.x = element_blank()), pe_abx,
                   ncol=1, align = "v", rel_heights = c(1,0.8))


plot_grid(pa_all, pp_all, pe_all, get_legend(pa),
          labels = c("a)", "b)", "c)", ""), hjust=0)
ggsave(here::here("Figures", "fig3.png"), bg="white", height = 8, width = 9)

all_runa %>%
  filter(time == max(time)) %>%
  filter(variable %in% c("IRa", "ISa")) %>%
  dcast(vaccination~variable) %>%
  mutate(IRa=IRa*100000, ISa=ISa*100000) %>%
  mutate(prop_IR = IRa/(IRa+ISa), prop_IS = ISa/(ISa+IRa))

all_runp %>%
  filter(time == max(time)) %>%
  filter(variable %in% c("IRa", "ISa")) %>%
  dcast(vaccination~variable) %>%
  mutate(IRa=IRa*100000, ISa=ISa*100000) %>%
  mutate(prop_IR = IRa/(IRa+ISa), prop_IS = ISa/(ISa+IRa))

all_rune %>%
  filter(time == max(time)) %>%
  filter(variable %in% c("IRa", "ISa")) %>%
  dcast(vaccination~variable) %>%
  mutate(IRa=IRa*100000, ISa=ISa*100000) %>%
  mutate(prop_IR = IRa/(IRa+ISa), prop_IS = ISa/(ISa+IRa))
