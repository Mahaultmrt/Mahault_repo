

p_aureus = ggplot() +
  geom_line(aes(seq(0:140), -log(1-I_vac_0_flu(0:140)*0.1)/(-log(1-I_vac_0_flu(0:140)*0.1)+0.0015)), linewidth=2) +
  geom_hline(yintercept=round(mean(-log(1-I_vac_0_flu(0:140)*0.1)/(-log(1-I_vac_0_flu(0:140)*0.1)+0.0015)),2), linewidth=2, linetype="dashed") +
  facet_grid(cols=vars("S. aureus - influenza - amoxicillin")) +
  annotate("text", x=17, y=0.2, label = round(mean(-log(1-I_vac_0_flu(0:140)*0.1)/(-log(1-I_vac_0_flu(0:140)*0.1)+0.0015)),2), size = 6) +
  labs(x = "Time (days)", y = "Proportion of daily antibiotic exposure\nattributable to viral infections") +
  scale_x_continuous(breaks = c(0,35,70,105,140)) +
  # scale_y_continuous(breaks = seq(0,0.003,0.001), limits=c(0,0.004)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text = element_text(size=11))

p_pneumo = ggplot() +
  geom_line(aes(seq(0:140), -log(1-I_vac_0_flu(0:140)*0.1)/(-log(1-I_vac_0_flu(0:140)*0.1)+0.0015)), linewidth=2) +
  geom_hline(yintercept=round(mean(-log(1-I_vac_0_flu(0:140)*0.1)/(-log(1-I_vac_0_flu(0:140)*0.1)+0.0015)),2), linewidth=2, linetype="dashed") +
  facet_grid(cols=vars("S. pneumoniae - influenza - amoxicillin")) +
  annotate("text", x=17, y=0.2, label = round(mean(-log(1-I_vac_0_flu(0:140)*0.1)/(-log(1-I_vac_0_flu(0:140)*0.1)+0.0015)),2), size = 6) +
  labs(x = "Time (days)", y = "Proportion of daily antibiotic exposure\nattributable to viral infections") +
  scale_x_continuous(breaks = c(0,35,70,105,140)) +
  # scale_y_continuous(breaks = seq(0,0.003,0.001), limits=c(0,0.004)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text = element_text(size=11))

p_coli = ggplot() +
  geom_line(aes(seq(0:140), -log(1-I_vac_0_rota(0:140)*0.1)/(-log(1-I_vac_0_rota(0:140)*0.1)+0.0015)), linewidth=2) +
  geom_hline(yintercept=round(mean(-log(1-I_vac_0_rota(0:140)*0.1)/(-log(1-I_vac_0_rota(0:140)*0.1)+0.0015)),2), linewidth=2, linetype="dashed") +
  facet_grid(cols=vars("E. coli - rotavirus - fluoroquinolone")) +
  annotate("text", x=17, y=0.5, label = round(mean(-log(1-I_vac_0_rota(0:140)*0.1)/(-log(1-I_vac_0_rota(0:140)*0.1)+0.0015)),2), size = 6) +
  labs(x = "Time (days)", y = "Proportion of daily antibiotic exposure\nattributable to viral infections") +
  scale_x_continuous(breaks = c(0,35,70,105,140)) +
  # scale_y_continuous(breaks = seq(0,0.003,0.001), limits=c(0,0.004)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text = element_text(size=11))




p_aureus_R = ggplot() +
  geom_line(aes(seq(0:140), (run2a$CRa+run2a$CR)/tail(run0a$CRa+run0a$CR, 1)),
            colour = "#BD5E00", linewidth=2) +
  facet_grid(cols=vars("S. aureus - influenza - amoxicillin")) +
  labs(x = "Time (days)", y = "Relative change in daily resistance\ncarriage prevalence over the season") +
  scale_x_continuous(breaks = c(0,35,70,105,140)) +
  # scale_y_continuous(breaks = seq(0,0.003,0.001), limits=c(0,0.004)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text = element_text(size=11))

p_pneumo_R = ggplot() +
  geom_line(aes(seq(0:140), (run2p$CRa+run2p$CR)/tail(run0p$CRa+run0p$CR, 1)),
            colour = "#BD5E00", linewidth=2) +
  facet_grid(cols=vars("S. pneumoniae - influenza - amoxicillin")) +
  labs(x = "Time (days)", y = "Relative change in daily resistance\ncarriage prevalence over the season") +
  scale_x_continuous(breaks = c(0,35,70,105,140)) +
  # scale_y_continuous(breaks = seq(0,0.003,0.001), limits=c(0,0.004)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text = element_text(size=11))

p_coli_R = ggplot() +
  geom_line(aes(seq(0:140), (run2e$CRa+run2e$CR)/tail(run0e$CRa+run0e$CR, 1)),
            colour = "#BD5E00", linewidth=2) +
  facet_grid(cols=vars("E. coli - rotavirus - fluoroquinolone")) +
  labs(x = "Time (days)", y = "Relative change in daily resistance\ncarriage prevalence over the season") +
  scale_x_continuous(breaks = c(0,35,70,105,140)) +
  # scale_y_continuous(breaks = seq(0,0.003,0.001), limits=c(0,0.004)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        axis.text = element_text(size=11))

plot_grid(p_aureus, p_pneumo, p_coli,
          p_aureus_R, p_pneumo_R, p_coli_R,
          labels = c("a)", "b)", "c)", "d)", "e)", "f)"), nrow = 2, hjust = 0)

ggsave(here::here("Figures", "figx.png"), bg="white", height=7)

