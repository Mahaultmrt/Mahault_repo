
source(here::here("Scripts", "SIR_Influenza.R"))

p1=ggplot() +
  geom_line(aes(seq(0:140), I_vac_0_flu(0:140)), linewidth=2, colour="grey70") +
  labs(x = "Time (days)", y = "Virus incidence") +#bquote({I^V}(t)~~+~~{I^NV}(t))) +
  scale_x_continuous(breaks = c(0,30,60,90,120)) +
  # scale_y_continuous(breaks = seq(0,0.03,0.01)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12))

p2=ggplot() +
  geom_line(aes(seq(0:140), I_vac_0_flu(0:140)*0.1), linewidth=2, colour="grey50") +
  labs(x = "Time (days)", y = "Abx exposure\nincidence") +#bquote({P^V}(t))) +
  scale_x_continuous(breaks = c(0,30,60,90,120)) +
  # scale_y_continuous(breaks = seq(0,0.003,0.001)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12))

p3=ggplot() +
  geom_line(aes(seq(0:140), rep(0.0015,141)), linetype="dashed", linewidth=2, colour="grey50") +
  labs(x = "Time (days)", y = "Baseline rate of\nabx exp.") +#bquote({theta^0}(t))) +
  scale_x_continuous(breaks = c(0,30,60,90,120)) +
  scale_y_continuous(breaks = seq(0.0009,0.0018,0.0003), limits=c(0.001,0.002)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12))

p4=ggplot() +
  geom_line(aes(seq(0:140), (-log(1-I_vac_0_flu(0:140)*0.1)+0.0015)), linewidth=2) +
  labs(x = "Time (days)", y = "Total rate of abx exp.") +#bquote({theta^V}(t))) +
  scale_x_continuous(breaks = c(0,30,60,90,120)) +
  # scale_y_continuous(breaks = seq(0,0.003,0.001), limits=c(0,0.004)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12))

mean(-log(1-I_vac_0_flu(0:140)*0.1)/(-log(1-I_vac_0_flu(0:140)*0.1)+0.0015))

plot_grid(plot_grid(p1,NULL,NULL,nrow=1,rel_widths = c(1,1,0.5)),
          plot_grid(NULL,p2,NULL,nrow=1,rel_widths = c(0.5,1,1)),
          plot_grid(p3,NULL,p4,nrow=1,rel_widths=c(1,0.5,1)),
          nrow=3)

ggsave(here::here("Figures", "fig1c.png"), bg="white", height=5.2, width = 8)
