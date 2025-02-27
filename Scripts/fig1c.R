
source(here::here("Impact of vaccination on ABR", "SIR_Influenza.R"))

p1=ggplot() +
  geom_line(aes(seq(0:100), I_vac_0(0:100)), linewidth=2, colour="grey70") +
  labs(x = "Time (days)", y = bquote({I^V}(t)~~+~~{I^NV}(t))) +
  scale_x_continuous(breaks = c(0,50,100)) +
  scale_y_continuous(breaks = seq(0,0.03,0.01)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12))

p2=ggplot() +
  geom_line(aes(seq(0:100), I_vac_0(0:100)*0.1), linewidth=2, colour="grey50") +
  labs(x = "Time (days)", y = bquote({P^V}(t))) +
  scale_x_continuous(breaks = c(0,50,100)) +
  scale_y_continuous(breaks = seq(0,0.003,0.001)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12))

p3=ggplot() +
  geom_line(aes(seq(0:100), rep(0.0014,101)), linetype="dashed", linewidth=2, colour="grey50") +
  labs(x = "Time (days)", y = bquote({theta^0}(t))) +
  scale_x_continuous(breaks = c(0,50,100)) +
  scale_y_continuous(breaks = seq(0,0.003,0.001), limits=c(0,0.0026)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12))

p4=ggplot() +
  geom_line(aes(seq(0:100), -log(1-I_vac_0(0:100)*0.1)+0.0014), linewidth=2) +
  labs(x = "Time (days)", y = bquote({theta^V}(t))) +
  scale_x_continuous(breaks = c(0,50,100)) +
  scale_y_continuous(breaks = seq(0,0.003,0.001), limits=c(0,0.004)) +
  theme_bw() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.text = element_text(size=12))

plot_grid(plot_grid(p1,NULL,NULL,nrow=1),
          plot_grid(NULL,p2,NULL,NULL,nrow=1,rel_widths = c(0.8,1.5,1,1)),
          plot_grid(p3,NULL,p4,nrow=1),
          nrow=3)

ggsave(here::here("Figures", "fig1c.png"), bg="white", height=4.65)
