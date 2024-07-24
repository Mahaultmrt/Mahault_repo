legend <- get_legend(run0_A)
run0_P <- run0_P + theme(legend.position = "none")
run0_P2<-ggdraw() +
  draw_plot(run0_P, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 1.A: S.Pneumoniae Colonization dynamics (without virus epidemic)", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)
run0_A<- run0_A + theme(legend.position = "none")
run0_A2<-ggdraw() +
  draw_plot(run0_A, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 1.B: S.Aureus Colonization dynamics (without virus epidemic)", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)
run0_E<- run0_E + theme(legend.position = "none")
run0_E2<-ggdraw() +
  draw_plot(run0_E, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 1.C: E.Coli Colonization dynamics without virus epidemics", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)



plot_grid(run0_P2,run0_A2,run0_E2,legend,ncol=2,rel_widths = c(1, 1, 1, 0.3))

legend2 <- get_legend(I_flu)
I_flu <- I_flu + theme(legend.position = "none")
I_flu2<-ggdraw() +
  draw_plot(I_flu, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 2.A: Incidence of infected people by influenza", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)
I_rota<- I_rota + theme(legend.position = "none")
I_rota2<-ggdraw() +
  draw_plot(I_rota, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 2.B: Incidence of infected people by rotavirus", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)

plot_grid(I_flu2,I_rota2,legend2,ncol=3,rel_widths = c(1, 1, 0.3))

legend3 <- get_legend(res_graphs_P)
res_graphs_P <- res_graphs_P + theme(legend.position = "none")
res_graphs_P2<-ggdraw() +
  draw_plot(res_graphs_P, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 3.A: SCI model of S.Pneumoniae depending on vaccine coverage", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)
res_graphs_A<- res_graphs_A + theme(legend.position = "none")
res_graphs_A2<-ggdraw() +
  draw_plot(res_graphs_A, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 3.B: SCI model of S.Aureus depending on vaccine coverage", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)
res_graphs_E<- res_graphs_E+ theme(legend.position = "none")
res_graphs_E2<-ggdraw() +
  draw_plot(res_graphs_E, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 3.C: SCI model of E.Coli depending on vaccine coverage", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)

plot_grid(res_graphs_P2,res_graphs_A2,res_graphs_E2,legend3,ncol=2,rel_widths = c(1, 1, 0.3))

legend4 <- get_legend(Cumulative_incidence_P)
Cumulative_incidence_P <- Cumulative_incidence_P + theme(legend.position = "none")
Cumulative_incidence_P2<-ggdraw() +
  draw_plot(Cumulative_incidence_P, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 4.A: S. Pneumonia cumulative incidence of infection depending on vaccine coverage", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)
Cumulative_incidence_A<- Cumulative_incidence_A + theme(legend.position = "none")
Cumulative_incidence_A2<-ggdraw() +
  draw_plot(Cumulative_incidence_A, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 4.B: S.Aureus cumulative incidence of infection depending on vaccine coverage", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)
Cumulative_incidence_E<- Cumulative_incidence_E+ theme(legend.position = "none")
Cumulative_incidence_E2<-ggdraw() +
  draw_plot(Cumulative_incidence_E, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 4.C: E.Coli cumulative incidence of infection depending on vaccine coverage", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)

plot_grid(Cumulative_incidence_P2,Cumulative_incidence_A2,Cumulative_incidence_E2,legend4,ncol=2,rel_widths = c(1, 1, 0.3))

legend5<- get_legend(h1_P)
h1_P <- h1_P + theme(legend.position = "none")
h1_P2<-ggdraw() +
  draw_plot(h1_P, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 5.A: S. Pneumonia cumulative incidence of infection depending on vaccine coverage and antibiotics", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)
h1_A<- h1_A + theme(legend.position = "none")
h1_A2<-ggdraw() +
  draw_plot(h1_A, 0, 0.05, 1, 0.95) +  
  draw_label("Figure 5.B: S.Aureus cumulative incidence of infection depending on vaccine coverage and antibiotics", fontface = 'italic', x = 0.5, y = 0.02, hjust = 0.5, size = 10)

plot_grid(h1_P2,h1_A2,legend5,ncol=3,rel_widths = c(1, 1, 0.5))

