legend <- get_legend(run0_A)
run0_P <- run0_P + theme(legend.position = "none")
run0_P2<-ggdraw() +
  draw_plot(run0_P) +  
  draw_plot_label("A")
run0_A<- run0_A + theme(legend.position = "none")
run0_A2<-ggdraw() +
  draw_plot(run0_A) +  
  draw_plot_label("B")
run0_E<- run0_E + theme(legend.position = "none")
run0_E2<-ggdraw() +
  draw_plot(run0_E) +  
  draw_plot_label("C")

plot_grid(run0_P2,run0_A2,run0_E2,legend,ncol=2,rel_widths = c(1, 1, 1, 0.3))

legend2 <- get_legend(I_flu)
I_flu <- I_flu + theme(legend.position = "none")
I_flu2<-ggdraw() +
  draw_plot(I_flu) +  
  draw_plot_label("A")
I_rota<- I_rota + theme(legend.position = "none")
I_rota2<-ggdraw() +
  draw_plot(I_rota) +  
  draw_plot_label("B")

plot_grid(I_flu2,I_rota2,legend2,ncol=3,rel_widths = c(1, 1, 0.3))

legend3 <- get_legend(res_graphs_P)
res_graphs_P <- res_graphs_P + theme(legend.position = "none")
res_graphs_P2<-ggdraw() +
  draw_plot(res_graphs_P) +  
  draw_plot_label("A")
res_graphs_A<- res_graphs_A + theme(legend.position = "none")
res_graphs_A2<-ggdraw() +
  draw_plot(res_graphs_A) +  
  draw_plot_label("B")
res_graphs_E<- res_graphs_E+ theme(legend.position = "none")
res_graphs_E2<-ggdraw() +
  draw_plot(res_graphs_E) +  
  draw_plot_label("C")

plot_grid(res_graphs_P2,res_graphs_A2,res_graphs_E2,legend3,ncol=2,rel_widths = c(1, 1, 0.3))

legend4 <- get_legend(Cumulative_incidence_P)
Cumulative_incidence_P <- Cumulative_incidence_P + theme(legend.position = "none")
Cumulative_incidence_P2<-ggdraw() +
  draw_plot(Cumulative_incidence_P) +  
  draw_plot_label("A")
Cumulative_incidence_A<- Cumulative_incidence_A + theme(legend.position = "none")
Cumulative_incidence_A2<-ggdraw() +
  draw_plot(Cumulative_incidence_A) +  
  draw_plot_label("B")
Cumulative_incidence_E<- Cumulative_incidence_E+ theme(legend.position = "none")
Cumulative_incidence_E2<-ggdraw() +
  draw_plot(Cumulative_incidence_E) +  
  draw_plot_label("C")

plot_grid(Cumulative_incidence_P2,Cumulative_incidence_A2,Cumulative_incidence_E2,legend4,ncol=2,rel_widths = c(1, 1, 0.3))

legend5<- get_legend(h1_A)
h1_P <- h1_P + theme(legend.position = "none")
h1_P2<-ggdraw() +
  draw_plot(h1_P) +  
  draw_plot_label("A")
h1_A<- h1_A + theme(legend.position = "none")
h1_A2<-ggdraw() +
  draw_plot(h1_A) +  
  draw_plot_label("B")

plot_grid(h1_P2,h1_A2,legend5,ncol=3,rel_widths = c(1, 1, 0.5))

legend6 <- get_legend(diff_P)
diff_P <- diff_P + theme(legend.position = "none")
diff_P2<-ggdraw() +
  draw_plot(diff_P) +  
  draw_plot_label("A")
diff_A<- diff_A + theme(legend.position = "none")
diff_A2<-ggdraw() +
  draw_plot(diff_A) +  
  draw_plot_label("B")
diff_E<- diff_E+ theme(legend.position = "none")
diff_E2<-ggdraw() +
  draw_plot(diff_E) +  
  draw_plot_label("C")

plot_grid(diff_P2,diff_A2,diff_E2,legend6,ncol=2,rel_widths = c(1, 1, 0.3))

legend7 <- get_legend(diff_sim_P)
diff_sim_P <- diff_sim_P + theme(legend.position = "none")
diff_sim_P2<-ggdraw() +
  draw_plot(diff_sim_P) +  
  draw_plot_label("A")
diff_sim_A<- diff_sim_A + theme(legend.position = "none")
diff_sim_A2<-ggdraw() +
  draw_plot(diff_sim_A) +  
  draw_plot_label("B")
diff_sim_E<- diff_sim_E+ theme(legend.position = "none")
diff_sim_E2<-ggdraw() +
  draw_plot(diff_sim_E) +  
  draw_plot_label("C")

plot_grid(diff_sim_P2,diff_sim_A2,diff_sim_E2,legend7,ncol=2,rel_widths = c(1, 1, 0.3))
