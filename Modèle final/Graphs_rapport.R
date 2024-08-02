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


h1_P2<-ggdraw() +
  draw_plot(h1_P) +  
  draw_plot_label("A")
h1_A2<-ggdraw() +
  draw_plot(h1_A) +  
  draw_plot_label("B")

plot_grid(h1_P2,h1_A2,ncol=1,rel_widths = c(1, 1, 0.5))

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

legend8 <- get_legend(pcorR_P)
pcorR_P <- pcorR_P + theme(legend.position = "none")
pcorR_P2<-ggdraw() +
  draw_plot(pcorR_P) +  
  draw_plot_label("A")
pcorR_A<- pcorR_A + theme(legend.position = "none")
pcorR_A2<-ggdraw() +
  draw_plot(pcorR_A) +  
  draw_plot_label("B")
pcorR_E<- pcorR_E+ theme(legend.position = "none")
pcorR_E2<-ggdraw() +
  draw_plot(pcorR_E) +  
  draw_plot_label("C")

plot_grid(pcorR_P2,pcorR_A2,pcorR_E2,legend8,ncol=2,rel_widths = c(1, 1, 0.3))

legend9 <- get_legend(pcorS_P)
pcorS_P <- pcorS_P + theme(legend.position = "none")
pcorS_P2<-ggdraw() +
  draw_plot(pcorS_P) +  
  draw_plot_label("A")
pcorS_A<- pcorS_A + theme(legend.position = "none")
pcorS_A2<-ggdraw() +
  draw_plot(pcorS_A) +  
  draw_plot_label("B")
pcorS_E<- pcorS_E+ theme(legend.position = "none")
pcorS_E2<-ggdraw() +
  draw_plot(pcorS_E) +  
  draw_plot_label("C")

plot_grid(pcorS_P2,pcorS_A2,pcorS_E2,legend9,ncol=2,rel_widths = c(1, 1, 0.3))


legend10 <- get_legend(pcorSR_P)
pcorSR_P <- pcorSR_P + theme(legend.position = "none")
pcorSR_P2<-ggdraw() +
  draw_plot(pcorSR_P) +  
  draw_plot_label("A")
pcorSR_A<- pcorSR_A + theme(legend.position = "none")
pcorSR_A2<-ggdraw() +
  draw_plot(pcorSR_A) +  
  draw_plot_label("B")
pcorSR_E<- pcorSR_E+ theme(legend.position = "none")
pcorSR_E2<-ggdraw() +
  draw_plot(pcorSR_E) +  
  draw_plot_label("C")

plot_grid(pcorSR_P2,pcorSR_A2,pcorSR_E2,legend10,ncol=2,rel_widths = c(1, 1, 0.3))


legend11 <- get_legend(psa_graph_P)
psa_graph_P <- psa_graph_P + theme(legend.position = "none")
psa_graph_P2<-ggdraw() +
  draw_plot(psa_graph_P) +  
  draw_plot_label("A")
psa_graph_A<- psa_graph_A + theme(legend.position = "none")
psa_graph_A2<-ggdraw() +
  draw_plot(psa_graph_A) +  
  draw_plot_label("B")
psa_graph_E<- psa_graph_E+ theme(legend.position = "none")
psa_graph_E2<-ggdraw() +
  draw_plot(psa_graph_E) +  
  draw_plot_label("C")

plot_grid(psa_graph_P2,psa_graph_A2,psa_graph_E2,legend11,ncol=2,rel_widths = c(1, 1, 0.3))


legend12 <- get_legend(psa_graph_P)
graph_diff_expP <- graph_diff_expP + theme(legend.position = "none")
graph_diff_expP2<-ggdraw() +
  draw_plot(graph_diff_expP) +  
  draw_plot_label("A")
graph_diff_expA<- graph_diff_expA + theme(legend.position = "none")
graph_diff_expA2<-ggdraw() +
  draw_plot(graph_diff_expA) +  
  draw_plot_label("B")
graph_diff_expE<- graph_diff_expE+ theme(legend.position = "none")
graph_diff_expE2<-ggdraw() +
  draw_plot(graph_diff_expE) +  
  draw_plot_label("C")

plot_grid(graph_diff_expP2,graph_diff_expA2,graph_diff_expE2,legend12,ncol=2,rel_widths = c(1, 1, 0.3))
