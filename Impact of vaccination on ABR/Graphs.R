library(scales)  

# Code avec les fonctions pour afficher les graphiques de sortie des modèles + heatmap

graph<- function(data,filter_values,title){
  color_values <- c("Sa" = "#499124", "CRa" = "#DE6F00", "CSa" = "#2072BA", "IRa" = "#BD5E00", "ISa" = "#163F9E", 
                    "S" = "#68CF33", "CR"="#FC7E00", "CS"="#2B9CFF")
  

  
  if(!is.null(filter_values))
  {
    
    p<-data %>%
      melt(id = "time") %>%
      filter(variable %in% filter_values) %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12,hjust = 0.5)) +
      labs(x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")+
      scale_color_manual(values = color_values)+
      ylim(0,0.8)
    
  }
  
  else{
    p<-data %>%
      melt(id = "time") %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12,hjust = 0.5)) +
      labs(x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")+
      scale_color_manual(values = color_values)+
      ylim(0,0.8)
    
    
  }
  
  return(p)
}
graph2<- function(data,filter_values,title){

  
  if(!is.null(filter_values))
  {
    p<-data %>%
      melt(id = "time") %>%
      filter(variable %in% filter_values) %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12,hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")
  }
  else{
    p<-data %>%
      melt(id = "time") %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")
  }
  
  return(p)
}


graph3<- function(data,filter_values,title){

  
  if(!is.null(filter_values))
  {
    p<-data %>%
      melt(id = "time") %>%
      filter(variable %in% filter_values) %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "value", colour = "Population:")+
      scale_y_continuous(breaks = function(x) pretty(x, n = 10))    
    
  }
  else{
    p<-data %>%
      melt(id = "time") %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")+
      scale_y_continuous(breaks = function(x) pretty(x, n = 10))    
    
  }
  
  return(p)
}

graph4<- function(data,filter_values,title){
  
  if(!is.null(filter_values))
  {
    p<-data %>%
      melt(id = "time") %>%
      filter(variable %in% filter_values) %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12,hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of newly infected individuals", colour = "Vaccination coverage:")+
      scale_color_discrete(labels = c("no_vaccination" = "0%", "vaccination_50" = "50%", "vaccination_80" = "80%"))
  }
  else{
    p<-data %>%
      melt(id = "time") %>%
      ggplot() +
      geom_line(aes(time, value, colour = variable), linewidth = 0.8) +
      theme_bw() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of newly infected individuals", colour = "Vaccination coverage:")+
      scale_color_discrete(labels = c("no_vaccination" = "0%", "vaccination_50" = "50%", "vaccination_80" = "80%"))
  }
  
  return(p)
}

#Graph pour afficher l'incidence cumulée et les personnes infectés au pic épidémic 
graph_I_R<-function(data){
  ggplot() +   
    geom_point(data=data, aes(x=vacc,y=max_propI,colour="Infected people at Epidemic peak"))+
    geom_line(data=data, aes(x=vacc,y=max_propI, colour="Infected people at Epidemic peak"))+
    geom_point(data=data, aes(x=vacc,y=last_propR, colour="Cumulative incidence"))+
    geom_line(data=data, aes(x=vacc,y=last_propR, colour="Cumulative incidence"))+
    labs(title = "Infected people at Epidemic peak \nand cumulative incidence according to vaccination", y = "Proportion of Individuals",
         x = "Vaccine coverage",size=6) +
    scale_colour_manual(name = "Legend", values = c("Infected people at Epidemic peak" = "purple", "Cumulative incidence" = "orange")) +
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold",hjust = 0.5))
}

heatmap <- function(data, x_var, y_var, fill_var, x_text = NULL, y_text = NULL, fill_text = NULL, title = NULL, low_col = "#D1C4E9", high_col = "#311B92", values = FALSE, var_text = NULL) {
  graph <- ggplot(data, aes_string(x = x_var, y = y_var, fill = fill_var)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = low_col, high = high_col) +
    labs(title = title,
         x = x_text,
         y = y_text,
         fill = fill_text) +
    theme_minimal()+
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) 
  
  if (values & !is.null(var_text)) {
    graph <- graph + geom_text(aes_string(label = var_text), color = "#FF8C00", fontface="bold", size = 4)
  }
  
  return(graph)
}

#Barplot de l'incidence cumulée des infections en fonction de la couverture vaccinale
graph_barplot<-function(data){
  ggplot(data, aes(fill=Strain, y=Value, x=vacc)) + 
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(name=" ",labels = c("LastIS" = "Cumulative incidence of infection (senstive strain)", 
                                          "LastIR" = "Cumulative incidence of infection (resistant strain)"),
                      values = c("LastIS" = "#163F9E", 
                                 "LastIR" = "#BD5E00")) +
    labs(x = "Vaccine coverage", y = "Cumulative incidence of infection (per 100,000)") +
    theme_minimal()+
    theme(axis.text = element_text(size = 12),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) 
    
}

# Graphique pour avoir les sorties de modèles en fonction des trois couvertures vaccinales
color_values <- c("Sa" = "#499124", "CRa" = "#DE6F00", "CSa" = "#2072BA", "IRa" = "#BD5E00", "ISa" = "#163F9E", 
                  "S" = "#68CF33", "CR" = "#FC7E00", "CS" = "#2B9CFF")

all_graph <-function(data,title){
  ggplot(data, aes(x = time, y = value, colour = variable)) +
    geom_line(linewidth = 0.8) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+
    labs(title = title, x = "Time (days)", y = "Proportion of Individuals", colour = "Population:") +
    facet_grid(. ~ vaccination) + 
    scale_color_manual(values = color_values)+
    theme_bw() 
  
}

# Fonction annexe pour avoir le pourcentage d'individus dans les compartiments
percentage_final<-function(data){
  if(!("S" %in% colnames(data)) || !("Sa" %in% colnames(data))){
    
    val_CR=tail(data$CRa, n = 1)*100+tail(data$CR, n = 1)*100
    val_CS=tail(data$CSa, n = 1)*100+tail(data$CS, n = 1)*100
    val_exp=tail(data$CRa, n = 1)*100+tail(data$CSa, n = 1)*100+tail(data$ISa, n = 1)*100+tail(data$IRa, n = 1)*100
    val_non_exp=tail(data$CR, n = 1)*100+tail(data$CS, n = 1)*100
    tab<-data.frame(compart=c("CR","CS","exp","non_exp"), Value=c(val_CR,val_CS,val_exp,val_non_exp))
  }
  else
  {
    val_S=tail(data$Sa, n = 1)*100+tail(data$S, n = 1)*100
    val_CR=tail(data$CRa, n = 1)*100+tail(data$CR, n = 1)*100
    val_CS=tail(data$CSa, n = 1)*100+tail(data$CS, n = 1)*100
    val_exp=tail(data$Sa, n = 1)*100+tail(data$CRa, n = 1)*100+tail(data$CSa, n = 1)*100+tail(data$ISa, n = 1)*100+tail(data$IRa, n = 1)*100
    val_non_exp=tail(data$S, n = 1)*100+tail(data$CR, n = 1)*100+tail(data$CS, n = 1)*100
    tab<-data.frame(compart=c("S","CR","CS","exp","non_exp"), Value=c(val_S,val_CR,val_CS,val_exp,val_non_exp))
    
  }
 
  return(tab)
}


diff_graph<- function(data,data2,data3,data4){
  
  if(is.null(data2) || is.null(data3) || is.null(data4)){
    ggplot() +   
      geom_point(data=data, aes(x=vacc,y=diffIR,colour="resistant strain"))+
      geom_line(data=data, aes(x=vacc,y=diffIR,colour="resistant strain"))+
      geom_point(data=data, aes(x=vacc,y=diffIS, colour="sensitive strain"))+
      geom_line(data=data, aes(x=vacc,y=diffIS, colour="sensitive strain"))+
      geom_point(data=data, aes(x=vacc,y=diffISIR, colour="sensitive strain + resistant strain"))+
      geom_line(data=data, aes(x=vacc,y=diffISIR, colour="sensitive strain + resistant strain"))+
      geom_hline(yintercept=0,linetype="dashed",alpha=0.5)+
      labs(y = "Relative difference in cumulative infection",
           x = "Vaccine coverage",size=6) +
      scale_colour_manual(name = "Difference in cumulative infection for:", values = c("resistant strain" = "#BD5E00", "sensitive strain" = "#163F9E",
                                                      "sensitive strain + resistant strain"="#4B0082")) +
      theme_bw()+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5))
  }
  else{
    ggplot() +   
      geom_point(data=data, aes(x=vacc,y=diffIR,colour="resistant strain"))+
      geom_line(data=data, aes(x=vacc,y=diffIR,colour="resistant strain"))+
      geom_point(data=data, aes(x=vacc,y=diffIS, colour="sensitive strain"))+
      geom_line(data=data, aes(x=vacc,y=diffIS, colour="sensitive strain"))+
      geom_point(data=data, aes(x=vacc,y=diffISIR, colour="sensitive strain + resistant strain"))+
      geom_line(data=data, aes(x=vacc,y=diffISIR, colour="sensitive strain + resistant strain"))+
      geom_point(data=data2, aes(x=vacc,y=incidence, colour="resistant strain"), size=0.3)+
      geom_point(data=data3, aes(x=vacc,y=incidence, colour="sensitive strain"), size=0.3)+
      geom_point(data=data4, aes(x=vacc,y=incidence, colour="sensitive strain + resistant strain"), size=0.3)+
      geom_hline(yintercept=0,linetype="dashed",alpha=0.5)+
      labs(y = "Relative difference in cumulative infection",
           x = "Vaccine coverage",size=6) +
      scale_colour_manual(name = "Difference in cumulative infection for:", values = c("resistant strain" = "#BD5E00", "sensitive strain" = "#163F9E",
                                                      "sensitive strain + resistant strain"="#4B0082")) +
      theme_bw()+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5))
  }
  
}

diff_graph_sim<- function(data){
  ggplot() +   
    geom_point(data=data, aes(x=vacc,y=mean_incidence_IR,colour="resistant strain"))+
    geom_line(data=data, aes(x=vacc,y=mean_incidence_IR,colour="resistant strain"))+
    geom_errorbar(data=data, aes(x=vacc,ymin=IR_ic_l,ymax=IR_ic_u,colour="resistant strain"),width=0.05)+
    geom_point(data=data, aes(x=vacc,y=mean_incidence_IS, colour="sensitive strain"))+
    geom_line(data=data, aes(x=vacc,y=mean_incidence_IS, colour="sensitive strain"))+
    geom_errorbar(data=data, aes(x=vacc,ymin=IS_ic_l,ymax=IS_ic_u,colour="sensitive strain"),width=0.05)+
    geom_point(data=data, aes(x=vacc,y=mean_incidence_ISR, colour="sensitive strain + resistant strain"))+
    geom_line(data=data, aes(x=vacc,y=mean_incidence_ISR, colour="sensitive strain + resistant strain"))+
    geom_errorbar(data=data, aes(x=vacc,ymin=ISR_ic_l,ymax=ISR_ic_u,colour="sensitive strain + resistant strain"),width=0.05)+
    geom_hline(yintercept=0,linetype="dashed",alpha=0.5)+
    labs(y = "Relative difference in cumulative infection",
         x = "Vaccine coverage",size=6) +
    scale_colour_manual(name = "Difference in cumulative infection for:", values = c("resistant strain" = "#BD5E00", "sensitive strain" = "#163F9E",
                                                    "sensitive strain + resistant strain"="#4B0082")) +
    theme_bw()+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold",hjust = 0.5))
  
}

graph_exp<-function(data,data2){
  if(!is.null(data2)){
    ggplot() +   
      geom_point(data=data, aes(x=vacc,y=diffexp,colour="exposure to antibiotics"))+
      geom_line(data=data, aes(x=vacc,y=diffexp,colour="exposure to antibiotics"))+
      geom_point(data=data2, aes(x=vacc,y=diffIR,colour="infection of resistant strain"))+
      geom_line(data=data2, aes(x=vacc,y=diffIR,colour="infection of resistant strain"))+
      labs(y = "Relative difference in cumulative incidence",
           x = "Vaccine coverage",size=6) +
      scale_colour_manual(name = "Difference in cumulative incidence of:", values = c("exposure to antibiotics" = "#FFDB58","infection of resistant strain"= "#BD5E00")) +
      theme_bw()+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5))
    
  }
  else{
    ggplot() +   
      geom_point(data=data, aes(x=vacc,y=diffexp,colour="exposure to antibiotics"))+
      geom_line(data=data, aes(x=vacc,y=diffexp,colour="exposure to antibiotics"))+
      labs(y = "Relative difference in cumulative incidence of exposure to antibiotics",
           x = "Vaccine coverage",size=6) +
      scale_colour_manual(name = "Difference in cumulative infection for:", values = c("exposure to antibiotics" = "#FFDB58")) +
      theme_bw()+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5))
  }
  
}
density_graph<- function(data){
  ggplot() +   
    geom_density(data=data, aes(x=incidenceR,y = ..count../sum(..count..),colour="Incidence of infection (resistant strain)",fill="Incidence of infection (resistant strain)"),alpha=0.5)+
    geom_density(data=data, aes(x=incidenceS,y = ..count../sum(..count..),colour="Incidence of infection (sensitive strain)",fill="Incidence of infection (sensitive strain)"),alpha=0.5)+
    #geom_density(data=data, aes(x=incidenceS+incidenceR, y = ..count../sum(..count..),colour="Incidence of infection",fill="Incidence of infection"),alpha=0.5)+
    labs(title = "Density of incidence of infection", y = "Density",
         x = "Incidence",size=6) +
    scale_colour_manual(name = "Legend", values = c("Incidence of infection (resistant strain)" = "#BD5E00", "Incidence of infection (sensitive strain)" = "#163F9E","Incidence of infection"="#4B0082")) +
    scale_fill_manual(name = "Legend", values = c("Incidence of infection (resistant strain)" = "#BD5E00", "Incidence of infection (sensitive strain)" = "#163F9E","Incidence of infection"="#4B0082"))+ 
    theme_bw()+
    theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold",hjust = 0.5))
  
}

#Corrélation partielle
graph_pcor<-function(data){
  ggplot(data) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 1, linetype = "solid") +
    geom_hline(yintercept = -1, linetype = "solid") +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_hline(yintercept = -0.5, linetype = "dashed") +
    geom_pointrange(aes(x = param, y = est, ymin = lower, ymax = upper), size = 0.8) +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=8)) +
    labs(colour = "", x = "Parameters", y = "Correlation coefficient with cumulative incidence")+
    scale_x_discrete(labels = c(bquote(alpha),
                                bquote(ATB),
                                bquote(beta),
                                bquote(fcost),
                                bquote(gamma),
                                bquote(omega),
                                bquote(theta)))
}

