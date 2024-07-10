# Code avec les fonctions pour afficher les graphiques de sortie des modèles + heatmap

graph<- function(data,filter_values,title){
  #data_name<-as.character(substitute(data))
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
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")+
      scale_color_manual(values = color_values)
    
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
      scale_color_manual(values = color_values)
    
    
  }
  
  return(p)
}
graph2<- function(data,filter_values,title){
  #data_name<-as.character(substitute(data))
  
  
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
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")
    
    
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
      labs(title=title,x = "Time (days)", y = "Proportion of Individuals", colour = "Population:")
    
    
  }
  
  return(p)
}

heatmap <- function(data, x_var, y_var, fill_var, x_text = NULL, y_text = NULL, fill_text = NULL, title = NULL, low_col = "#B39DDB", high_col = "#4B0082", values = FALSE, var_text = NULL) {
  graph <- ggplot(data, aes_string(x = x_var, y = y_var, fill = fill_var)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = low_col, high = high_col) +
    labs(title = title,
         x = x_text,
         y = y_text,
         fill = fill_text) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
    theme_minimal()
  
  if (values & !is.null(var_text)) {
    graph <- graph + geom_text(aes_string(label = var_text), color = "#FFB347", fontface="bold", size = 3)
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
    labs(title = "Cumulative incidence of infection depending on the vaccine coverage", x = "Vaccine coverage", y = "Cumulative incidence of infection (per 100,000)") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold",hjust = 0.5)) +
    theme_minimal()
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
