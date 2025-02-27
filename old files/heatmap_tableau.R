library(ggplot2)
tableau<-data.frame(
  Bacteria_virus=c("S.Pneumoniae/Inluenza","S.Aureus/Influenza","E.Coli/Rotavirus"),
  Overall_carrying_frenquency=c(1,2,4),
  Resistance_frequency=c(2,1,2),
  Transmission_rate=c(4,2,1),
  R0=c(1,1,5),
  Vaccine_efficacy=c(2,2,2))

tableau<-melt(tableau,id.vars="Bacteria_virus")

ggplot(tableau, aes(x = variable, y = Bacteria_virus, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#228B22", mid = "#FFDD57", high = "#E31A1C", midpoint = median(tableau$value)) +
  theme_minimal() +
  labs(title = "Heatmap des données bactériennes et virales",
       x = "Catégorie",
       y = "Virus/Bactérie",
       fill = "Valeur")
