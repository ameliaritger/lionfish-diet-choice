setwd("~/Desktop")
data<-read.csv('Curacao fish surveys.csv')
colnames(data)[colnames(data)=="ï..Species"] <- "Species"
number_ticks <- function(n) {function(limits) pretty(limits, n)}
library(ggplot2)

#density
ggplot(data,aes(x = Location, y = Density, fill = Species)) + 
  geom_bar(position = "fill",stat = "identity") +
  ylab(bquote("Density" ~ (fish ~m^-2)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black")) +       
  scale_y_continuous(expand = c(0, 0), breaks = number_ticks(10), labels = scales::percent) +
  theme(axis.title.x=element_blank(), axis.title.y = element_text(size = 14), axis.ticks.x=element_blank(), legend.title = element_blank(), text = element_text(size=15))

#biomass
ggplot(data,aes(x = Location, y = Biomass, fill = Species)) + 
  geom_bar(position = "fill",stat = "identity") +
  ylab(bquote("Biomass" ~ (g ~m^-2)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +       
  scale_y_continuous(expand = c(0, 0), breaks = number_ticks(10), labels = scales::percent) +
  theme(axis.title.x=element_blank(), axis.title.y = element_text(size = 14), axis.ticks.x=element_blank(), legend.title = element_blank(), text = element_text(size=15))
