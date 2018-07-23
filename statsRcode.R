##############################################################################
# 
# ELISA G. DE LOPE for PHARAMACELERA S.L.
# ATOM TYPES - MOPAC PROJECT (MSc Bioinformatics)
# 
# File: statsRcode.R
#
# Created on 26/06/2018
# 
# README: Script meant to generate plots and graphs from statistical data 
# related to electrostatic and cavitational hydrophobic descriptors obtained 
# from the parametrization of LogPelectrostatic and LogP cavitational with 
# Pharmscreen AT version (All rights reserved) in a training set (library of compounds).
#
# ----------------------------------------------------------------------------
# Set working directory, load libraries
setwd("~/AT_MOPAC_PY/dataset_results/finalsdf_parametrization/")
library("ggplot2")
library(lattice)
theme_set(theme_bw())

#-----------------------------------------------------------------------------
# Load LogP electrostatic data from file
mystats<-read.table(file = "Stats_parametrizationele.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str(mystats)

# sort data according to Standard Deviation
mystats <- mystats[order(mystats$SDev), ]  
mystats$AtomType <- factor(mystats$AtomType, levels=mystats$AtomType)

#-----------------------------------------------------------------------------
# Plot data: Mean, histogram-mode ('mode'), Standard Deviation, Outliers %
s<- ggplot(mystats, aes(mystats$AtomType, group = 1))+
  geom_line(aes(y = mystats$Mean, colour = "Mean"), size =0.7) +
  geom_line(aes(y=mystats$Mode, colour = "Parametrized value"), size = 0.8)+
  geom_line(aes(y=mystats$SDev, colour = "Standard Deviation"), alpha="0.3")+
  geom_line(aes(y=mystats$outliers_per100, colour = "Outliers (%)"), alpha="0.3")+
  geom_ribbon(aes(ymin=0,ymax=mystats$SDev), fill="yellow", alpha="0.7")+
  geom_ribbon(aes(ymin=0,ymax=mystats$outliers_per100), fill="gray65", alpha="0.3")+
  theme(axis.title=element_text(size=14), axis.text.y = element_text(size=11), axis.text.x = element_text(size=11, angle = 90),legend.position = c(0.15,0.2), legend.title = element_text(size=13), legend.text = element_text(size=12),plot.title = element_text(size=18),plot.subtitle = element_text(size=14))+
  labs(title = "Atomic logP electrostatic in training set", subtitle = "LogP electrostatic contribution per atom-type", x = "Atom-types", y= "")+
  scale_y_continuous(breaks=seq(-4,1,1))+
  scale_color_manual(name = "Parametrization statistics", 
                     values = c('Mean' = 'coral1',  'Parametrized value' = 'steelblue3', 'Standard Deviation' = 'goldenrod1', 'Outliers (%)' = 'gray65')) 
s

# Export to pdf 
pdf('Stats_ele.pdf', width=12, height=7)
s
dev.off()

# Generate occurrence graph
occurrence <- ggplot(mystats, aes(mystats$AtomType, group = 1))+
  geom_line(aes(y=mystats$Total_ats), colour= "skyblue")+
  theme(axis.title=element_text(size=14), axis.text.y = element_text(size=11), axis.text.x = element_text(size=11, angle = 90),legend.position = c(0.15,0.2), legend.title = element_text(size=13), legend.text = element_text(size=12),plot.title = element_text(size=18),plot.subtitle = element_text(size=14))+
  labs(title = "Occurrence of atom-types in training set", subtitle = "Number of atoms identified per atom-type", x = "Atom-types", y= "Occurrence")+
  geom_point(data = subset(mystats, Total_ats < 1000), 
           aes(x=data$AtomType, y=data$Total_ats), size = 2, color = "orange") + 
  scale_y_log10() 

occurrence

# Export to pdf and tiff
pdf('Occurrence_atomtypes.pdf', width=12, height=7)
occurrence
dev.off()

tiff('Occurrence_atomtypes.tiff', units="in", width=12, height=6, res=400)
occurrence
dev.off()



#-----------------------------------------------------------------------------
# Load LogP cavitational data from file
mystatscav<-read.table(file = "Stats_parametrizationcav.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str(mystatscav)

# sort data according to Standard Deviation
mystatscav <- mystatscav[order(mystatscav$SDev), ]  
mystatscav$AtomType <- factor(mystatscav$AtomType, levels=mystatscav$AtomType)

# Plot data: Mean, histogram-mode ('mode'), Standard Deviation, Outliers %
r<- ggplot(mystatscav, aes(mystatscav$AtomType, group = 1))+
  geom_line(aes(y=mystatscav$Mode, colour = "Parametrized value"), size = 0.8,linetype = 1)+
  geom_line(aes(y = mystatscav$Mean, colour = "Mean"), size =0.75, linetype = "dashed") +
  geom_line(aes(y=mystatscav$SDev, colour = "Standard Deviation"), alpha = "0.3")+
  geom_line(aes(y=mystatscav$outliers_per100, colour = "Outliers (%)"))+
  geom_ribbon(aes(ymin=0,ymax=mystatscav$SDev), fill="yellow", alpha="0.7")+
  geom_ribbon(aes(ymin=0,ymax=mystatscav$outliers_per100), fill="gray65", alpha="0.3")+
  theme(axis.title=element_text(size=14), axis.text.y = element_text(size=11), axis.text.x = element_text(size=11, angle = 90),legend.position = c(0.15,0.13), legend.title = element_text(size=13), legend.text = element_text(size=12),plot.title = element_text(size=18),plot.subtitle = element_text(size=14))+
  labs(title = "Atomic logP cavitational in training set", subtitle = "LogP cavitational contribution per atom-type", x = "Atom-types", y= "")+
  scale_y_continuous(breaks=seq(-0.4,0.1,0.1))+
  scale_color_manual(name = "Parametrization statistics", 
                     values = c('Mean' = 'coral1',  'Parametrized value' = 'olivedrab3', 'Standard Deviation' = 'goldenrod1', 'Outliers (%)' = 'gray65')) 

r

# Export to pdf
pdf('Stats_cav.pdf', width=12, height=7)
r
dev.off()
