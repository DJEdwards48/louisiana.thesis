#####NMDS for publication - Louisiana#####

##################
# Clear any variables from the R environment
##################

rm(list=ls())

###Loading Libraries
library(vegan)
library(ggplot2)
library(ggpubr)
library(rysgran)    ########Shepard Diagram
library(dplyr)
library(MASS)
library(RColorBrewer) ###### for coloring graph
library(tidyverse)    ##for data wrangling
library(ggrepel) #Label making in plot to prevent over lap

##################
# Set WD()
##################

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory
getwd()

#LoadingData
data.taxa <- read.csv("data.taxa_1.0.csv", header = T, sep = ",") #read in stream invertebrate data collated by TA
data.wq <- read.csv(file = "data.wq_1.1_wpoolID.csv", header = TRUE, stringsAsFactors = TRUE)
data.div <- read.csv(file = "RareAdultColonist2022_finalRao.csv", header = TRUE, stringsAsFactors = TRUE)
data.undermesh.taxa <- read.csv(file = "Undermesh.nmds.data.csv", header = TRUE, stringsAsFactors = TRUE)


####Prepping dfs
##turning NAs into ZEROs
data.taxa <- data.taxa %>% replace(is.na(.), 0)
summary(data.taxa)
#Filtering out control pools
data.undermesh.taxa <- data.undermesh.taxa %>% 
  replace(is.na(.), 0) %>%
  filter(!(Pool. == "7"))%>%
  filter(!(Pool. == "19"))%>%
  filter(!(Pool. == "31"))%>%
  filter(!(Pool. == "35"))%>%
  filter(!(Pool. == "45"))


##Selecting only columns of interest in WQ
data.wq.all <- data.wq %>%
  dplyr::select(Rao,
         FDiv,
         Richness,
         Water.Tannin,
         DO, 
         pH,
         Turbidity, 
         Temp, 
         Conductivity,
         Water.Color,
         Leaf_tannin,
         LMA,
         CN,
         CP,
         N.P)


data.wq.all2 <- data.wq.all %>%
  mutate(Rao_z = (Rao - mean(Rao, na.rm = TRUE))/sd(Rao, na.rm = TRUE))%>%  
  mutate(FDiv_z = (FDiv - mean(FDiv, na.rm = TRUE))/sd(FDiv, na.rm = TRUE))%>%  
  mutate(Richness_z = (Richness - mean(Richness, na.rm = TRUE))/sd(Richness, na.rm = TRUE))%>%  
  mutate(Water.Tannin_z = (Water.Tannin - mean(Water.Tannin, na.rm = TRUE))/sd(Water.Tannin, na.rm = TRUE))%>%  
  mutate(DO_z = (DO - mean(DO, na.rm = TRUE))/sd(DO, na.rm = TRUE))%>%  
  mutate(pH_z = (pH - mean(pH, na.rm = TRUE))/sd(pH, na.rm = TRUE))%>%  
  mutate(Turbidity_z = (Turbidity - mean(Turbidity, na.rm = TRUE))/sd(Turbidity, na.rm = TRUE))%>%  
  mutate(Temp_z = (Temp - mean(Temp, na.rm = TRUE))/sd(Temp, na.rm = TRUE))%>%  
  mutate(Conductivity_z = (Conductivity - mean(Conductivity, na.rm = TRUE))/sd(Conductivity, na.rm = TRUE))%>%  
  mutate(Water.Color_z = (Water.Color - mean(Water.Color, na.rm = TRUE))/sd(Water.Color, na.rm = TRUE))%>%  
  mutate(Leaf_tannin_z = (Leaf_tannin - mean(Leaf_tannin, na.rm = TRUE))/sd(Leaf_tannin, na.rm = TRUE))%>%  
  mutate(LMA_z = (LMA - mean(LMA, na.rm = TRUE))/sd(LMA, na.rm = TRUE))%>%  
  mutate(CN_z = (CN - mean(CN, na.rm = TRUE))/sd(CN, na.rm = TRUE))%>%  
  mutate(CP_z = (CP - mean(CP, na.rm = TRUE))/sd(CP, na.rm = TRUE))%>%  
  mutate(N.P_z = (N.P - mean(N.P, na.rm = TRUE))/sd(N.P, na.rm = TRUE))

data.wq.all2 <- data.wq.all2 %>%
  dplyr::select(Rao_z,
                FDiv_z,
                Richness_z,
                Water.Tannin_z,
                DO_z, 
                pH_z,
                Turbidity_z, 
                Temp_z, 
                Conductivity_z,
                Water.Color_z,
                Leaf_tannin_z,
                LMA_z,
                CN_z,
                CP_z,
                N.P_z)

####Runing NMDS, dimensions is set with "k" set to 3 to get stress below 0.2
nmds1<-metaMDS(data.taxa, k=2, autotransform = F, trymax = 1000)
nmds1
stressplot(nmds1)

#################EnvironmentalFIt###########################

###Which Sp. drive site distribution
nmds.spp.fit <- envfit(nmds1, data.taxa, permutations = 999)
head(nmds.spp.fit)

#setting up WQ fit
nmds.env <- envfit(nmds1, data.wq.all2, permutations = 999)
head(nmds.env)

###Plotting nMDS
ordiplot(nmds1, type = "n", main = "nMDS of Pool Communites", xlim = c(-1.25,1), ylim = c(-1,1))
orditorp(nmds1, display = "sites", labels = T, cex = 1)
plot(nmds.spp.fit, p.max = 0.01, col = "Red", cex = 0.8) # change the significance level of species shown with p.max
plot(nmds.env, col = "Blue", p.max = 0.05, cex = 0.7)


# Extract NMDS site coordinates and prepare them for ggplot
nmds_points <- as.data.frame(scores(nmds1, display = "sites"))
nmds_points$Site <- rownames(nmds_points)

# Convert fitted species and environmental vectors to data frames
spp_scores <- as.data.frame(scores(nmds.spp.fit, display = "vectors"))
spp_scores$Species <- rownames(spp_scores)
env_scores <- as.data.frame(scores(nmds.env, display = "vectors"))
env_scores$Variable <- rownames(env_scores)

# Applying filters for p values
spp_scores <- subset(spp_scores, nmds.spp.fit$vectors$pvals < 0.01)
env_scores <- subset(env_scores, nmds.env$vectors$pvals < 0.05)

library(ggrepel)

# Plot of colonizers (weekly)
GGnmds <- ggplot(data = nmds_points, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color = "Black") +
  geom_segment(data = spp_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  geom_text_repel(data = spp_scores, aes(x = NMDS1, y = NMDS2, label = Species),
            color = "red", size = 4, hjust = 1.1, box.padding = .5) +
  geom_segment(data = env_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue") +
  geom_text_repel(data = env_scores, aes(x = NMDS1, y = NMDS2, label = Variable),
            color = "blue", size = 4, hjust = -0.1, box.padding = .5) +
  labs(title = "nMDS of Colonizing Communities") +
  xlim(-1.25, 1) + ylim(-1, 1) +
  theme_classic()
GGnmds

ggsave("Colonizer.nmds.png", GGnmds, #Save the file to the current working directory, specifying file type (png) within the quotation marks
       width = 7, #Set plot width
       height = 7, #Set plot height
       units = "in") #Set plot width/height units

#############NMDS of Undermesh sample

####Runing NMDS, dimensions is set with "k" set to 3 to get stress below 0.2
####Transform taxa using hellinger transformation
taxa.hell<-decostand(data.undermesh.taxa[,-1], method = "hellinger")


nmds2<-metaMDS(taxa.hell, k=2, autotransform = F, trymax = 1000)
nmds2
stressplot(nmds2)

###Which Sp. drive site distribution
nmds.spp.fit <- envfit(nmds2, taxa.hell, permutations = 999)
head(nmds.spp.fit)

#setting up WQ fit
nmds.env <- envfit(nmds2, data.wq.all, permutations = 999)
head(nmds.env)

###Plotting nMDS
ordiplot(nmds2, type = "n", main = "nMDS of Undermesh Communites", xlim = c(-1.25,1), ylim = c(-1,1))
orditorp(nmds2, display = "sites", labels = T, cex = 1)
plot(nmds.spp.fit, p.max = 0.01, col = "Red", cex = 0.8) # change the significance level of species shown with p.max
plot(nmds.env, col = "Blue", p.max = 0.05, cex = 0.7)


# Extract NMDS site coordinates and prepare them for ggplot
nmds_points <- as.data.frame(scores(nmds2, display = "sites"))
nmds_points$Site <- rownames(nmds_points)

# Convert fitted species and environmental vectors to data frames
spp_scores <- as.data.frame(scores(nmds.spp.fit, display = "vectors"))
spp_scores$Species <- rownames(spp_scores)
env_scores <- as.data.frame(scores(nmds.env, display = "vectors"))
env_scores$Variable <- rownames(env_scores)

# Applying filters for p values
spp_scores <- subset(spp_scores, nmds.spp.fit$vectors$pvals < 0.01)
env_scores <- subset(env_scores, nmds.env$vectors$pvals < 0.05)




# Plot of Undermesh samples
ggplot(data = nmds_points, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color = "Black") +
  geom_segment(data = spp_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  geom_text_repel(data = spp_scores, aes(x = NMDS1, y = NMDS2, label = Species),
                  color = "red", size = 4, hjust = 1.1, box.padding = .5) +
  geom_segment(data = env_scores, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue") +
  geom_text_repel(data = env_scores, aes(x = NMDS1, y = NMDS2, label = Variable),
                  color = "blue", size = 4, hjust = -0.1, box.padding = .5) +
  labs(title = "nMDS of undermesh Communities") +
  xlim(-1.25, 1) + ylim(-1, 1) +
  theme_classic()

ggsave("Undermesh.nmds.png", GGnmds, #Save the file to the current working directory, specifying file type (png) within the quotation marks
       width = 7, #Set plot width
       height = 7, #Set plot height
       units = "in") #Set plot width/height units


