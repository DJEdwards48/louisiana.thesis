###Leaves VS No Leaf Code

## Reading in Data 
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory
getwd()
data <- read.csv("Leaf~NoLeaf.csv",header=TRUE)

###Importing data
data <- `Leaf~NoLeaf`
summary(data)

##Leaf_NoLeaf numerical>>>categorical
data$Leaf_NoLeaf2 <- as.factor(data$Leaf_NoLeaf)
summary(data)

#############Hemiptera###############

###Normality Test
shapiro.test(data$Hemiptera)
hist(data$Hemiptera)
###p-value = 1.427e-06

###Transformation Station
data$tHemiptera <- (log((data$Hemiptera)+1))
shapiro.test(data$tHemiptera)
hist(data$tHemiptera, breaks = 10)
summary(data$tHemiptera)
###P=0.007
###COULD NOT TRANSFORM  

aov(Hemiptera~Leaf_NoLeaf2, data=data)
model <- aov(Hemiptera~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.02480.02

plot(Hemiptera~Leaf_NoLeaf2, data=data)


###Coleoptera##########
###Normality Test
shapiro.test(data$Coleoptera)
hist(data$Coleoptera)
###p-value = 3.816e-06

###Transformation Station
data$tColeoptera <- ((data$Coleoptera)^0.5)
shapiro.test(data$tColeoptera)
hist(data$tColeoptera, breaks = 20)
###P=0.06397

aov(tColeoptera~Leaf_NoLeaf2, data=data)
model <- aov(tColeoptera~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.0799

plot(Coleoptera~Leaf_NoLeaf2, data=data)

###Hydrophilidae##########
###Normality Test
shapiro.test(data$Hydrophilid)
hist(data$Hydrophilid)
###p-value = 0.005031

###Transformation Station
data$tHydrophilid <- (log(((data$Hydrophilid)+1)))
shapiro.test(data$tHydrophilid)
hist(data$tHydrophilid, breaks = 7)
###P=0.06397



aov(Hydrophilid~Leaf_NoLeaf2, data=data)
model <- aov(Hydrophilid~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.0799

plot(Coleoptera~Leaf_NoLeaf2, data=data)

###Dytiscidae##########
###Normality Test
shapiro.test(data$Dytiscidae)
hist(data$Dytiscidae)
###p-value = 0.08847


aov(Dytiscidae~Leaf_NoLeaf2, data=data)
model <- aov(Dytiscidae~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.151

plot(Dytiscidae~Leaf_NoLeaf2, data=data)

#######BoxPlot
data$Leaf_NoLeaf2 <- factor(data$Leaf_NoLeaf2,levels = c(0,1), labels = c("No Leaves Added", "With Leaves Added"))
plot(Dytiscidae~Leaf_NoLeaf2, data=data, xlab="", ylab="Dytiscidae Abundence", main="Preference for Pools with Leaves by Dytiscidae" )

###Colonists##########
###Normality Test
shapiro.test(data$Colonists)
hist(data$Colonists)
###p-value = 0.04596

###Transformation Station
data$tColonists <- (log(data$Colonists))
shapiro.test(data$tColonists)
hist(data$tColonists, breaks = 20)
###P=0.4315

aov(tColonists~Leaf_NoLeaf2, data=data)
model <- aov(tColonists~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.00137


data$Leaf_NoLeaf2 <- factor(data$Leaf_NoLeaf2,levels = c(0,1), labels = c("No Leaves Added", "With Leaves Added"))

plot(Colonists~Leaf_NoLeaf2, data=data, xlab="", ylab="Colonizer Abundence", main="Preference for Pools with Leaves by Total Abundance" )
plot(Hydrophilid~Leaf_NoLeaf2, data=data )
#####PLOTS

######Glyphicus#########

###Normality Test
shapiro.test(data$Glyphicus)
hist(data$Glyphicus)
###p-value = 0.2013

aov(Glyphicus~Leaf_NoLeaf2, data=data)
model <- aov(Glyphicus~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.31

plot(Glyphicus~Leaf_NoLeaf2, data=data)


#####Fasciatus###########

###Normality Test
shapiro.test(data$Fasciatus)
hist(data$Fasciatus)
###p-value = 0.0008725

###Transformation Station
data$tFasciatus <- (((data$Fasciatus)+1)^0.5)
shapiro.test(data$tFasciatus)
hist(data$tFasciatus, breaks = 20)
###P=0.1422

aov(tFasciatus~Leaf_NoLeaf2, data=data)
model <- aov(tFasciatus~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.00137

#############Richness###############

###Normality Test
shapiro.test(data$Shannon)
hist(data$Shannon)
###p-value = 0.5709

aov(Shannon~Leaf_NoLeaf2, data=data)
model <- aov(Shannon~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.02480.02

plot(Shannon~Leaf_NoLeaf2, data=data)



###Hyla eggs##########
###Normality Test
shapiro.test(data$Hyla.eggs)
hist(data$Hyla.eggs)
###p-value = 9.447e-07

###Transformation Station
data$tHyla.eggs <- (log(data$Hyla.eggs))
shapiro.test(data$tHyla.eggs)
hist(data$tHyla.eggs, breaks = 10)
###P=0.4517

aov(tHyla.eggs~Leaf_NoLeaf2, data=data)
model <- aov(tHyla.eggs~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.824


data$Leaf_NoLeaf2 <- factor(data$Leaf_NoLeaf2,levels = c(0,1), labels = c("No Leaves Added", "With Leaves Added"))

plot(Hyla.eggs~Leaf_NoLeaf2, data=data, xlab="", ylab="Colonizer Abundence", main="Preference for Pools with Leaves of Hyla Eggs" )

###Hyla eggs##########
###Normality Test
shapiro.test(data$Hyla.eggs)
hist(data$Hyla.eggs)
###p-value = 7.401e-07

###Transformation Station
data$tHyla.eggs <- (sqrt(data$Hyla.eggs))
shapiro.test(data$tHyla.eggs)
hist(data$tHyla.eggs, breaks = 10)
###P=0.4517

aov(tHyla.eggs~Leaf_NoLeaf2, data=data)
model <- aov(tHyla.eggs~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.54


data$Leaf_NoLeaf2 <- factor(data$Leaf_NoLeaf2,levels = c(0,1), labels = c("No Leaves", "With Leaves"))

plot(Hyla.eggs~Leaf_NoLeaf2, data=data, xlab="", ylab="Colonizer Abundence", main="Hyla Eggs" )


###Chironomidae##########
###Normality Test
shapiro.test(data$Chironimidae)
hist(data$Chironimidae)
###p-value = 6.674e-07

###Transformation Station
data$tChironimidae <- (log(data$Chironimidae))
shapiro.test(data$tChironimidae)
hist(data$tChironimidae, breaks = 10)
###P=0.07146

aov(tChironimidae~Leaf_NoLeaf2, data=data)
model <- aov(tChironimidae~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.824


data$Leaf_NoLeaf2 <- factor(data$Leaf_NoLeaf2,levels = c(0,1), labels = c("No Leaves", "With Leaves"))

plot(Chironimidae~Leaf_NoLeaf2, data=data, xlab="", ylab="Colonizer Abundence", main="Chironomidae" )


###Libellulidae##########
###Normality Test
shapiro.test(data$Libellulidae)
hist(data$Libellulidae)
###p-value = 7.468e-06

###Transformation Station
data$tLibellulidae <- ((sqrt(data$Libellulidae)))
shapiro.test(data$tLibellulidae)
hist(data$tLibellulidae, breaks = 10)
###P=0.07146

aov(tLibellulidae~Leaf_NoLeaf2, data=data)
model <- aov(tLibellulidae~Leaf_NoLeaf2, data=data)
summary(model)
###P=0.0289


data$Leaf_NoLeaf2 <- factor(data$Leaf_NoLeaf2,levels = c(0,1), labels = c("No Leaves", "With Leaves"))

plot(Libellulidae~Leaf_NoLeaf2, data=data, xlab="", ylab="Colonizer Abundence", main="Libellulidae" )



