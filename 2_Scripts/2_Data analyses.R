#  
#  "Neo-tropical felid activity patterns in relation to potential prey and
#  intra-guild competitors in the Calakmul Biosphere Reserve, Mexico"                                                                                                #                                                                                                 #                                                                          #
#  Authors: Cristina Argudin-Violante, Owen Middleton, Kathy Slater, 
#           Esteban Dominguez-Bonilla & Patrick Doncaster                                                  
#  Published in: Biotropica   
#  Date (Accepted): April 2023

# Script 2: Analysing tidy camera trap data.
#
# Written by: Owen Middleton
# Contact: o.middleton@sussex.ac.uk
# Last updated: 01/04/2023

# The purpose of this code is to perform all of the analysis and produce figures for the above manuscript
#
#         1. Complete all circular statistics analyses:
#               - Rayleigh Tests of Unimodality
#               - Mardia Watson-Wheeler Tests (MWW)
#               - Pair-wise comparisons of species' proportions of activity overlap (& bootstrapped)
#               - Calculate confidence intervals of overlap estimates through bootstrapping
#         
#         2. Produce all figures used in the analysis
#               - Figure 2: predator-predator overlapping kernal density plots
#               - Figure 3: All species pair-wise plots of proportion of overlap and CIs
#               - Figure 4: 9 predator-prey overlapping kernal density plots
#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Housekeeping -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
if(1){
#  Load libraries ----
library(pacman)
pacman::p_load(lubridate, overlap, circular, tidyverse,
               reshape2, stringr, grDevices, scales, grid)
# Load data ----
data <- read.csv("./1_Data/Final 2019 Calakmul dataset.csv", header = TRUE)
data <- data[!is.na(data$Species),]
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Tidying data - getting ready for analysis -----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Get data for species of interest
# Reformatting time to radians
# Getting dawn and dusk times
if(1){

# Get data for species of interest ----
# Interested in wild non-human mammals and large birds that have >20 independent camera trap 
# records for activity pattern analysis.
if(1){
# 1. Mammal species
# Remove mammals with < 10 records
t <- as.data.frame(table(data$Mammal_Species))
t <- t %>% filter(Freq > 20)
species <- as.character(t$Var1)
# Remove empties, bats, humans, taira, grey fox (not including non-felid
# carnivores - out of scope)
mammal.species <- species[-which(species %in% c("","HUMANS",
                                                "Urocyon cinereoargenteus"))]
# 2. Large bird species
bird.species <- c("Meleagris ocellata", "Crax rubra")
# 3. Get all species of interest in the study.
all.species <- c(mammal.species, bird.species)
# 4. Subset data to include only those species of interest:
data <- data %>% filter(Species %in% all.species)
speciesnames <- unique(as.character(data$Species))
}

# Convert time to Radians ----
if(1){
# Get the value of the minute in the day...
data$NewTime <- as.character(factor(data$Time))
data$Hour <- substr(data$NewTime,1,2)
data$Mins <- substr(data$NewTime,4,5)
data$MinutesInDay <- as.numeric(data$Mins) +(as.numeric(data$Hour)*60)
# Convert to Radians
data$Radians <- (2*pi/1440)*data$MinutesInDay
}

# Get dawn and dusk ----
if(1){
dawn_minutes <- 60*6.5
dusk_minutes <- 60*19
}

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Perform all circular statistics ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----

# 1. Activity pattern categories (diurnal, nocturnal, crespuscular) ----
# Get data for Table 1 that gives percentage of records for each species that 
# fall into typical activity pattern categories.
if(1){
#Create empty matrix for storing outputs
activity_matrix <- matrix(ncol = 4, nrow = length(speciesnames))
# Assign categories to every activity record and get the % records in each period 
# for each species.
for (i in 1:length(speciesnames)){
  x <- data[data$Species == speciesnames[[i]],]$MinutesInDay
  nocturnal <- x[x < dawn_minutes - 60 |
                   x > dusk_minutes + 60]
  diurnal <- x[x > dawn_minutes + 60 &
                 x < dusk_minutes - 60]
  dawn <- x[x >= dawn_minutes - 60 &
                     x <= dawn_minutes+60]
  dusk <- x[x >= dusk_minutes-60 &
            x <= dusk_minutes+60]
  activity_matrix[i,1] <- round(length(diurnal)/length(x),digits = 2)
  activity_matrix[i,2] <- round(length(nocturnal)/length(x), digits = 2)
  activity_matrix[i,3] <- round((length(dawn)+length(dusk))/length(x), digits = 2)
  activity_matrix[i,4] <- speciesnames[i]
}

# Tidy matrix
activity_matrix <- as.data.frame(activity_matrix)
colnames(activity_matrix)[1:4] <- c("diurnal","nocturnal", "crepuscular", "species")
show(activity_matrix)
}

# 2. Mardia-Watson-Wheeler tests ----
# This analysis identifies whether species have matching
# or different activity patterns with data presented in Table S2.
if(1){

# Need to set up to store pair-wise comparisons of species proportion of overlap
predator <- c("Panthera onca", "Puma concolor", "Leopardus pardalis")

# Create empty matrix (3x13) to store output of pairwise comparisons
MWWStatmatrix <- matrix(ncol = length(predator), nrow = length(speciesnames))
MWWPmatrix <- matrix(ncol = length(predator), nrow = length(speciesnames))
colnames(MWWStatmatrix) <- predator
rownames(MWWStatmatrix) <- speciesnames
colnames(MWWPmatrix) <- predator
rownames(MWWPmatrix) <- speciesnames

for (i in 1:length(predator)) {
  for (j in 1:length(speciesnames)){
    x <- data[data$Species == speciesnames[[j]],]$Radians
    y <- data[data$Species == predator[[i]],]$Radians
    output <- watson.wheeler.test(list(x,y))
    output1 <- output$statistic
    output2 <- output$p.value
    MWWStatmatrix[j,i] <- output1
    MWWPmatrix[j,i] <- output2
}}

show(MWWPmatrix) 
show(MWWStatmatrix)

# A p-value > 0.05 shows no difference in pair-wise activity
# record distributions i.e. jaguars and ocelots
}

# 3. Pairwise bootstrapped activity pattern overlap AND confidence intervals ----
if(1){
# Create matrix to add results to
BOOTSTRAPPEDoverlapmatrix <- matrix(ncol = length(predator), nrow = length(speciesnames))
colnames(BOOTSTRAPPEDoverlapmatrix) <- names(predator)
rownames(BOOTSTRAPPEDoverlapmatrix) <- names(speciesnames)
CImatrix <- matrix(ncol = length(predator), nrow = length(speciesnames))
colnames(CImatrix) <- names(predator)
rownames(CImatrix) <- names(speciesnames)

# Run for-loop to get pairwise overlap and confidence intervals
# If one of the samples has a sample size < 50, Dhat 4 is used instead
for (i in 1:length(predator)) {
  for (j in 1:length(speciesnames)){
    x <- data[data$Species == speciesnames[[j]],]$Radians
    y <- data[data$Species == predator[[i]],]$Radians
    # Separate depending on the sample size
    if(length(x) < 50 | length(y) < 50) {
    xboot <- resample(x,1000)
    yboot <- resample(y,1000)
    xyest <- bootEst(xboot,yboot, type = "Dhat1")
    BSmean <- mean(xyest)
    BOOTSTRAPPEDoverlapmatrix[j,i] <- as.numeric(BSmean)
    print(i)
    Dhat <- overlapEst(x,y, type = "Dhat1")
    #Get confidence intervals
    tmp <- xyest
    #help <- bootCI(Dhat, tmp)
    help <- bootCIlogit(BSmean,tmp)
    lowerCI <- help[2,1]
    lowerCI <- as.character(lowerCI) 
    upperCI <- help[2,2]
    upperCI <- as.character(upperCI) 
    CImatrix[j,i] <- paste(lowerCI,upperCI)
    print(j)
    } else {
      # Otherwise use Dhat4
      xboot <- resample(x,1000)
      yboot <- resample(y,1000)
      xyest <- bootEst(xboot,yboot, type = "Dhat4")
      BSmean <- mean(xyest)
      BOOTSTRAPPEDoverlapmatrix[j,i] <- as.numeric(BSmean)
      print(i)
      Dhat <- overlapEst(x,y, type = "Dhat4")
      #Get confidence intervals
      tmp <- xyest
      #help <- bootCI(Dhat, tmp)
      help <- bootCIlogit(BSmean,tmp)
      lowerCI <- help[2,1]
      lowerCI <- as.character(lowerCI) 
      upperCI <- help[2,2]
      upperCI <- as.character(upperCI) 
      CImatrix[j,i] <- paste(lowerCI,upperCI)
      print(j)
    }
  }}

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Tidy results dataframes for plots ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Merge all dataframes ----
if(1){

# 1. Normal sample overlap estimates
tidyoverlap <- melt(overlapmatrix)

# 2. Bootstrapped overlaps
colnames(BOOTSTRAPPEDoverlapmatrix) <- predator
row.names(BOOTSTRAPPEDoverlapmatrix) <- speciesnames
tidyBOOTSTRAPPEDoverlap <- melt(BOOTSTRAPPEDoverlapmatrix)

# 3. Confidence intervals
CImatrix <- as.matrix(CImatrix)
colnames(CImatrix) <- predator
row.names(CImatrix) <- speciesnames
tidyCImatrix <- melt(CImatrix)
CIsplit <- str_split_fixed(tidyCImatrix$value," ",2)
CIsplit <- as.data.frame(CIsplit)
colnames(CIsplit)[1:2] <- paste(c("lower", "upper"))

# Create tidy dataframe for boostrapped estimates and confidence intervals
masteroverlap <- cbind(tidyBOOTSTRAPPEDoverlap,CIsplit)
str(masteroverlap)
masteroverlap$lower <- as.character(masteroverlap$lower)
masteroverlap$upper <- as.character(masteroverlap$upper)
masteroverlap$lower <- as.numeric(masteroverlap$lower)
masteroverlap$upper <- as.numeric(masteroverlap$upper)

# Remove same species comparisons
remove <- which(as.character(masteroverlap$Var1) == as.character(masteroverlap$Var2))
masteroverlap <- masteroverlap[-remove,]

str(masteroverlap)
masteroverlap$Var1 <- as.factor(masteroverlap$Var1)
masteroverlap$Var2 <- as.factor(masteroverlap$Var2)

levels(masteroverlap$Var1)
}

# Tidy dataframe ----
if(1){
# Rename some factor levels
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Crax rubra"] <- "Great Curassow"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Cuniculus paca"] <- "Paca"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Dasyprocta punctata"] <- "Agouti"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Mazama sp"] <- "B. Deer"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Meleagris ocellata"] <- "Turkey"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Nasua narica"] <- "Coati"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Odocoileus virginianus"] <- "WT. Deer"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Panthera onca"] <- "Jaguar"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Pecari tajacu"] <- "C. Peccary"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Puma concolor"] <- "Puma"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Leopardus pardalis"] <- "Ocelot"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Didelphis marsupialis"] <- "C. Opossum"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Tapirus bairdii"] <- "Tapir"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "HUMANS"] <- "Humans"


levels(masteroverlap$Var2)[levels(masteroverlap$Var2) == "Leopardus pardalis"] <- "Ocelot"
levels(masteroverlap$Var2)[levels(masteroverlap$Var2) == "Panthera onca"] <- "Jaguar"
levels(masteroverlap$Var2)[levels(masteroverlap$Var2) == "Puma concolor"] <- "Puma"

# Change orders of factors so that it is correct order for ggplot
masteroverlap$Var1 <- factor(masteroverlap$Var1, levels = c("Turkey","Great Curassow",
                                                            "C. Opossum","Agouti","Paca","Coati",
                                                            "B. Deer","C. Peccary","WT. Deer", "Tapir", "Humans",
                                                            "Ocelot","Puma","Jaguar"))

masteroverlap$Var2 <- factor(masteroverlap$Var2, levels = c("Jaguar","Puma","Ocelot"))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Generate manuscript figures ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# FIGURE 2: Temporal overlap overview ----
# Note, manual editing occurred post production to increase transparency of 
# pairwise species comparisons that are not likely to form predator-prey
# interactions. Further, number were placed on the figure to show whether
# interactions are likely to be of primary or secondary importance, using 
# information in Table S1.
if(1){
(p1 <- ggplot(masteroverlap, aes(x = Var1, y = value, colour = Var2)) +
  geom_point(stat = "identity", na.rm = TRUE, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, na.rm = TRUE, position = position_dodge(width = 0.5)) +
  geom_vline(aes(xintercept = 10.5), linetype = 1, size = 0.5, alpha = 0.5) +
  geom_vline(aes(xintercept = 6.5), linetype = 1, size = 0.5, alpha = 0.5) +
  geom_vline(aes(xintercept = 2.5), linetype = 1, size = 0.5, alpha = 0.5) +
  annotate(geom="text", x=13.3, y=0.13, label="Predators",color="black", size= 3) +
  annotate(geom="text", x=10.3, y=0.15, label="Prey (>15kg)",color="black", size= 3) +
  annotate(geom="text", x=6.3, y=0.15, label="Prey (<15kg)",color="black", size= 3) +
  annotate(geom="text", x=2.3, y=0.13, label="Bird prey",color="black", size= 3) +
  scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8,1), limits = c(0.1,1.1)) +
  scale_colour_manual(name = "Predator", values = c("green4","blue2", "grey2")) +
  #scale_colour_manual(name = "Predator", values = c("#330033","#FF6600", "#99CC00")) +
  xlab("Species") + ylab("Proportion of temporal overlap") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = NA, colour = NA, size = 0.25),
        legend.text=element_text(size = 12),  legend.title = element_text(size = 12),
        axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5)) +
  coord_flip())

# Save figure
ggsave("./3_Figures/Figure 1_Pairwise activity overlap.pdf", p1, 
       width = 7, height = 5.5, units = "in")
}

# FIGURE S1: All species individual activity patterns ----
if(1){
# Make density plots and add species IDs
a <- overlap::densityPlot(data[data$Species == "Panthera onca",]$Radians,rug = TRUE)
a$ID.pred <- "Jaguar"
b <- densityPlot(data[data$Species == "Puma concolor",]$Radians)
b$ID.pred <- "Puma"
c <- densityPlot(data[data$Species == "Leopardus pardalis",]$Radians)
c$ID.pred <- "Ocelot"
d <- densityPlot(data[data$Species == "Tapirus bairdii",]$Radians)
d$ID.pred <- "Tapir"
e <- densityPlot(data[data$Species == "Pecari tajacu",]$Radians)
e$ID.pred <- "Collared peccary"
f <- densityPlot(data[data$Species == "Mazama sp",]$Radians)
f$ID.pred <- "Brocket deer"
g <- densityPlot(data[data$Species == "Odocoileus virginianus",]$Radians)
g$ID.pred <- "White-tailed deer"
h <- densityPlot(data[data$Species == "Cuniculus paca",]$Radians, adjust = 0.8, rug = TRUE)
h$ID.pred <- "Paca"
i <- densityPlot(data[data$Species == "Nasua narica",]$Radians)
i$ID.pred <- "Coati"
j <- densityPlot(data[data$Species == "Dasyprocta punctata",]$Radians)
j$ID.pred <- "Agouti"
k <- densityPlot(data[data$Species == "Didelphis marsupialis",]$Radians, adjust = 0.8, rug = TRUE)
k$ID.pred <- "Common opossum"
l <- densityPlot(data[data$Species == "Crax rubra",]$Radians)
l$ID.pred <- "Great curassow"
m <- densityPlot(data[data$Species == "Meleagris ocellata",]$Radians)
m$ID.pred <- "Ocellated turkey"

all.densities <- rbind(a,b,c,d,e,f,g,h,i,j,k,l,m)

all.densities$ID.pred <- factor(all.densities$ID.pred, levels = c("Jaguar","Puma","Ocelot",
                                                                  "Tapir","Collared peccary", "Brocket deer","White-tailed deer",
                                                                  "Coati","Paca","Agouti","Common opossum",
                                                                  "Great curassow", "Ocellated turkey"))

# Generate figure:
(density.plots <- ggplot() +
    geom_line(data = all.densities, aes(x = x, y = y), size = 1, linetype = 1, colour = "black", alpha = 0.5) +
    geom_vline(xintercept = c(6.5,19), alpha = 0.5) +
    geom_vline(xintercept = c(5.5,7.5,18,20), linetype = 2, alpha = 0.5) +
    labs(y = "Proportion of activity schedule", x = "Time of day") +
    scale_x_continuous(breaks = c(0,6,12,18,24), labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
    facet_wrap(.~ID.pred, ncol = 3) +
    theme_bw() + theme(panel.grid = element_blank(),
                            strip.text.x = element_text(colour = "black", size = 16)))
# Save figure:
ggsave("./3_Figures/Figure S2_All species activity patterns.pdf", density.plots,
       width = 7, height = 9, units = "in")
}

# FIGURE 3: (Intra-guild) Felid activity pattern overlap ----
# Note, manual editing occurred post-production to remove the bottom right hand
# facet for this figure.
if(1){

# 1. Jaguar and puma
# Get density plot
p <- overlapPlot(data[data$Species == "Panthera onca",]$Radians,
                 data[data$Species == "Puma concolor",]$Radians)
# Extract values for ggplot
p$smaller <- p$densityA > p$densityB
p$shaded <- NA
p[p$smaller == TRUE,]$shaded <- p[p$smaller == TRUE,]$densityB
p[p$smaller != TRUE,]$shaded <- p[p$smaller != TRUE,]$densityA
p$ID.pred1 <- "Jaguar"
p$ID.pred2 <- "Puma"

# 2. Jaguar and ocelot
a <- overlapPlot(data[data$Species == "Panthera onca",]$Radians,
                 data[data$Species == "Leopardus pardalis",]$Radians)
# Extract values for ggplot
a$smaller <- a$densityA > a$densityB
a$shaded <- NA
a[a$smaller == TRUE,]$shaded <- a[a$smaller == TRUE,]$densityB
a[a$smaller != TRUE,]$shaded <- a[a$smaller != TRUE,]$densityA
a$ID.pred1 <- "Jaguar"
a$ID.pred2 <- "Ocelot"

# 3. Puma and ocelot
b <- overlapPlot(data[data$Species == "Puma concolor",]$Radians,
                 data[data$Species == "Leopardus pardalis",]$Radians)
# Extract values for ggplot
b$smaller <- b$densityA > b$densityB
b$shaded <- NA
b[b$smaller == TRUE,]$shaded <- b[b$smaller == TRUE,]$densityB
b[b$smaller != TRUE,]$shaded <- b[b$smaller != TRUE,]$densityA
b$ID.pred1 <- "Puma"
b$ID.pred2 <- "Ocelot"

# Bind the dataframes:
master <- rbind(p,a,b)

# Plot main figure
(p2 <- ggplot() +
  geom_line(data = master, aes(x = x, y = densityA), size = 1, linetype = 1, colour = "black", alpha = 0.5) +
  geom_line(data = master, aes(x = x, y = densityB), linetype = 1, colour = "blue", size = 1, alpha = 0.5) +
  geom_ribbon(data = master, aes(x = x, ymin = 0, ymax = shaded), alpha = 0.1) +
  geom_vline(xintercept = c(6.5,19), alpha = 0.5) +
  geom_vline(xintercept = c(5.5,7.5,18,20), linetype = 2, alpha = 0.5) +
  labs(y = "Proportion of activity schedule", x = "Time of day") +
  scale_x_continuous(breaks = c(0,6,12,18,24), labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  facet_grid(ID.pred2~ID.pred1) +
  theme_minimal() + theme(panel.grid = element_blank(),
                          strip.text.y = element_text(colour = "blue3", size = 16),
                          strip.text.x = element_text(colour = "black", size = 16)))

# Save the figure:
ggsave("./3_Figures/Figure 2_Intraguild activity patterns overlap.pdf", p2,
       width = 7, height = 5, units = "in")
}

# FIGURE 4: Predator-prey overlap ----
if(1){

  # 1. All predators with brocket deer:
if(1){
# Jaguar and brocket deer
p1 <- overlapPlot(data[data$Species == "Panthera onca",]$Radians,
                 data[data$Species == "Mazama sp",]$Radians)
# Extract values for ggplot
p1$smaller <- p1$densityA > p$densityB
p1$shaded <- NA
p1[p1$smaller == TRUE,]$shaded <- p1[p1$smaller == TRUE,]$densityB
p1[p1$smaller != TRUE,]$shaded <- p1[p1$smaller != TRUE,]$densityA
p1$ID.pred <- "Jaguar"
p1$ID.prey <- "Brocket deer"


# Puma and brocket deer
p2 <- overlapPlot(data[data$Species == "Puma concolor",]$Radians,
                 data[data$Species == "Mazama sp",]$Radians)
# Extract values for ggplot
p2$smaller <- p2$densityA > p2$densityB
p2$shaded <- NA
p2[p2$smaller == TRUE,]$shaded <- p2[p2$smaller == TRUE,]$densityB
p2[p2$smaller != TRUE,]$shaded <- p2[p2$smaller != TRUE,]$densityA
p2$ID.pred <- "Puma"
p2$ID.prey <- "Brocket deer"


# Ocelot and brocket deer
p3 <- overlapPlot(data[data$Species == "Leopardus pardalis",]$Radians,
                 data[data$Species == "Mazama sp",]$Radians)
# Extract values for ggplot
p3$smaller <- p3$densityA > p3$densityB
p3$shaded <- NA
p3[p3$smaller == TRUE,]$shaded <- p3[p3$smaller == TRUE,]$densityB
p3[p3$smaller != TRUE,]$shaded <- p3[p3$smaller != TRUE,]$densityA
p3$ID.pred <- "Ocelot"
p3$ID.prey <- "Brocket deer"
}

# 2. All predators with collared peccary:
if(1){
# Jaguar and Collared peccary
p4 <- overlapPlot(data[data$Species == "Panthera onca",]$Radians,
                 data[data$Species == "Pecari tajacu",]$Radians)
# Extract values for ggplot
p4$smaller <- p4$densityA > p4$densityB
p4$shaded <- NA
p4[p4$smaller == TRUE,]$shaded <- p4[p4$smaller == TRUE,]$densityB
p4[p4$smaller != TRUE,]$shaded <- p4[p4$smaller != TRUE,]$densityA
p4$ID.pred <- "Jaguar"
p4$ID.prey <- "Collared peccary"

# Puma and collared peccary
p5 <- overlapPlot(data[data$Species == "Puma concolor",]$Radians,
                 data[data$Species == "Pecari tajacu",]$Radians)
# Extract values for ggplot
p5$smaller <- p5$densityA > p5$densityB
p5$shaded <- NA
p5[p5$smaller == TRUE,]$shaded <- p5[p5$smaller == TRUE,]$densityB
p5[p5$smaller != TRUE,]$shaded <- p5[p5$smaller != TRUE,]$densityA
p5$ID.pred <- "Puma"
p5$ID.prey <- "Collared peccary"

# Ocelot and collared peccary
p6 <- overlapPlot(data[data$Species == "Leopardus pardalis",]$Radians,
                 data[data$Species == "Pecari tajacu",]$Radians)
# Extract values for ggplot
p6$smaller <- p6$densityA > p6$densityB
p6$shaded <- NA
p6[p6$smaller == TRUE,]$shaded <- p6[p6$smaller == TRUE,]$densityB
p6[p6$smaller != TRUE,]$shaded <- p6[p6$smaller != TRUE,]$densityA
p6$ID.pred <- "Ocelot"
p6$ID.prey <- "Collared peccary"
}

# 3. All predators and common opossum
if(1){
# Jaguar and common opossum
p7 <- overlapPlot(data[data$Species == "Panthera onca",]$Radians,
                 data[data$Species == "Didelphis marsupialis",]$Radians)
# Extract values for ggplot
p7$smaller <- p7$densityA > p7$densityB
p7$shaded <- NA
p7[p7$smaller == TRUE,]$shaded <- p7[p7$smaller == TRUE,]$densityB
p7[p7$smaller != TRUE,]$shaded <- p7[p7$smaller != TRUE,]$densityA
p7$ID.pred <- "Jaguar"
p7$ID.prey <- "Common opossum"

# Puma and common opossum
p8 <- overlapPlot(data[data$Species == "Puma concolor",]$Radians,
                 data[data$Species == "Didelphis marsupialis",]$Radians)
# Extract values for ggplot
p8$smaller <- p8$densityA > p8$densityB
p8$shaded <- NA
p8[p8$smaller == TRUE,]$shaded <- p8[p8$smaller == TRUE,]$densityB
p8[p8$smaller != TRUE,]$shaded <- p8[p8$smaller != TRUE,]$densityA
p8$ID.pred <- "Puma"
p8$ID.prey <- "Common opossum"

# Ocelot and common opossum
p9 <- overlapPlot(data[data$Species == "Leopardus pardalis",]$Radians,
                 data[data$Species == "Didelphis marsupialis",]$Radians)
# Extract values for ggplot
p9$smaller <- p9$densityA > p9$densityB
p9$shaded <- NA
p9[p9$smaller == TRUE,]$shaded <- p9[p9$smaller == TRUE,]$densityB
p9[p9$smaller != TRUE,]$shaded <- p9[p9$smaller != TRUE,]$densityA
p9$ID.pred <- "Ocelot"
p9$ID.prey <- "Common opossum"
}

master <- rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9)

(p3 <- ggplot() +
  geom_line(data = master, aes(x = x, y = densityA), size = 1, linetype = 1, colour = "black", alpha = 0.5) +
  geom_line(data = master, aes(x = x, y = densityB), linetype = 1, colour = "blue", size = 1, alpha = 0.5) +
  geom_ribbon(data = master, aes(x = x, ymin = 0, ymax = shaded), alpha = 0.1) +
  geom_vline(xintercept = c(6.5,19), alpha = 0.5) +
  geom_vline(xintercept = c(5.5,7.5,18,20), linetype = 2, alpha = 0.5) +
  labs(y = "Proportion of activity schedule", x = "Time of day") +
  scale_x_continuous(breaks = c(0,6,12,18,24), labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  facet_grid(ID.prey~ID.pred) +
  theme_minimal() + theme(panel.grid = element_blank(),
                          strip.text.y = element_text(colour = "blue3", size = 14),
                          strip.text.x = element_text(colour = "black", size = 14)))

# Save figure:
ggsave("./3_Figures/Figure 3_Predator prey overlap.pdf", p3,
       width = 9.5, height = 6.5, units = "in")

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# End of Script ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# Graveyard ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----

# Rayleigh Test of unimodality ----
# Note, this is no longer in the manuscript.
# Reviewers suggested removing this for the above
# descriptive analysis of activity patterns.
if(0){
  #Species names again
  speciesnames <- levels(factor(data$Species))
  #Create empty matrix for storing outputs
  rayleighmatrix <- matrix(ncol = 2, nrow = length(speciesnames))
  for (i in 1:length(speciesnames)){
    x <- data[data$Species == speciesnames[[i]],]$Radians
    y <- rayleigh.test(x, mu = NULL)
    rayleighmatrix[i,1] <- y$p.value
    rayleighmatrix[i,2] <- y$statistic
    print(i)
  }
  rayleighmatrix <- as.data.frame(rayleighmatrix)
  rayleighmatrix$Species <- speciesnames
  colnames(rayleighmatrix)[1:2] <- c("p-value","test.stat")
  show(rayleighmatrix) 
  # Testing the H0 that there is no unimodal distribution - 
  # A p < 0.05 shows a unimodal activity pattern
}

# Example of pairwise bootstap for single pairwise comparison ----
if(1){
#Example: Jaguar and brocket deer

# Calculate estimates of overlap:
Dhats <- overlapEst(data[data$Species == "Panthera onca",]$Radians,
                    data[data$Species == "Mazama sp",]$Radians, type = "Dhat4")

# Do smoothed bootstrap
# Generate 999 data sets:
jagboot <- resample(data[data$Species == "Panthera onca",]$Radians, 1000)
brocketboot <- resample(data[data$Species == "Mazama sp",]$Radians, 1000)
bs <- bootEst(jagboot, brocketboot, type="Dhat4") # Takes a while
mean(bs)

# Get confidence intervals
bootCI(Dhats[1], bs)['norm0', ]
bootCI(Dhats[1], bs)['basic0', ]
}

# 3. Pair-wise species activity pattern overlap ----
# Note, this is not used in the manuscript. The bootstrapped calculations of
# activity overlaps is used.
# Compute estimates of activity overlap, using Dhat4, based upon recommendations
# (Linkie & Rideout, 2014)
if(1){
  # Prepare output matrix for pairwise proportions of overlap
  overlapmatrix <- matrix(ncol = length(predator), nrow = length(speciesnames))
  colnames(overlapmatrix) <- names(predator)
  rownames(overlapmatrix) <- names(speciesnames)
  
  # Run loop to get all pair-wise overlaps
  for (i in 1:length(predator)) {
    for (j in 1:length(speciesnames)){
      x <- data[data$Species == speciesnames[[j]],]$Radians
      y <- data[data$Species == predator[[i]],]$Radians
      
      # If one of the samples has <50 samples, use Dhat 1, otherwise use Dhat 4
      if(length(x) < 50 | length(y) < 50) {
        
        overlapestimate <- overlapEst(x,y, type = "Dhat1")
        overlapmatrix[j,i] <- as.numeric(overlapestimate)
        
      } else {
        
        overlapestimate <- overlapEst(x,y, type = "Dhat4")
        overlapmatrix[j,i] <- as.numeric(overlapestimate)
        
      }
    }}
  
  colnames(overlapmatrix) <- predator
  row.names(overlapmatrix) <- speciesnames
  
  show(overlapmatrix)
}