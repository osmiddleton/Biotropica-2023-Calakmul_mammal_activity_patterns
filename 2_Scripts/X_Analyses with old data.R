##########################################################################################################
#                                                                                                        #
#  Neo-tropical felid activity patterns in relation to inter-specific interactions and prey              #
#       behaviour in the Calakmul Biosphere Reserve, Mexico'                                             #                                                         #
#                                                                                                        #
#  Written by: Owen Middleton    
#                                                                          #
#  Paper authors: OWen Middleton, Kathy Slater, Luis Pratas-Santiago, David Sima,                        #
#                 Antonio Lopez-Zen & Patrick Doncaster                                                  #
#                                                                                                        #    
#  In submission to: Biotropica                                                                          #                                                                                                        #
#                                                                                                        #


# The purpose of this code is to perform all of the analysis and produce figures for the above manuscript
#
#     The code aims to do the following:
#         1. Complete all analysis using circular statistics, including:
#               - Rayleigh Tests of Unimodality
#               - Mardia Watson-Wheeler Tests (MWW)
#               - Pair-wise comparisons of species' proportions of activity overlap (& bootstrapped)
#               - Calculate confidence intervals of overlap estimates through data bootstrapping
#         
#         2. Produce all figures used in the analysis
#               - All species pair-wise plots of proportion of overlap and CIs (Fig. 2)
#               - 9 predator-prey overlapping kernal density plots
#               - 3 predator-predator overlapping kernal density plots
#


#  ---- Load libraries ----

library(pacman)
pacman::p_load(lubridate, overlap, circular, tidyverse, reshape2, stringr, grDevices, scales, grid)


#  ---- Load data ----
data <- read.csv("./1_Data/Old/Biotropica2019 - Activity Pattern Data_UPDATED.csv", header = TRUE)
data <- data[!is.na(data$Species),]

# ---- 1. Annually consistent activity patterns: test using Mardia-Whatson Wheeler test ----

#   2016 data
species2016 <- data[data$Year == "2016",]
species2016 <- species2016[!is.na(species2016$Species),]

#   2017 data
species2017 <- data[data$Year == "2017",]
species2017 <- species2017[!is.na(species2017$Species),]

# MWW TEST FOR SPECIES YEAR PAIR-WISE RECORD DISTRIBUTION COMPARISON
# If they are identical record distributions between years, it is definitely safe to use

#create empty matrix to save output
species.comparison <- matrix(ncol = 3, nrow = 11)
# Get species names
speciesnames <- levels(data$Species)

# only data for one year for mice and white-tailed peccary, so delete
speciesnames
speciesnames <- speciesnames[-c(7,13)] 

for (i in 1:length(speciesnames)) {
  x <- species2016[species2016$Species == speciesnames[[i]],]
  y <- species2017[species2017$Species == speciesnames[[i]],]
  x <- x[c(3)]
  colnames(x) <- c("Species A")
  y <- y[c(3)] 
  colnames(y) <- c("Species B")
  output <- watson.wheeler.test(list(x$`Species A`,y$`Species B`))
  p <- output$p.value
  species.comparison[i,1] <- length(x$`Species A`)
  species.comparison[i,2] <- length(y$`Species B`)
  species.comparison[i,3] <- p
  print(i)
}

species.comparison <- as.data.frame(species.comparison)
species.comparison$Species <- speciesnames
colnames(species.comparison) <- c("2016.numbers", "2017.numbers", "MWW", "Species")
species.comparison$Sample.Size <- species.comparison$`2016.numbers`+ species.comparison$`2017.numbers`
show(species.comparison)

# The only species which do not have the same record distribution are: paca, great curassow and turkey, 
# however this is likely due to sampling issues with big differences between records between years

# Continue but worth bearing in mind...
#-------------------------------------------------------------------------------------------
#   STEP 2: RAYLEIGH TEST OF UNIMODALITY
#-------------------------------------------------------------------------------------------
#Species names again
speciesnames <- levels(data$Species)

#Create empty matrix for storing outputs
rayleighmatrix <- matrix(ncol = 2, nrow = length(speciesnames))

for (i in 1:length(speciesnames)){
  x <- data[data$Species == speciesnames[[i]],]$Radians
  y <- rayleigh.test(x, mu = NULL)
  rayleighmatrix[i,1] <- y$p.value
  rayleighmatrix[i,2] <- y$statistic
  print(i)
}

y$call
rayleighmatrix <- as.data.frame(rayleighmatrix)
rayleighmatrix$Species <- speciesnames
colnames(rayleighmatrix)[1:2] <- c("p-value","test.stat")
show(rayleighmatrix)

# Testing the H0 that there is no unimodal distribution - 
# A p-value greater than 0.05 shows a non-unimodal activity pattern

#-------------------------------------------------------------------------------------------
#   STEP 3: MARDIA WATSON WHEELER TESTS
#-------------------------------------------------------------------------------------------

# Need to set up to store pair-wise comparisons of species proportion of overlap
predator <- speciesnames[c(4,10,12)]

# Create empty matrix (3x11) to store output of pairwise comparisons
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

# A p-value < 0.05 shows difference in pair-wise activity record distributions i.e. jaguars and ocelots

#-------------------------------------------------------------------------------------------
#   STEP 4: ACTIVITY OVERLAP FOR PAIR-WISE SPECIES
#-------------------------------------------------------------------------------------------

# Compute estimate of activity overlap, using Dhat4, based upon recommendations (Linkie $ Rideout, 2014)
jagpecest<-overlapEst(data[data$Species == "Panthera onca",]$Radians,
                      data[data$Species == "Pecari tajacu",]$Radians, 
                      type = "Dhat4")
as.numeric(jagpecest)
par(mar = c(5,4,4,4)+0.1)

#Prepare output matrix for pairwise proportions of overlap
overlapmatrix <- matrix(ncol = length(predator), nrow = length(speciesnames))
colnames(overlapmatrix) <- names(predator)
rownames(overlapmatrix) <- names(speciesnames)

for (i in 1:length(predator)) {
  for (j in 1:length(speciesnames)){
    x <- data[data$Species == speciesnames[[j]],]$Radians
    y <- data[data$Species == predator[[i]],]$Radians
    overlapestimate <- overlapEst(x,y, type = "Dhat4")
    overlapmatrix[j,i] <- as.numeric(overlapestimate)
  }}

colnames(overlapmatrix) <- predator
row.names(overlapmatrix) <- speciesnames

show(overlapmatrix)


#-------------------------------------------------------------------------------------------
#   STEP 5: GET ACTIVITY OVERLAP CONFIDENCE INTERVALS
#-------------------------------------------------------------------------------------------

#Example: Jaguar and brocket deer

# Bootstrap data 1000 times
jagboot <- resample(data[data$Species == "Panthera onca",]$Radians, 1000)
dim(jagboot)
brocketboot <- resample(data[data$Species == "Manzama sp.",]$Radians, 1000)
dim(brocketboot)

# Using bootstrapped estimates, calculate bootstrapped proportions of overlap
jagbrocest2 <- bootEst(jagboot,brocketboot,adjust = c(NA, 1, NA))
dim(jagbrocest2)
BSmean <- colMeans(jagbrocest2)
as.numeric(BSmean[2])
tmp <- jagbrocest2[, 2]

# Calculate confidence intervals
bootCI(jagbrocest2[2], tmp)
bootCIlogit(jagbrocest2[2],tmp)
#### After this, in the same loop, all pairwise confidence intervals will be calculated...

# Create matrix to add results to
BOOTSTRAPPEDoverlapmatrix <- matrix(ncol = length(predator), nrow = length(speciesnames))
colnames(BOOTSTRAPPEDoverlapmatrix) <- names(predator)
rownames(BOOTSTRAPPEDoverlapmatrix) <- names(speciesnames)
CImatrix <- matrix(ncol = length(predator), nrow = length(speciesnames))
colnames(CImatrix) <- names(predator)
rownames(CImatrix) <- names(speciesnames)

for (i in 1:length(predator)) {
  for (j in 1:length(speciesnames)){
    x <- data[data$Species == speciesnames[[j]],]$Radians
    y <- data[data$Species == predator[[i]],]$Radians
    xboot <- resample(x,1000)
    yboot <- resample(y,1000)
    xyest <- bootEst(xboot,yboot,adjust = c(NA, 1, NA))
    BSmean <- colMeans(xyest)
    BOOTSTRAPPEDoverlapmatrix[j,i] <- as.numeric(BSmean[2])
    print(i)
    
    #Get confidence intervals
    tmp <- xyest[, 2]
    bootCI(xyest[2], tmp)
    help <- bootCIlogit(xyest[2],tmp)
    lowerCI <- help[5,1]
    lowerCI <- as.character(lowerCI) 
    upperCI <- help[5,2]
    upperCI <- as.character(upperCI) 
    CImatrix[j,i] <- paste(lowerCI,upperCI)
    print(j)
  }}

# Normal sample overlap estimates
tidyoverlap <- melt(overlapmatrix)

# Now, for bootstrapped overlaps
colnames(BOOTSTRAPPEDoverlapmatrix) <- predator
row.names(BOOTSTRAPPEDoverlapmatrix) <- speciesnames
tidyBOOTSTRAPPEDoverlap <- melt(BOOTSTRAPPEDoverlapmatrix)

# Now to work on the confidence intervals
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
masteroverlap <- masteroverlap[-c(4,22,36),]

str(masteroverlap)
masteroverlap$Var1 <- as.factor(masteroverlap$Var1)
masteroverlap$Var2 <- as.factor(masteroverlap$Var2)

levels(masteroverlap$Var1)

# Rename some factor levels
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Crax rubra"] <- "Great Curassow"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Cuniculus paca"] <- "Paca"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Dasyprocta punctata"] <- "Agouti"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Manzama sp."] <- "B. Deer"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Meleagris ocellata"] <- "Turkey"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Mouse"] <- "Mouse"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Nasua narica nelsoni"] <- "Coati"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Odocoileus virginianus"] <- "WT. Deer"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Panthera onca"] <- "Jaguar"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Pecari tajacu"] <- "C. Peccary"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Puma concolor"] <- "Puma"
levels(masteroverlap$Var1)[levels(masteroverlap$Var1) == "Leopardus pardalis"] <- "Ocelot"

levels(masteroverlap$Var2)[levels(masteroverlap$Var2) == "Leopardus pardalis"] <- "Ocelot"
levels(masteroverlap$Var2)[levels(masteroverlap$Var2) == "Panthera onca"] <- "Jaguar"
levels(masteroverlap$Var2)[levels(masteroverlap$Var2) == "Puma concolor"] <- "Puma"

# Change orders of factors so that it is correct order for ggplot
masteroverlap$Var1 <- factor(masteroverlap$Var1, levels = c("Mouse","Agouti","Paca","Coati",
                                                            "Turkey","Great Curassow",
                                                            "B. Deer","C. Peccary","WT. Deer",
                                                            "Ocelot","Puma","Jaguar"))

masteroverlap$Var2 <- factor(masteroverlap$Var2, levels = c("Jaguar","Puma","Ocelot"))

#------------------------------------------------------------------------------------------------------
#   STEP 6: PLOT FIGURES
#------------------------------------------------------------------------------------------------------

# Silhouettes were added in after
ggplot(masteroverlap, aes(x = Var1, y = value, colour = Var2)) +
  geom_point(stat = "identity", na.rm = TRUE, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, na.rm = TRUE, position = position_dodge(width = 0.5)) +
  geom_vline(aes(xintercept = 9.5), linetype = "dashed", size = 0.5) +
  scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8,1), limits = c(0.2,1.1)) +
  scale_colour_manual(name = "Predator", values = c("#330033","#FF6600", "#99CC00")) +
  xlab("Species") + ylab("Proportion of temporal overlap") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        legend.key = element_rect(fill = NA, colour = NA, size = 0.25),
        legend.text=element_text(size = 12),  legend.title = element_text(size = 12),
        axis.text=element_text(size=12), axis.title=element_text(size=12)) +
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5)) +
  coord_flip()


# Plot 9 predator-prey figures
# Assign correct format
par(mfrow = c(3,3),
    oma = c(1.5,1.5,0,0) + 0.5,
    mar = c(1,1,1,1) + 0.5)
# par(mfrow = c(1,1))

# Jaguar and white-tailed deer
overlapPlot(data[data$Species == "Panthera onca",]$Radians,
            data[data$Species == "Odocoileus virginianus",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.12),
            cex.axis = 1.5, yaxp = c(0,0.12,4))
#legend('top', c("Jaguar", "wT Deer"), lty=c(1,2), col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Puma and white-tailed deer
overlapPlot(data[data$Species == "Puma concolor",]$Radians,
            data[data$Species == "Odocoileus virginianus",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.12),
            cex.axis = 1.5, yaxp = c(0,0.12,4))
#legend('top', c("Puma", "WT Deer"), lty=c(1,2), col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Ocelot and white-tailed deer
overlapPlot(data[data$Species == "Leopardus pardalis",]$Radians,
            data[data$Species == "Odocoileus virginianus",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.12),
            cex.axis = 1.5, yaxp = c(0,0.12,4))
#legend('top', c("Ocelot", "WT Deer"), lty=c(1,2), col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Jaguar and Collared peccary
overlapPlot(data[data$Species == "Panthera onca",]$Radians,
            data[data$Species == "Pecari tajacu",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.12),
            cex.axis = 1.5, yaxp = c(0,0.12,4))
#legend('top', c("Jaguar", "C. Peccary"), lty=c(1,2), col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Puma and collared peccary
overlapPlot(data[data$Species == "Puma concolor",]$Radians,
            data[data$Species == "Pecari tajacu",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.12),
            cex.axis = 1.5, yaxp = c(0,0.12,4))
#legend('top', c("Puma", "C. Peccary"), lty=c(1,2), col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Ocelot and collared peccary
overlapPlot(data[data$Species == "Leopardus pardalis",]$Radians,
            data[data$Species == "Odocoileus virginianus",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.12),
            cex.axis = 1.5, yaxp = c(0,0.12,4))
#legend('top', c("Ocelot", "C. Peccary"), lty=c(1,2), col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Jaguar and paca
overlapPlot(data[data$Species == "Panthera onca",]$Radians,
            data[data$Species == "Cuniculus paca",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.12),
            cex.axis = 1.5, yaxp = c(0,0.12,4))
#legend('top', c("Jaguar", "Paca"), lty=c(1,2), col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Puma and paca
overlapPlot(data[data$Species == "Puma concolor",]$Radians,
            data[data$Species == "Cuniculus paca",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.12),
            cex.axis = 1.5, yaxp = c(0,0.12,4))
#legend('top', c("Puma", "Paca"), lty=c(1,2), col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Ocelot and paca
overlapPlot(data[data$Species == "Leopardus pardalis",]$Radians,
            data[data$Species == "Cuniculus paca",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.12),
            cex.axis = 1.5, yaxp = c(0,0.12,4))
#legend('top', c("Puma", "Paca"), lty=c(1,2), col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# SILHOUETTES AND AXES LABELS ERE ADDED IN AFTERWARDS

# Plot 3 predator-predator plots
# Get correct format
par(mfrow = c(2,2),
    oma = c(1.5,1.5,0,0) + 0.5,
    mar = c(1,1,1,1) + 0.5)

# Jaguar and puma
overlapPlot(data[data$Species == "Panthera onca",]$Radians,
            data[data$Species == "Puma concolor",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.1),
            cex.axis = 1.25, yaxp = c(0,0.1,4))
legend('top', c("Jaguar", "Puma"), lty=c(1,2), cex = 1.5,  col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)


# Jaguar and ocelot
overlapPlot(data[data$Species == "Panthera onca",]$Radians,
            data[data$Species == "Leopardus pardalis",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.1),
            cex.axis = 1.25, yaxp = c(0,0.1,4))
legend('top', c("Jaguar", "Ocelot"), lty=c(1,2), cex = 1.5,  col=c(1,4), bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Ocelot and puma
overlapPlot(data[data$Species == "Puma concolor",]$Radians,
            data[data$Species == "Leopardus pardalis",]$Radians,
            adjust=1, main=NA,xlab="",
            ylab="",lwd = 3, bty = "n", ylim = c(0,0.1),
            cex.axis = 1.25, yaxp = c(0,0.1,4))
legend('top', c("Puma", "Ocelot"), lty=c(1,2), col=c(1,4), cex = 1.5,  bty='n')
abline(v=c(6.5, (18+47/60) - 24), lty=1)
abline(v=c(5.5, (18+47/60) - 24), lty=3)
abline(v=c(7.5, (18+47/60) - 24), lty=3)
abline(v=c(19, (18+47/60) - 24), lty=1)
abline(v=c(18, (18+47/60) - 24), lty=3)
abline(v=c(20, (18+47/60) - 24), lty=3)

# Add in axes labels afterwards

#------------------------------------------------------------------------------------------------------
#     END OF CODE.
#------------------------------------------------------------------------------------------------------














