########################################################################################################
########################## Early sarcopterygian morphological disparity ################################
#############################                 Mai 2025                  ################################
########################################################################################################


getwd()
setwd("C:/Users/olivi/Documents/UQAR/DOC/Redaction-Articles-These/Article 2 - Sarcopterygian disparity/Analyses-files/Final-files")


library(geomorph)
library(ggplot2)
library(Momocs)
library(ggfortify)
library(ggthemes)
library(ggdark)
library(Morpho)
library(readxl)
library(dplyr)
library(vegan)
library(divDyn)
library(dispRity)
library(scales)
library(ggpubr)
library(forcats)

########################################################################################################
############################################# COLOR CODES ##############################################
########################################################################################################

# Color for epochs #
cols_epoch <- c("Lower Devonian" = "#E5AC4D", "Middle Devonian" = "#F1C868", "Upper Devonian" = "#F1E19D", "Mississippian" = "#678F66", "Pennsylvanian" = "#7EBCC6") #international chronostratigraphic chart colors

# Color for ages #
cols_age <- c("Emsian" = "#e9ba6a", "Eifelian" = "#F1D576", "Givetian" = "#F1E185", "Frasnian" = "#F2EDAD", "Famennian" = "#F2EDB3" , "Tournaisian" = "#8CB06C", "Visean"= "#A6B96C", "Serpukhovian" = "#BFC26B", "Bashkirian" = "#99C2B5", "Moscovian" = "#B3CBB9", "Gzhelian" = "#CCD4C7") #color codes following the international chronostratigraphic chart 

# Color for aquatic habitats #
cols_habitats <- c("Freshwater" = "#00c5ff", "Marine" = "#0080ff", "Estuary" = "#00ffb9")
cols_habitats2 <- c("Estuary" = "#00ffb9", "Bay" = "#ff9e00", "Oxbow lake/meandering river" = "#5a3801", "Calm  freshwater lake" = "#8407fa", "Dynamic freshwater lake" = "#470289", "Lagoon" = "#fd03f5", "Alluvial plain" = "#fdf503", "Delta" = "#fd0330","Coastal marine" = "#0d0299", "Reef" = "#049902")

# Color for groups #
cols_group <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Dipnoi" = "#FA8072", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F", "Elpistostegalia" = "#ffd496", "Tetrapoda" = "#6f4202") #Assign each group a color


########################################################################################################
################################# CONVEX HULL FUNCTION -- LAURENT HOULE ################################
########################################################################################################

# Creating convex hull: if there is only one observation for one category, the point is left alone in the graph
convex.hulls <- function(data, name.x, name.y, name.fill){
  
  library(tidyverse)
  library(grDevices)
  library(swaRm)
  
  data$X <- pca.class[[name.x]]
  data$Y <- data[[name.y]]
  data$FILL <- factor(data[[name.fill]])
  
  hull1 <- data %>%
    dplyr::slice(chull(X, Y))
  hull1 <- hull1[NULL,]
  SA <- c()
  ind <- 0
  for(j in 1:length(levels(data$FILL))){
    new <- filter(data, FILL == levels(data$FILL)[j])
    
    if(length(new$X) > 2){
      ind <- ind + 1
      hull.new <- new %>%
        dplyr::slice(chull(X, Y))
      
      Surf_area <- chull_area(hull.new$X,hull.new$Y)
      SA <- c(SA, Surf_area)
      names(SA)[ind] <- levels(data$FILL)[j]
      hull1 <- rbind(hull1,hull.new)
    }
    
  }
  
  return(list(table = hull1, surface_area = SA))
}
# end of the function

##############################################################################################################
##################################### SPECIES RICHNESS THROUGH TIME ##########################################
##############################################################################################################

data(stages) # loading dataset from divDyn package containing information about chronostratigraphic chart
stages_upd <- read_excel(file.choose(),1) #open stages-update.xlsx : the time for each stage from Ordovician, Silurian, Devonian and Carboniferous periods has been updated following the chronostratigraphic time chart 2024.

# diversity species in morphological matrix

list_sarco <- read_excel(file.choose(),1) #open Matrix-all-species.xlsx 
list_sarco<- list_sarco %>%
  mutate_if(is.character, as.factor)

list_sarco %>% mutate_at(c("max_ma", "min_ma"), as.numeric) #transform LAD/FAD column from character to numeric
list_sarco$me_ma <- apply(list_sarco[, c("max_ma", "min_ma")], 1, mean) # calculate the median age of each species

flDual <- fadlad(list_sarco, tax = "Species", age = c("max_ma", "min_ma")) # create a species x FAD/LAD matrix

# Occurrence sarcopterygian species #
Sarco_occu <- divDyn(list_sarco, tax = "Species", age = "me_ma", breaks = c(425.6, 407.5,389.4,371.3,353.2,335.1,317,298.9)) #calculate occurrence for the timebins
plot(Sarco_occu$me_ma, Sarco_occu$divSIB, xlim = rev(range(Sarco_occu$me_ma)))

list_sarco$mid <- stages$mid[list_sarco$Time_bin] #create a new column in the dataframe with the median age of the timebin
tsplot(stages_upd, shading="stage", boxes=c("short","series"),xlim=c(440,298.9), labels.args=list(cex=0.7), boxes.col=c("seriesCol", "systemCol")) #create background with the timechart
divDyn::ranges(list_sarco, tax="Species", bin=c("max_ma","min_ma"), labs=T, labels.args=list(cex=0.6, font = 3), occs=TRUE, filt="orig") #add duration of each species; font = 3 is italic

# Diveristy through time #
Sarco_div <- divDyn(list_sarco, bin = "Time_bin", tax = "Species")

tsplot(stages_upd, shading="stage", boxes=c("short", "series"), xlim=13:29, ylab="Richness (diversity)", ylim=c(0,25), labels.args = list(cex=0.7), boxes.col = c("seriesCol", "systemCol"))
lines(stages$mid[27:42], Sarco_div$divRT[27:42], col="black", lwd=2)

# diversity 279 early sarcopterygian species

list_all <- read_excel(file.choose(),1) #open Matrix-all-species-diversity.xlsx 
list_all<- list_all %>%
  mutate_if(is.character, as.factor)

list_all %>% mutate_at(c("max_ma", "min_ma"), as.numeric) #transform LAD/FAD column from character to numeric
list_all$me_ma <- apply(list_all[, c("max_ma", "min_ma")], 1, mean) # calculate the median age of each species

flDual_all <- fadlad(list_allo, tax = "Species", age = c("max_ma", "min_ma")) # create a species x FAD/LAD matrix

# Occurrence sarcopterygian species #
Sarco_occu_all <- divDyn(list_all, tax = "Species", age = "me_ma", breaks = c(425.6, 407.5,389.4,371.3,353.2,335.1,317,298.9)) #calculate occurrence for the timebins
plot(Sarco_occu_all$me_ma, Sarco_occu_all$divSIB, xlim = rev(range(Sarco_occu_all$me_ma)))

list_all$mid <- stages$mid[list_all$Time_bin] #create a new column in the dataframe with the median age of the timebin
tsplot(stages_upd, shading="stage", boxes=c("short","series"),xlim=c(440,298.9), labels.args=list(cex=0.7), boxes.col=c("seriesCol", "systemCol")) #create background with the timechart
divDyn::ranges(list_all, tax="Species", bin=c("max_ma","min_ma"), labs=T, labels.args=list(cex=0.6, font = 3), occs=TRUE, filt="orig") #add duration of each species; font = 3 is italic

# Diveristy through time #
Sarco_div_all <- divDyn(list_all, bin = "Time_bin", tax = "Species")

tsplot(stages_upd, shading="stage", boxes=c("short", "series"), xlim=13:29, ylab="Richness (diversity)", ylim=c(0,70), labels.args = list(cex=0.7), boxes.col = c("seriesCol", "systemCol"))
lines(stages$mid[27:42], Sarco_div_all$divRT[27:42], col="black", lwd=2)
lines(stages$mid[27:42], Sarco_div$divRT[27:42], col="blue", lty= 2, lwd=2)

##############################################################################################################
########################################## HABITAT OCCUPATION ################################################
##############################################################################################################

data_habitat <- read_excel(file.choose(),1)

data_habitat <- data_habitat %>%
  mutate_if(is.character, as.factor)

fq_habitat_age <- data_habitat %>%
  group_by(Paleoenvironments, Age) %>%
  summarise(count = n())
print(fq_habitat_age)

level_order <- factor(fq_habitat_age$Age, level = c("Ludfordian", "Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Viséan", "Serpukhovian", "Bashkirian", "Moscovian", "Kasimovian", "Gzhelian"))
 
ggplot(fq_habitat_age, aes(fill=Paleoenvironments, y=count, x=level_order)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = cols_habitats2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")

fq_habitat_gp <- data_habitat %>%
  group_by(`Aquatic habitat`, Group) %>%
  summarise(count = n())
print(fq_habitat_gp)

ggplot(fq_habitat_gp, aes(fill=Group, y=count, x=`Aquatic habitat`)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = cols_group) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right")


##############################################################################################################
########################################### FULLBODY DISPARITY ############################################
##############################################################################################################

######################################### Import datasets ################################################

# TPS files with landmarks and semilandmarks
Sarco_PC <- readland.tps("Postcranial-final.TPS",specID="imageID",negNA = TRUE,readcurves = TRUE,warnmsg = TRUE)
dim(Sarco_PC)

name_PC <- dimnames(Sarco_PC)[[3]]
name_PC

# Excel file with the list of species and their ages
Period_sarco_PC <- read_excel(file.choose(), 1) #Open List-species-PC.xlsx
Period_sarco_PC

Period_sarco_PC<- Period_sarco_PC %>%
  mutate_if(is.character, as.factor)

# Convert curves into semi-landmarks
matrice_PC <- rbind(define.sliders(12:41), define.sliders(42:61), define.sliders(62:71), define.sliders(72:107), define.sliders(108:207))

########################################## GPA and PCA ####################################################

# Procrustes superimposition #
data.super_PC <- gpagen(Sarco_PC, curves = matrice_PC, ProcD = FALSE)
attributes(data.super_PC)

plot(data.super_PC) 

# Principal component analyses #
pca_PC <- gm.prcomp(data.super_PC$coords)
pca_PC

plot_PC <- plot(pca_PC, axis1 = 1, axis2 = 2)
text(pca_PC[["x"]][,1], pca_PC[["x"]][,2], labels = name_PC) # PCA 1 vs 2

plot_PC2_PC <- plot(pca_PC, axis1=2, axis2=3)
text(pca_PC[["x"]][,2], pca_PC[["x"]][,3], labels = name_PC) # PCA 2 vs 3

# Saving PC scores #

PC.scores_PC <- pca_PC$x 
as.data.frame(PC.scores_PC) # Save PC scores as a data frame object

write.csv(PC.scores_PC,"PC.scores_PC.csv",row.names=TRUE) # Save PC scores as a csv file

# Extracting eigenvectors #

ev_PC <- pca_PC$rotation
write.csv(ev_PC, "eigenvectors-PC.csv")


##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_PC), aes(x=PC.scores_PC[,1], y=PC.scores_PC[,2], label = name_PC, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_PC$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 49.20 %", y = "PC2 = 14.94 %" ) 

#PC2 vs PC3 -- color epoch
ggplot(as.data.frame(PC.scores_PC), aes(x=PC.scores_PC[,2], y=PC.scores_PC[,3], label = name_PC, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_PC$Epoque), color = "black", size = 6, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 14.94 %", y = "PC3 = 10.28 %" )

#PC1 vs PC2 -- color Aquatic habitats 
ggplot(as.data.frame(PC.scores_PC), aes(x=PC.scores_PC[,1], y=PC.scores_PC[,2], label = name_PC, fontface = "italic")) +
  labs(fill="Aquatic habitats") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_PC$`Aquatic habitat`), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 49.20 %", y = "PC2 = 14.94 %" ) 

#PC2 vs PC3 -- color Aquatic habitats  
ggplot(as.data.frame(PC.scores_PC), aes(x=PC.scores_PC[,2], y=PC.scores_PC[,3], label = name_PC, fontface = "italic")) +
  labs(fill="Aquatic habitat") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_PC$`Aquatic habitat`), color = "black", size = 6, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 14.94 %", y = "PC3 = 10.28 %" )

#PC1 vs PC2 -- color groups

ggplot(as.data.frame(PC.scores_PC), aes(x=PC.scores_PC[,1], y=PC.scores_PC[,2], label = name_PC, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_PC$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 49.20 %", y = "PC2 = 14.94 %" ) 

#PC2 vs PC3 -- color groups
ggplot(as.data.frame(PC.scores_PC), aes(x=PC.scores_PC[,2], y=PC.scores_PC[,3], label = name_PC, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_PC$Group), color = "black", size = 6, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 14.94 %", y = "PC3 = 10.28 %" )

################### Calculating convex hulls for the PC1 vs PC2 #######################

# Add classifier variables (group, age, habitats) by hand to .csv file
# loading .csv file with PC scores and classifiers:
pca.class <- read.csv(file.choose()) #Open PC.scores_PC-class.csv

# for the epochs
CHull_epoch_PC <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Epoch")

CHull_epoch_PC.table <- CHull_epoch_PC$table
CHull_epoch_PC$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Epoch, shape = Epoch)) +
  labs(fill="Epoch") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Epoch), shape = 21, size = 3, stroke = 0.10) + 
  scale_fill_manual(values = cols_epoch) +
  scale_color_manual(values= cols_age)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 49.20 %", y = "PC2 = 14.94 %" ) +
  geom_polygon(data = CHull_epoch_PC.table, aes(x = X, y = Y), alpha = 0.5)

# for the Aquatic habitats 
CHull_habitat_PC <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Habitat")

CHull_habitat_PC.table <- CHull_habitat_PC$table
CHull_habitat_PC$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Habitat, shape = Habitat)) +
  labs(fill="Habitat") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Habitat), shape = 21, size = 3, stroke = 0.10) + 
  scale_fill_manual(values = cols_habitats) +
  scale_color_manual(values= cols_habitats)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 49.20 %", y = "PC2 = 14.94 %" ) +
  geom_polygon(data = CHull_habitat_PC.table, aes(x = X, y = Y), alpha = 0.5)

# for the groups
CHull_group_PC <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Group")

CHull_group_PC.table <- CHull_group_PC$table
CHull_group_PC$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Group, shape = Group)) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Group), shape = 21, size = 3, stroke = 0.10) + 
  scale_fill_manual(values = cols_group) +
  scale_color_manual(values= cols_group)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 49.20 %", y = "PC2 = 14.94 %" ) +
  geom_polygon(data = CHull_group_PC.table, aes(x = X, y = Y), alpha = 0.5)

################################# Morphological disparity analyses ###########################################

gdf_PC_epoch <- geomorph.data.frame(data.super_PC, species = Period_sarco_PC$Species, Time = Period_sarco_PC$Epoque) #gdf for the epochs

gdf_PC_age <- geomorph.data.frame(data.super_PC, species = Period_sarco_PC$Species, Time = Period_sarco_PC$Age) #gdf for the ages

gdf_PC_habitat <- geomorph.data.frame(data.super_PC, species = Period_sarco_PC$Species, Habitat = Period_sarco_PC$`Aquatic habitat`) #gdf for the aquatic habitats

gdf_PC_gp <- geomorph.data.frame(data.super_PC, species = Period_sarco_PC$Species, Group = Period_sarco_PC$Group) #gdf for the groups

# Sum of Procrustes variances : age ######################

SOV_mean_PC_age <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_PC_age, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_PC_age)

SOV_PC_age <- morphol.disparity(coords~Time,groups=~Time, data = gdf_PC_age, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_PC_age)

# Graphical representation #

cols_age_PC <- c("Pragian" = "#e9ba6a", "Emsian" = "#e9ba6a", "Eifelian" = "#F1D576", "Givetian" = "#F1E185", "Frasnian" = "#F2EDAD", "Famennian" = "#F2EDB3" , "Serpukhovian" = "#BFC26B", "Moscovian" = "#B3CBB9") #color codes following the international chronostratigraphic chart 

# figure of the variation of the disparity through time
SoV_PC.age <- c("0.000000000", "0.000000000", "0.01880830", "0.03004748", "0.02827263", "0.03001024", "0.03176942", "0.000000000") # copy the SOV_PC_age results
Age <- c("Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Serpukhovian", "Moscovian")
SoV_PC.age <- as.numeric(as.character(SoV_PC.age))
Age <- as.factor(as.character(Age))

DF_SoV.PC.age <- data.frame(Age = Age, SoV = SoV_PC.age)
DF_SoV.PC.age # Create data frame Age x SoV

level_order <- factor(DF_SoV.PC.age$Age, level = c("Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Serpukhovian", "Moscovian"))

ggplot(DF_SoV.PC.age, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_age_PC, size=5)+
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.04))+
  xlab("Age") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))

# Contribution of each time bin : Pie chart

PPV_PC.age <- c("2.912969", "1.661333", "9.021609", "9.635993", "17.553428", "20.444546", "35.528856", "3.241266") # copy the SOV_mean_PC_age results
Age <- c("Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Serpukhovian", "Moscovian")
PPV_PC.age <- as.numeric(as.character(PPV_PC.age))
Age <- as.factor(as.character(Age))

DF_PPV.age <- data.frame(Age = Age, PPV = PPV_PC.age)
DF_PPV.age # Create data frame Age x Proportion of variance

level_order <- factor(DF_PPV.age$Age, level = c("Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Serpukhovian", "Moscovian"))

# Compute percentages
DF_PPV.age$fraction = DF_PPV.age$PPV / sum(DF_PPV.age$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.age$ymax = cumsum(DF_PPV.age$fraction)

# Compute the bottom of each rectangle
DF_PPV.age$ymin = c(0, head(DF_PPV.age$ymax, n=-1))

DF_PPV.age$labelPosition <- (DF_PPV.age$ymax + DF_PPV.age$ymin) / 2

# Compute a good label
DF_PPV.age$label <- paste0(DF_PPV.age$PPV)

# Make the plot
ggplot(DF_PPV.age, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Age)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_age_PC) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.age, aes(x="", y=PPV, fill=Age)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_age_PC) +
  theme_void() #Pie chart

# Sum of Procrustes variances : Epochs #########################

SOV_mean_PC_epoch <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_PC_epoch, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_PC_epoch)

SOV_PC_epoch <- morphol.disparity(coords~Time,groups=~Time, data = gdf_PC_epoch, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_PC_epoch)

# figure of the variation of the disparity through time

SoV_PC.epoch <- c("0.0122252", "0.02508597", "0.03504269", "0.03176942", "0.00000000") # copy the SOV_PC_epoch results
Epoch <- c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
SoV_PC.epoch <- as.numeric(as.character(SoV_PC.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_SoV.PC.epoch <- data.frame(Epoch = Epoch, SoV = SoV_PC.epoch)
DF_SoV.PC.epoch # Create data frame Epoch x SoV 

level_order <- factor(DF_SoV.PC.epoch$Epoch, level = c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

ggplot(DF_SoV.PC.epoch, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_epoch, size=5)+
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.04))+
  xlab("Epoch") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))

# Contribution of each time bin : Pie chart

PPV_PC.epoch <- c("4.6","18.7","38","35.5","3.2") # copy the SOV_mean_PC_epoch results
Epoch <- c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
PPV_PC.epoch <- as.numeric(as.character(PPV_PC.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_PPV.epoch <- data.frame(Epoch = Epoch, PPV = PPV_PC.epoch)
DF_PPV.epoch # Create data frame Epoch x Proportion of variance

level_order <- factor(DF_PPV.epoch$Epoch, level = c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

# Compute percentages
DF_PPV.epoch$fraction = DF_PPV.epoch$PPV / sum(DF_PPV.epoch$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.epoch$ymax = cumsum(DF_PPV.epoch$fraction)

# Compute the bottom of each rectangle
DF_PPV.epoch$ymin = c(0, head(DF_PPV.epoch$ymax, n=-1))

DF_PPV.epoch$labelPosition <- (DF_PPV.epoch$ymax + DF_PPV.epoch$ymin) / 2

# Compute a good label
DF_PPV.epoch$label <- paste0(DF_PPV.epoch$PPV)

# Make the plot
ggplot(DF_PPV.epoch, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Epoch)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_epoch) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.epoch, aes(x="", y=PPV, fill=Epoch)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_epoch) +
  theme_void() #Pie chart

# Sum of Procrustes variances : Aquatic habitats #########################

SOV_mean_PC_ah <- morphol.disparity(coords ~ 1, groups= ~ Habitat, partial = TRUE, 
                                    data = gdf_PC_habitat, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_PC_ah)

SOV_PC_ah <- morphol.disparity(coords~Habitat,groups=~Habitat, data = gdf_PC_habitat, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_PC_ah)

# Graphical representation #

# Contribution of each time bin : Pie chart

PPV_PC.ah <- c("14.67022", "49.04465", "36.28513") # copy the SOV_mean_PC_ah results
Habitat <- c("Estuary", "Freshwater", "Marine")
PPV_PC.ah <- as.numeric(as.character(PPV_PC.ah))
Habitat <- as.factor(as.character(Habitat))

DF_PPV.ah <- data.frame(Habitat = Habitat, PPV = PPV_PC.ah)
DF_PPV.ah # Create data frame Habitat x Proportion of variance

# Compute percentages
DF_PPV.ah$fraction = DF_PPV.ah$PPV / sum(DF_PPV.ah$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.ah$ymax = cumsum(DF_PPV.ah$fraction)

# Compute the bottom of each rectangle
DF_PPV.ah$ymin = c(0, head(DF_PPV.ah$ymax, n=-1))

DF_PPV.ah$labelPosition <- (DF_PPV.ah$ymax + DF_PPV.ah$ymin) / 2

# Compute a good label
DF_PPV.ah$label <- paste0(DF_PPV.ah$PPV)

# Make the plot
ggplot(DF_PPV.ah, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Habitat)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_habitats) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.ah, aes(x="", y=PPV, fill=Habitat)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_habitats) +
  theme_void() #Pie chart

# Sum of Procrustes variances :  groups #########################

SOV_mean_PC_gp <- morphol.disparity(coords ~ 1, groups= ~ Group, partial = TRUE, 
                                    data = gdf_PC_gp, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_PC_gp)

SOV_PC_gp <- morphol.disparity(coords~Group,groups=~Group, data = gdf_PC_gp, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_PC_gp)

# Graphical representation #

cols <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Dipnoi" = "#FA8072", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F")

# Contribution of each time bin : Pie chart

PPV_PC.gp <- c("43.890", "22.945", "2.518", "16.966", "9.710", "3.971") # copy the SOV_mean_PC_gp results
Gp <- c("Actinistia", "Dipnoi", "Onychodontida", "Osteolepiform", "Porolepiform", "Rhizodontida")
PPV_PC.gp <- as.numeric(as.character(PPV_PC.gp))
Gp <- as.factor(as.character(Gp))

DF_PPV.gp <- data.frame(Group = Gp, PPV = PPV_PC.gp)
DF_PPV.gp # Create data frame Group x Proportion of variance

# Compute percentages
DF_PPV.gp$fraction = DF_PPV.gp$PPV / sum(DF_PPV.gp$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.gp$ymax = cumsum(DF_PPV.gp$fraction)

# Compute the bottom of each rectangle
DF_PPV.gp$ymin = c(0, head(DF_PPV.gp$ymax, n=-1))

DF_PPV.gp$labelPosition <- (DF_PPV.gp$ymax + DF_PPV.gp$ymin) / 2

# Compute a good label
DF_PPV.gp$label <- paste0(DF_PPV.gp$PPV)

# Make the plot
ggplot(DF_PPV.gp, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Gp)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.gp, aes(x="", y=PPV, fill=Gp)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols) +
  theme_void() #Pie chart


##############################################################################################################
############################################## CHEEK DISPARITY ###############################################
##############################################################################################################

######################################### Import datasets ################################################

Sarco_Cheek <- readland.tps("Cheek-final.TPS",specID="imageID",negNA = TRUE,warnmsg = TRUE)

dim(Sarco_Cheek)

name_Cheek <- dimnames(Sarco_Cheek)[[3]]
name_Cheek

Period_sarco_Cheek <- read_excel(file.choose(), 1) # Open List-species-cheek.xlsx
Period_sarco_Cheek

Period_sarco_Cheek<- Period_sarco_Cheek %>%
  mutate_if(is.character, as.factor)

################################ GPA and PCA ##################################

# Procrustes superimposition #
data.super_Cheek <- gpagen(Sarco_Cheek, ProcD = FALSE)
attributes(data.super_Cheek)

plot(data.super_Cheek) 

# Principal component analyses #
pca_Cheek <- gm.prcomp(data.super_Cheek$coords)
pca_Cheek

plot_Cheek <- plot(pca_Cheek, axis1 = 1, axis2 = 2)
text(pca_Cheek[["x"]][,1], pca_Cheek[["x"]][,2], labels = name_Cheek) # PCA 1 vs 2

plot_PC2_Cheek <- plot(pca_Cheek, axis1=2, axis2=3)
text(pca_Cheek[["x"]][,2], pca_Cheek[["x"]][,3], labels = name_Cheek) # PCA 2 vs 3

# Saving PC scores #

PC.scores_Cheek <- pca_Cheek$x 
as.data.frame(PC.scores_Cheek) # Save PC scores as a data frame object

write.csv(PC.scores_Cheek,"PC.scores_Cheek.csv",row.names=TRUE) # Save PC scores as a csv file

# Extracting eigenvectors #

ev_Cheek <- pca_Cheek$rotation
write.csv(ev_Cheek, "eigenvectors-Cheek.csv")

##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_Cheek), aes(x=PC.scores_Cheek[,1], y=PC.scores_Cheek[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_Cheek$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 40.25 %", y = "PC2 = 15.48 %" ) 

#PC2 vs PC3 -- color epoch
ggplot(as.data.frame(PC.scores_Cheek), aes(x=PC.scores_Cheek[,2], y=PC.scores_Cheek[,3], label = name_Cheek, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_Cheek$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.6,0.5), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 15.48 %", y = "PC3 = 9.40 %" )

#PC1 vs PC2 -- color Aquatic habitats 
ggplot(as.data.frame(PC.scores_Cheek), aes(x=PC.scores_Cheek[,1], y=PC.scores_Cheek[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Aquatic habitats") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_Cheek$`Aquatic habitat`), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 40.25 %", y = "PC2 = 15.48 %" ) 

#PC2 vs PC3 -- color Aquatic habitats  
ggplot(as.data.frame(PC.scores_Cheek), aes(x=PC.scores_Cheek[,2], y=PC.scores_Cheek[,3], label = name_Cheek, fontface = "italic")) +
  labs(fill="Aquatic habitat") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_Cheek$`Aquatic habitat`), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.6,0.5), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 15.48 %", y = "PC3 = 9.40 %" )

#PC1 vs PC2 -- color groups
ggplot(as.data.frame(PC.scores_Cheek), aes(x=PC.scores_Cheek[,1], y=PC.scores_Cheek[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_Cheek$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 40.25 %", y = "PC2 = 15.48 %" ) 

#PC2 vs PC3 -- color groups
ggplot(as.data.frame(PC.scores_Cheek), aes(x=PC.scores_Cheek[,2], y=PC.scores_Cheek[,3], label = name_Cheek, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_Cheek$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.6,0.5), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 15.48 %", y = "PC3 = 9.40 %" )

################### Calculating convex hulls for the PC1 vs PC2 #######################

# Add classifier variables (group, age, habitats) by hand to .csv file
# loading .csv file with PC scores and classifiers:
pca.class <- read.csv(file.choose()) #Open PC.scores_Cheek-class.csv

# for the epochs
CHull_epoch_Cheek <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Epoch")

CHull_epoch_Cheek.table <- CHull_epoch_Cheek$table
CHull_epoch_Cheek$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Epoch)) +
  labs(fill="Epoch") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Epoch), shape = 21, size = 5, stroke = 0) + 
  scale_fill_manual(values = cols_epoch) +
  scale_color_manual(values= cols_age)+
  lims(x=c(-0.7,0.6), y = c(-0.5,0.4)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 40.25 %", y = "PC2 = 15.48 %" ) +
  geom_polygon(data = CHull_epoch_Cheek.table, aes(x = X, y = Y), alpha = 0.5)

# for the Aquatic habitats 
CHull_habitat_Cheek <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Habitat")

CHull_habitat_Cheek.table <- CHull_habitat_Cheek$table
CHull_habitat_Cheek$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Habitat)) +
  labs(fill="Habitat") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Habitat), shape = 21, size = 5, stroke = 0.10) + 
  scale_fill_manual(values = cols_habitats) +
  scale_color_manual(values= cols_habitats)+
  lims(x=c(-0.7,0.6), y = c(-0.5,0.4)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 40.25 %", y = "PC2 = 15.48 %" ) +
  geom_polygon(data = CHull_habitat_Cheek.table, aes(x = X, y = Y), alpha = 0.5)

# for the groups
CHull_group_Cheek <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Group")

CHull_group_Cheek.table <- CHull_group_Cheek$table
CHull_group_Cheek$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Group, shape = Group)) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Group), shape = 21, size = 5, stroke = 0.10) + 
  scale_fill_manual(values = cols_group) +
  scale_color_manual(values= cols_group)+
  lims(x=c(-0.7,0.6), y = c(-0.5,0.4)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 40.25 %", y = "PC2 = 15.48 %" ) +
  geom_polygon(data = CHull_group_Cheek.table, aes(x = X, y = Y), alpha = 0.5)


################################# Morphological disparity analyses ###########################################

gdf_Cheek_epoch <- geomorph.data.frame(data.super_Cheek, species = Period_sarco_Cheek$Species, Time = Period_sarco_Cheek$Epoque) #gdf for the epochs

gdf_Cheek_age <- geomorph.data.frame(data.super_Cheek, species = Period_sarco_Cheek$Species, Time = Period_sarco_Cheek$Age) #gdf for the ages

gdf_Cheek_habitat <- geomorph.data.frame(data.super_Cheek, species = Period_sarco_Cheek$Species, Habitat = Period_sarco_Cheek$`Aquatic habitat`) #gdf for the aquatic habitats

gdf_Cheek_gp <- geomorph.data.frame(data.super_Cheek, species = Period_sarco_Cheek$Species, Group = Period_sarco_Cheek$Group) #gdf for the groups

# Sum of Procrustes variances : age ######################

SOV_mean_Cheek_age <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_Cheek_age, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_Cheek_age)

SOV_Cheek_age <- morphol.disparity(coords~Time,groups=~Time, data = gdf_Cheek_age, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_Cheek_age)

# figure of the variation of the disparity through time
SoV_Cheek.age <- c("0.00000000", "0.00000000", "0.02637913", "0.09050480", "0.06943954", "0.05843217", "0.03526274", "0.02398637", "0.07960483", "0.00000000") # copy the SOV_Cheek_age results
Age <- c("Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Visean", "Serpukhovian","Bashkirian")
SoV_Cheek.age <- as.numeric(as.character(SoV_Cheek.age))
Age <- as.factor(as.character(Age))

DF_SoV.Cheek.age <- data.frame(Age = Age, SoV = SoV_Cheek.age)
DF_SoV.Cheek.age # Create data frame Age x SoV

level_order <- factor(DF_SoV.Cheek.age$Age, level = c("Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Visean", "Serpukhovian","Bashkirian"))

cols_Cheek.age <- c("Pragian" = "#E5C468","Emsian" = "#e9ba6a", "Eifelian" = "#F1D576", "Givetian" = "#F1E185", "Frasnian" = "#F2EDAD", "Famennian" = "#F2EDB3" , "Tournaisian" = "#8CB06C", "Visean"= "#A6B96C", "Serpukhovian" = "#BFC26B", "Bashkirian" = "#99C2B5") #color codes following the international chronostratigraphic chart 

ggplot(DF_SoV.Cheek.age, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_Cheek.age, size=5)+
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.1))+
  xlab("Age") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))

# Contribution of each time bin: Pie chart

PPV_Cheek.age <- c("1.21", "1.45", "4.52", "20.56", "22.24", "14.87", "2.77", "2.40", "27.05", "2.93")
Age <- c("Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Visean", "Serpukhovian", "Bashkirian") # copy the SoV_mean_Cheek_age results
PPV_Cheek.age <- as.numeric(as.character(PPV_Cheek.age))
Age <- as.factor(as.character(Age))

DF_PPV.Cheek.age <- data.frame(Age = Age, PPV = PPV_Cheek.age)
DF_PPV.Cheek.age # Create data frame Age x Proportion of variance

level_order <- factor(DF_PPV.Cheek.age$Age, level = c("Pragian", "Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Visean", "Serpukhovian", "Bashkirian"))

# Compute percentages
DF_PPV.Cheek.age$fraction = DF_PPV.Cheek.age$PPV / sum(DF_PPV.Cheek.age$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.Cheek.age$ymax = cumsum(DF_PPV.Cheek.age$fraction)

# Compute the bottom of each rectangle
DF_PPV.Cheek.age$ymin = c(0, head(DF_PPV.Cheek.age$ymax, n=-1))

DF_PPV.Cheek.age$labelPosition <- (DF_PPV.Cheek.age$ymax + DF_PPV.Cheek.age$ymin) / 2

# Compute a good label
DF_PPV.Cheek.age$label <- paste0(DF_PPV.Cheek.age$PPV)

# Make the plot
ggplot(DF_PPV.Cheek.age, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Age)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_Cheek.age) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.Cheek.age, aes(x="", y=PPV, fill=Age)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_Cheek.age) +
  theme_void() #Pie chart

# Sum of Procrustes variances : Epochs #########################

SOV_mean_Cheek_epoch <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_Cheek_epoch, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_Cheek_epoch)

SOV_Cheek_epoch <- morphol.disparity(coords~Time,groups=~Time, data = gdf_Cheek_epoch, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells which time bin is most/least disparate
summary(SOV_Cheek_epoch)

# figure of the variation of the disparity through time

SoV_Cheek.epoch <- c("0.008206276", "0.07137672", "0.07021827", "0.1003253", "0.00000000") # copy the SOV_Cheek_epoch results
Epoch <- c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
SoV_Cheek.epoch <- as.numeric(as.character(SoV_Cheek.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_SoV.Cheek.epoch <- data.frame(Epoch = Epoch, SoV = SoV_Cheek.epoch)
DF_SoV.Cheek.epoch # Create data frame Epoch x SoV 

level_order <- factor(DF_SoV.Cheek.epoch$Epoch, level = c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

ggplot(DF_SoV.Cheek.epoch, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_epoch, size=5) +
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.11))+
  xlab("Epoch") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))

# Contribution of each time bin : Pie chart

PPV_Cheek.epoch <- c("2.661","25.081","37.107","32.222","2.923") # copy the SOV_mean_Cheek_epoch results
Epoch <- c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
PPV_Cheek.epoch <- as.numeric(as.character(PPV_Cheek.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_PPV.epoch <- data.frame(Epoch = Epoch, PPV = PPV_Cheek.epoch)
DF_PPV.epoch # Create data frame Epoch x Proportion of variance

level_order <- factor(DF_PPV.epoch$Epoch, level = c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

# Compute percentages
DF_PPV.epoch$fraction = DF_PPV.epoch$PPV / sum(DF_PPV.epoch$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.epoch$ymax = cumsum(DF_PPV.epoch$fraction)

# Compute the bottom of each rectangle
DF_PPV.epoch$ymin = c(0, head(DF_PPV.epoch$ymax, n=-1))

DF_PPV.epoch$labelPosition <- (DF_PPV.epoch$ymax + DF_PPV.epoch$ymin) / 2

# Compute a good label
DF_PPV.epoch$label <- paste0(DF_PPV.epoch$PPV)

# Make the plot
ggplot(DF_PPV.epoch, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Epoch)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_epoch) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.epoch, aes(x="", y=PPV, fill=Epoch)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_epoch) +
  theme_void() #Pie chart

# Sum of Procrustes variances : Aquatic habitats #########################

SOV_mean_Cheek_ah <- morphol.disparity(coords ~ 1, groups= ~ Habitat, partial = TRUE, 
                                    data = gdf_Cheek_habitat, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_Cheek_ah)

SOV_Cheek_ah <- morphol.disparity(coords~Habitat,groups=~Habitat, data = gdf_Cheek_habitat, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_Cheek_ah)

# Graphical representation #

# Contribution of each time bin : Pie chart

PPV_Cheek.ah <- c("15.68", "31.59", "52.73") # copy the SOV_mean_Cheek_ah results
Habitat <- c("Estuary", "Freshwater", "Marine")
PPV_Cheek.ah <- as.numeric(as.character(PPV_Cheek.ah))
Habitat <- as.factor(as.character(Habitat))

DF_PPV.ah <- data.frame(Habitat = Habitat, PPV = PPV_Cheek.ah)
DF_PPV.ah # Create data frame Habitat x Proportion of variance

# Compute percentages
DF_PPV.ah$fraction = DF_PPV.ah$PPV / sum(DF_PPV.ah$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.ah$ymax = cumsum(DF_PPV.ah$fraction)

# Compute the bottom of each rectangle
DF_PPV.ah$ymin = c(0, head(DF_PPV.ah$ymax, n=-1))

DF_PPV.ah$labelPosition <- (DF_PPV.ah$ymax + DF_PPV.ah$ymin) / 2

# Compute a good label
DF_PPV.ah$label <- paste0(DF_PPV.ah$PPV)

# Make the plot
ggplot(DF_PPV.ah, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Habitat)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_habitats) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.ah, aes(x="", y=PPV, fill=Habitat)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_habitats) +
  theme_void() #Pie chart

# Sum of Procrustes variances :  groups #########################

SOV_mean_Cheek_gp <- morphol.disparity(coords ~ 1, groups= ~ Group, partial = TRUE, 
                                    data = gdf_Cheek_gp, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_Cheek_gp)

SOV_Cheek_gp <- morphol.disparity(coords~Group,groups=~Group, data = gdf_Cheek_gp, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_Cheek_gp)

# Graphical representation #

cols <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F", "Elpistostegalia" = "#ffd496")

# Contribution of each time bin : Pie chart

PPV_Cheek.gp <- c("47.93", "3.25", "11.35", "24.26", "10.37", "2.84") # copy the SOV_mean_Cheek_gp results
Gp <- c("Actinistia", "Elpistostegalia", "Onychodontida", "Osteolepiform", "Porolepiform", "Rhizodontida")
PPV_Cheek.gp <- as.numeric(as.character(PPV_Cheek.gp))
Gp <- as.factor(as.character(Gp))

DF_PPV.gp <- data.frame(Group = Gp, PPV = PPV_Cheek.gp)
DF_PPV.gp # Create data frame Group x Proportion of variance

# Compute percentages
DF_PPV.gp$fraction = DF_PPV.gp$PPV / sum(DF_PPV.gp$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.gp$ymax = cumsum(DF_PPV.gp$fraction)

# Compute the bottom of each rectangle
DF_PPV.gp$ymin = c(0, head(DF_PPV.gp$ymax, n=-1))

DF_PPV.gp$labelPosition <- (DF_PPV.gp$ymax + DF_PPV.gp$ymin) / 2

# Compute a good label
DF_PPV.gp$label <- paste0(DF_PPV.gp$PPV)

# Make the plot
ggplot(DF_PPV.gp, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Gp)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.gp, aes(x="", y=PPV, fill=Gp)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols) +
  theme_void() #Pie chart


##############################################################################################################
############################################## SKULL ROOF DISPARITY ##########################################
##############################################################################################################

######################################### Import datasets ################################################

Sarco_SR <- readland.tps("Skullroof-final.TPS",specID="imageID",negNA = TRUE,warnmsg = TRUE)

dim(Sarco_SR)

name_SR <- dimnames(Sarco_SR)[[3]]
name_SR

Period_sarco_SR <- read_excel(file.choose(), 1) # Open List-species-SR.xlsx
Period_sarco_SR

Period_sarco_SR<- Period_sarco_SR %>%
  mutate_if(is.character, as.factor)

################################### GPA and PCA ##################################

# Procrustes superimposition #

data.super_SR <- gpagen(Sarco_SR, ProcD = FALSE)
attributes(data.super_SR)

plot(data.super_SR) 

# Asymmetry #

nbb <- as.character(c(1:47))
data.super_SR$ind=nbb # adding an ind vector to gpagen

n_pairs <- c(1,2,3,4,5,6,7,8,9,10,11,12,14,16,15,17,18,21,19,22,20,23)
pairs_matrix <- matrix(n_pairs, ncol=2, byrow = TRUE) # match paired landmarks
pairs_matrix


sym <- bilat.symmetry(A = data.super_SR$coords, ind=name_SR, object.sym = TRUE, land.pairs = pairs_matrix, iter = 149)
summary(sym)

plot(sym$symm.shape[,c(1,2),8])

ppca<-gm.prcomp(sym$symm.shape)
mmshape<-mshape(sym$symm.shape)

sym$symm.shape


# Principal component analyses #

pca_sym.SR <- gm.prcomp(sym$symm.shape) #with symmetrized coordinates
pca_sym.SR
plot_sym <- plot(pca_sym.SR)

plot_SR <- plot(pca_sym.SR, axis1 = 1, axis2 = 2)
text(pca_sym.SR[["x"]][,1], pca_sym.SR[["x"]][,2], labels = name_SR) # PCA 1 vs 2

plot_PC2_SR <- plot(pca_sym.SR, axis1=2, axis2=3)
text(pca_sym.SR[["x"]][,2], pca_sym.SR[["x"]][,3], labels = name_SR) # PCA 2 vs 3

# Saving PC scores #

PC.scores_SR <- pca_sym.SR$x 
as.data.frame(PC.scores_SR) # Save PC scores as a data frame object

write.csv(PC.scores_SR,"PC.scores_SR.csv",row.names=TRUE) # Save PC scores as a csv file

# Extracting eigenvectors #

ev_SR <- pca_sym.SR$rotation
write.csv(ev_SR, "eigenvectors-SR.csv")

##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_SR), aes(x=PC.scores_SR[,1], y=PC.scores_SR[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_SR$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 52.46 %", y = "PC2 = 17.59 %" ) 

#PC2 vs PC3 -- color epoch
ggplot(as.data.frame(PC.scores_SR), aes(x=PC.scores_SR[,2], y=PC.scores_SR[,3], label = name_SR, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_SR$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.3,0.3), y = c(-0.2,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 17.59 %", y = "PC3 = 10.59 %" )

#PC1 vs PC2 -- color Aquatic habitats 
ggplot(as.data.frame(PC.scores_SR), aes(x=PC.scores_SR[,1], y=PC.scores_SR[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Aquatic habitats") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_SR$`Aquatic habitat`), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 52.46 %", y = "PC2 = 17.59 %" ) 

#PC2 vs PC3 -- color Aquatic habitats  
ggplot(as.data.frame(PC.scores_SR), aes(x=PC.scores_SR[,2], y=PC.scores_SR[,3], label = name_SR, fontface = "italic")) +
  labs(fill="Aquatic habitat") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_SR$`Aquatic habitat`), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.3,0.3), y = c(-0.2,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 17.59 %", y = "PC3 = 10.59 %" )

#PC1 vs PC2 -- color groups
ggplot(as.data.frame(PC.scores_SR), aes(x=PC.scores_SR[,1], y=PC.scores_SR[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_SR$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 52.46 %", y = "PC2 = 17.59 %" ) 

#PC2 vs PC3 -- color groups
ggplot(as.data.frame(PC.scores_SR), aes(x=PC.scores_SR[,2], y=PC.scores_SR[,3], label = name_SR, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sarco_SR$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.3,0.3), y = c(-0.2,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 17.59 %", y = "PC3 = 10.59 %" )

################### Calculating convex hulls for the PC1 vs PC2 #######################

# Add classifier variables (group, age, habitats) by hand to .csv file
# loading .csv file with PC scores and classifiers:
pca.class <- read.csv(file.choose()) #Open PC.scores_SR-class.csv

# for the epochs
CHull_epoch_SR <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Epoch")

CHull_epoch_SR.table <- CHull_epoch_SR$table
CHull_epoch_SR$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Epoch, shape = Epoch)) +
  labs(fill="Epoch") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Epoch), shape = 21, size = 5, stroke = 0.10) + 
  scale_fill_manual(values = cols_epoch) +
  scale_color_manual(values= cols_age)+
  lims(x=c(-0.7,0.7), y = c(-0.45,0.35))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 52.46 %", y = "PC2 = 17.59 %" ) +
  geom_polygon(data = CHull_epoch_SR.table, aes(x = X, y = Y), alpha = 0.5)

# for the Aquatic habitats 
CHull_habitat_SR <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Habitat")

CHull_habitat_SR.table <- CHull_habitat_SR$table
CHull_habitat_SR$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Habitat, shape = Habitat)) +
  labs(fill="Habitat") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Habitat), shape = 21, size = 5, stroke = 0.10) + 
  scale_fill_manual(values = cols_habitats) +
  scale_color_manual(values= cols_habitats)+
  lims(x=c(-0.7,0.7), y = c(-0.45,0.35))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 52.46 %", y = "PC2 = 17.59 %" ) +
  geom_polygon(data = CHull_habitat_SR.table, aes(x = X, y = Y), alpha = 0.5)

# for the groups
CHull_group_SR <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Group")

CHull_group_SR.table <- CHull_group_SR$table
CHull_group_SR$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Group)) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Group), shape = 21, size = 5, stroke = 0.10) + 
  scale_fill_manual(values = cols_group) +
  scale_color_manual(values= cols_group)+
  lims(x=c(-0.7,0.7), y = c(-0.45,0.35))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 52.46 %", y = "PC2 = 17.59 %" ) +
  geom_polygon(data = CHull_group_SR.table, aes(x = X, y = Y), alpha = 0.5)


################################# Morphological disparity analyses ###########################################

gdf_SR_epoch <- geomorph.data.frame(sym$symm.shape, species = Period_sarco_SR$Species, Time = Period_sarco_SR$Epoque) #gdf for the epochs
names(gdf_SR_epoch)<-c("coords", "species", "Time")

gdf_SR_age <- geomorph.data.frame(sym$symm.shape, species = Period_sarco_SR$Species, Time = Period_sarco_SR$Age) #gdf for the ages
names(gdf_SR_age)<-c("coords", "species", "Time")

gdf_SR_habitat <- geomorph.data.frame(sym$symm.shape, species = Period_sarco_SR$Species, Habitat = Period_sarco_SR$`Aquatic habitat`) #gdf for the aquatic habitats
names(gdf_SR_habitat)<-c("coords", "species", "Habitat")

gdf_SR_gp <- geomorph.data.frame(sym$symm.shape, species = Period_sarco_SR$Species, Group = Period_sarco_SR$Group) #gdf for the groups
names(gdf_SR_gp)<-c("coords", "species", "Group")

# Sum of Procrustes variances : age ######################

SOV_mean_SR_age <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_SR_age, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_SR_age)

SOV_SR_age <- morphol.disparity(coords~Time,groups=~Time, data = gdf_SR_age, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_SR_age)

# figure of the variation of the disparity through time
SoV_SR.age <- c("0.05657880", "0.04834485", "0.06864643", "0.06752854", "0.04274543", "0.00007766136", "0.06004026", "0.03116380", "0.00000000", "0.00000000", "0.00000000") # copy the SOV_SR_age results
Age <- c("Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Visean", "Serpukhovian", "Bashkirian", "Moscovian", "Gzhelian")
SoV_SR.age <- as.numeric(as.character(SoV_SR.age))
Age <- as.factor(as.character(Age))

DF_SoV.SR.age <- data.frame(Age = Age, SoV = SoV_SR.age)
DF_SoV.SR.age # Create data frame Age x SoV

level_order <- factor(DF_SoV.SR.age$Age, level = c("Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Visean", "Serpukhovian", "Moscovian", "Bashkirian", "Gzhelian"))

cols_SR.age <- c("Emsian" = "#e9ba6a", "Eifelian" = "#F1D576", "Givetian" = "#F1E185", "Frasnian" = "#F2EDAD", "Famennian" = "#F2EDB3" , "Tournaisian" = "#8CB06C", "Visean"= "#A6B96C", "Serpukhovian" = "#BFC26B", "Bashkirian" = "#99C2B5","Moscovian" = "#B3CBB9", "Gzhelian" = "#CCD4C7") #color codes following the international chronostratigraphic chart 

ggplot(DF_SoV.SR.age, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_SR.age, size=5)+
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.1))+
  xlab("Age") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))

# Contribution of each time bin: Pie chart

PPV_SR.age <- c("5.09", "5.11", "14.65", "34.13", "14.5", "3.49", "4.39", "9.04", "3.53", "3.17", "2.9")
Age <- c("Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Visean", "Serpukhovian", "Moscovian", "Bashkirian", "Gzhelian") # copy the SoV_mean_SR_age results
PPV_SR.age <- as.numeric(as.character(PPV_SR.age))
Age <- as.factor(as.character(Age))

DF_PPV.SR.age <- data.frame(Age = Age, PPV = PPV_SR.age)
DF_PPV.SR.age # Create data frame Age x Proportion of variance

level_order <- factor(DF_PPV.SR.age$Age, level = c("Emsian", "Eifelian", "Givetian", "Frasnian", "Famennian", "Tournaisian", "Visean", "Serpukhovian", "Moscovian", "Bashkirian", "Gzhelian"))

# Compute percentages
DF_PPV.SR.age$fraction = DF_PPV.SR.age$PPV / sum(DF_PPV.SR.age$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.SR.age$ymax = cumsum(DF_PPV.SR.age$fraction)

# Compute the bottom of each rectangle
DF_PPV.SR.age$ymin = c(0, head(DF_PPV.SR.age$ymax, n=-1))

DF_PPV.SR.age$labelPosition <- (DF_PPV.SR.age$ymax + DF_PPV.SR.age$ymin) / 2

# Compute a good label
DF_PPV.SR.age$label <- paste0(DF_PPV.SR.age$PPV)

# Make the plot
ggplot(DF_PPV.SR.age, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Age)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_SR.age) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.SR.age, aes(x="", y=PPV, fill=Age)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_SR.age) +
  theme_void() #Pie chart

# Sum of Procrustes variances : Epochs #########################

SOV_mean_SR_epoch <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_SR_epoch, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_SR_epoch)

SOV_SR_epoch <- morphol.disparity(coords~Time,groups=~Time, data = gdf_SR_epoch, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells which time bin is most/least disparate
summary(SOV_SR_epoch)

# figure of the variation of the disparity through time

SoV_SR.epoch <- c("0.05657880", "0.06329561", "0.06069786", "0.04564711", "0.08464444") # copy the SOV_SR_epoch results
Epoch <- c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
SoV_SR.epoch <- as.numeric(as.character(SoV_SR.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_SoV.SR.epoch <- data.frame(Epoch = Epoch, SoV = SoV_SR.epoch)
DF_SoV.SR.epoch # Create data frame Epoch x SoV 

level_order <- factor(DF_SoV.SR.epoch$Epoch, level = c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

ggplot(DF_SoV.SR.epoch, aes(x=level_order, y=SoV, group = 1)) +
  geom_line(color = "black", linewidth=1) +
  geom_point(color = cols_epoch, size=5) +
  scale_y_continuous(labels = label_number(accuracy = 0.001),limits = c(0,0.11))+
  xlab("Epoch") +
  ylab("Procrustes variance")  +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(color="black", size=20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(angle = 90, hjust=0.5, size = 14))

# Contribution of each time bin : Pie chart

PPV_SR.epoch <- c("5.09","19.76","48.63","16.92","9.6") # copy the SOV_mean_SR_epoch results
Epoch <- c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian")
PPV_SR.epoch <- as.numeric(as.character(PPV_SR.epoch))
Epoch <- as.factor(as.character(Epoch))

DF_PPV.epoch <- data.frame(Epoch = Epoch, PPV = PPV_SR.epoch)
DF_PPV.epoch # Create data frame Epoch x Proportion of variance

level_order <- factor(DF_PPV.epoch$Epoch, level = c("Lower Devonian", "Middle Devonian", "Upper Devonian", "Mississippian", "Pennsylvanian"))

# Compute percentages
DF_PPV.epoch$fraction = DF_PPV.epoch$PPV / sum(DF_PPV.epoch$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.epoch$ymax = cumsum(DF_PPV.epoch$fraction)

# Compute the bottom of each rectangle
DF_PPV.epoch$ymin = c(0, head(DF_PPV.epoch$ymax, n=-1))

DF_PPV.epoch$labelPosition <- (DF_PPV.epoch$ymax + DF_PPV.epoch$ymin) / 2

# Compute a good label
DF_PPV.epoch$label <- paste0(DF_PPV.epoch$PPV)

# Make the plot
ggplot(DF_PPV.epoch, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Epoch)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_epoch) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.epoch, aes(x="", y=PPV, fill=Epoch)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_epoch) +
  theme_void() #Pie chart

# Sum of Procrustes variances : Aquatic habitats #########################

SOV_mean_SR_ah <- morphol.disparity(coords ~ 1, groups= ~ Habitat, partial = TRUE, 
                                       data = gdf_SR_habitat, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_SR_ah)

SOV_SR_ah <- morphol.disparity(coords~Habitat,groups=~Habitat, data = gdf_SR_habitat, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_SR_ah)

# Graphical representation #

# Contribution of each time bin : Pie chart

PPV_SR.ah <- c("15.60", "57.56", "26.84") # copy the SOV_mean_SR_ah results
Habitat <- c("Estuary", "Freshwater", "Marine")
PPV_SR.ah <- as.numeric(as.character(PPV_SR.ah))
Habitat <- as.factor(as.character(Habitat))

DF_PPV.ah <- data.frame(Habitat = Habitat, PPV = PPV_SR.ah)
DF_PPV.ah # Create data frame Habitat x Proportion of variance

# Compute percentages
DF_PPV.ah$fraction = DF_PPV.ah$PPV / sum(DF_PPV.ah$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.ah$ymax = cumsum(DF_PPV.ah$fraction)

# Compute the bottom of each rectangle
DF_PPV.ah$ymin = c(0, head(DF_PPV.ah$ymax, n=-1))

DF_PPV.ah$labelPosition <- (DF_PPV.ah$ymax + DF_PPV.ah$ymin) / 2

# Compute a good label
DF_PPV.ah$label <- paste0(DF_PPV.ah$PPV)

# Make the plot
ggplot(DF_PPV.ah, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Habitat)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_habitats) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.ah, aes(x="", y=PPV, fill=Habitat)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_habitats) +
  theme_void() #Pie chart

# Sum of Procrustes variances :  groups #########################

SOV_mean_SR_gp <- morphol.disparity(coords ~ 1, groups= ~ Group, partial = TRUE, 
                                       data = gdf_SR_gp, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_SR_gp)

SOV_SR_gp <- morphol.disparity(coords~Group,groups=~Group, data = gdf_SR_gp, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_SR_gp)

# Graphical representation #
cols.SR <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Dipnoi" = "#FA8072", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F", "Elpistostegalia" = "#FFEBCD", "Tetrapoda" = "#D2B48C") #Assign each group a color


# Contribution of each time bin : Pie chart

PPV_SR.gp <- c("12.33", "28.85", "3.57", "8.50", "14.44", "8.66", "5.28", "18.37") # copy the SOV_mean_SR_gp results
Gp <- c("Actinistia", "Dipnoi", "Elpistostegalia", "Onychodontida", "Osteolepiform", "Porolepiform", "Rhizodontida", "Tetrapoda")
PPV_SR.gp <- as.numeric(as.character(PPV_SR.gp))
Gp <- as.factor(as.character(Gp))

DF_PPV.gp <- data.frame(Group = Gp, PPV = PPV_SR.gp)
DF_PPV.gp # Create data frame Group x Proportion of variance

# Compute percentages
DF_PPV.gp$fraction = DF_PPV.gp$PPV / sum(DF_PPV.gp$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.gp$ymax = cumsum(DF_PPV.gp$fraction)

# Compute the bottom of each rectangle
DF_PPV.gp$ymin = c(0, head(DF_PPV.gp$ymax, n=-1))

DF_PPV.gp$labelPosition <- (DF_PPV.gp$ymax + DF_PPV.gp$ymin) / 2

# Compute a good label
DF_PPV.gp$label <- paste0(DF_PPV.gp$PPV)

# Make the plot
ggplot(DF_PPV.gp, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Gp)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols.SR) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.gp, aes(x="", y=PPV, fill=Gp)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols.SR) +
  theme_void() #Pie chart

##############################################################################################
############################## PRECISE PALEOENVIRONMENTS #####################################
##############################################################################################

# we removed species for which no precise information about paleoenvironment is available

##############################################################################################################
########################################### FULLBODY DISPARITY ############################################
##############################################################################################################

######################################### Import datasets ################################################

# TPS files with landmarks and semilandmarks
Sarco_PC <- readland.tps("Postcranial-final.TPS",specID="imageID",negNA = TRUE,readcurves = TRUE,warnmsg = TRUE)
dim(Sarco_PC)
dimnames(Sarco_PC)[[3]] 
Sarco_PC_pal <- Sarco_PC[,,-2][,,-17][,,-23]

name_PC_pal <- dimnames(Sarco_PC_pal)[[3]]
name_PC_pal

# Excel file with the list of species and their ages
Period_sarco_PC <- read_excel(file.choose(), 1) #Open List-species-PC.xlsx
Period_sarco_PC

Period_sarco_PC<- Period_sarco_PC %>%
  mutate_if(is.character, as.factor)

Period_PC_pal <- Period_sarco_PC %>%  filter(!row_number() %in% c(2,18,25))

# Convert curves into semi-landmarks
matrice_PC_pal <- rbind(define.sliders(12:41), define.sliders(42:61), define.sliders(62:71), define.sliders(72:107), define.sliders(108:207))

########################################## GPA and PCA ####################################################

# Procrustes superimposition #
data.super_PC_pal <- gpagen(Sarco_PC_pal, curves = matrice_PC_pal, ProcD = FALSE)
attributes(data.super_PC_pal)

plot(data.super_PC_pal) 

# Principal component analyses #
pca_PC_pal <- gm.prcomp(data.super_PC_pal$coords)
pca_PC_pal

plot_PC_pal <- plot(pca_PC_pal, axis1 = 1, axis2 = 2)
text(pca_PC_pal[["x"]][,1], pca_PC_pal[["x"]][,2], labels = name_PC_pal) # PCA 1 vs 2

plot_PC2_PC_pal <- plot(pca_PC_pal, axis1=2, axis2=3)
text(pca_PC_pal[["x"]][,2], pca_PC_pal[["x"]][,3], labels = name_PC_pal) # PCA 2 vs 3

# Saving PC scores #

PC.scores_PC_pal <- pca_PC_pal$x 
as.data.frame(PC.scores_PC_pal) # Save PC scores as a data frame object

write.csv(PC.scores_PC_pal,"PC.scores_PC_pal.csv",row.names=TRUE) # Save PC scores as a csv file

##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_PC_pal), aes(x=PC.scores_PC_pal[,1], y=PC.scores_PC_pal[,2], label = name_PC_pal, fontface = "italic")) +
  labs(fill="Paleoenvironments") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_PC_pal$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 52.32 %", y = "PC2 = 13.19 %" ) 

#PC2 vs PC3 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_PC_pal), aes(x=PC.scores_PC_pal[,2], y=PC.scores_PC_pal[,3], label = name_PC_pal, fontface = "italic")) +
  labs(fill="Paleoenvironments") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_PC_pal$Paleoenvironments), color = "black", size = 6, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 13.19 %", y = "PC3 = 8.84 %" )


################### Calculating convex hulls for the PC1 vs PC2 #######################

# Add classifier variables (group, age, habitats) by hand to .csv file
# loading .csv file with PC scores and classifiers:
pca.class <- read.csv(file.choose()) #Open PC.scores_PC_pal.csv

# for the more precise paleoenvironments
CHull_habitat2_PC <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Habitat2")

CHull_habitat2_PC.table <- CHull_habitat2_PC$table
CHull_habitat2_PC$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Habitat2, shape = Habitat2)) +
  labs(fill="Habitat2") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Habitat2), shape = 21, size = 3, stroke = 0.10) + 
  scale_fill_manual(values = cols_habitats2) +
  scale_color_manual(values= cols_habitats2)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 52.32 %", y = "PC2 = 13.19 %") +
  geom_polygon(data = CHull_habitat2_PC.table, aes(x = X, y = Y), alpha = 0.5)


################################# Morphological disparity analyses ###########################################

gdf_PC_habitat2 <- geomorph.data.frame(data.super_PC_pal, species = Period_PC_pal$Species, Habitat = Period_PC_pal$Paleoenvironments) #gdf for the more precise paleoenvironments

# Sum of Procrustes variances : More precise aquatic habitats #########################

SOV_mean_PC_pal <- morphol.disparity(coords ~ 1, groups= ~ Habitat, partial = TRUE, 
                                     data = gdf_PC_habitat2, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_PC_pal)

SOV_PC_pal <- morphol.disparity(coords~Habitat,groups=~Habitat, data = gdf_PC_habitat2, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_PC_pal)

# Graphical representation #

# Contribution of each time bin : Pie chart

PPV_PC.pal <- c("5.184","33.539","13.879","3.355","8.415","16.212","5.170","14.246") # copy the SOV_mean_PC_pal results
Habitat2 <- c("Alluvial plain", "Bay", "Calm  freshwater lake", "Delta", "Dynamic freshwater lake", "Estuary", "Lagoon", "Oxbow lake/meandering river")
PPV_PC.pal <- as.numeric(as.character(PPV_PC.pal))
Habitat2 <- as.factor(as.character(Habitat2))

DF_PPV.pal <- data.frame(Habitats = Habitat2, PPV = PPV_PC.pal)
DF_PPV.pal # Create data frame Paleoenvironment x Proportion of variance

# Compute percentages
DF_PPV.pal$fraction = DF_PPV.pal$PPV / sum(DF_PPV.pal$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.pal$ymax = cumsum(DF_PPV.pal$fraction)

# Compute the bottom of each rectangle
DF_PPV.pal$ymin = c(0, head(DF_PPV.pal$ymax, n=-1))

DF_PPV.pal$labelPosition <- (DF_PPV.pal$ymax + DF_PPV.pal$ymin) / 2

# Compute a good label
DF_PPV.pal$label <- paste0(DF_PPV.pal$PPV)

# Make the plot
ggplot(DF_PPV.pal, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Habitat2)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_habitats2) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.pal, aes(x="", y=PPV, fill=Habitat2)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_habitats2) +
  theme_void() #Pie chart


##############################################################################################################
############################################## CHEEK DISPARITY ###############################################
##############################################################################################################

######################################### Import datasets ################################################

Sarco_Cheek <- readland.tps("Cheek-final.TPS",specID="imageID",negNA = TRUE,warnmsg = TRUE)
dim(Sarco_Cheek)

Sarco_Cheek_pal <- Sarco_Cheek[,,-2][,,-2][,,-2][,,-8][,,-10][,,-20][,,-28][,,-28]
name_Cheek_pal <- dimnames(Sarco_Cheek_pal)[[3]]
name_Cheek_pal

Period_sarco_Cheek <- read_excel(file.choose(), 1) # Open List-species-cheek.xlsx
Period_sarco_Cheek

Period_sarco_Cheek<- Period_sarco_Cheek %>%
  mutate_if(is.character, as.factor)
Period_Cheek_pal <- Period_sarco_Cheek %>%  filter(!row_number() %in% c(2,3,4,11,14,25,34,35))

################################ GPA and PCA ##################################

# Procrustes superimposition #
data.super_Cheek_pal <- gpagen(Sarco_Cheek_pal, ProcD = FALSE)
attributes(data.super_Cheek_pal)

plot(data.super_Cheek_pal) 

# Principal component analyses #
pca_Cheek_pal <- gm.prcomp(data.super_Cheek_pal$coords)
pca_Cheek_pal

plot_Cheek_pal <- plot(pca_Cheek_pal, axis1 = 1, axis2 = 2)
text(pca_Cheek_pal[["x"]][,1], pca_Cheek_pal[["x"]][,2], labels = name_Cheek_pal) # PCA 1 vs 2

plot_PC2_Cheek_pal <- plot(pca_Cheek_pal, axis1=2, axis2=3)
text(pca_Cheek_pal[["x"]][,2], pca_Cheek_pal[["x"]][,3], labels = name_Cheek_pal) # PCA 2 vs 3

# Saving PC scores #

PC.scores_Cheek_pal <- pca_Cheek_pal$x 
as.data.frame(PC.scores_Cheek_pal) # Save PC scores as a data frame object

write.csv(PC.scores_Cheek_pal,"PC.scores_Cheek_pal.csv",row.names=TRUE) # Save PC scores as a csv file

##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_Cheek_pal), aes(x=PC.scores_Cheek_pal[,1], y=PC.scores_Cheek_pal[,2], label = name_Cheek_pal, fontface = "italic")) +
  labs(fill="Paleoenvironments") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Cheek_pal$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 45.11 %", y = "PC2 = 16.04 %" ) 

#PC2 vs PC3 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_Cheek_pal), aes(x=PC.scores_Cheek_pal[,2], y=PC.scores_Cheek_pal[,3], label = name_Cheek_pal, fontface = "italic")) +
  labs(fill="Paleoenvironments") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Cheek_pal$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.6,0.5), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 16.04 %", y = "PC3 = 9.68 %" )

################### Calculating convex hulls for the PC1 vs PC2 #######################

# Add classifier variables (group, age, habitats) by hand to .csv file
# loading .csv file with PC scores and classifiers:
pca.class <- read.csv(file.choose()) #Open PC.scores_Cheek_pal.csv

# for the more precise paleoenvironments
CHull_habitat2_Cheek <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Habitat2")

CHull_habitat2_Cheek.table <- CHull_habitat2_Cheek$table
CHull_habitat2_Cheek$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Habitat2)) +
  labs(fill="Habitat2") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Habitat2), shape = 21, size = 5, stroke = 0.10) + 
  scale_fill_manual(values = cols_habitats2) +
  scale_color_manual(values= cols_habitats2)+
  lims(x=c(-0.7,0.6), y = c(-0.5,0.4)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 45.11 %", y = "PC2 = 16.04 %" ) +
  geom_polygon(data = CHull_habitat2_Cheek.table, aes(x = X, y = Y), alpha = 0.5)

################################# Morphological disparity analyses ###########################################

gdf_Cheek_habitat2 <- geomorph.data.frame(data.super_Cheek_pal, species = Period_Cheek_pal$Species, Habitat = Period_Cheek_pal$Paleoenvironments) #gdf for the more precise paleoenvironments

# Sum of Procrustes variances : More precise aquatic habitats #########################

SOV_mean_Cheek_pal <- morphol.disparity(coords ~ 1, groups= ~ Habitat, partial = TRUE, 
                                        data = gdf_Cheek_habitat2, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_Cheek_pal)

SOV_Cheek_pal <- morphol.disparity(coords~Habitat,groups=~Habitat, data = gdf_Cheek_habitat2, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells which time bin is most/least disparate

summary(SOV_Cheek_pal)

# Graphical representation #

# Contribution of each time bin : Pie chart

PPV_Cheek.pal <- c("7.063","30.340","5.651","1.675","1.141","18.406","21.131","4.925","9.668") # copy the SOV_mean_Cheek_pal results
Habitat2 <- c("Alluvial plain", "Bay", "Calm  freshwater lake", "Coastal marine","Dynamic freshwater lake", "Estuary", "Lagoon", "Oxbow lake/meandering river", "Reef")
PPV_Cheek.pal <- as.numeric(as.character(PPV_Cheek.pal))
Habitat2 <- as.factor(as.character(Habitat2))

DF_PPV.pal <- data.frame(Habitats = Habitat2, PPV = PPV_Cheek.pal)
DF_PPV.pal # Create data frame Paleoenvironment x Proportion of variance

# Compute percentages
DF_PPV.pal$fraction = DF_PPV.pal$PPV / sum(DF_PPV.pal$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.pal$ymax = cumsum(DF_PPV.pal$fraction)

# Compute the bottom of each rectangle
DF_PPV.pal$ymin = c(0, head(DF_PPV.pal$ymax, n=-1))

DF_PPV.pal$labelPosition <- (DF_PPV.pal$ymax + DF_PPV.pal$ymin) / 2

# Compute a good label
DF_PPV.pal$label <- paste0(DF_PPV.pal$PPV)

# Make the plot
ggplot(DF_PPV.pal, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Habitat2)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_habitats2) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.pal, aes(x="", y=PPV, fill=Habitat2)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_habitats2) +
  theme_void() #Pie chart

##############################################################################################################
############################################## SKULL ROOF DISPARITY ##########################################
##############################################################################################################

######################################### Import datasets ################################################

Sarco_SR <- readland.tps("Skullroof-final.TPS",specID="imageID",negNA = TRUE,warnmsg = TRUE)
dim(Sarco_SR)

Sarco_SR_pal <- Sarco_SR[,,-6][,,-6][,,-6][,,-13][,,-24][,,-27][,,-37]
names_SR_pal <- dimnames(Sarco_SR_pal)[[3]]
names_SR_pal

Period_sarco_SR <- read_excel(file.choose(), 1) # Open List-species-SR.xlsx
Period_sarco_SR

Period_sarco_SR<- Period_sarco_SR %>%
  mutate_if(is.character, as.factor)
Period_SR_pal <- Period_sarco_SR %>% filter(!row_number() %in% c(6,7,8,16,28,32,43))

################################### GPA and PCA ##################################

# Procrustes superimposition #

data.super_SR_pal <- gpagen(Sarco_SR_pal, ProcD = FALSE)
attributes(data.super_SR_pal)

plot(data.super_SR_pal) 

# Asymmetry #

nbb <- as.character(c(1:47))
data.super_SR_pal$ind=nbb # adding an ind vector to gpagen

n_pairs <- c(1,2,3,4,5,6,7,8,9,10,11,12,14,16,15,17,18,21,19,22,20,23)
pairs_matrix <- matrix(n_pairs, ncol=2, byrow = TRUE) # match paired landmarks
pairs_matrix


sym_pal <- bilat.symmetry(A = data.super_SR_pal$coords, ind=names_SR_pal, object.sym = TRUE, land.pairs = pairs_matrix, iter = 149)
summary(sym_pal)

plot(sym_pal$symm.shape[,c(1,2),8])

ppca_pal<-gm.prcomp(sym_pal$symm.shape)
mmshape_pal<-mshape(sym_pal$symm.shape)

sym_pal$symm.shape


# Principal component analyses #

pca_sym.SR_pal <- gm.prcomp(sym_pal$symm.shape) #with symmetrized coordinates
pca_sym.SR_pal
plot_sym_pal <- plot(pca_sym.SR_pal)

plot_SR_pal <- plot(pca_sym.SR_pal, axis1 = 1, axis2 = 2)
text(pca_sym.SR_pal[["x"]][,1], pca_sym.SR_pal[["x"]][,2], labels = names_SR_pal) # PCA 1 vs 2

plot_PC2_SR_pal <- plot(pca_sym.SR_pal, axis1=2, axis2=3)
text(pca_sym.SR_pal[["x"]][,2], pca_sym.SR_pal[["x"]][,3], labels = names_SR_pal) # PCA 2 vs 3

# Saving PC scores #

PC.scores_SR_pal <- pca_sym.SR_pal$x 
as.data.frame(PC.scores_SR_pal) # Save PC scores as a data frame object

write.csv(PC.scores_SR_pal,"PC.scores_SR_pal.csv",row.names=TRUE) # Save PC scores as a csv file

##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_SR_pal), aes(x=PC.scores_SR_pal[,1], y=PC.scores_SR_pal[,2], label = names_SR_pal, fontface = "italic")) +
  labs(fill="Paleoenvironments") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_SR_pal$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 49.95 %", y = "PC2 = 18.58 %" )  

#PC2 vs PC3 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_SR_pal), aes(x=PC.scores_SR_pal[,2], y=PC.scores_SR_pal[,3], label = names_SR_pal, fontface = "italic")) +
  labs(fill="Paleoenvironments") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_SR_pal$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.3,0.3), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 25, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 29, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 29, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC2 = 18.58 %", y = "PC3 = 11.67 %" )

################### Calculating convex hulls for the PC1 vs PC2 #######################

# Add classifier variables (group, age, habitats) by hand to .csv file
# loading .csv file with PC scores and classifiers:
pca.class <- read.csv(file.choose()) #Open PC.scores_SR_pal.csv

# for the more precise paleoenvironments
CHull_habitat2_SR <- convex.hulls(data = pca.class, name.x = "Comp1", name.y = "Comp2", name.fill = "Habitat2")

CHull_habitat2_SR.table <- CHull_habitat2_SR$table
CHull_habitat2_SR$surface_area

ggplot(pca.class, aes(x = Comp1,y = Comp2, fill = Habitat2)) +
  labs(fill="Habitat2") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Habitat2), shape = 21, size = 5, stroke = 0.10) + 
  scale_fill_manual(values = cols_habitats2) +
  scale_color_manual(values= cols_habitats2)+
  lims(x=c(-0.7,0.7), y = c(-0.45,0.35))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 49.95 %", y = "PC2 = 18.58 %" ) +
  geom_polygon(data = CHull_habitat2_SR.table, aes(x = X, y = Y), alpha = 0.5)

################################# Morphological disparity analyses ###########################################

gdf_SR_habitat2 <- geomorph.data.frame(sym_pal$symm.shape, species = Period_SR_pal$Species, Habitat = Period_SR_pal$Paleoenvironments) #gdf for the more precise paleoenvironments
names(gdf_SR_habitat2)<-c("coords", "species", "Habitat")

# Sum of Procrustes variances : More precise aquatic habitats #########################

SOV_mean_SR_pal <- morphol.disparity(coords ~ 1, groups= ~ Habitat, partial = TRUE, 
                                     data = gdf_SR_habitat2, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_SR_pal)

SOV_SR_pal <- morphol.disparity(coords~Habitat,groups=~Habitat, data = gdf_SR_habitat2, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells which time bin is most/least disparate

summary(SOV_SR_pal)

# Graphical representation #

# Contribution of each time bin : Pie chart

PPV_SR.pal <- c("5.861", "7.774", "19.761", "2.107", "18.222", "18.067", "12.270", "15.938") # copy the SOV_mean_SR_pal results
Habitat2 <- c("Alluvial plain", "Bay", "Calm  freshwater lake", "Coastal marine", "Estuary", "Lagoon", "Oxbow lake/meandering river", "Reef")
PPV_SR.pal <- as.numeric(as.character(PPV_SR.pal))
Habitat2 <- as.factor(as.character(Habitat2))

DF_PPV.pal <- data.frame(Habitats = Habitat2, PPV = PPV_SR.pal)
DF_PPV.pal # Create data frame Paleoenvironment x Proportion of variance

# Compute percentages
DF_PPV.pal$fraction = DF_PPV.pal$PPV / sum(DF_PPV.pal$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.pal$ymax = cumsum(DF_PPV.pal$fraction)

# Compute the bottom of each rectangle
DF_PPV.pal$ymin = c(0, head(DF_PPV.pal$ymax, n=-1))

DF_PPV.pal$labelPosition <- (DF_PPV.pal$ymax + DF_PPV.pal$ymin) / 2

# Compute a good label
DF_PPV.pal$label <- paste0(DF_PPV.pal$PPV)

# Make the plot
ggplot(DF_PPV.pal, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Habitat2)) +
  geom_rect() + 
  coord_polar(theta="y") +
  scale_fill_manual(values = cols_habitats2) +
  xlim(c(2, 4)) +
  geom_label( x=3.5, aes(y=labelPosition, label=label, size=6)) +
  theme_void() +
  theme(legend.position = "none") #donut chart

ggplot(DF_PPV.pal, aes(x="", y=PPV, fill=Habitat2)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols_habitats2) +
  theme_void() #Pie chart

################################################################################
########################## Figure SOV through ages #############################
################################################################################

SOV <- data.frame(Time_bin = c(30:42), 
                  SOV_PC = c(0.000000000, 0.000000000, 0.01880830, 0.03004748, 0.02827263, 0.03001024, NA, NA, 0.03176942, NA, 0.000000000, NA, NA), 
                  SOV_Cheek = c(0.00000000, 0.00000000, 0.02637913, 0.09050480, 0.06943954, 0.05843217, 0.03526274, 0.02398637, 0.07960483, 0.00000000, NA, NA, NA),
                  SOV_SR = c(NA, 0.05657880, 0.04834485, 0.06864643, 0.06752854, 0.04274543, 0.00007766136, 0.06004026, 0.03116380, 0.00000000, 0.00000000, NA, 0.00000000))



tsplot(stages_upd, shading="stage", boxes=c("short", "series"), xlim=16:29, ylab="SOV", ylim=c(0,0.1), labels.args = list(cex=1), boxes.col = c("seriesCol", "systemCol"))

lines(stages_upd$mid[17:29], SOV$SOV_PC, col="blue", lwd = 4)
points(stages_upd$mid[17:29], SOV$SOV_PC, col="blue", pch = 15, cex = 2)
lines(stages_upd$mid[17:29], SOV$SOV_Cheek, col="red", lwd = 4)
points(stages_upd$mid[17:29], SOV$SOV_Cheek, col="red", pch = 16, cex= 2)
lines(stages_upd$mid[17:29], SOV$SOV_SR, col="#03be36", lwd = 4)
points(stages_upd$mid[17:29], SOV$SOV_SR, col="#03be36", pch = 17, cex = 2)


################################################################################
########################## Warp Grids for each PC ##############################
################################################################################

# Use this to find which specimen is closest to the mean shape:
PC <- findMeanSpec(data.super_PC$coords) # turned out to be "Osteolepis macrolepidotus"
Cheek <- findMeanSpec(data.super_Cheek$coords) # turned out to be "Thursius macrolepidotus"
SR <- findMeanSpec(sym$symm.shape) # turned out to be "Edenopteron keithcrooki"

# Estimate the mean shape for a set of aligned specimens
msh_PC <- mshape(data.super_PC$coords)
msh_Cheek <- mshape(data.super_Cheek$coords)
msh_SR <- mshape(sym$symm.shape)

#Compare the minimum and maximum values to the global consensus:

# Full body 
plotRefToTarget(pca_PC$shapes$shapes.comp1$min, msh_PC) 
plotRefToTarget(msh_PC, pca_PC$shapes$shapes.comp1$max) # PC 1
plotRefToTarget(msh_PC, data.super_PC$coords[,,18], mag = 0.8) # target one specimen

plotRefToTarget(pca_PC$shapes$shapes.comp2$min, msh_PC) 
plotRefToTarget(msh_PC, pca_PC$shapes$shapes.comp2$max) # PC 2

# Cheek
plotRefToTarget(pca_Cheek$shapes$shapes.comp1$min, msh_Cheek) 
plotRefToTarget(msh_Cheek, pca_Cheek$shapes$shapes.comp1$max) # PC 1
plotRefToTarget(msh_Cheek, data.super_Cheek$coords[,,37], mag = 0.8) # target one specimen

plotRefToTarget(pca_Cheek$shapes$shapes.comp2$min, msh_Cheek) 
plotRefToTarget(msh_Cheek, pca_Cheek$shapes$shapes.comp2$max) # PC 2

# Skull roof

plotRefToTarget(pca_sym.SR$shapes$shapes.comp1$min, msh_SR) 
plotRefToTarget(msh_SR, pca_sym.SR$shapes$shapes.comp1$max) # PC 1
plotRefToTarget(msh_SR, data.super_SR$coords[,,35], mag = 0.8) # target one specimen

plotRefToTarget(pca_sym.SR$shapes$shapes.comp2$min, msh_SR) 
plotRefToTarget(msh_SR, pca_sym.SR$shapes$shapes.comp2$max) # PC 2


###############################################################################
######################## LAGERSTATTE EFFECT ###################################
###############################################################################

# to test the Lagerstatte effect we are removing Bear Gulch (Serpukhovian Lagerstatte) from the dataset

################################################################################
################################## FULLBODY SHAPE ##############################
################################################################################

# TPS files with landmarks and semilandmarks
Sarco_PC <- readland.tps("Postcranial-final.TPS",specID="imageID",negNA = TRUE,readcurves = TRUE,warnmsg = TRUE)
dim(Sarco_PC)

Sarco_PC_fd<- Sarco_PC[,,-1][,,-4][,,-13][,,-19][,,-26] # removing species from Bear Gulch
name_PC_fd <- dimnames(Sarco_PC_fd)[[3]]


# Excel file with the list of species and their ages
Period_sarco_PC <- read_excel(file.choose(), 1) #Open List-species-PC.xlsx
Period_sarco_PC<- Period_sarco_PC %>%
  mutate_if(is.character, as.factor)

Period_PC_fd <- Period_sarco_PC %>%  filter(!row_number() %in% c(1,4,13,19,26))

# Convert curves into semi-landmarks
matrice_PC_fd <- rbind(define.sliders(12:41), define.sliders(42:61), define.sliders(62:71), define.sliders(72:107), define.sliders(108:207))

########################################## GPA and PCA ####################################################

# Procrustes superimposition #
data.super_PC_fd <- gpagen(Sarco_PC_fd, curves = matrice_PC_fd, ProcD = FALSE)
attributes(data.super_PC_fd)

plot(data.super_PC_fd) 

# Principal component analyses #
pca_PC_fd <- gm.prcomp(data.super_PC_fd$coords)
pca_PC_fd

plot_PCA_fd <- plot(pca_PC_fd, axis1 = 1, axis2 = 2)
text(pca_PC_fd[["x"]][,1], pca_PC_fd[["x"]][,2], labels = name_PC_fd) # PCA 1 vs 2

# Saving PC scores #

PC.scores_PC_fd <- pca_PC_fd$x 
as.data.frame(PC.scores_PC_fd) # Save PC scores as a data frame object


##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_PC_fd), aes(x=PC.scores_PC_fd[,1], y=PC.scores_PC_fd[,2], label = name_PC_fd, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_PC_fd$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 43.54 %", y = "PC2 = 19.73 %" ) 

#PC1 vs PC2 -- color Aquatic habitats 
ggplot(as.data.frame(PC.scores_PC_fd), aes(x=PC.scores_PC_fd[,1], y=PC.scores_PC_fd[,2], label = name_PC_fd, fontface = "italic")) +
  labs(fill="Aquatic habitat") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_PC_fd$'Aquatic habitat'), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 43.54 %", y = "PC2 = 19.73 %" ) 

#PC1 vs PC2 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_PC_fd), aes(x=PC.scores_PC_fd[,1], y=PC.scores_PC_fd[,2], label = name_PC_fd, fontface = "italic")) +
  labs(fill="Paleoenvironment") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_PC_fd$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 43.54 %", y = "PC2 = 19.73 %" ) 

#PC1 vs PC2 -- color groups

ggplot(as.data.frame(PC.scores_PC_fd), aes(x=PC.scores_PC_fd[,1], y=PC.scores_PC_fd[,2], label = name_PC_fd, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_PC_fd$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 43.54 %", y = "PC2 = 19.73 %" ) 

################################# Morphological disparity analyses ###########################################

gdf_PC_age_fd <- geomorph.data.frame(data.super_PC_fd, species = Period_PC_fd$Species, Time = Period_PC_fd$Age) #gdf for the epochs

gdf_PC_gp_fd <- geomorph.data.frame(data.super_PC_fd, species = Period_PC_fd$Species, Group = Period_PC_fd$Group) #gdf for the groups

# Sum of Procrustes variances : age ######################

SOV_mean_PC_age_fd <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_PC_age_fd, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_PC_age_fd)

SOV_PC_age_fd <- morphol.disparity(coords~Time,groups=~Time, data = gdf_PC_age_fd, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_PC_age_fd)


# Sum of Procrustes variances :  groups #########################

SOV_mean_PC_gp_fd <- morphol.disparity(coords ~ 1, groups= ~ Group, partial = TRUE, 
                                    data = gdf_PC_gp_fd, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_PC_gp_fd)

SOV_PC_gp_fd <- morphol.disparity(coords~Group,groups=~Group, data = gdf_PC_gp_fd, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_PC_gp_fd)

# Graphical representation #

cols <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Dipnoi" = "#FA8072", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F")

# Contribution of each time bin : Pie chart

PPV_PC.gp <- c("15.61", "31.47", "3.91", "28.87", "15.19", "4.94") # copy the SOV_mean_PC_gp results
Gp <- c("Actinistia", "Dipnoi", "Onychodontida", "Osteolepiform", "Porolepiform", "Rhizodontida")
PPV_PC.gp <- as.numeric(as.character(PPV_PC.gp))
Gp <- as.factor(as.character(Gp))

DF_PPV.gp <- data.frame(Group = Gp, PPV = PPV_PC.gp)
DF_PPV.gp # Create data frame Group x Proportion of variance

# Compute percentages
DF_PPV.gp$fraction = DF_PPV.gp$PPV / sum(DF_PPV.gp$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.gp$ymax = cumsum(DF_PPV.gp$fraction)

# Compute the bottom of each rectangle
DF_PPV.gp$ymin = c(0, head(DF_PPV.gp$ymax, n=-1))

DF_PPV.gp$labelPosition <- (DF_PPV.gp$ymax + DF_PPV.gp$ymin) / 2

# Compute a good label
DF_PPV.gp$label <- paste0(DF_PPV.gp$PPV)

# Make the plot
ggplot(DF_PPV.gp, aes(x="", y=PPV, fill=Gp)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols) +
  theme_void() #Pie chart


###############################################################################
################################### CHEEK SHAPE ###############################
###############################################################################

Sarco_Cheek <- readland.tps("Cheek-final.TPS",specID="imageID",negNA = TRUE,warnmsg = TRUE)
Sarco_Cheek_fd <- Sarco_Cheek[,,-1][,,-7][,,-17][,,-31]

name_Cheek <- dimnames(Sarco_Cheek_fd)[[3]]

Period_sarco_Cheek <- read_excel(file.choose(), 1) # Open List-species-cheek.xlsx
Period_sarco_Cheek<- Period_sarco_Cheek %>%
  mutate_if(is.character, as.factor)

Period_Cheek_fd <- Period_sarco_Cheek %>%  filter(!row_number() %in% c(1,7,17,31))

################################ GPA and PCA ##################################

# Procrustes superimposition #
data.super_Cheek_fd <- gpagen(Sarco_Cheek_fd, ProcD = FALSE)
attributes(data.super_Cheek_fd)

plot(data.super_Cheek_fd) 

# Principal component analyses #
pca_Cheek_fd <- gm.prcomp(data.super_Cheek_fd$coords)

plot_Cheek_fd <- plot(pca_Cheek_fd, axis1 = 1, axis2 = 2)
text(pca_Cheek_fd[["x"]][,1], pca_Cheek_fd[["x"]][,2], labels = name_Cheek) # PCA 1 vs 2

# Saving PC scores #

PC.scores_Cheek_fd <- pca_Cheek_fd$x 
as.data.frame(PC.scores_Cheek_fd) # Save PC scores as a data frame object


##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_Cheek_fd), aes(x=PC.scores_Cheek_fd[,1], y=PC.scores_Cheek_fd[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Cheek_fd$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 37.76 %", y = "PC2 = 16.81 %" ) 

#PC1 vs PC2 -- color Aquatic habitats 
ggplot(as.data.frame(PC.scores_Cheek_fd), aes(x=PC.scores_Cheek_fd[,1], y=PC.scores_Cheek_fd[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Aquatic habitat") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Cheek_fd$`Aquatic habitat`), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 37.76 %", y = "PC2 = 16.81 %" ) 

#PC1 vs PC2 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_Cheek_fd), aes(x=PC.scores_Cheek_fd[,1], y=PC.scores_Cheek_fd[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Paleoenvrionment") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Cheek_fd$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 37.76 %", y = "PC2 = 16.81 %" ) 

#PC1 vs PC2 -- color groups
ggplot(as.data.frame(PC.scores_Cheek_fd), aes(x=PC.scores_Cheek_fd[,1], y=PC.scores_Cheek_fd[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_Cheek_fd$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 37.76 %", y = "PC2 = 16.81 %" )

################################# Morphological disparity analyses ###########################################

gdf_Cheek_age_fd <- geomorph.data.frame(data.super_Cheek_fd, species = Period_Cheek_fd$Species, Time = Period_Cheek_fd$Age) #gdf for the ages

gdf_Cheek_gp_fd <- geomorph.data.frame(data.super_Cheek_fd, species = Period_Cheek_fd$Species, Group = Period_Cheek_fd$Group) #gdf for the groups

# Sum of Procrustes variances : age ######################

SOV_mean_Cheek_age_fd <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_Cheek_age_fd, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_Cheek_age_fd)

SOV_Cheek_age_fd <- morphol.disparity(coords~Time,groups=~Time, data = gdf_Cheek_age_fd, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_Cheek_age_fd)

# Sum of Procrustes variances :  groups #########################

SOV_mean_Cheek_gp_fd <- morphol.disparity(coords ~ 1, groups= ~ Group, partial = TRUE, 
                                             data = gdf_Cheek_gp_fd, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_Cheek_gp_fd)

SOV_Cheek_gp_fd <- morphol.disparity(coords~Group,groups=~Group, data = gdf_Cheek_gp_fd, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_Cheek_gp_fd)

# Graphical representation #

cols <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F", "Elpistostegalia" = "#ffd496")

# Contribution of each time bin : Pie chart

PPV_Cheek.gp <- c("26.67", "3.78", "14.26", "32.18", "19.48", "3.64") # copy the SOV_mean_Cheek_gp results
Gp <- c("Actinistia", "Elpistostegalia", "Onychodontida", "Osteolepiform", "Porolepiform", "Rhizodontida")
PPV_Cheek.gp <- as.numeric(as.character(PPV_Cheek.gp))
Gp <- as.factor(as.character(Gp))

DF_PPV.gp <- data.frame(Group = Gp, PPV = PPV_Cheek.gp)
DF_PPV.gp # Create data frame Group x Proportion of variance

# Compute percentages
DF_PPV.gp$fraction = DF_PPV.gp$PPV / sum(DF_PPV.gp$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.gp$ymax = cumsum(DF_PPV.gp$fraction)

# Compute the bottom of each rectangle
DF_PPV.gp$ymin = c(0, head(DF_PPV.gp$ymax, n=-1))

DF_PPV.gp$labelPosition <- (DF_PPV.gp$ymax + DF_PPV.gp$ymin) / 2

# Compute a good label
DF_PPV.gp$label <- paste0(DF_PPV.gp$PPV)

# Make the plot
ggplot(DF_PPV.gp, aes(x="", y=PPV, fill=Gp)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols) +
  theme_void() #Pie chart

###############################################################################
############################# SKULL ROOF SHAPE ################################
###############################################################################

Sarco_SR <- readland.tps("Skullroof-final.TPS",specID="imageID",negNA = TRUE,warnmsg = TRUE)
Sarco_SR_fd <- Sarco_SR[,,-2][,,-11][,,-26][,,-40]

name_SR <- dimnames(Sarco_SR_fd)[[3]]

Period_sarco_SR <- read_excel(file.choose(), 1) # Open List-species-SR.xlsx
Period_sarco_SR<- Period_sarco_SR %>%
  mutate_if(is.character, as.factor)

Period_SR_fd<- Period_sarco_SR %>%  filter(!row_number() %in% c(2,11,26,40))

################################### GPA and PCA ##################################

# Procrustes superimposition #

data.super_SR_fd <- gpagen(Sarco_SR_fd, ProcD = FALSE)
attributes(data.super_SR_fd)

plot(data.super_SR_fd) 

# Asymmetry #

nbb <- as.character(c(1:47))
data.super_SR_fd$ind=nbb # adding an ind vector to gpagen

n_pairs <- c(1,2,3,4,5,6,7,8,9,10,11,12,14,16,15,17,18,21,19,22,20,23)
pairs_matrix <- matrix(n_pairs, ncol=2, byrow = TRUE) # match paired landmarks
pairs_matrix


sym_fd <- bilat.symmetry(A = data.super_SR_fd$coords, ind=name_SR, object.sym = TRUE, land.pairs = pairs_matrix, iter = 149)
summary(sym_fd)

plot(sym_fd$symm.shape[,c(1,2),8])

# Principal component analyses #

pca_sym.SR_fd <- gm.prcomp(sym_fd$symm.shape) #with symmetrized coordinates
plot_sym_fd <- plot(pca_sym.SR_fd)

plot_SR_fd <- plot(pca_sym.SR_fd, axis1 = 1, axis2 = 2)
text(pca_sym.SR_fd[["x"]][,1], pca_sym.SR_fd[["x"]][,2], labels = name_SR) # PCA 1 vs 2

# Saving PC scores #

PC.scores_SR_fd <- pca_sym.SR_fd$x 
as.data.frame(PC.scores_SR_fd) # Save PC scores as a data frame object


##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_SR_fd), aes(x=PC.scores_SR_fd[,1], y=PC.scores_SR_fd[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_SR_fd$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 51.63 %", y = "PC2 = 18.27 %" ) 

#PC1 vs PC2 -- color Aquatic habitats
ggplot(as.data.frame(PC.scores_SR_fd), aes(x=PC.scores_SR_fd[,1], y=PC.scores_SR_fd[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Aquatic habitats") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_SR_fd$'Aquatic habitat'), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 51.63 %", y = "PC2 = 18.27 %" ) 

#PC1 vs PC2 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_SR_fd), aes(x=PC.scores_SR_fd[,1], y=PC.scores_SR_fd[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Paleoenvironment") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_SR_fd$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 51.63 %", y = "PC2 = 18.27 %" ) 

#PC1 vs PC2 -- color groups
ggplot(as.data.frame(PC.scores_SR_fd), aes(x=PC.scores_SR_fd[,1], y=PC.scores_SR_fd[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_SR_fd$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 51.63 %", y = "PC2 = 18.27 %" ) 


################################# Morphological disparity analyses ###########################################

gdf_SR_age_fd <- geomorph.data.frame(sym_fd$symm.shape, species = Period_SR_fd$Species, Time = Period_SR_fd$Age) #gdf for the ages
names(gdf_SR_age_fd)<-c("coords", "species", "Time")

gdf_SR_gp_fd <- geomorph.data.frame(sym_fd$symm.shape, species = Period_SR_fd$Species, Group = Period_SR_fd$Group) #gdf for the groups
names(gdf_SR_gp_fd)<-c("coords", "species", "Group")

# Sum of Procrustes variances : age ######################

SOV_mean_SR_age_fd <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_SR_age_fd, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_SR_age_fd)

SOV_SR_age_fd <- morphol.disparity(coords~Time,groups=~Time, data = gdf_SR_age_fd, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_SR_age_fd)

# Sum of Procrustes variances :  groups #########################

SOV_mean_SR_gp_fd <- morphol.disparity(coords ~ 1, groups= ~ Group, partial = TRUE, 
                                          data = gdf_SR_gp_fd, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_SR_gp_fd)

SOV_SR_gp_fd <- morphol.disparity(coords~Group,groups=~Group, data = gdf_SR_gp_fd, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_SR_gp_fd)

# Graphical representation #
cols.SR <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Dipnoi" = "#FA8072", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F", "Elpistostegalia" = "#FFEBCD", "Tetrapoda" = "#D2B48C") #Assign each group a color


# Contribution of each time bin : Pie chart

PPV_SR.gp <- c("6.02", "31.65", "3.82", "9.72", "15.69", "7.00", "6.03", "20.07") # copy the SOV_mean_SR_gp results
Gp <- c("Actinistia", "Dipnoi", "Elpistostegalia", "Onychodontida", "Osteolepiform", "Porolepiform", "Rhizodontida", "Tetrapoda")
PPV_SR.gp <- as.numeric(as.character(PPV_SR.gp))
Gp <- as.factor(as.character(Gp))

DF_PPV.gp <- data.frame(Group = Gp, PPV = PPV_SR.gp)
DF_PPV.gp # Create data frame Group x Proportion of variance

# Compute percentages
DF_PPV.gp$fraction = DF_PPV.gp$PPV / sum(DF_PPV.gp$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.gp$ymax = cumsum(DF_PPV.gp$fraction)

# Compute the bottom of each rectangle
DF_PPV.gp$ymin = c(0, head(DF_PPV.gp$ymax, n=-1))

DF_PPV.gp$labelPosition <- (DF_PPV.gp$ymax + DF_PPV.gp$ymin) / 2

# Compute a good label
DF_PPV.gp$label <- paste0(DF_PPV.gp$PPV)

# Make the plot
ggplot(DF_PPV.gp, aes(x="", y=PPV, fill=Gp)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols.SR) +
  theme_void() #Pie chart


###################################### Variation of SOV through time

SOV_fd <- data.frame(Time_bin = c(30:42), 
                        SOV_PC = c(0.000000000, 0.000000000, 0.02901119, 0.02972264, 0.03290473, 0.04746851, NA, NA, 0.00000000, NA, 0.000000000, NA, NA), 
                        SOV_Cheek = c(0.00000000, 0.00000000, 0.02635731, 0.09035513, 0.06504816, 0.08450898, 0.03524049, 0.03101225, NA, 0.00000000, NA, NA, NA),
                        SOV_SR = c(NA, 0.05665448, 0.04833384, 0.06842920, 0.06190026, 0.04289169, 0.00007728467, 0.05913751, 0.00000000, 0.00000000, 0.00000000, NA, 0.00000000))

tsplot(stages_upd, shading="stage", boxes=c("short", "series"), xlim=16:29, ylab="SOV", ylim=c(0,0.1), labels.args = list(cex=1), boxes.col = c("seriesCol", "systemCol"))

lines(stages_upd$mid[17:29], SOV_fd$SOV_PC, col="blue", lwd = 4)
points(stages_upd$mid[17:29], SOV_fd$SOV_PC, col="blue", pch = 15, cex = 2)
lines(stages_upd$mid[17:29], SOV_fd$SOV_Cheek, col="red", lwd = 4)
points(stages_upd$mid[17:29], SOV_fd$SOV_Cheek, col="red", pch = 16, cex= 2)
lines(stages_upd$mid[17:29], SOV_fd$SOV_SR, col="#03be36", lwd = 4)
points(stages_upd$mid[17:29], SOV_fd$SOV_SR, col="#03be36", pch = 17, cex = 2)


################################################################################
######################## Allenypterus effect ###################################
################################################################################


################################################################################
################################## FULLBODY SHAPE ##############################
################################################################################

# TPS files with landmarks and semilandmarks
Sarco_PC <- readland.tps("Postcranial-final.TPS",specID="imageID",negNA = TRUE,readcurves = TRUE,warnmsg = TRUE)
dim(Sarco_PC)

Sarco_PC_sAlle <- Sarco_PC[,,-1]
name_PC <- dimnames(Sarco_PC_sAlle)[[3]]


# Excel file with the list of species and their ages
Period_sarco_PC <- read_excel(file.choose(), 1) #Open List-species-PC.xlsx
Period_sarco_PC<- Period_sarco_PC %>%
  mutate_if(is.character, as.factor)

Period_sAlle_PC <- Period_sarco_PC %>%  filter(!row_number() %in% 1)


# Convert curves into semi-landmarks
matrice_PC <- rbind(define.sliders(12:41), define.sliders(42:61), define.sliders(62:71), define.sliders(72:107), define.sliders(108:207))

########################################## GPA and PCA ####################################################

# Procrustes superimposition #
data.super_PC_sAlle <- gpagen(Sarco_PC_sAlle, curves = matrice_PC, ProcD = FALSE)
attributes(data.super_PC_sAlle)

plot(data.super_PC_sAlle) 

# Principal component analyses #
pca_PC_sAlle <- gm.prcomp(data.super_PC_sAlle$coords)
pca_PC_sAlle

plot_PCA_sAlle <- plot(pca_PC_sAlle, axis1 = 1, axis2 = 2)
text(pca_PC_sAlle[["x"]][,1], pca_PC_sAlle[["x"]][,2], labels = name_PC) # PCA 1 vs 2

# Saving PC scores #

PC.scores_PC_sAlle <- pca_PC_sAlle$x 
as.data.frame(PC.scores_PC_sAlle) # Save PC scores as a data frame object


##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_PC_sAlle), aes(x=PC.scores_PC_sAlle[,1], y=PC.scores_PC_sAlle[,2], label = name_PC, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_PC$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 44.13 %", y = "PC2 = 18.22 %" ) 


#PC1 vs PC2 -- color Aquatic habitats 
ggplot(as.data.frame(PC.scores_PC_sAlle), aes(x=PC.scores_PC_sAlle[,1], y=PC.scores_PC_sAlle[,2], label = name_PC, fontface = "italic")) +
  labs(fill="Aquatic habitats") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_PC$`Aquatic habitat`), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 44.13 %", y = "PC2 = 18.22 %" ) 


#PC1 vs PC2 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_PC_sAlle), aes(x=PC.scores_PC_sAlle[,1], y=PC.scores_PC_sAlle[,2], label = name_PC, fontface = "italic")) +
  labs(fill="Paleoenvironments") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_PC$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 44.13 %", y = "PC2 = 18.22 %" ) 

#PC1 vs PC2 -- color groups

ggplot(as.data.frame(PC.scores_PC_sAlle), aes(x=PC.scores_PC_sAlle[,1], y=PC.scores_PC_sAlle[,2], label = name_PC, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_PC$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.5,0.3), y = c(-0.3,0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 44.13 %", y = "PC2 = 18.22 %" ) 


################################# Morphological disparity analyses ###########################################

gdf_PC_age_sAlle <- geomorph.data.frame(data.super_PC_sAlle, species = Period_sAlle_PC$Species, Time = Period_sAlle_PC$Age) #gdf for the age

gdf_PC_gp_sAlle <- geomorph.data.frame(data.super_PC_sAlle, species = Period_sAlle_PC$Species, Group = Period_sAlle_PC$Group) #gdf for the groups

# Sum of Procrustes variances : Epochs #########################

SOV_mean_PC_age_sAlle <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_PC_age_sAlle, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_PC_age_sAlle)

SOV_PC_age_sAlle <- morphol.disparity(coords~Time,groups=~Time, data = gdf_PC_age_sAlle, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_PC_age_sAlle)

# Sum of Procrustes variances :  groups #########################

SOV_mean_PC_gp_sAlle <- morphol.disparity(coords ~ 1, groups= ~ Group, partial = TRUE, 
                                    data = gdf_PC_gp_sAlle, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_PC_gp_sAlle)

SOV_PC_gp <- morphol.disparity(coords~Group,groups=~Group, data = gdf_PC_gp, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_PC_gp)

# Graphical representation #

cols <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Dipnoi" = "#FA8072", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F")

# Contribution of each time bin : Pie chart

PPV_PC.gp <- c("34.39", "26.80", "3.51", "19.53", "11.24", "4.53") # copy the SOV_mean_PC_gp results
Gp <- c("Actinistia", "Dipnoi", "Onychodontida", "Osteolepiform", "Porolepiform", "Rhizodontida")
PPV_PC.gp <- as.numeric(as.character(PPV_PC.gp))
Gp <- as.factor(as.character(Gp))

DF_PPV.gp <- data.frame(Group = Gp, PPV = PPV_PC.gp)
DF_PPV.gp # Create data frame Group x Proportion of variance

# Compute percentages
DF_PPV.gp$fraction = DF_PPV.gp$PPV / sum(DF_PPV.gp$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.gp$ymax = cumsum(DF_PPV.gp$fraction)

# Compute the bottom of each rectangle
DF_PPV.gp$ymin = c(0, head(DF_PPV.gp$ymax, n=-1))

DF_PPV.gp$labelPosition <- (DF_PPV.gp$ymax + DF_PPV.gp$ymin) / 2

# Compute a good label
DF_PPV.gp$label <- paste0(DF_PPV.gp$PPV)

# Make the plot

ggplot(DF_PPV.gp, aes(x="", y=PPV, fill=Gp)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols) +
  theme_void() #Pie chart

###############################################################################
################################### CHEEK SHAPE ###############################
###############################################################################

Sarco_Cheek <- readland.tps("Cheek-final.TPS",specID="imageID",negNA = TRUE,warnmsg = TRUE)
Sarco_Cheek_sAlle <- Sarco_Cheek[,,-1]

name_Cheek <- dimnames(Sarco_Cheek_sAlle)[[3]]

Period_sarco_Cheek <- read_excel(file.choose(), 1) # Open List-species-cheek.xlsx
Period_sarco_Cheek<- Period_sarco_Cheek %>%
  mutate_if(is.character, as.factor)

Period_sAlle_Cheek <- Period_sarco_Cheek %>%  filter(!row_number() %in% 1)

################################ GPA and PCA ##################################

# Procrustes superimposition #
data.super_Cheek_sAlle <- gpagen(Sarco_Cheek_sAlle, ProcD = FALSE)
attributes(data.super_Cheek_sAlle)

plot(data.super_Cheek_sAlle) 

# Principal component analyses #
pca_Cheek_sAlle <- gm.prcomp(data.super_Cheek_sAlle$coords)

plot_Cheek_sAlle <- plot(pca_Cheek_sAlle, axis1 = 1, axis2 = 2)
text(pca_Cheek_sAlle[["x"]][,1], pca_Cheek_sAlle[["x"]][,2], labels = name_Cheek) # PCA 1 vs 2

# Saving PC scores #

PC.scores_Cheek_sAlle <- pca_Cheek_sAlle$x 
as.data.frame(PC.scores_Cheek_sAlle) # Save PC scores as a data frame object


##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_Cheek_sAlle), aes(x=PC.scores_Cheek_sAlle[,1], y=PC.scores_Cheek_sAlle[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_Cheek$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 34.90 %", y = "PC2 = 17.29 %" ) 

#PC1 vs PC2 -- color Aquatic habitats 
ggplot(as.data.frame(PC.scores_Cheek_sAlle), aes(x=PC.scores_Cheek_sAlle[,1], y=PC.scores_Cheek_sAlle[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Aquatic habitat") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_Cheek$`Aquatic habitat`), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 34.90 %", y = "PC2 = 17.29 %" ) 

#PC1 vs PC2 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_Cheek_sAlle), aes(x=PC.scores_Cheek_sAlle[,1], y=PC.scores_Cheek_sAlle[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Paleoenvrionment") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_Cheek$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 34.90 %", y = "PC2 = 17.29 %" ) 

#PC1 vs PC2 -- color groups
ggplot(as.data.frame(PC.scores_Cheek_sAlle), aes(x=PC.scores_Cheek_sAlle[,1], y=PC.scores_Cheek_sAlle[,2], label = name_Cheek, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_Cheek$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.8,0.6), y = c(-0.6,0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 34.90 %", y = "PC2 = 17.29 %" )

################################# Morphological disparity analyses ###########################################

gdf_Cheek_age_sAlle <- geomorph.data.frame(data.super_Cheek_sAlle, species = Period_sAlle_Cheek$Species, Time = Period_sAlle_Cheek$Age) #gdf for the ages

gdf_Cheek_gp_sAlle <- geomorph.data.frame(data.super_Cheek_sAlle, species = Period_sAlle_Cheek$Species, Group = Period_sAlle_Cheek$Group) #gdf for the groups

# Sum of Procrustes variances : age ######################

SOV_mean_Cheek_age_sAlle <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_Cheek_age_sAlle, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_Cheek_age_sAlle)

SOV_Cheek_age_sAlle <- morphol.disparity(coords~Time,groups=~Time, data = gdf_Cheek_age_sAlle, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_Cheek_age_sAlle)

# Sum of Procrustes variances :  groups #########################

SOV_mean_Cheek_gp_sAlle <- morphol.disparity(coords ~ 1, groups= ~ Group, partial = TRUE, 
                                       data = gdf_Cheek_gp_sAlle, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_Cheek_gp_sAlle)

SOV_Cheek_gp_sAlle <- morphol.disparity(coords~Group,groups=~Group, data = gdf_Cheek_gp_sAlle, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_Cheek_gp_sAlle)

# Graphical representation #

cols <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F", "Elpistostegalia" = "#ffd496")

# Contribution of each time bin : Pie chart

PPV_Cheek.gp <- c("43.08", "3.51", "13.08", "25.62", "11.41", "3.3") # copy the SOV_mean_Cheek_gp results
Gp <- c("Actinistia", "Elpistostegalia", "Onychodontida", "Osteolepiform", "Porolepiform", "Rhizodontida")
PPV_Cheek.gp <- as.numeric(as.character(PPV_Cheek.gp))
Gp <- as.factor(as.character(Gp))

DF_PPV.gp <- data.frame(Group = Gp, PPV = PPV_Cheek.gp)
DF_PPV.gp # Create data frame Group x Proportion of variance

# Compute percentages
DF_PPV.gp$fraction = DF_PPV.gp$PPV / sum(DF_PPV.gp$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.gp$ymax = cumsum(DF_PPV.gp$fraction)

# Compute the bottom of each rectangle
DF_PPV.gp$ymin = c(0, head(DF_PPV.gp$ymax, n=-1))

DF_PPV.gp$labelPosition <- (DF_PPV.gp$ymax + DF_PPV.gp$ymin) / 2

# Compute a good label
DF_PPV.gp$label <- paste0(DF_PPV.gp$PPV)

# Make the plot
ggplot(DF_PPV.gp, aes(x="", y=PPV, fill=Gp)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols) +
  theme_void() #Pie chart

###############################################################################
############################# SKULL ROOF SHAPE ################################
###############################################################################

Sarco_SR <- readland.tps("Skullroof-final.TPS",specID="imageID",negNA = TRUE,warnmsg = TRUE)
Sarco_SR_sAlle <- Sarco_SR[,,-2]

name_SR <- dimnames(Sarco_SR_sAlle)[[3]]
name_SR

Period_sarco_SR <- read_excel(file.choose(), 1) # Open List-species-SR.xlsx
Period_sarco_SR<- Period_sarco_SR %>%
  mutate_if(is.character, as.factor)

Period_sAlle_SR <- Period_sarco_SR %>%  filter(!row_number() %in% 2)

################################### GPA and PCA ##################################

# Procrustes superimposition #

data.super_SR_sAlle <- gpagen(Sarco_SR_sAlle, ProcD = FALSE)
attributes(data.super_SR_sAlle)

plot(data.super_SR_sAlle) 

# Asymmetry #

nbb <- as.character(c(1:47))
data.super_SR_sAlle$ind=nbb # adding an ind vector to gpagen

n_pairs <- c(1,2,3,4,5,6,7,8,9,10,11,12,14,16,15,17,18,21,19,22,20,23)
pairs_matrix <- matrix(n_pairs, ncol=2, byrow = TRUE) # match paired landmarks
pairs_matrix


sym_sAlle <- bilat.symmetry(A = data.super_SR_sAlle$coords, ind=name_SR, object.sym = TRUE, land.pairs = pairs_matrix, iter = 149)
summary(sym_sAlle)

plot(sym_sAlle$symm.shape[,c(1,2),8])

# Principal component analyses #

pca_sym.SR_sAlle <- gm.prcomp(sym_sAlle$symm.shape) #with symmetrized coordinates
plot_sym_sAlle <- plot(pca_sym.SR_sAlle)

plot_SR_sAlle <- plot(pca_sym.SR_sAlle, axis1 = 1, axis2 = 2)
text(pca_sym.SR_sAlle[["x"]][,1], pca_sym.SR_sAlle[["x"]][,2], labels = name_SR) # PCA 1 vs 2


# Saving PC scores #

PC.scores_SR_sAlle <- pca_sym.SR_sAlle$x 
as.data.frame(PC.scores_SR_sAlle) # Save PC scores as a data frame object


##################################### Graphical representations ############################################

########### PCA ###########

#PC1 vs PC2 -- color epoch
ggplot(as.data.frame(PC.scores_SR_sAlle), aes(x=PC.scores_SR_sAlle[,1], y=PC.scores_SR_sAlle[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Epoque") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_SR$Epoque), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_epoch)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 53.50 %", y = "PC2 = 17.72 %" ) 

#PC1 vs PC2 -- color Aquatic habitats
ggplot(as.data.frame(PC.scores_SR_sAlle), aes(x=PC.scores_SR_sAlle[,1], y=PC.scores_SR_sAlle[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Aquatic habitats") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_SR$'Aquatic habitat'), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 53.50 %", y = "PC2 = 17.72 %" ) 

#PC1 vs PC2 -- color more precise paleoenvironments
ggplot(as.data.frame(PC.scores_SR_sAlle), aes(x=PC.scores_SR_sAlle[,1], y=PC.scores_SR_sAlle[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Paleoenvironment") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_SR$Paleoenvironments), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_habitats2)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 53.50 %", y = "PC2 = 17.72 %" ) 

#PC1 vs PC2 -- color groups
ggplot(as.data.frame(PC.scores_SR_sAlle), aes(x=PC.scores_SR_sAlle[,1], y=PC.scores_SR_sAlle[,2], label = name_SR, fontface = "italic")) +
  labs(fill="Group") +
  coord_fixed(ratio = 1) +
  geom_point(aes(fill = Period_sAlle_SR$Group), color = "black", size = 4, shape = 21, stroke = 0.10)  +
  scale_fill_manual(values = cols_group)+
  lims(x=c(-0.4,0.4), y = c(-0.3,0.3)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 15, color="black", margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.title.x = element_text(size = 15, margin = margin(t = 15, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +
  labs(x = "PC1 = 53.50 %", y = "PC2 = 17.72 %" ) 


################################# Morphological disparity analyses ###########################################

gdf_SR_age_sAlle <- geomorph.data.frame(sym_sAlle$symm.shape, species = Period_sAlle_SR$Species, Time = Period_sAlle_SR$Age) #gdf for the ages
names(gdf_SR_age_sAlle)<-c("coords", "species", "Time")

gdf_SR_gp_sAlle <- geomorph.data.frame(sym_sAlle$symm.shape, species = Period_sAlle_SR$Species, Group = Period_sAlle_SR$Group) #gdf for the groups
names(gdf_SR_gp_sAlle)<-c("coords", "species", "Group")

# Sum of Procrustes variances : age ######################

SOV_mean_SR_age_sAlle <- morphol.disparity(coords ~ 1, groups= ~ Time, partial = TRUE, data = gdf_SR_age_sAlle, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean
summary(SOV_mean_SR_age_sAlle)

SOV_SR_age_sAlle <- morphol.disparity(coords~Time,groups=~Time, data = gdf_SR_age_sAlle, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate
summary(SOV_SR_age_sAlle)

# Sum of Procrustes variances :  groups #########################

SOV_mean_SR_gp_sAlle <- morphol.disparity(coords ~ 1, groups= ~ Group, partial = TRUE, 
                                    data = gdf_SR_gp_sAlle, iter = 999, print.progress = FALSE) #calculate disparity of the time bins and comparing it to the grand mean

summary(SOV_mean_SR_gp_sAlle)

SOV_SR_gp_sAlle <- morphol.disparity(coords~Group,groups=~Group, data = gdf_SR_gp_sAlle, iter = 999) # calculate the disparity within each bin, and compare that to the disparity within other bins. This tells me which time bin is most/least disparate

summary(SOV_SR_gp_sAlle)

# Graphical representation #
cols.SR <- c("Actinistia" = "#87CEFA", "Porolepiform" = "#FFA500", "Dipnoi" = "#FA8072", "Onychodontida" = "#1E90FF", "Osteolepiform" = "#32CD32", "Rhizodontida" = "#ADFF2F", "Elpistostegalia" = "#FFEBCD", "Tetrapoda" = "#D2B48C") #Assign each group a color


# Contribution of each time bin : Pie chart

PPV_SR.gp <- c("10.14", "29.88", "3.62", "8.75", "14.45", "8.83", "5.34", "18.99") # copy the SOV_mean_SR_gp results
Gp <- c("Actinistia", "Dipnoi", "Elpistostegalia", "Onychodontida", "Osteolepiform", "Porolepiform", "Rhizodontida", "Tetrapoda")
PPV_SR.gp <- as.numeric(as.character(PPV_SR.gp))
Gp <- as.factor(as.character(Gp))

DF_PPV.gp <- data.frame(Group = Gp, PPV = PPV_SR.gp)
DF_PPV.gp # Create data frame Group x Proportion of variance

# Compute percentages
DF_PPV.gp$fraction = DF_PPV.gp$PPV / sum(DF_PPV.gp$PPV)

# Compute the cumulative percentages (top of each rectangle)
DF_PPV.gp$ymax = cumsum(DF_PPV.gp$fraction)

# Compute the bottom of each rectangle
DF_PPV.gp$ymin = c(0, head(DF_PPV.gp$ymax, n=-1))

DF_PPV.gp$labelPosition <- (DF_PPV.gp$ymax + DF_PPV.gp$ymin) / 2

# Compute a good label
DF_PPV.gp$label <- paste0(DF_PPV.gp$PPV)

# Make the plot
ggplot(DF_PPV.gp, aes(x="", y=PPV, fill=Gp)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = cols.SR) +
  theme_void() #Pie chart

###################################### Variation of SOV through time

SOV_sAlle <- data.frame(Time_bin = c(30:42), 
                        SOV_PC = c(0.000000000, 0.000000000, 0.018704538, 0.029802063, 0.028033571, 0.030018082, NA, NA, 0.008151817, NA, 0.000000000, NA, NA), 
                        SOV_Cheek = c(0.00000000, 0.00000000, 0.02636855, 0.09037388, 0.06939746, 0.05822921, 0.03525328, 0.02399639, 0.06931725, 0.00000000, NA, NA, NA),
                        SOV_SR = c(NA, 0.05651637, 0.04820343, 0.06847243, 0.06749470, 0.04273853, 0.00007770336, 0.05998263, 0.01851974, 0.00000000, 0.00000000, NA, 0.00000000))

tsplot(stages_upd, shading="stage", boxes=c("short", "series"), xlim=16:29, ylab="SOV", ylim=c(0,0.1), labels.args = list(cex=1), boxes.col = c("seriesCol", "systemCol"))

lines(stages_upd$mid[17:29], SOV_sAlle$SOV_PC, col="blue", lwd = 4)
points(stages_upd$mid[17:29], SOV_sAlle$SOV_PC, col="blue", pch = 15, cex = 2)
lines(stages_upd$mid[17:29], SOV_sAlle$SOV_Cheek, col="red", lwd = 4)
points(stages_upd$mid[17:29], SOV_sAlle$SOV_Cheek, col="red", pch = 16, cex= 2)
lines(stages_upd$mid[17:29], SOV_sAlle$SOV_SR, col="#03be36", lwd = 4)
points(stages_upd$mid[17:29], SOV_sAlle$SOV_SR, col="#03be36", pch = 17, cex = 2)
