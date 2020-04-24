## R code for Fleischer et al. Galapagos mockingbirds


###############################################################################################################

setwd("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS")


# load packages

#install.packages("vegan")

library(vegan)
library(expss)
library(ggplot2)
library(phyloseq)
library(gridExtra)
library(ecodist)
library(tidyr)
library(tidyverse)
library(metagMisc)
require(scales)
require(reshape2)
library(microbiome)
library(adegenet)
library(poppr)
library(BEDASSLE)
library(coda)

####################################FULL MICROBIOME ###########################################
####################################FULL MICROBIOME ###########################################
####################################FULL MICROBIOME ###########################################
####################################FULL MICROBIOME ###########################################
####################################FULL MICROBIOME ###########################################
####################################FULL MICROBIOME ###########################################
####################################FULL MICROBIOME ###########################################


mockingbird_microbiome<-readRDS("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\PHYLOSEQ OBJECTS\\mockingbird_unrarefied.rds") 


#################################################### Generate rarefraction curves ##########################################################
#################################################### Generate rarefraction curves ##########################################################
#################################################### Generate rarefraction curves ##########################################################
#################################################### Generate rarefraction curves ##########################################################
#################################################### Generate rarefraction curves ##########################################################
#################################################### Generate rarefraction curves ##########################################################
#################################################### Generate rarefraction curves ##########################################################


## function to generate rarefaction curves (copied from https://github.com/mahendra-mariadassou/phyloseq-extended/blob/master/R/graphical_methods.R)

## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}


###################################### generate rarefaction curve plot ############################
###################################### generate rarefaction curve plot ############################
###################################### generate rarefaction curve plot ############################
###################################### generate rarefaction curve plot ############################
###################################### generate rarefaction curve plot ############################
###################################### generate rarefaction curve plot ############################



p <- ggrare(mockingbird_microbiome, color= 'sample_Species', step = 1000, se= FALSE)+ xlim(c(0,15000))
p + xlim(c(0,15000))
p1<- p + facet_wrap(~sample_Species, scales="free")+geom_vline(xintercept=10000)+theme_bw()+theme(legend.position ="none")
p1


rarefaction_figure<-p1

table(sample_data(mockingbird_microbiome)$Island)



################################# RAREFY ###################################
################################# RAREFY ###################################
################################# RAREFY ###################################
################################# RAREFY ###################################
################################# RAREFY ###################################
################################# RAREFY ###################################

mockingbirds_rare<- rarefy_even_depth(mockingbird_microbiome, sample.size = 10000,rngsee = 100, replace = TRUE, trimOTUs=TRUE,verbose=TRUE)

mockingbirds_rare

###bar plot

bar_plot<-plot_bar(mockingbirds_rare, fill = 'Phylum')+
  facet_wrap(~sample_Species+Island, scales= "free")+
  theme(legend.position = "none")


################################# alpha diversity ##############################
################################# alpha diversity ##############################
################################# alpha diversity ##############################
################################# alpha diversity ##############################
################################# alpha diversity ##############################
################################# alpha diversity ##############################


Richness<-estimate_richness(mockingbirds_rare, measures=c("Observed","Simpson"))
#faiths_richness<-estimate_pd(phylo_object) #phylogenetic diversity #this uses btools
faiths_richness<-phyloseq_phylo_div(mockingbirds_rare, measures = c("PD")) #another method using metagMisc

Alpha<-sample_data(mockingbirds_rare)
names(Alpha)

Alpha$Observed<-Richness$Observed
Alpha$Simpson<-Richness$Simpson
Alpha$Faiths<-faiths_richness$PD

alpha_long<-gather(Alpha, Measure, Richness, Observed, Simpson, Faiths, factor_key = TRUE)


supp.labs1 <- c("Observed ASVs", "Simpson Index","Faiths PD Index")
names(supp.labs1) <- c("Observed", "Simpson", "Faiths")

supp.labs1

##  Plot
Alpha.sample_Species<-ggplot(alpha_long, aes(x =sample_Species, y = Richness, fill=sample_Species))+
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + 
  # stat_compare_means(method = "anova", label.y = c(0.3),label.x =2)+
  scale_fill_manual(values=c("royalblue","orange","darkolivegreen3","darkorchid3"))+
  facet_wrap(~Measure, scales = "free", ncol = 3, labeller = labeller(Measure = supp.labs1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background = element_rect(colour="black", fill="white"))+
  labs(x ="sample_Species")+
  theme(strip.text = element_text(size=11))+
  geom_jitter(aes(fill = sample_Species), size=1, colour = "black", alpha = 0.7, width =0.1)

Alpha.sample_Species

####################supplementary figures ####################
####################supplementary figures ####################
####################supplementary figures ####################
####################supplementary figures ####################



#rarefaction_figure
#bar_plot
#Alpha.sample_Species




####################################### BETA DIVERSITY ############################################
####################################### BETA DIVERSITY ############################################
####################################### BETA DIVERSITY ############################################
####################################### BETA DIVERSITY ############################################
####################################### BETA DIVERSITY ############################################

mockingbirds_rare1 <- core(mockingbirds_rare, detection = 0, prevalence = .05)



names(sample_data(mockingbirds_rare1))[8]<-'Species'

# take out and edit metadata
mb_metadata<-data.frame(sample_data(mockingbirds_rare1))
mb_metadata$Species <-factor(mb_metadata$Species , levels=c("M.parvulus", "M.trifasciatus", "M.macdonaldi", "M.melanotis"))
mb_metadata$Island<-factor(mb_metadata$Island, levels=c("Marchena","StCruz","Isabela","Santiago","Rabida","Fernandina","Pinta","StFe","Gardner","Champion","Espanola","SanCristobal"))

sample_data(mockingbirds_rare1)<- mb_metadata


# ordinate

set.seed(1)

mockingbirds_mb_uni_mds <- ordinate(
  physeq = mockingbirds_rare1, 
  method = "MDS", 
  distance = "unifrac")



#extract axis info from pcoa

mockingbirds_mb_uni_mds_vectors<-mockingbirds_mb_uni_mds$vectors
mockingbirds_mb_uni_mds_vectors<-as.data.frame(mockingbirds_mb_uni_mds_vectors)

# extract only axis 1+2 into new df

mb_df_axis<-mockingbirds_mb_uni_mds_vectors[,1:2]
mb_df_axis

# add species and island info
# check if Ids metadata file and axis df are ordered the same!

mb_metadata

mb_metadata$feature.id
mb_metadata$Island
mb_metadata$Species

mb_df_axis$Island<-mb_metadata$Island
mb_df_axis$Species<-mb_metadata$Species

plot_ordination(mockingbirds_rare1, mockingbirds_mb_uni_mds, type="Samples", color="Species", title="taxa")


# plot and make ellipse around species

mb_df_axis_full<-mb_df_axis

mbuni_mds<-ggplot(mb_df_axis_full)+
  geom_point(aes(x=Axis.1,y=Axis.2, color=Species, shape=Island), size = 3) + 
  stat_ellipse(aes(x=Axis.1,y=Axis.2,color=Species),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("royalblue3","orange","darkolivegreen3","darkorchid3")) + 
  scale_shape_manual(values=c(16,15,17,18,8,16,15,16,17)) +  
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Full microbiome") + 
  labs(x = "Axis 1 [14.8 %]", y ="Axis 2 [12.7 %]")+
  theme(legend.position = "none")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=18))


mbuni_mds


################################################### MHC Adonis test

mb_unifrac <- phyloseq::distance(mockingbirds_rare1, method = "unifrac")
mb_unifrac

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(mockingbirds_rare1))
#sampledf

#Test
#adonis2(mb_unifrac ~ Species+Species/Island,strata = sampledf$Species, data = sampledf)
adonis_mb<-adonis2(mb_unifrac ~ Species+Island, data = sampledf) #doesnt make a difference if use strata or not

adonis_mb

perm_mb<-permustats(adonis_mb)
mb_full_results<-summary(perm_mb) #ADONIS RESULTS MHC
mb_full_results





#############################################################


#################### DISTANCE MATRIX - MEAN DISTANCES BETWEEN CENTROIDS FORE CORE MICROBIOME (used in bedassle analysis)

mb_distance <- phyloseq::distance(mockingbirds_rare1, method = "unifrac")
mb_distance

# make a data frame from the sample_data
sample_data_mb <- data.frame(sample_data(mockingbirds_rare1))
sample_data_mb 

island_list_character<-as.character(sample_data_mb$Island)



mb_centroids_meandist<-meandist(mb_distance, island_list_character)
diag(mb_centroids_meandist)<-0
mb_meandist<-mb_centroids_meandist[,1:9]

#################################final list of object for full microbiome ######################
#################################final list of object for full microbiome ######################
#################################final list of object for full microbiome ######################
#################################final list of object for full microbiome ######################
#################################final list of object for full microbiome ######################

mbuni_mds
mb_full_results
mb_meandist


############################## weighted unifrac plot (Supplementary material

set.seed(1)

mockingbirds_mb_wuni_mds <- ordinate(
  physeq = mockingbirds_rare1, 
  method = "MDS", 
  distance = "wunifrac")



#extract axis info from pcoa

mockingbirds_mb_wuni_mds_vectors<-mockingbirds_mb_wuni_mds$vectors
mockingbirds_mb_wuni_mds_vectors<-as.data.frame(mockingbirds_mb_wuni_mds_vectors)

# extract only axis 1+2 into new df

mb_df_axis<-mockingbirds_mb_wuni_mds_vectors[,1:2]
mb_df_axis

# add species and island info
# check if Ids metadata file and axis df are ordered the same!

mb_metadata

mb_metadata$feature.id
mb_metadata$Island
mb_metadata$Species

mb_df_axis$Island<-mb_metadata$Island
mb_df_axis$Species<-mb_metadata$Species

plot_ordination(mockingbirds_rare1, mockingbirds_mb_wuni_mds, type="Samples", color="Species", title="taxa")


# plot and make ellipse around species

mb_df_axis_full_weighted<-mb_df_axis

mbuni_mds_weighted<-ggplot(mb_df_axis_full_weighted)+
  geom_point(aes(x=Axis.1,y=Axis.2, color=Species, shape=Island), size = 3) + 
  stat_ellipse(aes(x=Axis.1,y=Axis.2,color=Species),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("royalblue3","orange","darkolivegreen3","darkorchid3")) + 
  scale_shape_manual(values=c(16,15,17,18,8,16,15,16,17)) +  
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Weighted Unifrac: Full microbiome") + 
  labs(x = "Axis 1 [34.1 %]", y ="Axis 2 [13.5 %]")+
 # theme(legend.position = "none")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=18))


mbuni_mds_weighted

mb_wunifrac <- phyloseq::distance(mockingbirds_rare1, method = "wunifrac")
mb_wunifrac

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(mockingbirds_rare1))
#sampledf

#Test
adonis_mb_weighted<-adonis2(mb_wunifrac ~ Species+Island, data = sampledf) #doesnt make a difference if use strata or not

adonis_mb_weighted
perm_mb_weighted<-permustats(adonis_mb_weighted)
summary(perm_mb_weighted) #ADONIS RESULTS WEIGHTED




######################################### REPEAT FOR (SPECIES) CORE MICROBIOME ######################
######################################### REPEAT FOR (SPECIES) CORE MICROBIOME ######################
######################################### REPEAT FOR (SPECIES) CORE MICROBIOME ######################
######################################### REPEAT FOR (SPECIES) CORE MICROBIOME ######################
######################################### REPEAT FOR (SPECIES) CORE MICROBIOME ######################
######################################### REPEAT FOR (SPECIES) CORE MICROBIOME ######################

### core Microbiome UniFrac MDS 


mockingbird_core<-readRDS("mockingbird_island_core.rds")
mockingbird_core<-readRDS("mockingbird_island_core_css.rds")



names(sample_data(mockingbird_core))[8]<-'Species'

# take out and edit metadata
mb_metadata<-data.frame(sample_data(mockingbird_core))
mb_metadata$Species <-factor(mb_metadata$Species , levels=c("M.parvulus", "M.trifasciatus", "M.macdonaldi", "M.melanotis"))
mb_metadata$Island<-factor(mb_metadata$Island, levels=c("Marchena","StCruz","Isabela","Santiago","Rabida","Fernandina","Pinta","StFe","Gardner","Champion","Espanola","SanCristobal"))

sample_data(mockingbird_core)<- mb_metadata


# ordinate

set.seed(1)

mockingbirds_mb_uni_mds <- ordinate(
  physeq = mockingbird_core, 
  method = "MDS", 
  distance = "unifrac")



#extract axis info from pcoa

mockingbirds_mb_uni_mds_vectors<-mockingbirds_mb_uni_mds$vectors
mockingbirds_mb_uni_mds_vectors<-as.data.frame(mockingbirds_mb_uni_mds_vectors)

# extract only axis 1+2 into new df

mb_df_axis<-mockingbirds_mb_uni_mds_vectors[,1:2]
mb_df_axis

# add species and island info
# check if Ids metadata file and axis df are ordered the same!

mb_metadata

mb_metadata$feature.id
mb_metadata$Island
mb_metadata$Species

mb_df_axis$Island<-mb_metadata$Island
mb_df_axis$Species<-mb_metadata$Species

mb_df_axis_core<-mb_df_axis

plot_ordination(mockingbird_core, mockingbirds_mb_uni_mds, type="Samples", color="Species", title="taxa")

# plot and make ellipse around species

mbuni_mds_core<-ggplot(mb_df_axis_core)+
  geom_point(aes(x=Axis.1,y=Axis.2, color=Species, shape=Island), size=3) + 
  stat_ellipse(aes(x=Axis.1,y=Axis.2,color=Species),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("royalblue3","orange","darkolivegreen3","darkorchid3")) + 
  scale_shape_manual(values=c(16,15,17,18,8,16,15,16,17)) +  
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Core microbiome") + 
  labs(x = "Axis 1 [26.8 %]", y ="Axis 2 [16.0 %]")+
  theme(legend.position = "none")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=18))
  

mbuni_mds_core




###################################################  Adonis test

mb_unifrac <- phyloseq::distance(mockingbird_core, method = "unifrac")
mb_unifrac

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(mockingbird_core))
sampledf

#Test
#adonis2(mb_unifrac ~ Species+Species/Island,strata = sampledf$Species, data = sampledf)
adonis_mb<-adonis2(mb_unifrac ~ Species+Island, data = sampledf) #doesnt make a difference if use strata or not

perm_mb<-permustats(adonis_mb)
mb_core_results<-summary(perm_mb) #ADONIS RESULTS MHC
mb_core_results


#############################################################
#################### DISTANCE MATRIX - MEAN DISTANCES BETWEEN CENTROIDS FORE CORE MICROBIOME (used in bedassle analysis)

mb_core_distance <- phyloseq::distance(mockingbird_core, method = "unifrac")
mb_core_distance

# make a data frame from the sample_data
sample_data_mb <- data.frame(sample_data(mockingbird_core))
sample_data_mb 

island_list_character<-as.character(sample_data_mb$Island)

#library('usedist')

#alternative method 2

mb_core_centroids_meandist<-meandist(mb_core_distance, island_list_character)
diag(mb_core_centroids_meandist)<-0
mb_core_meandist<-mb_core_centroids_meandist[,1:9]



##################final list of object for core microbiome ######################
##################final list of object for core microbiome ######################
##################final list of object for core microbiome ######################
##################final list of object for core microbiome ######################
##################final list of object for core microbiome ######################

mbuni_mds_core
mb_core_results
mb_core_meandist




################################################### MHC #################################################
################################################### MHC #################################################
################################################### MHC #################################################
################################################### MHC #################################################
################################################### MHC #################################################
################################################### MHC #################################################
################################################### MHC #################################################
################################################### MHC #################################################
################################################### MHC #################################################
################################################### MHC #################################################


### MHC UniFrac MDS

mhc_phylo<-readRDS("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\mhc_phyloseq.rds")

mhc_phylo

mhc_tree<-phy_tree(mhc_phylo)

# exclude islands we don't need
mhc_phylo1<- subset_samples(mhc_phylo, Island!="Fernandina")
mhc_phylo1<- subset_samples(mhc_phylo1, Island!="Pinta")
mhc_phylo1<- subset_samples(mhc_phylo1, Island!="StFe")
mhc_phylo1

# take out and edit metadata file
metadata<-data.frame(sample_data(mhc_phylo1))
metadata$Species <-factor(metadata$Species , levels=c("M.parvulus", "M.trifasciatus", "M.macdonaldi", "M.melanotis"))
metadata$Island<-factor(metadata$Island, levels=c("Marchena","StCruz","Isabela","Santiago","Rabida","Fernandina","Pinta","StFe","Gardner","Champion","Espanola","SanCristobal"))

sample_data(mhc_phylo1)<- metadata

levels(sample_data(mhc_phylo1)$Species)
levels(sample_data(mhc_phylo1)$Island)


# ordination
set.seed(1)

mockingbirds_mhc_uni_mds <- ordinate(
  physeq = mhc_phylo1, 
  method = "MDS", 
  distance = "unifrac")

str(mockingbirds_mhc_uni_mds )
plot(mockingbirds_mhc_uni_mds$values$Eigenvalues)
plot(mockingbirds_mhc_uni_mds$values$Relative_eig)

#extract asix info from pcoa

mockingbirds_mhc_uni_mds_vectors<-mockingbirds_mhc_uni_mds$vectors
mockingbirds_mhc_uni_mds_vectors<-as.data.frame(mockingbirds_mhc_uni_mds_vectors)

# extract only axis 1+2 into new df

mhc_df_axis<-mockingbirds_mhc_uni_mds_vectors[,1:2]
mhc_df_axis

# add species and island info
# check if Ids in metadata file and axis df are ordered the same!

metadata$ID
metadata$Island
metadata$Species

mhc_df_axis$Island<-metadata$Island
mhc_df_axis$Species<-metadata$Species

plot_ordination(mhc_phylo1, mockingbirds_mhc_uni_mds, type="Samples", color="Species", title="MHC")


# plot and make ellipse around species
### MHC BETA PLOT

mhcuni_mds<-ggplot(mhc_df_axis)+
  geom_point(aes(x=Axis.1,y=Axis.2, color=Species, shape=Island), size=3) + 
  stat_ellipse(aes(x=Axis.1,y=Axis.2,color=Species),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("royalblue3","orange","darkolivegreen3","darkorchid3")) + 
  scale_shape_manual(values=c(16,15,17,18,8,16,15,16,17)) +  
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("MHC") + 
  labs(x = "Axis 1 [22 %]", y ="Axis 2 [16.4 %]")+
  theme(legend.position = "none")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=18))
  


mhcuni_mds


################################################### MHC Adonis test

mhc_unifrac <- phyloseq::distance(mhc_phylo1, method = "unifrac")
mhc_unifrac

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(mhc_phylo1))
sampledf

#Test
adonis2(mhc_unifrac ~ Species+Species/Island,strata = sampledf$Species, data = sampledf)
adonis_mhc<-adonis2(mhc_unifrac ~ Species+Island, data = sampledf) #doesnt make a difference if use strata or not

perm_mhc<-permustats(adonis_mhc)
mhc_results<-summary(perm_mhc) #ADONIS RESULTS MHC



#################### DISTANCE MATRIX - MEAN DISTANCES BETWEEN CENTROIDS

mhc_distance <- phyloseq::distance(mhc_phylo1, method = "unifrac") 
sample_data_mhc <- data.frame(sample_data(mhc_phylo1)) # make data frame of metadata

island_list_character<-as.character(sample_data_mhc$Island)

#library('usedist')


#pairwise mean distance

mhc_centroids_meandist<-meandist(mhc_distance, sample_data_mhc$Island)
diag(mhc_centroids_meandist)<-0
mhc_meandist<-mhc_centroids_meandist[,1:9]


island.order<-rownames(mb_meandist)

mhc_meandist<-mhc_meandist[,c(island.order)]
mhc_meandist<-mhc_meandist[c(island.order),]



############################ final objects for MHC ###################
############################ final objects for MHC ###################
############################ final objects for MHC ###################
############################ final objects for MHC ###################

mhcuni_mds
mhc_results
mhc_meandist 




############################################### MICROSATELLITES ############################
############################################### MICROSATELLITES ############################
############################################### MICROSATELLITES ############################
############################################### MICROSATELLITES ############################
############################################### MICROSATELLITES ############################
############################################### MICROSATELLITES ############################
############################################### MICROSATELLITES ############################
############################################### MICROSATELLITES ############################



### Microsatellite PCA

#install.packages("adegenet")

MicrosatTable<-read.csv("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\microsat_individual_alleles.csv", sep = ',', header = TRUE)

microsat<-as.data.frame(MicrosatTable)

microsat


## subset df and only use loci columns in df2genind
# we will add individual names and population info later to the genind object


microsat_genind<-df2genind(microsat, ploidy = 2,sep="_")

microsat_genind

## add individual names as rownames to object

row.names(microsat_genind@tab)<-c("237",
                                  "1055",
                                  "1109",
                                  "1116",
                                  "1123",
                                  "1125",
                                  "1179",
                                  "1344",
                                  "143001",
                                  "143005",
                                  "143008",
                                  "143011",
                                  "143026",
                                  "143050",
                                  "143052",
                                  "143053",
                                  "143054",
                                  "143055",
                                  "143090",
                                  "143096",
                                  "143105",
                                  "143110",
                                  "143112",
                                  "143113",
                                  "143114",
                                  "143115",
                                  "143116",
                                  "143118",
                                  "143119",
                                  "143120",
                                  "143121",
                                  "143122",
                                  "143123",
                                  "143124",
                                  "143128",
                                  "143130",
                                  "143134",
                                  "143138",
                                  "143144",
                                  "143153",
                                  "143154",
                                  "143157",
                                  "143165",
                                  "143166",
                                  "143168",
                                  "143169",
                                  "143171",
                                  "143172",
                                  "143173",
                                  "143175",
                                  "143176",
                                  "143181",
                                  "143182",
                                  "143184",
                                  "143185",
                                  "143187",
                                  "143190",
                                  "143192",
                                  "143193",
                                  "143196",
                                  "143197",
                                  "143201",
                                  "143203",
                                  "143204",
                                  "143209",
                                  "143215",
                                  "143216",
                                  "143220",
                                  "143223",
                                  "143225",
                                  "143231",
                                  "143236",
                                  "143238",
                                  "143239",
                                  "143243",
                                  "143245",
                                  "143248",
                                  "143255",
                                  "143256",
                                  "143257",
                                  "143260",
                                  "143264",
                                  "143268",
                                  "143269",
                                  "143273",
                                  "143278",
                                  "143324",
                                  "143326",
                                  "143328",
                                  "143330",
                                  "143331",
                                  "143334",
                                  "143335",
                                  "143336",
                                  "143338",
                                  "143339",
                                  "143340",
                                  "143341",
                                  "143342",
                                  "143343",
                                  "143345",
                                  "143368",
                                  "143369",
                                  "143371",
                                  "143372",
                                  "143375",
                                  "143376",
                                  "143379",
                                  "143380",
                                  "143383",
                                  "143386",
                                  "143392",
                                  "Cris04",
                                  "Cris06",
                                  "Cris08",
                                  "Cris13",
                                  "Cris14",
                                  "Cris15",
                                  "EM2",
                                  "GarF01",
                                  "GarF02",
                                  "GarF04",
                                  "GarF08",
                                  "GarF10",
                                  "GarF11",
                                  "GarF12",
                                  "GarF13",
                                  "GarF14",
                                  "GarF18",
                                  "GarF22",
                                  "IM10",
                                  "IM11",
                                  "IM12",
                                  "IM26",
                                  "IM7")

microsat_genind@tab # check it worked


## add population (island) info to the genind object

microsat_genind@pop<- as.factor(c("Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "StCruz",
                                  "Gardner",
                                  "StCruz",
                                  "StCruz",
                                  "StCruz",
                                  "StCruz",
                                  "Champion",
                                  "Champion",
                                  "Champion",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "Santiago",
                                  "StCruz",
                                  "StCruz",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "Espanola",
                                  "StCruz",
                                  "StCruz",
                                  "StCruz",
                                  "StCruz",
                                  "StCruz",
                                  "StCruz",
                                  "StCruz",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "Marchena",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "Champion",
                                  "Champion",
                                  "Champion",
                                  "Champion",
                                  "Champion",
                                  "Champion",
                                  "Champion",
                                  "Champion",
                                  "Champion",
                                  "Champion",
                                  "Gardner",
                                  "Gardner",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "Rabida",
                                  "StCruz",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "SanCristobal",
                                  "Espanola",
                                  "Gardner",
                                  "Gardner",
                                  "Gardner",
                                  "Gardner",
                                  "Gardner",
                                  "Gardner",
                                  "Gardner",
                                  "Gardner",
                                  "Gardner",
                                  "Gardner",
                                  "Gardner",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela",
                                  "Isabela"))

microsat_genind@pop

microsat_genind


# check some variables in the data to be sure everything is stored correct

head(pop(microsat_genind), 50)

barplot(table(pop(microsat_genind)), col=funky(17), las=3,xlab="Population", ylab="Sample size")

temp <- summary(microsat_genind)
temp

# Hexp = expected heterozygosity, Hobs = observed heterozygosity
plot(temp$Hexp, temp$Hobs, pch=20, cex=3, xlim=c(.4,1), ylim=c(.4,1))
abline(0,1,lty=2)




############# Multivariante analyses (PCA)

x.mock<- tab(microsat_genind, freq=TRUE, NA.method="mean")
x.mock




pca.mock <- dudi.pca(x.mock, center=TRUE, scale=FALSE)


## select number of axis: 2

# You can reproduce this result non-interactively with: 
# dudi.pca(df = x.mock, center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)

pca.mock1 <- dudi.pca(df = x.mock, center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
str(pca.mock1)


# some exploratory plots with adegenet
s.label(pca.mock1$li)
s.class(pca.mock1$li, fac=pop(microsat_genind), col=funky(15))

s.class(pca.mock1$li, fac=pop(microsat_genind),col=transp(funky(15),.6),axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.mock1$eig[1:50],3,1,2, ratio=.3)


# Extract axis 1+2 into new df to plot pca in ggplot
pcadf<-pca.mock1$li


# add island and species info to df

pcadf$Island<-microsat_genind@pop


pcadf$Island
pcadf$Island<-factor(pcadf$Island, levels=c("Marchena","StCruz","Isabela","Santiago","Rabida","Gardner","Champion","Espanola","SanCristobal"))
pcadf$Island

pcadf$Species<- c("M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.parvulus",
                  "M.trifasciatus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.macdonaldi",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.melanotis",
                  "M.melanotis",
                  "M.melanotis",
                  "M.melanotis",
                  "M.melanotis",
                  "M.melanotis",
                  "M.melanotis",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.melanotis",
                  "M.melanotis",
                  "M.melanotis",
                  "M.melanotis",
                  "M.melanotis",
                  "M.melanotis",
                  "M.macdonaldi",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.trifasciatus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus",
                  "M.parvulus")

# bring species into correct order for plot legend
pcadf$Species<-factor(pcadf$Species, levels=c("M.parvulus","M.trifasciatus","M.macdonaldi","M.melanotis"))
pcadf$Species


# plot

eigenvalues<-pca.mock1$eig

eigenvalues[1]/sum(eigenvalues) #percent explained by PC1
eigenvalues[2]/sum(eigenvalues) #percent explained by PC2

micromds<-ggplot(pcadf)+
  geom_point(aes(x=Axis1,y=Axis2, color=Species, shape=Island), size=3) + 
  stat_ellipse(aes(x=Axis1,y=Axis2,color=Species),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("royalblue3","orange","darkolivegreen3","darkorchid3")) + 
  scale_shape_manual(values=c(16,15,17,18,8,16,15,16,17)) +  
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Microsatellites") + 
  labs(x = "Axis 1 [24.7 %]", y ="Axis 2 [13.8 %]")+
  theme(legend.position = "none")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size = 18))

micromds 



############################################ generate mean distance matrix




microsat_nei.dist<-nei.dist(x.mock, warning=T)
head(pcadf)
microsat_nei.meandist<-meandist(microsat_nei.dist, pcadf$Island)
diag(microsat_nei.meandist)<-0
microsat_nei.meandist<-microsat_nei.meandist[,1:9]

island.order<-rownames(mb_meandist)

microsat_nei.meandist<-microsat_nei.meandist[,c(island.order)]
microsat_nei.meandist<-microsat_nei.meandist[c(island.order),]

microsat_meandist<-microsat_nei.meandist



############stats 

adonis_microsat<-adonis2(microsat_nei.dist ~ Species+Island, data = pcadf) #doesnt make a difference if use strata or not

perm_microsat<-permustats(adonis_microsat)
microsat_results<-summary(perm_microsat) #ADONIS RESULTS MHC





##################################

micromds 
microsat_results
microsat_meandist






########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################
########################### MANTEL TESTS ############################################

island.order<-rownames(mb_meandist)

##import geographic distance 

setwd("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS")

geographic_distance <- read.csv("geographic_mean_distance.csv", row.names=1)
geographic_distance<-as.matrix(geographic_distance)

geographic_distance<-geographic_distance[,c(island.order)]
geographic_distance<-geographic_distance[c(island.order),]




###############lets look at our mean pairwise distance matrices

mb_meandist
mb_core_meandist
mhc_meandist
microsat_meandist
geographic_distance


### standardize by standard deviation

geographic_distance1<-geographic_distance/sd(geographic_distance)

mb_core_meandist1<-mb_core_meandist/sd(mb_core_meandist)
mb_meandist1<-mb_meandist/sd(mb_meandist)

mhc_meandist1<-mhc_meandist/sd(mhc_meandist)
microsat_meandist1<-microsat_meandist/sd(microsat_meandist)



#library(ecodist)

##correlations between variables

vegan::mantel(as.dist(mhc_meandist1),as.dist(microsat_meandist1), permutations = 9999)
vegan::mantel(as.dist(mhc_meandist1),as.dist(geographic_distance1), permutations = 9999)
vegan::mantel(as.dist(microsat_meandist1),as.dist(geographic_distance1), permutations = 9999)



################################################################################################
################################################################################################
################################################################################################
################################################################################################



#from now on lets use mhc_centroids_meandist1 and microsat_meandist1 and geographic_distance1

#import matrix coding whether comparison is within or between species

ComparisonType <- read.csv("C:/Users/risel/Dropbox/Sommer postdoc/Ramona mockingbirds/Mockingbirds analysis/FINAL ANALYSIS/ComparisonType_matrix.csv", row.names=1)

#########CORE MICROBIOME
#########CORE MICROBIOME
#########CORE MICROBIOME
#########CORE MICROBIOME


# one variable
model1<- MRM(as.dist(mb_core_meandist1) ~ as.dist(mhc_meandist1), nperm = 9999)
model2<-MRM(as.dist(mb_core_meandist1) ~ as.dist(microsat_meandist1), nperm = 9999)
model3<-MRM(as.dist(mb_core_meandist1) ~ as.dist(geographic_distance1), nperm = 9999)
model4<-MRM(as.dist(mb_core_meandist1) ~ as.dist(ComparisonType), nperm = 9999)

# two variables

model5<-MRM(as.dist(mb_core_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(microsat_meandist1), nperm = 9999)
model6<-MRM(as.dist(mb_core_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(ComparisonType), nperm = 9999)
model7<-MRM(as.dist(mb_core_meandist1) ~ as.dist(microsat_meandist1)+ as.dist(ComparisonType), nperm = 9999)
model8<-MRM(as.dist(mb_core_meandist1) ~ as.dist(geographic_distance1)+ as.dist(ComparisonType), nperm = 9999)
model9<-MRM(as.dist(mb_core_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(geographic_distance1), nperm = 9999)
model10<-MRM(as.dist(mb_core_meandist1) ~ as.dist(microsat_meandist1)+ as.dist(geographic_distance1), nperm = 9999)

#3 variables

model11<-MRM(as.dist(mb_core_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(microsat_meandist1)+ as.dist(ComparisonType), nperm = 9999)
model12<-MRM(as.dist(mb_core_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(microsat_meandist1)+ as.dist(geographic_distance1), nperm = 9999)
model13<-MRM(as.dist(mb_core_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(ComparisonType)+ as.dist(geographic_distance1), nperm = 9999)
model14<-MRM(as.dist(mb_core_meandist1) ~ as.dist(microsat_meandist1)+ as.dist(ComparisonType)+ as.dist(geographic_distance1), nperm = 9999)

#4 variables
model15<-MRM(as.dist(mb_core_meandist1) ~ as.dist(mhc_meandist1)+as.dist(microsat_meandist1)+ as.dist(ComparisonType)+ as.dist(geographic_distance1), nperm = 9999)


# make table 28 rows and 9 cold

core_stats<-data.frame(matrix(ncol = 9, nrow = 30))
x <- c("name", "age", "gender")
colnames(core_stats) <- c('Model', 'Parameter','MHC', 'Microsatellite', 'Geograhic distance', 'Comparison type', 'Model R^2', 'Model F', 'Model p')

core_stats$Model<-c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14, 15, 15)
core_stats$Parameter<-c('F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p')

core_stats[1,3]<-model1$coef[2,1]
core_stats[2,3]<-model1$coef[2,2]
core_stats[1,7]<-model1$r.squared[1]
core_stats[1,8]<-model1$F.test[1]
core_stats[1,9]<-model1$F.test[2]

core_stats[3,4]<-model2$coef[2,1]
core_stats[4,4]<-model2$coef[2,2]
core_stats[3,7]<-model2$r.squared[1]
core_stats[3,8]<-model2$F.test[1]
core_stats[3,9]<-model2$F.test[2]

core_stats[5,5]<-model3$coef[2,1]
core_stats[6,5]<-model3$coef[2,2]
core_stats[5,7]<-model3$r.squared[1]
core_stats[5,8]<-model3$F.test[1]
core_stats[5,9]<-model3$F.test[2]

core_stats[7,6]<-model4$coef[2,1]
core_stats[8,6]<-model4$coef[2,2]
core_stats[7,7]<-model4$r.squared[1]
core_stats[7,8]<-model4$F.test[1]
core_stats[7,9]<-model4$F.test[2]

core_stats[9,3]<-model5$coef[2,1]
core_stats[10,3]<-model5$coef[2,2]
core_stats[9,4]<-model5$coef[3,1]
core_stats[10,4]<-model5$coef[3,2]
core_stats[9,7]<-model5$r.squared[1]
core_stats[9,8]<-model5$F.test[1]
core_stats[9,9]<-model5$F.test[2]

core_stats[11,3]<-model6$coef[2,1]
core_stats[12,3]<-model6$coef[2,2]
core_stats[11,6]<-model6$coef[3,1]
core_stats[12,6]<-model6$coef[3,2]
core_stats[11,7]<-model5$r.squared[1]
core_stats[11,8]<-model5$F.test[1]
core_stats[11,9]<-model5$F.test[2]

core_stats[13,4]<-model7$coef[2,1]
core_stats[14,4]<-model7$coef[2,2]
core_stats[13,6]<-model7$coef[3,1]
core_stats[14,6]<-model7$coef[3,2]
core_stats[13,7]<-model7$r.squared[1]
core_stats[13,8]<-model7$F.test[1]
core_stats[13,9]<-model7$F.test[2]

core_stats[15,5]<-model8$coef[2,1]
core_stats[16,5]<-model8$coef[2,2]
core_stats[15,6]<-model8$coef[3,1]
core_stats[16,6]<-model8$coef[3,2]
core_stats[15,7]<-model8$r.squared[1]
core_stats[15,8]<-model8$F.test[1]
core_stats[15,9]<-model8$F.test[2]

core_stats[17,3]<-model9$coef[2,1]
core_stats[18,3]<-model9$coef[2,2]
core_stats[17,5]<-model9$coef[3,1]
core_stats[18,5]<-model9$coef[3,2]
core_stats[17,7]<-model9$r.squared[1]
core_stats[17,8]<-model9$F.test[1]
core_stats[17,9]<-model9$F.test[2]

core_stats[19,4]<-model10$coef[2,1]
core_stats[20,4]<-model10$coef[2,2]
core_stats[19,5]<-model10$coef[3,1]
core_stats[20,5]<-model10$coef[3,2]
core_stats[19,7]<-model10$r.squared[1]
core_stats[19,8]<-model10$F.test[1]
core_stats[19,9]<-model10$F.test[2]

core_stats[21,3]<-model11$coef[2,1]
core_stats[22,3]<-model11$coef[2,2]
core_stats[21,4]<-model11$coef[3,1]
core_stats[22,4]<-model11$coef[3,2]
core_stats[21,6]<-model11$coef[4,1]
core_stats[22,6]<-model11$coef[4,2]
core_stats[21,7]<-model11$r.squared[1]
core_stats[21,8]<-model11$F.test[1]
core_stats[21,9]<-model11$F.test[2]


core_stats[23,3]<-model12$coef[2,1]
core_stats[24,3]<-model12$coef[2,2]
core_stats[23,4]<-model12$coef[3,1]
core_stats[24,4]<-model12$coef[3,2]
core_stats[23,5]<-model12$coef[4,1]
core_stats[24,5]<-model12$coef[4,2]
core_stats[23,7]<-model12$r.squared[1]
core_stats[23,8]<-model12$F.test[1]
core_stats[23,9]<-model12$F.test[2]


core_stats[25,3]<-model13$coef[2,1]
core_stats[26,3]<-model13$coef[2,2]
core_stats[25,6]<-model13$coef[3,1]
core_stats[26,6]<-model13$coef[3,2]
core_stats[25,5]<-model13$coef[4,1]
core_stats[26,5]<-model13$coef[4,2]
core_stats[25,7]<-model13$r.squared[1]
core_stats[25,8]<-model13$F.test[1]
core_stats[25,9]<-model13$F.test[2]

core_stats[27,4]<-model14$coef[2,1]
core_stats[28,4]<-model14$coef[2,2]
core_stats[27,6]<-model14$coef[3,1]
core_stats[28,6]<-model14$coef[3,2]
core_stats[27,5]<-model14$coef[4,1]
core_stats[28,5]<-model14$coef[4,2]
core_stats[27,7]<-model14$r.squared[1]
core_stats[27,8]<-model14$F.test[1]
core_stats[27,9]<-model14$F.test[2]

core_stats[29,3]<-model15$coef[2,1]
core_stats[30,3]<-model15$coef[2,2]
core_stats[29,4]<-model15$coef[3,1]
core_stats[30,4]<-model15$coef[3,2]
core_stats[29,6]<-model15$coef[4,1]
core_stats[30,6]<-model15$coef[4,2]
core_stats[29,5]<-model15$coef[5,1]
core_stats[30,5]<-model15$coef[5,2]
core_stats[29,7]<-model15$r.squared[1]
core_stats[29,8]<-model15$F.test[1]
core_stats[29,9]<-model15$F.test[2]


#Table S2

#write.csv(core_stats, 'core_mantel_tests.csv')

#########################################################
#########################################################
#########################################################
#########################################################


# one variable
model1<- MRM(as.dist(mb_meandist1) ~ as.dist(mhc_meandist1), nperm = 9999)
model2<-MRM(as.dist(mb_meandist1) ~ as.dist(microsat_meandist1), nperm = 9999)
model3<-MRM(as.dist(mb_meandist1) ~ as.dist(geographic_distance1), nperm = 9999)
model4<-MRM(as.dist(mb_meandist1) ~ as.dist(ComparisonType), nperm = 9999)

# two variables

model5<-MRM(as.dist(mb_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(microsat_meandist1), nperm = 9999)
model6<-MRM(as.dist(mb_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(ComparisonType), nperm = 9999)
model7<-MRM(as.dist(mb_meandist1) ~ as.dist(microsat_meandist1)+ as.dist(ComparisonType), nperm = 9999)
model8<-MRM(as.dist(mb_meandist1) ~ as.dist(geographic_distance1)+ as.dist(ComparisonType), nperm = 9999)
model9<-MRM(as.dist(mb_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(geographic_distance1), nperm = 9999)
model10<-MRM(as.dist(mb_meandist1) ~ as.dist(microsat_meandist1)+ as.dist(geographic_distance1), nperm = 9999)

#3 variables

model11<-MRM(as.dist(mb_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(microsat_meandist1)+ as.dist(ComparisonType), nperm = 9999)
model12<-MRM(as.dist(mb_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(microsat_meandist1)+ as.dist(geographic_distance1), nperm = 9999)
model13<-MRM(as.dist(mb_meandist1) ~ as.dist(mhc_meandist1)+ as.dist(ComparisonType)+ as.dist(geographic_distance1), nperm = 9999)
model14<-MRM(as.dist(mb_meandist1) ~ as.dist(microsat_meandist1)+ as.dist(ComparisonType)+ as.dist(geographic_distance1), nperm = 9999)

#4 variables
model15<-MRM(as.dist(mb_meandist1) ~ as.dist(mhc_meandist1)+as.dist(microsat_meandist1)+ as.dist(ComparisonType)+ as.dist(geographic_distance1), nperm = 9999)


# make table 28 rows and 9 cold

microbiome_stats<-data.frame(matrix(ncol = 9, nrow = 30))
x <- c("name", "age", "gender")
colnames(microbiome_stats) <- c('Model', 'Parameter','MHC', 'Microsatellite', 'Geograhic distance', 'Comparison type', 'Model R^2', 'Model F', 'Model p')

microbiome_stats$Model<-c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14, 15, 15)
microbiome_stats$Parameter<-c('F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p','F','p')

microbiome_stats[1,3]<-model1$coef[2,1]
microbiome_stats[2,3]<-model1$coef[2,2]
microbiome_stats[1,7]<-model1$r.squared[1]
microbiome_stats[1,8]<-model1$F.test[1]
microbiome_stats[1,9]<-model1$F.test[2]

microbiome_stats[3,4]<-model2$coef[2,1]
microbiome_stats[4,4]<-model2$coef[2,2]
microbiome_stats[3,7]<-model2$r.squared[1]
microbiome_stats[3,8]<-model2$F.test[1]
microbiome_stats[3,9]<-model2$F.test[2]

microbiome_stats[5,5]<-model3$coef[2,1]
microbiome_stats[6,5]<-model3$coef[2,2]
microbiome_stats[5,7]<-model3$r.squared[1]
microbiome_stats[5,8]<-model3$F.test[1]
microbiome_stats[5,9]<-model3$F.test[2]

microbiome_stats[7,6]<-model4$coef[2,1]
microbiome_stats[8,6]<-model4$coef[2,2]
microbiome_stats[7,7]<-model4$r.squared[1]
microbiome_stats[7,8]<-model4$F.test[1]
microbiome_stats[7,9]<-model4$F.test[2]

microbiome_stats[9,3]<-model5$coef[2,1]
microbiome_stats[10,3]<-model5$coef[2,2]
microbiome_stats[9,4]<-model5$coef[3,1]
microbiome_stats[10,4]<-model5$coef[3,2]
microbiome_stats[9,7]<-model5$r.squared[1]
microbiome_stats[9,8]<-model5$F.test[1]
microbiome_stats[9,9]<-model5$F.test[2]

microbiome_stats[11,3]<-model6$coef[2,1]
microbiome_stats[12,3]<-model6$coef[2,2]
microbiome_stats[11,6]<-model6$coef[3,1]
microbiome_stats[12,6]<-model6$coef[3,2]
microbiome_stats[11,7]<-model5$r.squared[1]
microbiome_stats[11,8]<-model5$F.test[1]
microbiome_stats[11,9]<-model5$F.test[2]

microbiome_stats[13,4]<-model7$coef[2,1]
microbiome_stats[14,4]<-model7$coef[2,2]
microbiome_stats[13,6]<-model7$coef[3,1]
microbiome_stats[14,6]<-model7$coef[3,2]
microbiome_stats[13,7]<-model7$r.squared[1]
microbiome_stats[13,8]<-model7$F.test[1]
microbiome_stats[13,9]<-model7$F.test[2]

microbiome_stats[15,5]<-model8$coef[2,1]
microbiome_stats[16,5]<-model8$coef[2,2]
microbiome_stats[15,6]<-model8$coef[3,1]
microbiome_stats[16,6]<-model8$coef[3,2]
microbiome_stats[15,7]<-model8$r.squared[1]
microbiome_stats[15,8]<-model8$F.test[1]
microbiome_stats[15,9]<-model8$F.test[2]

microbiome_stats[17,3]<-model9$coef[2,1]
microbiome_stats[18,3]<-model9$coef[2,2]
microbiome_stats[17,5]<-model9$coef[3,1]
microbiome_stats[18,5]<-model9$coef[3,2]
microbiome_stats[17,7]<-model9$r.squared[1]
microbiome_stats[17,8]<-model9$F.test[1]
microbiome_stats[17,9]<-model9$F.test[2]

microbiome_stats[19,4]<-model10$coef[2,1]
microbiome_stats[20,4]<-model10$coef[2,2]
microbiome_stats[19,5]<-model10$coef[3,1]
microbiome_stats[20,5]<-model10$coef[3,2]
microbiome_stats[19,7]<-model10$r.squared[1]
microbiome_stats[19,8]<-model10$F.test[1]
microbiome_stats[19,9]<-model10$F.test[2]

microbiome_stats[21,3]<-model11$coef[2,1]
microbiome_stats[22,3]<-model11$coef[2,2]
microbiome_stats[21,4]<-model11$coef[3,1]
microbiome_stats[22,4]<-model11$coef[3,2]
microbiome_stats[21,6]<-model11$coef[4,1]
microbiome_stats[22,6]<-model11$coef[4,2]
microbiome_stats[21,7]<-model11$r.squared[1]
microbiome_stats[21,8]<-model11$F.test[1]
microbiome_stats[21,9]<-model11$F.test[2]


microbiome_stats[23,3]<-model12$coef[2,1]
microbiome_stats[24,3]<-model12$coef[2,2]
microbiome_stats[23,4]<-model12$coef[3,1]
microbiome_stats[24,4]<-model12$coef[3,2]
microbiome_stats[23,5]<-model12$coef[4,1]
microbiome_stats[24,5]<-model12$coef[4,2]
microbiome_stats[23,7]<-model12$r.squared[1]
microbiome_stats[23,8]<-model12$F.test[1]
microbiome_stats[23,9]<-model12$F.test[2]


microbiome_stats[25,3]<-model13$coef[2,1]
microbiome_stats[26,3]<-model13$coef[2,2]
microbiome_stats[25,6]<-model13$coef[3,1]
microbiome_stats[26,6]<-model13$coef[3,2]
microbiome_stats[25,5]<-model13$coef[4,1]
microbiome_stats[26,5]<-model13$coef[4,2]
microbiome_stats[25,7]<-model13$r.squared[1]
microbiome_stats[25,8]<-model13$F.test[1]
microbiome_stats[25,9]<-model13$F.test[2]

microbiome_stats[27,4]<-model14$coef[2,1]
microbiome_stats[28,4]<-model14$coef[2,2]
microbiome_stats[27,6]<-model14$coef[3,1]
microbiome_stats[28,6]<-model14$coef[3,2]
microbiome_stats[27,5]<-model14$coef[4,1]
microbiome_stats[28,5]<-model14$coef[4,2]
microbiome_stats[27,7]<-model14$r.squared[1]
microbiome_stats[27,8]<-model14$F.test[1]
microbiome_stats[27,9]<-model14$F.test[2]

microbiome_stats[29,3]<-model15$coef[2,1]
microbiome_stats[30,3]<-model15$coef[2,2]
microbiome_stats[29,4]<-model15$coef[3,1]
microbiome_stats[30,4]<-model15$coef[3,2]
microbiome_stats[29,6]<-model15$coef[4,1]
microbiome_stats[30,6]<-model15$coef[4,2]
microbiome_stats[29,5]<-model15$coef[5,1]
microbiome_stats[30,5]<-model15$coef[5,2]
microbiome_stats[29,7]<-model15$r.squared[1]
microbiome_stats[29,8]<-model15$F.test[1]
microbiome_stats[29,9]<-model15$F.test[2]



#Table S1
#write.csv(microbiome_stats, 'microbiome_mantel_tests.csv')


################################## FIGRUES ####################################################################################
################################## FIGRUES ####################################################################################
################################## FIGRUES ####################################################################################
################################## FIGRUES ####################################################################################
################################## FIGRUES ####################################################################################
################################## FIGRUES ####################################################################################
################################## FIGRUES ####################################################################################



pairwise_distances_mb<-data.frame(mb_meandist)
pairwise_distances_mb$Island_from<-rownames(pairwise_distances_mb)
pairwise_distances_mb <- gather(pairwise_distances_mb, Island_to, Distance_mb, Champion:StCruz, factor_key=TRUE)


pairwise_distances_core_mb<-data.frame(mb_core_meandist)
pairwise_distances_core_mb$Island_from<-rownames(pairwise_distances_core_mb)
pairwise_distances_core_mb <- gather(pairwise_distances_core_mb, Island_to, Distance_core_mb, Champion:StCruz, factor_key=TRUE)


pairwise_distances_microsat<-data.frame(microsat_meandist)
pairwise_distances_microsat$Island_from<-rownames(pairwise_distances_microsat)
pairwise_distances_microsat <- gather(pairwise_distances_microsat, Island_to, Distance_microsat, Champion:StCruz, factor_key=TRUE)


pairwise_distances_mhc<-data.frame(mhc_meandist)
pairwise_distances_mhc$Island_from<-rownames(pairwise_distances_mhc)
pairwise_distances_mhc <- gather(pairwise_distances_mhc, Island_to, Distance_mhc, Champion:StCruz, factor_key=TRUE)


pairwise_distances_km<-data.frame(geographic_distance)
pairwise_distances_km$Island_from<-rownames(pairwise_distances_km)
pairwise_distances_km <- gather(pairwise_distances_km, Island_to, Distance_km, Champion:StCruz, factor_key=TRUE)

pairwise_distances<-cbind(pairwise_distances_mb, pairwise_distances_core_mb, pairwise_distances_mhc, pairwise_distances_microsat, pairwise_distances_km)

pairwise_distances<-pairwise_distances[,c(1,2,3,6,9,12,15)]

metadata<-data.frame(sample_data(mockingbird_core))

metadata<-metadata[,c(8,4)]
metadata<-distinct(metadata)

library(expss)
pairwise_distances$Species_from<-vlookup(pairwise_distances$Island_from, metadata, lookup_column = "Island", result_column = "Species")
pairwise_distances$Species_to<-vlookup(pairwise_distances$Island_to, metadata, lookup_column = "Island", result_column = "Species")

pairwise_distances$ComparisonType<-ifelse(pairwise_distances$Species_from == pairwise_distances$Species_to, "Within species", "Between species")

str(pairwise_distances)

pairwise_distances<-subset(pairwise_distances, Distance_mb > 0)

pairwise_distances1 <- gather(pairwise_distances, Distance_type, Distance, Distance_mhc:Distance_km, factor_key=TRUE)
data_long

#pairwise_distances2<-subset(pairwise_distances1, Distance_type != "Distance_km")


corr_plot_mb<-ggplot(pairwise_distances1, aes( y = Distance_mb, x = Distance))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, aes(group=ComparisonType, col = ComparisonType))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, linetype = "dashed", colour = "grey")+
  geom_point(aes(col = ComparisonType), size = 2, alpha = 0.7)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(2, "lines"))+
 facet_wrap(~Distance_type, scales = "free", ncol = 4)+
  scale_color_manual(values =c("yellow3","darkgreen"))+
  labs(y= "Full microbiome dissimilarity \n (Unweighted Unifrac)")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=14))

corr_plot_mb


corr_plot_core<-ggplot(pairwise_distances1, aes( y = Distance_core_mb, x = Distance))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, aes(group=ComparisonType, col = ComparisonType))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, linetype = "dashed", colour = "grey")+
  geom_point(aes(col = ComparisonType), size = 2, alpha = 0.7)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(2, "lines"))+
  facet_wrap(~Distance_type, scales = "free", ncol = 4)+
  # theme(strip.text.x = element_text(size=12),strip.background = element_blank()) +
  scale_color_manual(values =c("yellow3","darkgreen"))+
  labs(y= "Core microbiome dissimilarity \n (Unweighted Unifrac)")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=14))

corr_plot_core



grid.arrange(corr_plot_mb, corr_plot_core)

################################################### fig S4 ###################
################################################### fig S4 ###################
################################################### fig S4 ###################
################################################### fig S4 ###################

pairwise_distances

MHC_micro<-ggplot(pairwise_distances, aes( y = Distance_mhc, x = Distance_microsat))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, aes(group=ComparisonType, col = ComparisonType))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, linetype = "dashed", colour = "grey")+
  geom_point(aes(col = ComparisonType), size = 2, alpha = 0.7)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(2, "lines"))+
 # facet_wrap(~Distance_type, scales = "free", ncol = 4)+
  # theme(strip.text.x = element_text(size=12),strip.background = element_blank()) +
  scale_color_manual(values =c("yellow3","darkgreen"))+
 labs(y= "MHC distance", x = "Microsatellite distance")+
 # theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=14))

MHC_km<-ggplot(pairwise_distances, aes( y = Distance_mhc, x = Distance_km))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, aes(group=ComparisonType, col = ComparisonType))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, linetype = "dashed", colour = "grey")+
  geom_point(aes(col = ComparisonType), size = 2, alpha = 0.7)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(2, "lines"))+
  # facet_wrap(~Distance_type, scales = "free", ncol = 4)+
  # theme(strip.text.x = element_text(size=12),strip.background = element_blank()) +
  scale_color_manual(values =c("yellow3","darkgreen"))+
  labs(y= "MHC distance", x = "Geographic distance (km)")+
  # theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=14))

micro_km<-ggplot(pairwise_distances, aes( y = Distance_microsat, x = Distance_km))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, aes(group=ComparisonType, col = ComparisonType))+
  stat_smooth(method = "lm", alpha=0.7, se=FALSE, linetype = "dashed", colour = "grey")+
  geom_point(aes(col = ComparisonType), size = 2, alpha = 0.7)+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(2, "lines"))+
  # facet_wrap(~Distance_type, scales = "free", ncol = 4)+
  # theme(strip.text.x = element_text(size=12),strip.background = element_blank()) +
  scale_color_manual(values =c("yellow3","darkgreen"))+
  labs(y= "Microsatellite distance", x = "Geographic distance (km)")+
  # theme(axis.title.x = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=14))

grid.arrange(MHC_micro, MHC_km, micro_km, ncol = 3)


####################################################
####################################################
####################################################


mbuni_mds<-ggplot(mb_df_axis_full)+
  geom_point(aes(x=Axis.1,y=Axis.2, color=Species, shape=Island), size = 2.5) + 
  stat_ellipse(aes(x=Axis.1,y=Axis.2,color=Species),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("royalblue3","orange","darkolivegreen3","darkorchid3")) + 
  scale_shape_manual(values=c(16,15,17,18,8,16,15,16,17)) +  
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Full microbiome (Unweighted Unifrac)") + 
  labs(x = "Axis 1 [10.4 %]", y ="Axis 2 [8.3 %]")+
  theme(legend.position = "none")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=18))+
  ylim(c(-0.5, 0.6))


mbuni_mds

mbuni_mds_core<-ggplot(mb_df_axis_core)+
  geom_point(aes(x=Axis.1,y=Axis.2, color=Species, shape=Island), size=2.5) + 
  stat_ellipse(aes(x=Axis.1,y=Axis.2,color=Species),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("royalblue3","orange","darkolivegreen3","darkorchid3")) + 
  scale_shape_manual(values=c(16,15,17,18,8,16,15,16,17)) +  
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Island core microbiome (Unweighted Unifrac)") + 
  labs(x = "Axis 1 [26.8 %]", y ="Axis 2 [16.0 %]")+
  theme(legend.position = "none")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=18))+
  ylim(c(-0.8, 0.4))


mbuni_mds_core


mhcuni_mds<-ggplot(mhc_df_axis)+
  geom_point(aes(x=Axis.1,y=Axis.2, color=Species, shape=Island), size=2.5) + 
  stat_ellipse(aes(x=Axis.1,y=Axis.2,color=Species),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("royalblue3","orange","darkolivegreen3","darkorchid3")) + 
  scale_shape_manual(values=c(16,15,17,18,8,16,15,16,17)) +  
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("MHC") + 
  labs(x = "Axis 1 [22 %]", y ="Axis 2 [16.4 %]")+
  theme(legend.position = "none")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size=18))+
  ylim(c(-0.65, 0.75))



mhcuni_mds


micromds<-ggplot(pcadf)+
  geom_point(aes(x=Axis1,y=Axis2, color=Species, shape=Island), size=2.5) + 
  stat_ellipse(aes(x=Axis1,y=Axis2,color=Species),type = "norm", size = 1) +
  theme_bw() + 
  scale_color_manual(values=c("royalblue3","orange","darkolivegreen3","darkorchid3")) + 
  scale_shape_manual(values=c(16,15,17,18,8,16,15,16,17)) +  
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) +
  ggtitle("Microsatellites") + 
  labs(x = "Axis 1 [24.7 %]", y ="Axis 2 [13.8 %]")+
  theme(legend.position = "none")+
  theme(axis.text = element_text( size=12), axis.title = element_text( size=14))+
  theme(plot.title = element_text(size = 18))

micromds 

grid.arrange(mbuni_mds, mbuni_mds_core, mhcuni_mds, micromds)


##################################



############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################
############################################################ BEDASSLE ANALYSIS #######################################



setwd("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS")


microbiome<-readRDS("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\PHYLOSEQ OBJECTS\\mockingbird_unrarefied.rds")

microbiome<- rarefy_even_depth(microbiome, sample.size = 10000,rngsee = 100, replace = TRUE, trimOTUs=TRUE,verbose=TRUE)

summary(taxa_sums(microbiome))

head(sample_data(microbiome))


#remove singletons

microbiome_1 <- filter_taxa(microbiome, function (x) {sum(x > 0) > 1}, prune=TRUE)

microbiome_2 <- core(microbiome_1, detection = 0, prevalence = .05)

#covert to presence absence

microbiome_pa<-phyloseq_standardize_otu_abundance(microbiome_2, method = "pa")

################# merge by island

otu_island <- merge_samples(microbiome_pa, "Island")

otu_table(otu_island)

otu1<-data.frame(otu_table(otu_island)) #this is t3

otu1<-as.matrix(otu1)

#################

metadata<-data.frame(sample_data(microbiome_pa))

sample_sizes<-data.frame(table(metadata$Island))



#t5 = sample sizes (individuals per population)

sample_size_matrix <- as.data.frame(matrix(ncol = ncol(otu1), nrow = nrow(otu1)))
rownames(sample_size_matrix)  <- sample_sizes$Var1
colnames(sample_size_matrix) <- colnames(otu1)

#rownames(t5)

sample_size_matrix[,1:length(sample_size_matrix)] <- sample_sizes$Freq 

sample_size_matrix<-as.matrix(sample_size_matrix)

str(otu1) #1228 ASVs
str(sample_size_matrix)

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

####################### island core microbiome

mockingbird_core<-readRDS("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\mockingbird_island_core_css.rds")

mockingbird_core

#covert to presence absence

microbiome_pa<-phyloseq_standardize_otu_abundance(mockingbird_core, method = "pa")

################# merge by island

otu_island <- merge_samples(microbiome_pa, "Island")

otu_table(otu_island)

otu_core<-data.frame(otu_table(otu_island)) #this is t3

otu_core<-as.matrix(otu_core)

#################

metadata<-data.frame(sample_data(microbiome_pa))

sample_sizes_core<-data.frame(table(metadata$Island))



#t5 = sample sizes (individuals per population)

sample_size_matrix_core <- as.data.frame(matrix(ncol = ncol(otu_core), nrow = nrow(otu_core)))
rownames(sample_size_matrix_core)  <- sample_sizes_core$Var1
colnames(sample_size_matrix_core) <- colnames(otu_core)

#rownames(t5)

sample_size_matrix_core[,1:length(sample_size_matrix_core)] <- sample_sizes_core$Freq 

sample_size_matrix_core<-as.matrix(sample_size_matrix_core)

str(otu_core) #267 ASVs island core 
str(sample_size_matrix_core)

######################################################################################



######################################################################################

## D: Pairwise geographic distance between populations (in km)

geographic_distance1

################################################

## E: Pairwise ecological distance(s). Make them into a list.

#we have two genetic variables - mhc and microsat

##MICROSAT 


microsat_meandist1

# or microsat_centroids_meandist1

#################### MHC

mhc_meandist1

##make into list
genetics_list <- list(microsat_meandist1, mhc_meandist1)

names(genetics_list)<-c('MICROSAT', 'MHC')


################################################################
################################################################
################################################################
################################################################

#MB_FINAL is for <2% prevalence 15 mill iterations
#CORE_FINAL is for 30% species core and 15 mill iterations

#Call the Markov chain Monte Carlo for the beta binomial model
# FULL MICROBIOME AND USING MICROSAT DISTANCE MATRIX FROM TREE

MCMC_BB(
  counts = otu1,
  sample_sizes = sample_size_matrix,
  D = geographic_distance1,
  E = genetics_list,
  k = 9,
  loci = 1228,
  delta = 0.000001, 
  aD_stp = 0.01, 
  aE_stp = 0.1, 
  a2_stp = 0.03, 
  thetas_stp = 0.1,
  phi_stp = 30, 
  mu_stp = 0.3, 
  ngen = 15000000,  
  printfreq = 10000,
  savefreq = 10000,
  samplefreq = 100,
  directory =  "C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\BEDASSLE_OUTPUT",
  prefix = "MB_1228_15mil_",
  continue = FALSE,
  continuing.params = NULL
)



######################################################################################
######################################################################################
######################################################################################

############################################################################

setwd("C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\BEDASSLE_OUTPUT")

show(load("MB_1228_15mil_MCMC_output1.Robj"))



## run below code on local machine
plot_all_acceptance_rates('C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\bedassle_output\\
                          MB_1228_15mil_MCMC_output1.Robj') #shown in viewer

#make a vector popuilation names
pop_vector <- as.vector(rownames(sample_size_matrix))
plot_all_trace('MB_1228_15mil_MCMC_output1.Robj', percent.burnin = 0, thinning = 1, population.names = pop_vector)

posterior.predictive.sample(MCMC.output = 
                              "C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\bedassle_output\\MB_1228_15mil_MCMC_output1.Robj", 
                            posterior.predictive.sample.size = 1000, output.file = "post1000", prefix = "MB_1228_15mil_")

plot_posterior_predictive_samples("MB_1228_15mil_post1000.Robj", save.figure = TRUE, figure.name = "MB_1228_15mil_post1000_1.jpg")

#use a 60% burn-in 
plot_all_trace('MB_1228_15mil_post1000.Robj', percent.burnin = 60, thinning = 1, population.names = pop_vector)

################################################
################################################
################################################


#correct the aE : aD ratios
#1 is geology
microsat_ratio <- ((as.vector(aE[1,])* sd(geographic_distance))/(aD* sd(microsat_meandist)))
#2 is soil pc1
mhc_ratio <- ((as.vector(aE[2,])* sd(geographic_distance))/(aD * sd(mhc_meandist)))

#calculate output after 60% burn-in
#library(coda)
all.param = as.mcmc(cbind(microsat_ratio[-(1:90000)], mhc_ratio[-(1:90000)]))
s1 <- summary(all.param)
s2 <- as.data.frame(s1$statistics)
s2 <- s2[1:2,1:2]
s3 <- as.data.frame(HPDinterval(all.param, prob=0.95))
s3 <- cbind(s2, s3)
rownames(s3) <- c("microsat_ratio", "mhc_ratio")
s3 <- round(s3, digits = 5)
s3
str(s3)

quantile_microsat<-quantile(microsat_meandist)
quantile_mhc<-quantile(mhc_meandist)

s3$quantile_diff<-c((quantile_microsat[4]-quantile_microsat[2]),(quantile_mhc[4]-quantile_mhc[2]))

s3$normalized<-s3$Mean * s3$quantile_diff

MB_RESULTS<-s3
MB_RESULTS


############################################## core microbiome ##################
############################################## core microbiome ##################
############################################## core microbiome ##################
############################################## core microbiome ##################
############################################## core microbiome ##################

#Call the Markov chain Monte Carlo for the beta binomial model
# core island microbiome
MCMC_BB(
  counts = otu_core,
  sample_sizes = sample_size_matrix_core,
  D = geographic_distance1,
  E = genetics_list,
  k = 9,
  loci = 267,
  delta = 0.000001, 
  aD_stp = 0.01, 
  aE_stp = 0.1, 
  a2_stp = 0.03, 
  thetas_stp = 0.1,
  phi_stp = 30, 
  mu_stp = 0.3, 
  ngen = 20000000,  
  printfreq = 10000,
  savefreq = 10000,
  samplefreq = 100,
  directory =  "C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\BEDASSLE_OUTPUT",
  prefix = "ISLAND_CORE_CSS_20MIL_",
  continue = FALSE,
  continuing.params = NULL
)



##################################################


show(load("ISLAND_CORE_CSS_20MIL_MCMC_output1.Robj"))



## run below code on local machine
plot_all_acceptance_rates('C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\bedassle_output\\ISLAND_CORE_CSS_20MIL_MCMC_output1.Robj') #shown in viewer

#make a vector popuilation names
pop_vector <- as.vector(rownames(sample_size_matrix_core))
plot_all_trace('ISLAND_CORE_CSS_20MIL_MCMC_output1.Robj', percent.burnin = 0, thinning = 1, population.names = pop_vector)

posterior.predictive.sample(MCMC.output = 
                              "C:\\Users\\risel\\Dropbox\\Sommer postdoc\\Ramona mockingbirds\\Mockingbirds analysis\\FINAL ANALYSIS\\bedassle_output\\ISLAND_CORE_CSS_20MIL_MCMC_output1.Robj", 
                            posterior.predictive.sample.size = 1000, output.file = "post1000", prefix = "ISLAND_CORE_CSS_20MIL_")

plot_posterior_predictive_samples("ISLAND_CORE_CSS_20MIL_post1000.Robj", save.figure = TRUE, figure.name = "ISLAND_CORE_CSS_20MIL_post1000.jpg")

#use a 60% burn-in 
plot_all_trace('ISLAND_CORE_CSS_20MIL_MCMC_output1.Robj', percent.burnin = 60, thinning = 1, population.names = pop_vector)




geographic_distance

#correct the aE : aD ratios
#1 is geology
microsat_ratio <- ((as.vector(aE[1,])* sd(geographic_distance))/(aD* sd(microsat_meandist)))
#2 is soil pc1
mhc_ratio <- ((as.vector(aE[2,])* sd(geographic_distance))/(aD * sd(mhc_meandist)))

#calculate output after 60% burn-in
#library(coda)
all.param = as.mcmc(cbind(microsat_ratio[-(1:120000)], mhc_ratio[-(1:120000)]))
s1 <- summary(all.param)
s2 <- as.data.frame(s1$statistics)
s2 <- s2[1:2,1:2]
s3 <- as.data.frame(HPDinterval(all.param, prob=0.95))
s3 <- cbind(s2, s3)
rownames(s3) <- c("microsat_ratio", "mhc_ratio")
s3 <- round(s3, digits = 5)
s3
str(s3)

quantile_microsat<-quantile(microsat_meandist)
quantile_mhc<-quantile(mhc_meandist)

s3$quantile_diff<-c((quantile_microsat[4]-quantile_microsat[2]),(quantile_mhc[4]-quantile_mhc[2]))

s3$normalized<-s3$Mean * s3$quantile_diff

CORE_RESULTS<-s3

CORE_RESULTS


###################################

