

### LOAD PACKAGES ###
library(multcomp)
library(tidyverse)
library(codyn)
library(vegan)
library(abdiv)
library(lme4)
library(lmerTest)
library(emmeans)
# pairwiseAdonis needs R version 4.2.1 #
install.packages("devtools")
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)


# Read in data for treatments #
treatments <- read.csv("NutientLimitationLongleafPine2024/NutNet_treatments.csv")%>%
  mutate(treatment=paste(N,P,K,drought,sep="_")) %>%
  mutate(Nutrients=paste(N,P,K,sep="_"))

# Read in data for taxa lifeform #
lifeform <- read.csv("NutientLimitationLongleafPine2024/NutNet_taxa_descriptions.csv")

# Read in data for biomass #
biomass <- read.csv("NutientLimitationLongleafPine2024/NPKDNet_Paper_Biomass.csv")%>%
  spread(functionalgroup,mass)%>%
  rename(forbbiomass = Forb) %>%
  rename(graminoidbiomass = Graminoid) %>%
  rename(woodybiomass = Woody) %>%
  rename(legumebiomass = Legume) %>%
  rename(litterbiomass = Litter) %>%
  rename(pineneedlebiomass = Pine_Needles) %>%
  group_by(year, block, plot) %>%
  mutate(totalbiomass = sum(forbbiomass, graminoidbiomass, woodybiomass, legumebiomass, na.rm=T)) %>%
  mutate(legumebiomass=ifelse(is.na(legumebiomass), 0, legumebiomass)) %>%
  mutate(forbbiomass=ifelse(is.na(forbbiomass), 0, forbbiomass)) %>%
  mutate(woodybiomass=ifelse(is.na(woodybiomass), 0, woodybiomass)) %>%
  dplyr::select(-site,-Bryophyte,-PYD,-pineneedlebiomass,-litterbiomass)
biomass <- merge(biomass, treatments, by=c("block", "plot"),all=T)
biomass[biomass$year=="2022" & biomass$plot=="29", "forbbiomass"] <- NA
biomass[biomass$year=="2022" & biomass$plot=="29", "legumebiomass"] <- NA
biomass[biomass$year=="2022" & biomass$plot=="29", "woodybiomass"] <- NA
biomass[biomass$year=="2022" & biomass$plot=="29", "totalbiomass"] <- NA


# Read in Data Set for SpComp
speciescomposition <- read.csv("NutientLimitationLongleafPine2024/NPKDNet_SppComp.csv") %>%
  # remove mid-season 2022 data and remove NutNet taxa that aren't vascular plants #
  filter(timepoint!="Mid" & !taxa%in% c("Bare_ground","Overstory","Rock","Litter", "Unknown_mushroom", "Unknown_red_mushroom", "Unknown_moss","Unknown_lichen")) %>%
  # join with NutNet lifeform and lifespan df #
  full_join(lifeform, by="taxa") %>%
  select(-site, -date, -note_cover)

# Subset df to remove pre-treatment data (2019)
# wrangle species composition data for calculating metrics #
SpComp <- subset(speciescomposition, year!="2019") %>%
  # average cover in early and late sampling time points (keep cover for spp in only 1 season) # 
  group_by(block, plot, year, lifeform, PhotoGram, taxa) %>%
  reframe(maxcover = max(cover)) 

# Check for taxa typos #
#MaxUnique <- MaxSpComp %>%
 # reframe(Taxa = unique(taxa))

# calculate Berger-Parker dominance #
MaxBPdom <- SpComp %>%
  group_by(year, block, plot) %>%
  reframe(MaxBPdom = berger_parker_d(maxcover))

# calculate wiregrass cover #
Maxwiregrasscover <- SpComp %>%
  filter(taxa %in% c("Aristida_stricta"))%>%
  group_by(year, block, plot) %>%
  reframe(Maxwirecov = maxcover)

# calculate fern cover #
Maxferncover <- SpComp %>%
  filter(taxa %in% c("Pteridium_aquilinum"))%>%
  group_by(year, block, plot) %>%
  reframe(Maxferncov = maxcover)

# calculate C4 grass cover #
MaxC4cover <- SpComp %>%
  filter(PhotoGram %in% c("C4"))%>%
  group_by(year, block, plot) %>%
  reframe(MaxC4cov = sum(maxcover))

# calculate C3 grass cover #
MaxC3cover <- SpComp %>%
  filter(PhotoGram %in% c("C3"))%>%
  group_by(year, block, plot) %>%
  reframe(MaxC3cov = sum(maxcover))

# calculate functional group cover #
FunctionalCover <- SpComp %>%
  group_by(year, block, plot, lifeform)%>%
  summarise(maxcover = sum(maxcover))%>%
  group_by(year, block, plot)%>%
  spread(lifeform,maxcover)%>%
  mutate(Woody=sum(Shrub,Tree,Vine, na.rm=T))

# identify dominant species #
#Maxdominantspp <- MaxSpComp %>%
  #group_by(year, block, plot) %>%
  #mutate(dominant = max(maxcover)) %>%
  #reframe(Maxdomspp1=ifelse(dominant==maxcover, taxa, NA)) %>%
  #drop_na()
#write.csv(Maxdominantspp, "NutientLimitationLongleafPine2024\\MaxDominantSpecies.csv", row.names=FALSE)
# Read in data for dominant species #
Maxdominantspp <- read.csv("NutientLimitationLongleafPine2024/MaxDominantSpecies.csv",na.strings="")


# calculate biodiversity metrics #
MaxCommStructure <- SpComp%>%
  community_structure(time.var="year", abundance.var="maxcover",replicate.var="plot", metric=c("Evar")) %>%
  rename(MaxRichness = richness, MaxEvar=Evar)

#MaxCommStructure <- MaxSpComp%>%
  #group_by(year, block, plot)%>%
  #summarise(MaxRichness = length(unique(taxa)))

MaxCommDiv <- SpComp %>%
  community_diversity(time.var="year", abundance.var="maxcover", replicate.var="plot", metric=c("Shannon")) %>%
  rename(MaxShannon = Shannon)

# merge community structure and diversity dfs with all other metric dfs#
MaxCommMetrics <- full_join(MaxCommDiv,MaxCommStructure, by=c("year", "plot")) %>%
  full_join(MaxBPdom, by=c("year", "plot")) %>%
  full_join(Maxwiregrasscover, by=c("year", "block", "plot")) %>%
  full_join(Maxferncover, by=c("year", "block", "plot")) %>%
  full_join(MaxC4cover, by=c("year", "block", "plot")) %>%
  full_join(MaxC3cover, by=c("year", "block", "plot")) %>%
  full_join(FunctionalCover, by=c("year", "block", "plot")) %>%
  full_join(Maxdominantspp, by=c("year", "block", "plot")) %>%
  full_join(biomass, by=c("year", "block", "plot"))
MaxCommMetrics$year <- as.factor(MaxCommMetrics$year)
MaxCommMetrics$block <- as.factor(MaxCommMetrics$block)
MaxCommMetrics$plot <- as.factor(MaxCommMetrics$plot)
MaxCommMetrics$N <- as.factor(MaxCommMetrics$N)
MaxCommMetrics$P <- as.factor(MaxCommMetrics$P)
MaxCommMetrics$K <- as.factor(MaxCommMetrics$K)
MaxCommMetrics$drought <- as.factor(MaxCommMetrics$drought)
MaxCommMetrics$treatment <- as.factor(MaxCommMetrics$treatment)

# subset to keep just drought plots and controls #
DroughtMax <- filter(MaxCommMetrics, treatment=="C_C_C_C" | treatment=="C_C_C_Drought" | treatment=="Nitrogen_Phosphorus_Potassium_C" | 
                       treatment=="Nitrogen_Phosphorus_Potassium_Drought") %>%
  mutate(Nutrients=paste(N,P,K,sep="_"))
## DROP 2019 FROM ALL ANALYSES ##
DroughtMax <- subset(DroughtMax, year!="2019")   


# subset out all droughted plots and remove drought from model #
NoDroughtMax <- subset(MaxCommMetrics, drought!="Drought") %>%
  mutate(treatment=ifelse(treatment=="C_C_C_C","Control", ifelse(treatment=="Nitrogen_C_C_C", "N",
              ifelse(treatment=="C_Phosphorus_C_C", "P", ifelse(treatment=="C_C_Potassium_C", "K",
           ifelse(treatment=="Nitrogen_Phosphorus_C_C","NP",ifelse(treatment=="Nitrogen_C_Potassium_C", "NK", 
            ifelse(treatment=="C_Phosphorus_Potassium_C", "PK", ifelse(treatment=="Nitrogen_Phosphorus_Potassium_C", "NPK", NA)))))))))

## DROP 2019 FROM ALL ANALYSES ##
NoDroughtMax <- subset(NoDroughtMax, year!="2019") 


############################################################################################################


################################################################################################################################
## Plant Species Composition Ordinations and Stats ####

MaxlongFormat <- SpComp %>%
  select(-lifeform, -PhotoGram) %>%
  group_by(year, block, plot) %>%
  # convert df to wide format to account for pres/abs of all species for each plot #
  spread(taxa, maxcover) %>%
  # convert back to long format so spp recorded as NA that were not in a plot #
  gather(taxa, maxcover, 4:65)
# if spp was absent in plot, give a zero cover #
MaxlongFormat[is.na(MaxlongFormat)]<-0


wide <- merge(MaxlongFormat, treatments, by=c("block", "plot"))
wideFormat <- subset(wide, drought!="Drought" & Nutrients=="C_C_C" | Nutrients=="Nitrogen_C_C" | Nutrients=="C_Phosphorus_C" | Nutrients=="C_C_Potassium") %>%
  select(-treatment) %>%
  group_by(year, block, plot)%>%
  spread(taxa, maxcover)

#wideFormat <- subset(wideFormat, year=="2023")
  

Maxord <- wideFormat[,9:ncol(wideFormat)]
Maxord1<- wideFormat[,1:8]
Maxord2<- metaMDS(Maxord,distance="bray", k=2, trymax=30)
MaxNMDS <- data.frame(Maxord1,scores(Maxord2,display="sites")) %>%
  rename(MaxNMDS1 = NMDS1, MaxNMDS2 = NMDS2)

## PerMANOVA to test for significant differences between centroids (treatments) ##
# year*nutrients not significant #
PerMANOVA <- adonis2(formula = Maxord ~ Nutrients, 
                     data=Maxord1,permutations = 999, strata = Maxord1$block, method = "bray")
PerMANOVA

Posthoc_NPK<-pairwise.adonis(Maxord,factors=Maxord1$Nutrients, p.adjust.m = "BH")
Posthoc_NPK 

#Make a new dataframe and calculate the dissimilarity of the Species_Matrix dataframe
BC_Distance_Matrix <- vegdist(Maxord)
#Run a dissimilarity matrix (PermDisp) 
Dispersion <- betadisper(BC_Distance_Matrix,Maxord1$Nutrients)

anova(Dispersion)
#another way to do the anova #
adonis2(dist(Dispersion$distances) ~ Maxord1$Nutrients)
TukeyHSD(Dispersion)
# an alternative to TukeyHSD
permutest(Dispersion,pairwise = T, permutations = 999)


########################################################################################################################
          ## STATISTICS ####
# mixed effect models accounting for repeated measures #
##########################################################################
rich <- lmerTest::lmer(data=NoDroughtMax, MaxRichness ~ year * N * P * K + (1|block) + (1|plot))
anova(rich)
emmeans(rich, pairwise ~ P*K, adjust="BH")
####
# Drought x Nutrients (NPK together) #
####
rich2 <- lmerTest::lmer(data=DroughtMax, MaxRichness ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(rich2)


div <- lmerTest::lmer(data=NoDroughtMax, MaxShannon ~ year * N * P * K + (1|block) + (1|plot))
anova(div)
emmeans(div, pairwise ~ N|year, adjust="BH")
emmeans(div, pairwise ~ P*K, adjust="BH")
####
# Drought x Nutrients (NPK together) #
####
div2 <- lmerTest::lmer(data=DroughtMax, MaxShannon ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(div2)
emmeans(div2, pairwise ~ Nutrients|year, adjust="BH")


BPdom <- lmerTest::lmer(data=NoDroughtMax, MaxBPdom ~ year * N * P * K + (1|block) + (1|plot))
anova(BPdom)
emmeans(BPdom, pairwise ~ N|year, adjust="BH")
emmeans(BPdom, pairwise ~ N*P*K, adjust="BH")
####
# Drought x Nutrients (NPK together) #
#### 
BPdom2 <- lmerTest::lmer(data=DroughtMax, MaxBPdom ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(BPdom2)


TotalBio <- lmerTest::lmer(data=NoDroughtMax, totalbiomass ~ year * N * P * K + (1|block) + (1|plot))
anova(TotalBio)
emmeans(TotalBio, pairwise ~ N|year, adjust="BH")
emmeans(TotalBio, pairwise ~ K|year, adjust="BH")
emmeans(TotalBio, pairwise ~ N*P|year, adjust="BH")
emmeans(TotalBio, pairwise ~ N*K|year, adjust="BH")
####
# Drought x Nutrients (NPK together) #
####
TotalBio2 <- lmerTest::lmer(data=DroughtMax, totalbiomass ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(TotalBio2)


GramBio <- lmerTest::lmer(data=NoDroughtMax, graminoidbiomass ~ year * N * P * K + (1|block) + (1|plot))
anova(GramBio)
emmeans(GramBio, pairwise ~ N|year, adjust="BH")
emmeans(GramBio, pairwise ~ K|year, adjust="BH")
emmeans(GramBio, pairwise ~ N*K|year, adjust="BH")
emmeans(GramBio, pairwise ~ N*P|year, adjust="BH")
####
# Drought x Nutrients (NPK together) #
####
GramBio2 <- lmerTest::lmer(data=DroughtMax, graminoidbiomass ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(GramBio2)


ForbBio <- lmerTest::lmer(data=NoDroughtMax, forbbiomass ~ year * N * P * K + (1|block) + (1|plot))
anova(ForbBio)
emmeans(ForbBio, pairwise ~ N|year, adjust="BH")
emmeans(ForbBio, pairwise ~ P*K|year, adjust="BH")
####
# Drought x Nutrients (NPK together) #
####
ForbBio2 <- lmerTest::lmer(data=DroughtMax, forbbiomass ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(ForbBio2)




######################################################################################################################################

## NITROGEN EFFECTS ####
####

## Fig 2a ##
NRich <- NoDroughtMax %>%
  group_by(N) %>%
  reframe(mean = mean(MaxRichness),Sd = sd(MaxRichness),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(NRich, aes(x=N, y=mean, fill=N)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") + ylab("Richness") +
  scale_x_discrete(labels=c("C", "N")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 2b ##
Nshannon <- NoDroughtMax %>%
  group_by(year,N) %>%
  reframe(mean = mean(MaxShannon),Sd = sd(MaxShannon),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Nshannon, aes(x=year, y=mean, group=N,color=N, linetype=N)) +
  geom_line(size=3) + geom_point(aes(color=N),size=10) +
  geom_errorbar(aes(color=N,ymin=mean-Se,ymax=mean+Se),size=1.5, width=0.5)+
  annotate("text", x= 3, y = 2.41, label= "*", size = 20)+ 
  annotate("text", x= 4, y = 2.27, label= "*", size = 20)+
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  scale_x_discrete(labels=c("1", "2", "3", "4")) +
  ylab("Shannon Diversity") + xlab("Year") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(size = 15),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 40), 
        legend.title=element_blank(),legend.position="none",axis.text.y=element_text(size = 40))

## Fig 2c ##
Ndom <- NoDroughtMax %>%
  group_by(year,N) %>%
  reframe(mean = mean(MaxBPdom),Sd = sd(MaxBPdom),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Ndom, aes(x=year, y=mean, group=N,color=N, linetype=N)) +
  geom_line(size=3) + geom_point(aes(color=N),size=10) + 
  geom_errorbar(aes(color=N,ymin=mean-Se,ymax=mean+Se),size=1.5, width=0.5)+
  annotate("text", x= 3, y = 0.33, label= "**", size = 20)+
  annotate("text", x= 4, y = 0.35, label= "**", size = 20)+
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  scale_x_discrete(labels=c("1", "2", "3", "4")) +
  ylab("B-P Dominance") + xlab("Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(size = 15),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 40), 
        legend.title=element_blank(),legend.position="none",axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 4a ##
Ntotalbiomass <- NoDroughtMax %>%
  group_by(year,N) %>%
  reframe(mean = mean(totalbiomass, na.rm=T),Sd = sd(totalbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Ntotalbiomass, aes(x=year, y=mean, group=N,color=N, linetype=N)) +
  geom_line(size=3) + geom_point(aes(color=N),size=10) + 
  geom_errorbar(aes(color=N,ymin=mean-Se,ymax=mean+Se),size=1.5, width=0.5)+
  annotate("text", x= 2, y = 605, label= "***", size = 20)+ 
  annotate("text", x= 3, y = 425, label= "**", size = 20)+ 
  annotate("text", x= 4, y = 385, label= "**", size = 20)+
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_x_discrete(labels=c("1", "2", "3", "4")) +
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  ylab(bquote("Total Biomass"~(g/m^2))) + xlab ("Year") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size = 40),
        legend.position="none",axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 4b ##
Nforbbiomass <- NoDroughtMax %>%
  group_by(year,N) %>%
  reframe(mean = mean(forbbiomass, na.rm=T),Sd = sd(forbbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Nforbbiomass, aes(x=year, y=mean, group=N,color=N, linetype=N)) +
  geom_line(size=3) + geom_point(aes(color=N),size=10) + 
  geom_errorbar(aes(color=N,ymin=mean-Se,ymax=mean+Se),size=1.5, width=0.5)+
  annotate("text", x= 2, y = 120, label= "**", size = 20)+ 
  annotate("text", x= 3, y = 109, label= "**", size = 20)+ 
  annotate("text", x= 4, y = 111, label= "**", size = 20)+ 
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  scale_x_discrete(labels=c("1", "2", "3", "4")) +
  ylab(bquote("Forb Biomass"~(g/m^2))) + xlab("Year") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 40),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none",axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 4c ##
NGramBio <- NoDroughtMax %>%
  group_by(year,N) %>%
  reframe(mean = mean(graminoidbiomass, na.rm=T),Sd = sd(graminoidbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(NGramBio, aes(x=year, y=mean, group=N,color=N, linetype=N)) +
  geom_line(size=3) + geom_point(aes(color=N),size=10) + 
  geom_errorbar(aes(color=N,ymin=mean-Se,ymax=mean+Se),size=1.5, width=0.5)+
  annotate("text", x= 2, y = 470, label= "***", size = 20)+ ylim(50,500) +
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_color_manual(values=c("#56B4E9","#E69F00")) + ylim(50,500)+
  scale_x_discrete(labels=c("1", "2", "3", "4")) +
  ylab(bquote("Graminoid Biomass"~(g/m^2))) + xlab("Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 40),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none",axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #


####
## PHOSPHORUS EFFECTS ####
####

## Fig 2d ##
PRich <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(MaxRichness),Sd = sd(MaxRichness),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(PRich, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") + ylab("Richness") +
  scale_x_discrete(labels=c("C", "P")) +
  scale_fill_manual(values=c("#56B4E9","#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 2e ##
PShannon <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(MaxShannon),Sd = sd(MaxShannon),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(PShannon, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") + ylab("Shannon Diversity") +
  scale_x_discrete(labels=c("C", "P")) +
  ylim(0,2.5)+ scale_fill_manual(values=c("#56B4E9","#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 2f ##
Pdom <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(MaxBPdom),Sd = sd(MaxBPdom),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Pdom, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") + ylab("B-P Dominance") +
  scale_x_discrete(labels=c("C", "P")) + ylim(0,0.3)+
  scale_fill_manual(values=c("#56B4E9","#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

# Fig 4d ##
Ptotalbio <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(totalbiomass, na.rm=T),Sd = sd(totalbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Ptotalbio, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  annotate("text", x= 1.5, y = 375, label= "*", size = 20)+ 
  scale_x_discrete(labels=c("C", "P")) +
  scale_fill_manual(values=c("#56B4E9","#999999")) +
  xlab("Treatment") + ylab(bquote("Total Biomass"~(g/m^2))) + ylim(0,400) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 4e ##
Pforbbio <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(forbbiomass, na.rm=T),Sd = sd(forbbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Pforbbio, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab(bquote("Forb Biomass"~(g/m^2))) +
  annotate("text", x= 1.5, y = 80, label= "**", size = 20, color="black")+ 
  scale_x_discrete(labels=c("C", "P")) +
  scale_fill_manual(values=c("#56B4E9","#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 4f ##
Pgrambio <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(graminoidbiomass,na.rm=T),Sd = sd(graminoidbiomass,na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Pgrambio, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  ylab(bquote("Graminoid Biomass"~(g/m^2))) +  xlab("Treatment") +
  scale_x_discrete(labels=c("C","P")) +
  scale_fill_manual(values=c("#56B4E9","#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 40),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(size = 40), legend.position="none",axis.text.y=element_text(size = 40))
# 745 x 745 #



####
## POTASSIUM EFFECTS ####
####

## Fig 2g ##
KRich <- NoDroughtMax %>%
  group_by(K) %>%
  reframe(mean = mean(MaxRichness),Sd = sd(MaxRichness),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(KRich, aes(x=K, y=mean, fill=K)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") + ylab("Richness") +
  scale_x_discrete(labels=c("C", "K")) +
  scale_fill_manual(values=c("#56B4E9","palevioletred4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 2h ##
KShannon <- NoDroughtMax %>%
  group_by(K) %>%
  reframe(mean = mean(MaxShannon),Sd = sd(MaxShannon),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(KShannon, aes(x=K, y=mean, fill=K)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") + ylab("Shannon Diversity") +
  scale_x_discrete(labels=c("C", "K")) +
  scale_fill_manual(values=c("#56B4E9","palevioletred4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 2i ##
Kdom <- NoDroughtMax %>%
  group_by(K) %>%
  reframe(mean = mean(MaxBPdom),Sd = sd(MaxBPdom),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Kdom, aes(x=K, y=mean, fill=K)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") + ylab("B-P Dominance") +
  scale_x_discrete(labels=c("C", "K")) + ylim(0,0.3)+
  scale_fill_manual(values=c("#56B4E9","palevioletred4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

##Fig 4g ##
Ktotalbiomass <- NoDroughtMax %>%
  group_by(year,K) %>%
  reframe(mean = mean(totalbiomass, na.rm=T),Sd = sd(totalbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Ktotalbiomass, aes(x=year, y=mean, group=K,color=K, linetype=K)) +
  geom_line(size=3) + geom_point(aes(color=K),size=10) + 
  geom_errorbar(aes(color=K,ymin=mean-Se,ymax=mean+Se),size=1.5, width=0.5)+
  annotate("text", x= 2, y = 570, label= "**", size = 20)+ 
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_color_manual(values=c("#56B4E9","palevioletred4")) +
  scale_x_discrete(labels=c("1", "2", "3", "4")) +
  ylab(bquote("Total Biomass"~(g/m^2))) + xlab("Year") + ylim(175,600) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_text(size = 40),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",axis.text.x=element_text(size = 40))
# 745 x 745 #

## Fig 4h ##
Kforbbio <- NoDroughtMax %>%
  group_by(K) %>%
  reframe(mean = mean(forbbiomass, na.rm=T),Sd = sd(forbbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Kforbbio, aes(x=K, y=mean, fill=K)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab(bquote("Forb Biomass"~(g/m^2))) +
  scale_x_discrete(labels=c("C", "K")) +
  scale_fill_manual(values=c("#56B4E9","palevioletred4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

## Fig 4i ##
KGramBio <- NoDroughtMax %>%
  group_by(year,K) %>%
  reframe(mean = mean(graminoidbiomass, na.rm=T),Sd = sd(graminoidbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(KGramBio, aes(x=year, y=mean, group=K,color=K, linetype=K)) +
  geom_line(size=3) + geom_point(aes(color=K),size=10) + 
  geom_errorbar(aes(color=K,ymin=mean-Se,ymax=mean+Se),size=1.5, width=0.5)+
  annotate("text", x= 2, y = 445, label= "**", size = 20)+ 
  ylim(80,450) +
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_color_manual(values=c("#56B4E9","palevioletred4")) + 
  scale_x_discrete(labels=c("1", "2", "3", "4")) +
  ylab(bquote("Graminoid Biomass"~(g/m^2))) + xlab("Year")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 40),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none",axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #



####
## NITROGEN X PHOSPHORUS EFFECTS ####
####

NPtotalbiomass <- NoDroughtMax %>%
  group_by(year,N,P) %>%
  reframe(mean = mean(totalbiomass, na.rm=T),Sd = sd(totalbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))%>%
  mutate(trt = paste(N,P, sep=("_")))
NPtotalbiomass$trt = factor(NPtotalbiomass$trt, levels=c("C_C","Nitrogen_C", "C_Phosphorus", "Nitrogen_Phosphorus"))

## Fig 5a ##
ggplot(subset(NPtotalbiomass, year=="2020"),aes(x=trt, y=mean, fill=trt)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  ylab(bquote("Total Biomass"~(g/m^2))) + ylim(0,800) +
  scale_x_discrete(labels=c("C", "N", "P", "N x P")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00","#999999","indianred3")) +
  xlab("Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(size = 15),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 30),axis.text.x=element_text(size = 20), 
        legend.position="none",axis.text.y=element_text(size = 25))
# 745 x 745 #

## Fig 5b ##
ggplot(subset(NPtotalbiomass, year=="2021"),aes(x=trt, y=mean, fill=trt)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  annotate("text", x= 1, y = 410, label= "ab", size = 9)+ 
  annotate("text", x= 2, y = 575, label= "b", size = 9)+
  annotate("text", x= 3, y = 285, label= "a", size = 9)+ 
  annotate("text", x= 4, y = 800, label= "c", size = 9)+
  ylab(bquote("Total Biomass"~(g/m^2))) + ylim(0,800) +
  scale_x_discrete(labels=c("C", "N", "P", "N x P")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00","#999999","indianred3")) +
  xlab("Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(size = 15),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 30),axis.text.x=element_text(size = 20), 
        legend.position="none",axis.text.y=element_text(size = 25))
# 745 x 745 #

## Fig 5c ##
ggplot(subset(NPtotalbiomass, year=="2022"),aes(x=trt, y=mean, fill=trt)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  annotate("text", x= 1, y = 300, label= "a", size = 9)+ 
  annotate("text", x= 2, y = 410, label= "ab", size = 9)+
  annotate("text", x= 3, y = 305, label= "a", size = 9)+ 
  annotate("text", x= 4, y = 530, label= "b", size = 9)+
  ylab(bquote("Total Biomass"~(g/m^2))) + ylim(0,800) +
  scale_x_discrete(labels=c("C", "N", "P", "N x P")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00","#999999","indianred3")) +
  xlab("Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(size = 15),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 30),axis.text.x=element_text(size = 20), 
        legend.position="none",axis.text.y=element_text(size = 25))
# 745 x 745 #

## Fig 5d ##
ggplot(subset(NPtotalbiomass, year=="2023"),aes(x=trt, y=mean, fill=trt)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  annotate("text", x= 1, y = 220, label= "a", size = 9)+ 
  annotate("text", x= 2, y = 455, label= "b", size = 9)+
  annotate("text", x= 3, y = 370, label= "ab", size = 9)+ 
  annotate("text", x= 4, y = 420, label= "b", size = 9)+
  ylab(bquote("Total Biomass"~(g/m^2))) + ylim(0,800) +
  scale_x_discrete(labels=c("C", "N", "P", "N x P")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00","#999999","indianred3")) +
  xlab("Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(size = 15),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 30),axis.text.x=element_text(size = 20), 
        legend.position="none",axis.text.y=element_text(size = 25))
# 745 x 745 #

## Fig 6 ##
NPforbbio <- NoDroughtMax %>%
  group_by(N,P) %>%
  reframe(mean = mean(forbbiomass, na.rm=T),Sd = sd(forbbiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n)) %>%
  mutate(trt = paste(N,P, sep=("_")))
NPforbbio$trt = factor(NPforbbio$trt, levels=c("C_C","Nitrogen_C", "C_Phosphorus", "Nitrogen_Phosphorus"))

ggplot(NPforbbio, aes(x=trt, y=mean, fill=trt)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  ylab(bquote("Forb Biomass"~(g/m^2))) + xlab("Treatment") + ylim(0,175) +
  scale_x_discrete(labels=c("C", "N", "P", "N x P")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00","#999999","indianred3")) +
  annotate("text", x= 1, y = 70, label= "a", size = 9, color="black")+ 
  annotate("text", x= 2, y = 90, label= "a", size = 9, color="black")+ 
  annotate("text", x= 3, y = 75, label= "a", size = 9, color="black")+ 
  annotate("text", x= 4, y = 140, label= "b", size = 9, color="black")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(size = 15),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 30),axis.text.x=element_text(size = 20), 
        legend.position="none",axis.text.y=element_text(size = 25))
# 745 x 745 #


#######################################################################################################


############### RACs #############
SpCompfunc <- subset(SpComp, year!="2019") %>%
  group_by(year, block, plot, lifeform, taxa) %>%
  summarise(maxtaxacover = max(cover))
SpCompfunc2 <- SpCompfunc %>%
  group_by(year, block, plot, lifeform)%>%
  summarise(maxfunccover = sum(maxtaxacover))%>%
  group_by(year, block, plot)%>%
  spread(lifeform,maxfunccover)%>%
  mutate(Woody=sum(Shrub,Tree,Vine, na.rm=T))%>%
  gather(lifeform, maxfunccover,4:10)

MaxSpCompfunc <- merge(SpCompfunc2,treatments, by=c("block","plot"),all=TRUE)%>%
  filter(drought!="Drought") 

SpCompfunc1 <- SpCompfunc %>%
  group_by(year, block, plot) %>%
  mutate(totalcover = sum(maxtaxacover)) %>%
  group_by(year, block, plot, lifeform, taxa) %>%
  summarise(relativecover=(maxtaxacover/totalcover)*100)

SpCompfunc2 <- SpCompfunc1 %>%
  group_by(year, block, plot, lifeform)%>%
  summarise(maxfunccover = sum(relativecover))%>%
  group_by(year, block, plot)%>%
  spread(lifeform,maxfunccover)%>%
  mutate(Woody=sum(Shrub,Tree,Vine, na.rm=T))%>%
  gather(lifeform, maxfunccover,4:10)

MaxSpCompfunc <- merge(SpCompfunc2,treatments, by=c("block","plot"),all=TRUE)%>%
  filter(drought!="Drought") 

## For graph with summed woody cover ####
rac <- subset(MaxSpCompfunc, lifeform!="Shrub" & lifeform!="Tree" & lifeform!="Vine" & 
                Nutrients!="C_C_Potassium" & Nutrients!="C_Phosphorus_Potassium" 
              & Nutrients!="Nitrogen_C_Potassium" & Nutrients!="Nitrogen_Phosphorus_Potassium") %>%
  group_by(year, Nutrients, lifeform)%>%
  summarise(Abundance=mean(maxfunccover, na.rm=T)) 
rac2 <- rac %>%
  mutate(Rank=rank(-Abundance, ties.method="average")) %>%
  ungroup()
rac2$year <- as.factor(rac2$year)
rac2$Nutrients = factor(rac2$Nutrients, levels=c("C_C_C","Nitrogen_C_C", "C_Phosphorus_C", "Nitrogen_Phosphorus_C"))
rac2$Nutrients <- as.factor(rac2$Nutrients)


## Fig 3 ##
ggplot(data=rac2, aes(x=Rank, y=Abundance))+
  geom_line(size=1.5)+
  geom_point(aes(color=lifeform, shape =lifeform),size=12)+
  scale_shape_manual(name="Functional Group",breaks = c("Herb", "Graminoid","Legume", "Woody"),
                     labels = c("Forb", "Graminoid","Legume", "Woody"),
                     values = c(15, 16, 17, 18)) + 
  scale_color_manual(name="Functional Group",breaks = c("Herb", "Graminoid","Legume", "Woody"),
                       labels = c("Forb", "Graminoid","Legume", "Woody"),
                       values=c("#497381", "#548F01", "#CFA3EE","#D58A60"))+
  theme(panel.spacing = unit(2, "lines"), 
        plot.title = element_text(size = 30),axis.line = element_line(colour = "black"),
        text = element_text(size = 40),axis.text.x=element_text(size = 30),axis.text.y=element_text(size = 30))+
  ylim(0,150)+
  facet_grid(rows=vars(Nutrients), cols=vars(year), scales="fixed", 
    labeller = as_labeller(c(C_C_C="Control",Nitrogen_C_C="N",C_Phosphorus_C="P",
    Nitrogen_Phosphorus_C="N x P","2020"=" Year 1","2021"=" Year 2","2022"=" Year 3","2023"=" Year 4")))+
    ylab("Relative Abundance")+
    theme(strip.background=element_rect(fill="grey91", color="darkgrey"), strip.text = element_text(size = 30))
# 2400 x 1600 #

#######################################################################################################

####################################################################
## GRAPHS FOR MANUSCRIPT SUPPLEMENTAL DOC ####
##################################################################
ferns <- lmerTest::lmer(data=NoDroughtMax, Maxferncov ~ year * N * P * K + (1|block) + (1|plot))
anova(ferns)
emmeans(ferns, pairwise ~ N|year, adjust="BH")
####
# Drought x Nutrients (NPK together) #
####
ferns2 <- lmerTest::lmer(data=DroughtMax, Maxferncov ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(ferns2)




##############################################################################
##################################################################
S2 <- lmerTest::lmer(data=NoDroughtMax, Maxwirecov ~ year * N * P * K + (1|block) + (1|plot))
anova(S2)
emmeans(S2, pairwise ~ N|year, adjust="BH")
emmeans(S3, pairwise ~ N*K|year, adjust="BH")
####
# Drought x Nutrients (NPK together) #
####
a <- lmerTest::lmer(data=DroughtMax, Maxwirecov ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(a)


S3 <- lmerTest::lmer(data=NoDroughtMax, MaxEvar ~ year * N * P * K + (1|block) + (1|plot))
anova(S3)
emmeans(S3, pairwise ~ P*K|year, adjust="BH")
emmeans(S3, pairwise ~ N*K|year, adjust="BH")
####
# Drought x Nutrients (NPK together) #
####
b <- lmerTest::lmer(data=DroughtMax, MaxEvar ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(b)


S4 <- lmerTest::lmer(data=NoDroughtMax, woodybiomass ~ year * N * P * K + (1|block) + (1|plot))
anova(S4)
####
# Drought x Nutrients (NPK together) #
####
c <- lmerTest::lmer(data=DroughtMax, woodybiomass ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(c)


S5 <- lmerTest::lmer(data=NoDroughtMax, legumebiomass ~ year * N * P * K + (1|block) + (1|plot))
anova(S5)
####
# Drought x Nutrients (NPK together) #
####
a <- lmerTest::lmer(data=DroughtMax, legumebiomass ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(a)
emmeans(a, pairwise ~ drought|year, adjust="BH")



S6 <- lmerTest::lmer(data=NoDroughtMax, Maxferncov ~ year * N * P * K + (1|block) + (1|plot))
anova(S6)
emmeans(S6, pairwise ~ N|year, adjust="BH")
emmeans(S6, pairwise ~ N*P|year, adjust="BH")
####
# Drought x Nutrients (NPK together) #
####
a <- lmerTest::lmer(data=DroughtMax, legumebiomass ~ year * Nutrients * drought + (1|block) + (1|plot))
anova(a)
emmeans(a, pairwise ~ drought|year, adjust="BH")




####################################################

####
## NITROGEN EFFECTS ####
####

# Fig S7a #
NEvar <- NoDroughtMax %>%
  group_by(N) %>%
  reframe(mean = mean(MaxEvar),Sd = sd(MaxEvar),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(NEvar, aes(x=N, y=mean, fill=N)) +
  geom_bar(stat='identity', width=0.75) + ylim(0.0,0.5)+
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab("Evenness (Evar)") +
  scale_x_discrete(labels=c("C", "N")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #  


# Fig S7b #
Nwirecov <- NoDroughtMax %>%
  group_by(year,N) %>%
  reframe(mean = mean(Maxwirecov),Sd = sd(Maxwirecov),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Nwirecov, aes(x=year, y=mean, group=N,color=N, linetype=N)) +
  geom_line(size=3) + geom_point(aes(color=N),size=10) + 
  geom_errorbar(aes(color=N,ymin=mean-Se,ymax=mean+Se),size=1.5, width=0.5)+
  annotate("text", x= 3, y = 50, label= "**", size = 20)+ 
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  scale_x_discrete(labels=c("1", "2", "3", "4")) +
  ylab("Wiregrass Cover (%)") + xlab("Year") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 40),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none",axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #


# Fig S7c #
Nlegume <- NoDroughtMax %>%
  group_by(N) %>%
  reframe(mean = mean(legumebiomass, na.rm=T),Sd = sd(legumebiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Nlegume, aes(x=N, y=mean, fill=N)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab(bquote("Legume Biomass"~(g/m^2)))+ scale_x_discrete(labels=c("C", "N")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +ylim(0,35)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 # 


# Fig S7d #
NWoody <- NoDroughtMax %>%
  group_by(N) %>%
  reframe(mean = mean(woodybiomass, na.rm=T),Sd = sd(woodybiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(NWoody, aes(x=N, y=mean, fill=N)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab(bquote("Woody Biomass"~(g/m^2)))+ scale_x_discrete(labels=c("C", "N")) +
  scale_fill_manual(values=c("#56B4E9","#E69F00")) +
  annotate("text", x= 1.5, y = 40, label= "**", size = 20, color="black")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #

# Fig ??? #
Nferncov <- NoDroughtMax %>%
  group_by(year,N) %>%
  reframe(mean = mean(Maxferncov, na.rm=T),Sd = sd(Maxferncov, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Nferncov, aes(x=year, y=mean, group=N,color=N, linetype=N)) +
  geom_line(size=3) + geom_point(aes(color=N),size=10) + 
  geom_errorbar(aes(color=N,ymin=mean-Se,ymax=mean+Se),size=1.5, width=0.5)+
  annotate("text", x= 2, y = 41, label= "**", size = 20)+ 
  annotate("text", x= 4, y = 49, label= "*", size = 20)+ 
  scale_linetype_manual(values=c("solid", "longdash"))+
  scale_color_manual(values=c("#56B4E9","#E69F00")) +
  scale_x_discrete(labels=c("1", "2", "3", "4")) +
  ylab("Bracken Fern Cover (%)") + xlab("Year") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 40),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="none",axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #



####
## PHOSPHORUS EFFECTS ####
####

# Fig S7e #
PEvar <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(MaxEvar),Sd = sd(MaxEvar),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(PEvar, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab("Evenness (Evar)") +
  scale_x_discrete(labels=c("C", "P")) + ylim(0.0,0.5)+
  scale_fill_manual(values=c("#56B4E9","#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #


# Fig S7f #
Pwirecov <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(Maxwirecov),Sd = sd(Maxwirecov),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Pwirecov, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab("Wiregrass Cover (%)") +
  scale_x_discrete(labels=c("C", "P")) +
  scale_fill_manual(values=c("#56B4E9","#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #  


# Fig S7g #
Plegume <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(legumebiomass, na.rm=T),Sd = sd(legumebiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Plegume, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab(bquote("Legume Biomass"~(g/m^2)))+ scale_x_discrete(labels=c("C", "P")) +
  scale_fill_manual(values=c("#56B4E9","#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #


# Fig S7h #
PWoody <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(woodybiomass, na.rm=T),Sd = sd(woodybiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(PWoody, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab(bquote("Woody Biomass"~(g/m^2)))+ scale_x_discrete(labels=c("C", "P")) +
  scale_fill_manual(values=c("#56B4E9","#999999")) +ylim(0,40)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 # 


PFern <- NoDroughtMax %>%
  group_by(P) %>%
  reframe(mean = mean(Maxferncov, na.rm=T),Sd = sd(Maxferncov, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(PFern, aes(x=P, y=mean, fill=P)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab("Bracken Fern Cover (%)") +
  scale_x_discrete(labels=c("C", "P")) + ylim(0,35)+
  scale_fill_manual(values=c("#56B4E9","#999999")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #




####
## POTASSIUM EFFECTS ####
####

# Fig S7i #
KEvar <- NoDroughtMax %>%
  group_by(K) %>%
  reframe(mean = mean(MaxEvar),Sd = sd(MaxEvar),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(KEvar, aes(x=K, y=mean, fill=K)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab("Evenness (Evar)") +
  scale_x_discrete(labels=c("C", "K")) +
  annotate("text", x= 1.5, y = 0.5, label= "**", size = 20, color="black")+ 
  scale_fill_manual(values=c("#56B4E9","palevioletred4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #


# Fig S7j #
Kwirecov <- NoDroughtMax %>%
  group_by(K) %>%
  reframe(mean = mean(Maxwirecov),Sd = sd(Maxwirecov),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Kwirecov, aes(x=K, y=mean, fill=K)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab("Wiregrass Cover (%)") +
  scale_x_discrete(labels=c("C", "K")) +
  scale_fill_manual(values=c("#56B4E9","palevioletred4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #


# Fig S7k #
Klegume <- NoDroughtMax %>%
  group_by(K) %>%
  reframe(mean = mean(legumebiomass, na.rm=T),Sd = sd(legumebiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(Klegume, aes(x=K, y=mean, fill=K)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab(bquote("Legume Biomass"~(g/m^2)))+ scale_x_discrete(labels=c("C", "K")) +
  scale_fill_manual(values=c("#56B4E9","palevioletred4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 # 


# Fig S7l #
KWoody <- NoDroughtMax %>%
  group_by(K) %>%
  reframe(mean = mean(woodybiomass, na.rm=T),Sd = sd(woodybiomass, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(KWoody, aes(x=K, y=mean, fill=K)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab(bquote("Woody Biomass"~(g/m^2)))+ scale_x_discrete(labels=c("C", "K")) +
  scale_fill_manual(values=c("#56B4E9","palevioletred4")) + ylim(0,40)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 # 

# Fig S7i #
KFern <- NoDroughtMax %>%
  group_by(K) %>%
  reframe(mean = mean(Maxferncov, na.rm=T),Sd = sd(Maxferncov, na.rm=T),n=n()) %>%
  mutate(Se = Sd/sqrt(n))
ggplot(KFern, aes(x=K, y=mean, fill=K)) +
  geom_bar(stat='identity', width=0.75) +
  geom_errorbar(aes(ymin=mean-Se, ymax=mean+Se), width=0.3, size=1.5, position=position_dodge(0.9)) +
  xlab("Treatment") +
  ylab("Bracken Fern Cover (%)") + ylim(0,35)+
  scale_x_discrete(labels=c("C", "K")) +
  annotate("text", x= 1.5, y = 33, label= "*", size = 20, color="black")+ 
  scale_fill_manual(values=c("#56B4E9","palevioletred4")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 40),legend.position="none",
        axis.text.x=element_text(size = 40), axis.text.y=element_text(size = 40))
# 745 x 745 #