#Library packages
pkgs=c("tidyverse","ggplot2", "phyloseq","grid","scales","microbiomeSeq", "devtools","igraph",
       "vegan","cowplot")

for(pkg in pkgs){ 
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

#read bacteria RDS file
physeq_rs= readRDS("physeq-rs-bacteria.rds")
physeq_rs
#read fungi RDS file
physeq_rs_fungi= readRDS("physeq-rs-fungi.rds")
physeq_rs_fungi

#rarefy
set.seed(1)
min(sample_sums(physeq_rs))
physeq_rs_rare = rarefy_even_depth(physeq_rs, sample.size = min(sample_sums(physeq_rs)),
                                   rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
physeq_rs_rare

saveRDS(physeq_rs_rare, "physeq_rs_rare.rds")
physeq_rs_rare= readRDS("physeq_rs_rare.rds")

min(sample_sums(physeq_rs_fungi))
physeq_rs_fungi_rare = rarefy_even_depth(physeq_rs_fungi, sample.size = min(sample_sums(physeq_rs_fungi)),
                                         rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
physeq_rs_fungi_rare

saveRDS(physeq_rs_fungi_rare, "physeq_rs_fungi_rare.rds")
physeq_rs_fungi_rare= readRDS("physeq_rs_fungi_rare.rds")

theme_set(theme_light())
colors <- c(
  "green3","#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "red3",
  "#AD6F3B", "#673770","#D14285", "#8569D5","#652926"
)
#order factors
sample_data(physeq_rs_rare)$Time=factor(sample_data(physeq_rs_rare)$Time,levels = c('Mar','Apr','May'))
sample_data(physeq_rs_rare)$SoilType=factor(sample_data(physeq_rs_rare)$SoilType,levels = c('XL','KL','DM'))
sample_data(physeq_rs_rare)$Treatment=factor(sample_data(physeq_rs_rare)$Treatment)
sample_data(physeq_rs_fungi_rare)$Time=factor(sample_data(physeq_rs_fungi_rare)$Time,levels = c('Mar','Apr','May'))
sample_data(physeq_rs_fungi_rare)$SoilType=factor(sample_data(physeq_rs_fungi_rare)$SoilType,levels = c('XL','KL','DM'))
sample_data(physeq_rs_fungi_rare)$Treatment=factor(sample_data(physeq_rs_fungi_rare)$Treatment)

######################################################################################################################################
#1.whole database facter= Time
######################################################################################################################################
#alpha diversity
p1 <- plot_anova_diversity(physeq_rs_rare, method = c("shannon"), 
                           grouping_column="Time", pValueCutoff = 0.05)
p1

#beta diversity
ord_PCoA_wunifrac = ordinate(physeq_rs_rare, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
ord_PCoA_bray = ordinate(physeq_rs_rare, "PCoA", "bray")  # perform PCoA on bray-curtis distance
#PCoA_bray
p_PCoA_bary <- plot_ordination(physeq_rs_rare, ord_PCoA_bray, color = "Time"
) + 
  ggtitle("PCoA on Bray-Curtis distance")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = Time), alpha = 0.7, size = 1)+theme_light()  
p_PCoA_bary
#PCoA_wunifrac
p_PCoA_wunifrac <- plot_ordination(physeq_rs_rare, ord_PCoA_wunifrac, color = "Time"
) + 
  ggtitle("PCoA on Weighted-UniFrac distance")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = Time), alpha = 0.7, size = 1)+theme_light()  
p_PCoA_wunifrac
#Adonis test on whole dataset using compartment as factor (Part of Table 2) 
metadata<- data.frame(sample_data(physeq_rs_rare)[,])
metadata
#time
bray.dist = phyloseq::distance(physeq_rs_rare, method="bray")
adonis(bray.dist ~ Time, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare, method="wunifrac")
adonis(weihted.dist ~ Time, data=metadata, permutations = 1000)
#soiltype
bray.dist = phyloseq::distance(physeq_rs_rare, method="bray")
adonis(bray.dist ~ SoilType, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare, method="wunifrac")
adonis(weihted.dist ~ SoilType, data=metadata, permutations = 1000)
#treatment
bray.dist = phyloseq::distance(physeq_rs_rare, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

#barplot
physeq_rs_rare_top100 = prune_taxa(names(sort(taxa_sums(physeq_rs_rare), TRUE)[1:100]), physeq_rs_rare)
physeq_rs_rare_top100_class = tax_glom(physeq_rs_rare_top100, taxrank="Class", NArm=FALSE)%>% prune_taxa(taxa_sums(.) > 0, .)
physeq_rs_rare_top100_class_mean = merge_samples(physeq_rs_rare_top100_class,"Mean",fun=mean)
physeq_rs_rare_top100_class_mean_ra = transform_sample_counts(physeq_rs_rare_top100_class_mean, function(x){x / sum(x)})
p2=plot_bar(physeq_rs_rare_top100_class_mean_ra, fill="Class") + facet_wrap(~Time, scales= "free_x", nrow=1)+geom_bar(stat="identity", position="stack")
p2 = p2 + scale_fill_manual(values = colors)+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1)) 
p2

#alpha diversity
p1_fungi <- plot_anova_diversity(physeq_rs_fungi_rare, method = c("shannon"), # "simpson",
                                 grouping_column="Time", pValueCutoff = 0.05)
p1_fungi

#beta diversity
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare, "PCoA", "bray")  # perform PCoA on bray-curtis distance
#PCoA_bray
p_PCoA_bary_fungi <- plot_ordination(physeq_rs_fungi_rare, ord_PCoA_bray_fungi, color = "Time"
) + 
  ggtitle("Weighted - UniFrac")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = Time), alpha = 0.7, size = 1) +theme_light()  
p_PCoA_bary_fungi
#PCoA_wunifrac
p_PCoA_wunifrac_fungi <- plot_ordination(physeq_rs_fungi_rare, ord_PCoA_wunifrac_fungi, color = "Time"
) + 
  ggtitle("PCoA on Weighted-UniFrac distance")+
  scale_color_manual(values = c("black", "red","blue" )) +   
  geom_point(aes(color = Time), alpha = 0.7, size = 1) +theme_light()  
p_PCoA_wunifrac_fungi
#Adonis test on whole dataset using compartment as factor (Part of Table 2) 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare)[,])
metadata
#time
bray.dist = phyloseq::distance(physeq_rs_fungi_rare, method="bray")
adonis(bray.dist ~ Time, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare, method="wunifrac")
adonis(weihted.dist ~ Time, data=metadata, permutations = 1000)
#soiltype
bray.dist = phyloseq::distance(physeq_rs_fungi_rare, method="bray")
adonis(bray.dist ~ SoilType, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare, method="wunifrac")
adonis(weihted.dist ~ SoilType, data=metadata, permutations = 1000)
#treatment
bray.dist = phyloseq::distance(physeq_rs_fungi_rare, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

#barplot
physeq_rs_fungi_rare_top100 = prune_taxa(names(sort(taxa_sums(physeq_rs_fungi_rare), TRUE)[1:100]), physeq_rs_fungi_rare)
physeq_rs_fungi_rare_top100_class = tax_glom(physeq_rs_fungi_rare_top100, taxrank="Class", NArm=FALSE)%>% prune_taxa(taxa_sums(.) > 0, .)
physeq_rs_fungi_rare_top100_class_mean = merge_samples(physeq_rs_fungi_rare_top100_class,"Mean",fun=mean)
physeq_rs_fungi_rare_top100_class_mean_ra = transform_sample_counts(physeq_rs_fungi_rare_top100_class_mean, function(x){x / sum(x)})
p2_fungi=plot_bar(physeq_rs_fungi_rare_top100_class_mean_ra, fill="Class") + facet_wrap(~Time, scales= "free_x", nrow=1)+geom_bar(stat="identity", position="stack")
p2_fungi = p2_fungi + scale_fill_manual(values = colors) +theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1)) 
p2_fungi

library("gridExtra")
library("cowplot")
library("ggpubr")
gt= grid.arrange(arrangeGrob(p1, p_PCoA_bary, p_PCoA_wunifrac, ncol = 3), 
                 arrangeGrob(p1_fungi,p_PCoA_bary_fungi,p_PCoA_wunifrac_fungi, ncol = 3),
                 arrangeGrob(p2,p2_fungi, ncol = 2),
                 nrow = 3)
# changeto ggplot # add lables
p_all <- as_ggplot(gt) +                               
  draw_plot_label(label = c("a", "b", "c", "d","e","f","g","h"), size = 15,
                  x = c(0,0.35,0.65,0,0.35,0.65,0,0.5), y = c(1,1,1,0.65,0.65,0.65,0.35,0.35)) 
p_all

ggsave("Fig. 2.pdf", p_all,width = 30, height = 30, units = "cm")

######################################################################################################################################
#2.Mar   facter= SoilType
######################################################################################################################################

physeq_rs_rare_Mar = subset_samples(physeq_rs_rare, Time =="Mar") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_Mar, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_Mar, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary1 <- plot_ordination(physeq_rs_rare_Mar, ord_PCoA_bray, color = "SoilType") + 
  ggtitle("Mar_soiltype")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_bary1
#PCoA_wunifrac
p_PCoA_wunifrac1 <- plot_ordination(physeq_rs_rare_Mar, ord_PCoA_wunifrac, color = "SoilType") +
  ggtitle("Mar_soiltype")+
  scale_color_manual(values = c("black", "red","blue" )) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)  
p_PCoA_wunifrac1
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_Mar)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_Mar, method="bray")
adonis(bray.dist ~ SoilType, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_Mar, method="wunifrac")
adonis(weihted.dist ~ SoilType, data=metadata, permutations = 1000)

####################################
#2.1Mar-DM   facter= Treatment
####################################
physeq_rs_rare_Mar_DM = subset_samples(physeq_rs_rare, Time =="Mar") %>% subset_samples(., SoilType =="DM") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_Mar_DM, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_Mar_DM, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary2 <- plot_ordination(physeq_rs_rare_Mar_DM, ord_PCoA_bray,color = "Treatment") + 
  ggtitle("Mar_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_bary2
#PCoA_wunifrac
p_PCoA_wunifrac2 <- plot_ordination(physeq_rs_rare_Mar_DM, ord_PCoA_wunifrac, color = "Treatment") + 
  ggtitle("Mar_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac2
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_Mar_DM)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_Mar_DM, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_Mar_DM, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#2.2Mar-KL   facter= Treatment
####################################
physeq_rs_rare_Mar_KL = subset_samples(physeq_rs_rare, Time =="Mar") %>% subset_samples(., SoilType =="KL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_Mar_KL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_Mar_KL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary3 <- plot_ordination(physeq_rs_rare_Mar_KL, ord_PCoA_bray, color = "Treatment") + 
  ggtitle("Mar_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_bary3
#PCoA_wunifrac
p_PCoA_wunifrac3 <- plot_ordination(physeq_rs_rare_Mar_KL, ord_PCoA_wunifrac, color = "Treatment") + 
  ggtitle("Mar_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac3
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_Mar_KL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_Mar_KL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_Mar_KL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#2.3Mar-XL   facter= Treatment
####################################
physeq_rs_rare_Mar_XL = subset_samples(physeq_rs_rare, Time =="Mar") %>% subset_samples(., SoilType =="XL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_Mar_XL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_Mar_XL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary4 <- plot_ordination(physeq_rs_rare_Mar_XL, ord_PCoA_bray, color = "Treatment") + 
  ggtitle("Mar_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_bary4
#PCoA_wunifrac
p_PCoA_wunifrac4 <- plot_ordination(physeq_rs_rare_Mar_XL, ord_PCoA_wunifrac, color = "Treatment") + 
  ggtitle("Mar_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac4
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_Mar_XL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_Mar_XL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_Mar_XL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

######################################################################################################################################
#3.Apr   facter= SoilType
######################################################################################################################################

physeq_rs_rare_Apr = subset_samples(physeq_rs_rare, Time =="Apr") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_Apr, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_Apr, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary5 <- plot_ordination(physeq_rs_rare_Apr, ord_PCoA_bray, color = "SoilType") + 
  ggtitle("Apr_soiltype")+
  scale_color_manual(values = c("black", "red","blue" )) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_bary5
#PCoA_wunifrac
p_PCoA_wunifrac5 <- plot_ordination(physeq_rs_rare_Apr, ord_PCoA_wunifrac, color = "SoilType") +
  ggtitle("Apr_soiltype")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac5
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_Apr)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_Apr, method="bray")
adonis(bray.dist ~ SoilType, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_Apr, method="wunifrac")
adonis(weihted.dist ~ SoilType, data=metadata, permutations = 1000)

####################################
#3.1Apr-DM   facter= Treatment
####################################
physeq_rs_rare_Apr_DM = subset_samples(physeq_rs_rare, Time =="Apr") %>% subset_samples(., SoilType =="DM") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_Apr_DM, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_Apr_DM, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary6 <- plot_ordination(physeq_rs_rare_Apr_DM, ord_PCoA_bray, color = "Treatment") + 
  ggtitle("Apr_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_bary6
#PCoA_wunifrac
p_PCoA_wunifrac6 <- plot_ordination(physeq_rs_rare_Apr_DM, ord_PCoA_wunifrac, color = "Treatment") + 
  ggtitle("Apr_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac6
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_Apr_DM)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_Apr_DM, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_Apr_DM, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#3.2Apr-KL   facter= Treatment
####################################
physeq_rs_rare_Apr_KL = subset_samples(physeq_rs_rare, Time =="Apr") %>% subset_samples(., SoilType =="KL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_Apr_KL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_Apr_KL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary7 <- plot_ordination(physeq_rs_rare_Apr_KL, ord_PCoA_bray, color = "Treatment") + 
  ggtitle("Apr_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_bary7
#PCoA_wunifrac
p_PCoA_wunifrac7 <- plot_ordination(physeq_rs_rare_Apr_KL, ord_PCoA_wunifrac, color = "Treatment") + 
  ggtitle("Apr_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac7
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_Apr_KL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_Apr_KL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_Apr_KL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#3.3Apr-XL   facter= Treatment
####################################
physeq_rs_rare_Apr_XL = subset_samples(physeq_rs_rare, Time =="Apr") %>% subset_samples(., SoilType =="XL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_Apr_XL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_Apr_XL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary8 <- plot_ordination(physeq_rs_rare_Apr_XL, ord_PCoA_bray, color = "Treatment") + 
  ggtitle("Apr_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_bary8
#PCoA_wunifrac
p_PCoA_wunifrac8 <- plot_ordination(physeq_rs_rare_Apr_XL, ord_PCoA_wunifrac, color = "Treatment") + 
  ggtitle("Apr_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac8
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_Apr_XL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_Apr_XL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_Apr_XL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

######################################################################################################################################
#4.May   facter= SoilType
######################################################################################################################################

physeq_rs_rare_May = subset_samples(physeq_rs_rare, Time =="May") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_May, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_May, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary9 <- plot_ordination(physeq_rs_rare_May, ord_PCoA_bray, color = "SoilType") + 
  ggtitle("May_soiltype")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_bary9
#PCoA_wunifrac
p_PCoA_wunifrac9 <- plot_ordination(physeq_rs_rare_May, ord_PCoA_wunifrac, color = "SoilType") +
  ggtitle("May_soiltype")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac9
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_May)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_May, method="bray")
adonis(bray.dist ~ SoilType, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_May, method="wunifrac")
adonis(weihted.dist ~ SoilType, data=metadata, permutations = 1000)

####################################
#4.1May-DM   facter= Treatment
####################################
physeq_rs_rare_May_DM = subset_samples(physeq_rs_rare, Time =="May") %>% subset_samples(., SoilType =="DM") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_May_DM, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_May_DM, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary10 <- plot_ordination(physeq_rs_rare_May_DM, ord_PCoA_bray, color = "Treatment") + 
  ggtitle("May_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_bary10
#PCoA_wunifrac
p_PCoA_wunifrac10 <- plot_ordination(physeq_rs_rare_May_DM, ord_PCoA_wunifrac, color = "Treatment") + 
  ggtitle("May_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac10
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_May_DM)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_May_DM, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_May_DM, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#4.2May-KL   facter= Treatment
####################################
physeq_rs_rare_May_KL = subset_samples(physeq_rs_rare, Time =="May") %>% subset_samples(., SoilType =="KL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_May_KL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_May_KL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary11 <- plot_ordination(physeq_rs_rare_May_KL, ord_PCoA_bray, color = "Treatment") + 
  ggtitle("May_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_bary11
#PCoA_wunifrac
p_PCoA_wunifrac11 <- plot_ordination(physeq_rs_rare_May_KL, ord_PCoA_wunifrac, color = "Treatment") + 
  ggtitle("May_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac11
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_May_KL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_May_KL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_May_KL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#4.3May-XL   facter= Treatment
####################################
physeq_rs_rare_May_XL = subset_samples(physeq_rs_rare, Time =="May") %>% subset_samples(., SoilType =="XL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray = ordinate(physeq_rs_rare_May_XL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac = ordinate(physeq_rs_rare_May_XL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary12 <- plot_ordination(physeq_rs_rare_May_XL, ord_PCoA_bray, color = "Treatment") + 
  ggtitle("May_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_bary12
#PCoA_wunifrac
p_PCoA_wunifrac12 <- plot_ordination(physeq_rs_rare_May_XL, ord_PCoA_wunifrac, color = "Treatment") + 
  ggtitle("May_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac12
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_rare_May_XL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_rare_May_XL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_rare_May_XL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)


#fungi
######################################################################################################################################
#2.Mar   facter= SoilType
######################################################################################################################################

physeq_rs_fungi_rare_Mar = subset_samples(physeq_rs_fungi_rare, Time =="Mar") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_Mar, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_Mar, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi1 <- plot_ordination(physeq_rs_fungi_rare_Mar, ord_PCoA_bray_fungi, color = "SoilType") + 
  ggtitle("Mar_soiltype")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi1
#PCoA_wunifrac
p_PCoA_wunifrac_fungi1 <- plot_ordination(physeq_rs_fungi_rare_Mar, ord_PCoA_wunifrac_fungi, color = "SoilType") +
  ggtitle("Mar_soiltype")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi1
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_Mar)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_Mar, method="bray")
adonis(bray.dist ~ SoilType, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_Mar, method="wunifrac")
adonis(weihted.dist ~ SoilType, data=metadata, permutations = 1000)

####################################
#2.1Mar-DM   facter= Treatment
####################################
physeq_rs_fungi_rare_Mar_DM = subset_samples(physeq_rs_fungi_rare, Time =="Mar") %>% subset_samples(., SoilType =="DM") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_Mar_DM, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_Mar_DM, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi2 <- plot_ordination(physeq_rs_fungi_rare_Mar_DM, ord_PCoA_bray_fungi,color = "Treatment") + 
  ggtitle("Mar_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi2
#PCoA_wunifrac
p_PCoA_wunifrac_fungi2 <- plot_ordination(physeq_rs_fungi_rare_Mar_DM, ord_PCoA_wunifrac_fungi, color = "Treatment") + 
  ggtitle("Mar_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi2
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_Mar_DM)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_Mar_DM, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_Mar_DM, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#2.2Mar-KL   facter= Treatment
####################################
physeq_rs_fungi_rare_Mar_KL = subset_samples(physeq_rs_fungi_rare, Time =="Mar") %>% subset_samples(., SoilType =="KL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_Mar_KL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_Mar_KL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi3 <- plot_ordination(physeq_rs_fungi_rare_Mar_KL, ord_PCoA_bray_fungi, color = "Treatment") + 
  ggtitle("Mar_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi3
#PCoA_wunifrac
p_PCoA_wunifrac_fungi3 <- plot_ordination(physeq_rs_fungi_rare_Mar_KL, ord_PCoA_wunifrac_fungi, color = "Treatment") + 
  ggtitle("Mar_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi3
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_Mar_KL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_Mar_KL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_Mar_KL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#2.3Mar-XL   facter= Treatment
####################################
physeq_rs_fungi_rare_Mar_XL = subset_samples(physeq_rs_fungi_rare, Time =="Mar") %>% subset_samples(., SoilType =="XL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_Mar_XL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_Mar_XL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi4 <- plot_ordination(physeq_rs_fungi_rare_Mar_XL, ord_PCoA_bray_fungi, color = "Treatment") + 
  ggtitle("Mar_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi4
#PCoA_wunifrac
p_PCoA_wunifrac_fungi4 <- plot_ordination(physeq_rs_fungi_rare_Mar_XL, ord_PCoA_wunifrac_fungi, color = "Treatment") + 
  ggtitle("Mar_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi4
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_Mar_XL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_Mar_XL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_Mar_XL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

######################################################################################################################################
#3.Apr   facter= SoilType
######################################################################################################################################

physeq_rs_fungi_rare_Apr = subset_samples(physeq_rs_fungi_rare, Time =="Apr") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_Apr, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_Apr, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi5 <- plot_ordination(physeq_rs_fungi_rare_Apr, ord_PCoA_bray_fungi, color = "SoilType") + 
  ggtitle("Apr_soiltype")+
  scale_color_manual(values = c("black", "red","blue" )) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi5
#PCoA_wunifrac
p_PCoA_wunifrac_fungi5 <- plot_ordination(physeq_rs_fungi_rare_Apr, ord_PCoA_wunifrac_fungi, color = "SoilType") +
  ggtitle("Apr_soiltype")+
  scale_color_manual(values = c("black", "red","blue")) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi5
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_Apr)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_Apr, method="bray")
adonis(bray.dist ~ SoilType, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_Apr, method="wunifrac")
adonis(weihted.dist ~ SoilType, data=metadata, permutations = 1000)

####################################
#3.1Apr-DM   facter= Treatment
####################################
physeq_rs_fungi_rare_Apr_DM = subset_samples(physeq_rs_fungi_rare, Time =="Apr") %>% subset_samples(., SoilType =="DM") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_Apr_DM, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_Apr_DM, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi6 <- plot_ordination(physeq_rs_fungi_rare_Apr_DM, ord_PCoA_bray_fungi, color = "Treatment") + 
  ggtitle("Apr_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi6
#PCoA_wunifrac
p_PCoA_wunifrac_fungi6 <- plot_ordination(physeq_rs_fungi_rare_Apr_DM, ord_PCoA_wunifrac_fungi, color = "Treatment") + 
  
  ggtitle("Apr_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi6
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_Apr_DM)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_Apr_DM, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_Apr_DM, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#3.2Apr-KL   facter= Treatment
####################################
physeq_rs_fungi_rare_Apr_KL = subset_samples(physeq_rs_fungi_rare, Time =="Apr") %>% subset_samples(., SoilType =="KL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_Apr_KL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_Apr_KL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi7 <- plot_ordination(physeq_rs_fungi_rare_Apr_KL, ord_PCoA_bray_fungi, color = "Treatment") + 
  
  ggtitle("Apr_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi7
#PCoA_wunifrac
p_PCoA_wunifrac_fungi7 <- plot_ordination(physeq_rs_fungi_rare_Apr_KL, ord_PCoA_wunifrac_fungi, color = "Treatment") + 
  
  ggtitle("Apr_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi7
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_Apr_KL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_Apr_KL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_Apr_KL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#3.3Apr-XL   facter= Treatment
####################################
physeq_rs_fungi_rare_Apr_XL = subset_samples(physeq_rs_fungi_rare, Time =="Apr") %>% subset_samples(., SoilType =="XL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_Apr_XL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_Apr_XL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi8 <- plot_ordination(physeq_rs_fungi_rare_Apr_XL, ord_PCoA_bray_fungi, color = "Treatment") + 
  
  ggtitle("Apr_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi8
#PCoA_wunifrac
p_PCoA_wunifrac_fungi8 <- plot_ordination(physeq_rs_fungi_rare_Apr_XL, ord_PCoA_wunifrac_fungi, color = "Treatment") + 
  
  ggtitle("Apr_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi8
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_Apr_XL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_Apr_XL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_Apr_XL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

######################################################################################################################################
#4.May   facter= SoilType
######################################################################################################################################

physeq_rs_fungi_rare_May = subset_samples(physeq_rs_fungi_rare, Time =="May") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_May, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_May, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi9 <- plot_ordination(physeq_rs_fungi_rare_May, ord_PCoA_bray_fungi, color = "SoilType") + 
  ggtitle("May_soiltype")+
  scale_color_manual(values = c("black", "red","blue" )) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi9
#PCoA_wunifrac
p_PCoA_wunifrac_fungi9 <- plot_ordination(physeq_rs_fungi_rare_May, ord_PCoA_wunifrac_fungi, color = "SoilType") +
  ggtitle("May_soiltype")+
  scale_color_manual(values = c("black", "red","blue" )) +   
  geom_point(aes(color = SoilType), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi9
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_May)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_May, method="bray")
adonis(bray.dist ~ SoilType, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_May, method="wunifrac")
adonis(weihted.dist ~ SoilType, data=metadata, permutations = 1000)

####################################
#4.1May-DM   facter= Treatment
####################################
physeq_rs_fungi_rare_May_DM = subset_samples(physeq_rs_fungi_rare, Time =="May") %>% subset_samples(., SoilType =="DM") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_May_DM, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_May_DM, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi10 <- plot_ordination(physeq_rs_fungi_rare_May_DM, ord_PCoA_bray_fungi, color = "Treatment") + 
  ggtitle("May_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi10
#PCoA_wunifrac
p_PCoA_wunifrac_fungi10 <- plot_ordination(physeq_rs_fungi_rare_May_DM, ord_PCoA_wunifrac_fungi, color = "Treatment") + 
  ggtitle("May_DM_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) + theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi10
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_May_DM)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_May_DM, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_May_DM, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#4.2May-KL   facter= Treatment
####################################
physeq_rs_fungi_rare_May_KL = subset_samples(physeq_rs_fungi_rare, Time =="May") %>% subset_samples(., SoilType =="KL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_May_KL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_May_KL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi11 <- plot_ordination(physeq_rs_fungi_rare_May_KL, ord_PCoA_bray_fungi, color = "Treatment") + 
  ggtitle("May_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi11
#PCoA_wunifrac
p_PCoA_wunifrac_fungi11 <- plot_ordination(physeq_rs_fungi_rare_May_KL, ord_PCoA_wunifrac_fungi, color = "Treatment") + 
  ggtitle("May_KL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1)+ theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi11
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_May_KL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_May_KL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_May_KL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)

####################################
#4.3May-XL   facter= Treatment
####################################
physeq_rs_fungi_rare_May_XL = subset_samples(physeq_rs_fungi_rare, Time =="May") %>% subset_samples(., SoilType =="XL") %>% prune_taxa(taxa_sums(.) > 0, .)
#beta diversity
ord_PCoA_bray_fungi = ordinate(physeq_rs_fungi_rare_May_XL, "PCoA", "bray")  # perform PCoA on bray-curtis distance
ord_PCoA_wunifrac_fungi = ordinate(physeq_rs_fungi_rare_May_XL, "PCoA", "unifrac", weighted = TRUE)  # weighted-UniFrac
#PCoA_bray
p_PCoA_bary_fungi12 <- plot_ordination(physeq_rs_fungi_rare_May_XL, ord_PCoA_bray_fungi, color = "Treatment") + 
  ggtitle("May_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_bary_fungi12
#PCoA_wunifrac
p_PCoA_wunifrac_fungi12 <- plot_ordination(physeq_rs_fungi_rare_May_XL, ord_PCoA_wunifrac_fungi, color = "Treatment") + 
  ggtitle("May_XL_treatment")+
  scale_color_manual(values = c("green3", "magenta", "#ffae19" ,"yellow",
                                "#4daf4a", "#1919ff")) +   
  geom_point(aes(color = Treatment), alpha = 0.7, size = 1) +
  theme_light()+theme(aspect.ratio=1)
p_PCoA_wunifrac_fungi12
#Adonis test on whole dataset using compartment as factor 
metadata<- data.frame(sample_data(physeq_rs_fungi_rare_May_XL)[,])
metadata
bray.dist = phyloseq::distance(physeq_rs_fungi_rare_May_XL, method="bray")
adonis(bray.dist ~ Treatment, data=metadata, permutations = 1000)
weihted.dist = phyloseq::distance(physeq_rs_fungi_rare_May_XL, method="wunifrac")
adonis(weihted.dist ~ Treatment, data=metadata, permutations = 1000)


#plot
(p_all = plot_grid(p_PCoA_bary1,p_PCoA_bary2,p_PCoA_bary3,p_PCoA_bary4,p_PCoA_bary5,p_PCoA_bary6,
                   p_PCoA_bary7,p_PCoA_bary8,p_PCoA_bary9,p_PCoA_bary10,p_PCoA_bary11,p_PCoA_bary12,
                   labels = c("a", "b", "c", "d","e","f","g","h","i","j","k","l"), ncol = 4))
p_all
ggsave("Fig. 3.pdf", p_all,width = 30, height = 30, units = "cm")

(p_all = plot_grid(p_PCoA_wunifrac1,p_PCoA_wunifrac2,p_PCoA_wunifrac3,p_PCoA_wunifrac4,p_PCoA_wunifrac5,p_PCoA_wunifrac6,
                   p_PCoA_wunifrac7,p_PCoA_wunifrac8,p_PCoA_wunifrac9,p_PCoA_wunifrac10,p_PCoA_wunifrac11,p_PCoA_wunifrac12,
                   labels = c("a", "b", "c", "d","e","f","g","h","i","j","k","l"), ncol = 4))
ggsave("Supplementary Fig. 1.pdf", p_all,width = 30, height = 15, units = "cm")

(p_all = plot_grid(p_PCoA_bary_fungi1,p_PCoA_bary_fungi2,p_PCoA_bary_fungi3,p_PCoA_bary_fungi4,p_PCoA_bary_fungi5,p_PCoA_bary_fungi6,
                   p_PCoA_bary_fungi7,p_PCoA_bary_fungi8,p_PCoA_bary_fungi9,p_PCoA_bary_fungi10,p_PCoA_bary_fungi11,p_PCoA_bary_fungi12,
                   labels = c("a", "b", "c", "d","e","f","g","h","i","j","k","l"), ncol = 4))
ggsave("Fig. 4.pdf", p_all,width = 30, height = 30, units = "cm")

(p_all = plot_grid(p_PCoA_wunifrac_fungi1,p_PCoA_wunifrac_fungi2,p_PCoA_wunifrac_fungi3,p_PCoA_wunifrac_fungi4,p_PCoA_wunifrac_fungi5,p_PCoA_wunifrac_fungi6,
                   p_PCoA_wunifrac_fungi7,p_PCoA_wunifrac_fungi8,p_PCoA_wunifrac_fungi9,p_PCoA_wunifrac_fungi10,p_PCoA_wunifrac_fungi11,p_PCoA_wunifrac_fungi12,
                   labels = c("a", "b", "c", "d","e","f","g","h","i","j","k","l"), ncol = 4))
ggsave("Supplementary Fig. 2.pdf", p_all,width = 30, height = 15, units = "cm")


######################################################################################################################################
#5.DESeq2 
######################################################################################################################################
###Mar-Apr

physeq_rs_rare_top100_genus <- taxa_level(physeq_rs_rare_top100, "Genus")%>% prune_taxa(taxa_sums(.) > 0, .)
physeq_rs_rare_top100_genus_Apr_Mar = subset_samples(physeq_rs_rare_top100_genus, Time !="May") %>% prune_taxa(taxa_sums(.) > 0, .)
deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Apr_Mar, grouping_column = "Time", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)
p5 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p5)
ggsave(paste0("Fig. 5.pdf"), p5, width = 30, height = 30, units = "cm")

p6 <- plot_MDA(deseq_sig$importance)
print(p6)

physeq_rs_rare_top100_genus_Apr_Mar_pseudomonas<- subset_taxa(physeq_rs_rare_top100_genus_Apr_Mar, Genus=="Pseudomonas")
physeq_rs_rare_top100_genus
###Mar-May

physeq_rs_rare_top100_genus_Mar_May = subset_samples(physeq_rs_rare_top100_genus, Time !="Apr") %>% prune_taxa(taxa_sums(.) > 0, .)
deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Mar_May, grouping_column = "Time", output_norm = "log-relative",
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

p5 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p5)

ggsave(paste0("Mar_Mar.pdf"), p5, width = 30, height = 30, units = "cm")

p6 <- plot_MDA(deseq_sig$importance)
print(p6)

########################################################
#5.1Mar XL-DM  
########################################################

physeq_rs_rare_top100_genus_Mar_XL_DM=subset_samples(physeq_rs_rare_top100_genus, Time =="Mar") %>% subset_samples(., SoilType !="KL") %>% prune_taxa(taxa_sums(.) > 0, .)
sample_data(physeq_rs_rare_top100_genus_Mar_XL_DM)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Mar_XL_DM, grouping_column = "SoilType", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)


p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Mar_XL_DM.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#5.2Mar KL_DM  
########################################################

physeq_rs_rare_top100_genus_Mar_KL_DM=subset_samples(physeq_rs_rare_top100_genus, Time =="Mar") %>% subset_samples(., SoilType !="XL") %>% prune_taxa(taxa_sums(.) > 0, .)
sample_data(physeq_rs_rare_top100_genus_Mar_KL_DM)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Mar_KL_DM, grouping_column = "SoilType", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Mar_KL_DM.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#5.3Mar_DM CK-T1
########################################################
physeq_rs_rare_top100_genus_Mar_DM =subset_samples(physeq_rs_rare_top100_genus, Time =="Mar")%>% 
  subset_samples(., SoilType =="DM")  %>% prune_taxa(taxa_sums(.) > 0, .)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Mar_DM, grouping_column = "Treatment", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)
p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Mar_DM_treatment.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#5.4Mar_KL CK-T1
########################################################
physeq_rs_rare_top100_genus_Mar_KL =subset_samples(physeq_rs_rare_top100_genus, Time =="Mar")%>% 
  subset_samples(., SoilType =="KL")  %>% prune_taxa(taxa_sums(.) > 0, .)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Mar_KL, grouping_column = "Treatment", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)


p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Mar_KL_treatment.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#5.5Mar_XL CK-T1
########################################################
physeq_rs_rare_top100_genus_Mar_XL =subset_samples(physeq_rs_rare_top100_genus, Time =="Mar")%>% 
  subset_samples(., SoilType =="XL")  %>% prune_taxa(taxa_sums(.) > 0, .)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Mar_XL, grouping_column = "Treatment", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)


p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Mar_XL_treatment.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

#To generate plot for multiple testing corrections.
#p1= plot_corrections(kw_sig$SignfeaturesTable)
table_to_write <- kw_sig$SignfeaturesTable
#Detailed information about significantly differentially expressed features is can be obtained by;
table_to_write
#To get a stand alone visual representation of important features as obtained by random forest classifer, The plot shows Taxa description and corresponding Mean Decrease Accuracy values. This shows a bit more details since in the significant features plot we show only the ranks of these features but in this, it is the measured values of Mean Decrease Accuracy.
p1 <- plot_MDA(kw_sig$importance)
print(p1)
########################################################
#6.1Apr XL-DM  
########################################################

physeq_rs_rare_top100_genus_Apr_XL_DM=subset_samples(physeq_rs_rare_top100_genus, Time =="Apr") %>% subset_samples(., SoilType !="KL") %>% prune_taxa(taxa_sums(.) > 0, .)
sample_data(physeq_rs_rare_top100_genus_Apr_XL_DM)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Apr_XL_DM, grouping_column = "SoilType", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)


p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Apr_XL_DM.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

#
########################################################
#6.2Apr KL_DM  
########################################################

physeq_rs_rare_top100_genus_Apr_KL_DM=subset_samples(physeq_rs_rare_top100_genus, Time =="Apr") %>% subset_samples(., SoilType !="XL") %>% prune_taxa(taxa_sums(.) > 0, .)
sample_data(physeq_rs_rare_top100_genus_Apr_KL_DM)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Apr_KL_DM, grouping_column = "SoilType", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)


p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Apr_KL_DM.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#6.3Apr_DM CK-T1
########################################################
physeq_rs_rare_top100_genus_Apr_DM =subset_samples(physeq_rs_rare_top100_genus, Time =="Apr")%>% 
  subset_samples(., SoilType =="DM")  %>% prune_taxa(taxa_sums(.) > 0, .)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Apr_DM, grouping_column = "Treatment", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)


p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Apr_DM_treatment.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#6.4Apr_KL CK-T1
########################################################
physeq_rs_rare_top100_genus_Apr_KL =subset_samples(physeq_rs_rare_top100_genus, Time =="Apr")%>% 
  subset_samples(., SoilType =="KL")  %>% prune_taxa(taxa_sums(.) > 0, .)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Apr_KL, grouping_column = "Treatment", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Apr_KL_treatment.pdf"), p1, width = 30, height = 30, units = "cm")

p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#6.5Apr_XL CK-T1
########################################################
physeq_rs_rare_top100_genus_Apr_XL =subset_samples(physeq_rs_rare_top100_genus, Time =="Apr")%>% 
  subset_samples(., SoilType =="XL")  %>% prune_taxa(taxa_sums(.) > 0, .)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_Apr_XL, grouping_column = "Treatment", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("Apr_XL_treatment.pdf"), p1, width = 30, height = 30, units = "cm")

p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#7.1May XL-DM  
########################################################

physeq_rs_rare_top100_genus_May_XL_DM=subset_samples(physeq_rs_rare_top100_genus, Time =="May") %>% subset_samples(., SoilType !="KL") %>% prune_taxa(taxa_sums(.) > 0, .)
sample_data(physeq_rs_rare_top100_genus_May_XL_DM)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_May_XL_DM, grouping_column = "SoilType", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("May_XL_DM.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#7.2May KL_DM  
########################################################

physeq_rs_rare_top100_genus_May_KL_DM=subset_samples(physeq_rs_rare_top100_genus, Time =="May") %>% subset_samples(., SoilType !="XL") %>% prune_taxa(taxa_sums(.) > 0, .)
sample_data(physeq_rs_rare_top100_genus_May_KL_DM)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_May_KL_DM, grouping_column = "SoilType", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("May_KL_DM.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#7.3May_DM CK-T1
########################################################
physeq_rs_rare_top100_genus_May_DM =subset_samples(physeq_rs_rare_top100_genus, Time =="May")%>% 
  subset_samples(., SoilType =="DM")  %>% prune_taxa(taxa_sums(.) > 0, .)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_May_DM, grouping_column = "Treatment", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("May_DM_treatment.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#7.4May_KL CK-T1
########################################################
physeq_rs_rare_top100_genus_May_KL =subset_samples(physeq_rs_rare_top100_genus, Time =="May")%>% 
  subset_samples(., SoilType =="KL")  %>% prune_taxa(taxa_sums(.) > 0, .)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_May_KL, grouping_column = "Treatment", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("May_KL_treatment.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

########################################################
#7.5May_XL CK-T1
########################################################
physeq_rs_rare_top100_genus_May_XL =subset_samples(physeq_rs_rare_top100_genus, Time =="May")%>% 
  subset_samples(., SoilType =="XL")  %>% prune_taxa(taxa_sums(.) > 0, .)

deseq_sig <- differential_abundance(physeq_rs_rare_top100_genus_May_XL, grouping_column = "Treatment", output_norm = "log-relative", 
                                    pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)

p1 <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
print(p1)
ggsave(paste0("May_XL_treatment.pdf"), p1, width = 30, height = 30, units = "cm")


p2 <- plot_MDA(deseq_sig$importance)
print(p2)

######################################################################################################################################
#6CAP
#####################################################################################################################################
#bacteria
bray_not_na <- phyloseq::distance(physeq = physeq_rs_rare, method = "bray")
# CAP ordinate
cap_ord <- ordinate(
  physeq = physeq_rs_rare, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ pH + ORG + N + P + K + Fe + Mn + Cu)
# CAP plot
cap_plot <- plot_ordination(
  physeq = physeq_rs_rare, 
  ordination = cap_ord, 
  color = "SoilType", 
  axes = c(1,2)
) + 
  aes(shape = Treatment) + 
  geom_point(aes(colour = SoilType), alpha = 0.6, size = 1) + 
  geom_point(colour = "grey90", size = 0.4) + 
  scale_color_manual(values = c("#a65628", "red","darkorchid3" ,
                                "#4daf4a", "#1919ff", "#ffae19", "magenta")
  )+theme_light()+theme(aspect.ratio=1)
# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))
# Make a new graphic
p = cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .3, 
    data = arrowdf, 
    color = "#1919ff", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 1.3,  
    data = arrowdf, 
    show.legend = FALSE
  )+theme_light()+theme(aspect.ratio=1)
p

#Do a permutational ANOVA on constrained axes used in ordination
d=anova(cap_ord)
d
###############################################################################################
#fungi

bray_not_na <- phyloseq::distance(physeq = physeq_rs_fungi_rare, method = "bray")
# CAP ordinate
cap_ord <- ordinate(
  physeq = physeq_rs_fungi_rare, 
  method = "CAP",
  distance = bray_not_na,
  formula = ~ pH + ORG + N + P + K + Fe + Mn + Cu)
# CAP plot
cap_plot <- plot_ordination(
  physeq = physeq_rs_fungi_rare, 
  ordination = cap_ord, 
  color = "SoilType", 
  axes = c(1,2)
) + 
  aes(shape = Treatment) + 
  geom_point(aes(colour = SoilType), alpha = 0.6, size = 1) + 
  geom_point(colour = "grey90", size = 0.4) + 
  scale_color_manual(values = c("#a65628", "red","darkorchid3" ,
                                "#4daf4a", "#1919ff", "#ffae19", "magenta")
  )
# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))
# Make a new graphic
p1 = cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .3, 
    data = arrowdf, 
    color = "#1919ff", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 1.3,  
    data = arrowdf, 
    show.legend = FALSE
  )+theme_light()+theme(aspect.ratio=1)
p1
#Do a permutational ANOVA on constrained axes used in ordination
d=anova(cap_ord)
d

(p_all = plot_grid(p,p1,labels = c("a", "b"), ncol = 2))
ggsave("Supplementary Fig. 3.pdf", p_all,width = 15, height = 15, units = "cm")

#########################################################################################################
#boxplot for pseudomonas
physeq_rs_rare_top100_pseudomonas<- subset_taxa(physeq_rs_rare_top100, Genus=="Pseudomonas")
physeq_rs_rare_top100_pseudomonas_genus<-tax_glom(physeq_rs_rare_top100_pseudomonas, taxrank="Genus", NArm=FALSE)%>% prune_taxa(taxa_sums(.) > 0, .)

#creat pseudo
physeq_rs_rare_top100_pseudomonas_genus_mean
otu_t<-as.data.frame(t(otu_table(physeq_rs_rare_top100_pseudomonas_genus)))
otu_t["name"]<-row.names(otu_t)
sam<-as.data.frame(sample_data(physeq_rs_rare_top100_pseudomonas_genus))
sam["name"]<-row.names(sam)
pseudo<-left_join(sam,otu_t)
pseudo

pseudo$Time=factor(pseudo$Time,levels = c('Mar','Apr','May'))
pseudo$SoilType=factor(pseudo$SoilType,levels = c('DM','KL','XL'))
pseudo$Treatment=factor(pseudo$Treatment)
#
pseudo_Apr_Mar<-pseudo %>% filter(Time!="May")
pseudo_Apr_Mar

c1 = ggplot(pseudo_Apr_Mar, aes(x=Time, y=ASV3, color=Time)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="", y=paste("Pseudomonas Abundance")) + theme_light() +theme(legend.position="none")+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_color_manual(values = c("black", "red","blue"))
c1
ggsave("c1.pdf", c1,width = 5, height = 5, units = "cm")

#
pseudo_May_Mar<-pseudo %>% filter(Time!="Apr")
pseudo_May_Mar
c2 = ggplot(pseudo_May_Mar, aes(x=Time, y=ASV3, color=Time)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="", y=paste("Pseudomonas Abundance")) + theme_light() +theme(legend.position="none")+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_color_manual(values = c("black", "red","blue"))
c2
ggsave("c2.pdf", c2,width = 5, height = 5, units = "cm")

#
pseudo_Mar_XL_DM<-pseudo %>% filter(Time=="Mar") %>% filter(SoilType!="KL")
pseudo_Mar_XL_DM
c3 = ggplot(pseudo_Mar_XL_DM, aes(x=SoilType, y=ASV3, color=SoilType)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="", y=paste("Pseudomonas Abundance")) + theme_light() +theme(legend.position="none")+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_color_manual(values = c("black", "red","blue"))
c3
ggsave("c3.pdf", c3,width = 5, height = 5, units = "cm")

#
pseudo_Apr_XL_DM<-pseudo %>% filter(Time=="Apr") %>% filter(SoilType!="KL")
pseudo_Apr_XL_DM
c4 = ggplot(pseudo_Apr_XL_DM, aes(x=SoilType, y=ASV3, color=SoilType)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="", y=paste("Pseudomonas Abundance")) + theme_light() +theme(legend.position="none")+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_color_manual(values = c("black", "red","blue"))
c4
ggsave("c4.pdf", c4,width = 5, height = 5, units = "cm")

#
pseudo_Apr_KL_DM<-pseudo %>% filter(Time=="Apr") %>% filter(SoilType!="XL")
pseudo_Apr_KL_DM
c5 = ggplot(pseudo_Apr_KL_DM, aes(x=SoilType, y=ASV3, color=SoilType)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="", y=paste("Pseudomonas Abundance")) + theme_light() +theme(legend.position="none")+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_color_manual(values = c("black", "red","blue"))
c5
ggsave("c5.pdf", c5,width = 5, height = 5, units = "cm")

#
pseudo_May_XL_DM<-pseudo %>% filter(Time=="May") %>% filter(SoilType!="KL")
pseudo_May_XL_DM
c6 = ggplot(pseudo_May_XL_DM, aes(x=SoilType, y=ASV3, color=SoilType)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="", y=paste("Pseudomonas Abundance")) + theme_light() +theme(legend.position="none")+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_color_manual(values = c("black", "red","blue"))
c6
ggsave("c6.pdf", c6,width = 5, height = 5, units = "cm")

#
pseudo_May_KL_DM<-pseudo %>% filter(Time=="May") %>% filter(SoilType!="XL")
pseudo_May_KL_DM
c7 = ggplot(pseudo_Apr_KL_DM, aes(x=SoilType, y=ASV3, color=SoilType)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="", y=paste("Pseudomonas Abundance")) + theme_light() +theme(legend.position="none")+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_color_manual(values = c("black", "red","blue"))
c7
ggsave("c7.pdf", c7,width = 5, height = 5, units = "cm")

#
pseudo_Apr_XL_T1_CK<-pseudo %>% filter(Time=="Apr") %>% filter(SoilType=="XL")
pseudo_Apr_XL_T1_CK
c8 = ggplot(pseudo_Apr_XL_T1_CK, aes(x=Treatment, y=ASV3, color=Treatment)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x="", y=paste("Pseudomonas Abundance")) + theme_light() +theme(legend.position="none")+
  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  scale_color_manual(values = c("black", "red","blue"))
c8
ggsave("c8.pdf", c8,width = 5, height = 5, units = "cm")

#The details that do not affect the result are modified by retouching software.

