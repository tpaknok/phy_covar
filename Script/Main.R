setwd("C:/Users/pakno/OneDrive - University of Toronto/NutNet")

Nutnet <- read.csv("Data/comb-by-plot_2024-05-31.csv")
Cover <- read.csv("Data/full-cover_2024-05-31.csv")
Env <- read.csv("Data/comb-by-plot-clim-soil_2024-05-31.csv")
Biomass <- read.csv("Data/full-biomass_2024-05-31.csv")

###
library(tidyverse)
library(vegan)
library(rtrees)
library(ape)
library(glmmTMB)
library(phytools)
library(picante)
library(adiv)
library(lefse)
library(stringr)

library(V.PhyloMaker2)
library(U.Taxonstand)
library(kewr)

source("Script/powo_distribution.R")
source("Script/COV.AB.SES.R")

Nutnet_year <- Nutnet %>%
  group_by(site_name,site_code,block,plot) %>%
  dplyr::summarize(year_range = max(year)-min(year))

plotID <- Nutnet_year %>% filter(year_range >= 10) %>% mutate(ID = paste0(site_code,"_",block,"_",plot)) %>% ungroup() %>% select(ID)

### Correction
Cover$Taxon[Cover$Taxon == "AMPELAMUS LAEVIS"] <- "CYNANCHUM LAEVE"
Cover$Taxon[Cover$Taxon == "CHAETACANTHUS BURCHELLII"] <- "DYSCHORISTE BURCHELLII"
Cover$Taxon[Cover$Taxon == "CHAMAESYCE SP."] <- "EUPHORBIA SP."
Cover$Taxon[Cover$Taxon == "PSORALIDIUM TENUIFLORUM"] <- "PEDIOMELUM TENUIFLORUM"
Cover$Taxon[Cover$Taxon == "SPARTINA DENSIFLORA"] <- "SPOROBOLUS DENSIFLORA"
Cover$Taxon[Cover$Taxon == "STENOSIPHON LINIFOLIUS"] <- "OENOTHERA LINIFOLIUS"
Cover$Taxon <- trimws(Cover$Taxon)

Cover_year <- Cover %>%
  group_by(site_code,block,plot,subplot) %>%
  filter(trt == "Control") %>%
  dplyr::summarize(year_unique = length(unique(year))) %>%
  filter(year_unique >= 10) %>%
  mutate(ID = paste0(site_code,"_",block,"_",plot,"_",subplot)) %>%
  ungroup() %>%
  select(ID)

NAT_Status <- Cover %>%
  group_by(functional_group,local_provenance) %>%
  dplyr::summarize(n())

Cover_subset <- Cover %>%
  mutate(ID = paste0(site_code,"_",block,"_",plot,"_",subplot)) %>%
  filter(ID %in% Cover_year$ID) %>%
  filter(live == 1) %>%
  filter(Family != "NULL") %>%
  mutate(Genus = word(Taxon,sep=" ",1)) %>%
  filter(max_cover <= 100) %>%
  filter(Genus != "UNKNOWN") %>%
  filter(trt == "Control") %>%
  filter(trt != "Fence" & trt != "NPK+Fence" & trt != "Other") %>%
  filter(functional_group != "BRYOPHYTE") %>%
  filter(functional_group != "NULL") #%>%
 # filter(year_trt >= 1) #%>%
#filter(local_provenance == "NAT" | local_provenance == "INT")

Biomass_subset <- Biomass %>%
  mutate(ID = paste0(site_code,"_",block,"_",plot,"_",subplot)) %>%
  filter(ID %in% Cover_year$ID) %>%
  filter(live == 1) %>%
  filter(trt != "Fence" & trt != "NPK+Fence" & trt != "Other") %>%
  filter(category != "BRYOPHYTE") %>%
  filter(category != "NULL") %>%
  #filter(year_trt >= 1) %>%
  group_by(ID,year) %>%
  dplyr::summarize(sum_biomass = sum(mass))

#Cover_subset[grepl("_sp\\.",Cover_subset$Taxon),"Taxon"] <- paste0(Cover_subset[grepl("_sp\\.",Cover_subset$Taxon),"Taxon"],"_",Cover_subset[grepl("_sp\\.",Cover_subset$Taxon),"ID"])
word_length <- str_count(Cover_subset$Taxon,"\\w+")
Cover_subset[word_length == 1,"Taxon"] <- paste0(Cover_subset[word_length == 1,"Taxon"]," SP.")
Cover_subset$Taxon <- word(Cover_subset$Taxon,1,2,sep=" ")
Cover_subset$Taxon<- gsub(" ","_",Cover_subset$Taxon)
Cover_subset$Taxon <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", Cover_subset$Taxon, perl=TRUE)
Cover_subset$Genus <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", Cover_subset$Genus, perl=TRUE)

all_sp <- Cover_subset[,c("Family","Genus","Taxon")]
all_sp <- all_sp[!duplicated(all_sp$Taxon),]
colnames(all_sp)[c(1,2,3)] <- c("family","genus","species")
all_sp <- all_sp[,c(3,2,1)]

clean_taxon <- nameMatch(spList=gsub("_"," ",all_sp$species),spSource=database,author=F)
clean_taxon[clean_taxon$Submitted_Name == clean_taxon$Submitted_Genus,"Submitted_Name"] <- paste0(clean_taxon[clean_taxon$Submitted_Name == clean_taxon$Submitted_Genus,"Submitted_Name"]," sp.")
clean_taxon[clean_taxon$Accepted_SPNAME == clean_taxon$Submitted_Genus,"Accepted_SPNAME"] <- paste0(clean_taxon[clean_taxon$Accepted_SPNAME == clean_taxon$Submitted_Genus,"Accepted_SPNAME"]," sp.")

all_sp2 <- clean_taxon[,c("Accepted_SPNAME","Genus_in_database","Family")]
colnames(all_sp2) <- c("species","genus","family")
all_sp2[all_sp2$genus == "Thesium","family"] <- "Santalaceae" #in the database Thesiaceae was used. 

fixed_name <- clean_taxon[clean_taxon$Submitted_Name != clean_taxon$Accepted_SPNAME,]
Cover_subset <- Cover_subset %>%
  mutate(Taxon = gsub("_"," ",Cover_subset$Taxon)) %>%
  left_join(clean_taxon[,c("Submitted_Name","Accepted_SPNAME")],by=join_by(Taxon == Submitted_Name))

Cover_subset$Accepted_SPNAME <- gsub(" ","_",Cover_subset$Accepted_SPNAME)

Cover_subset$full_name <- paste0(Cover_subset$site_code,"_",Cover_subset$block,"_",Cover_subset$plot,"_",Cover_subset$subplot)
plot_result <- list()

for (i in 1:length(unique(Cover_subset$ID))) {
  message(i)
  
  Cover_df_subset <- Cover_subset %>%
    ungroup() %>%
    filter(ID == unique(Cover_subset$ID)[[i]]) %>%
    mutate(Epithet = word(Accepted_SPNAME,2,2,"_")) %>%
    filter(Epithet != "sp.")
  
  Cover_df_subset$Genus <- word(Cover_df_subset$Accepted_SPNAME,1,sep="_")
  
  foo_df <- Cover_df_subset %>%
    group_by(Accepted_SPNAME) %>%
    dplyr::count() %>%
    as.data.frame()
  
  Cover_df_subset <- Cover_df_subset %>%
    left_join(foo_df,join_by(Accepted_SPNAME))
  
  comp_df <- Cover_df_subset %>%
    select(year,Accepted_SPNAME,max_cover) %>%
    group_by(year,Accepted_SPNAME) %>%
    dplyr::summarize(max_cover = sum(max_cover)) %>%
    pivot_wider(id_cols=year,
                names_from= Accepted_SPNAME,
                values_from = max_cover,
                values_fill = 0)
  
  plot_info <- unique(Cover_df_subset[,c("site_code","block","plot","subplot","trt")])
  plot_info$year <- length(unique(comp_df$year))
  plot_info$gamma <- ncol(comp_df[,-1])
  plot_info$alpha <- mean(vegan::specnumber(comp_df[,-1]))
  
  sp_info <- unique(Cover_df_subset[,c("Accepted_SPNAME","Family","Genus","n")])
  sp_info$avg_cover <- colMeans(comp_df[,-1])
  
  plot_result[[i]] <- SES.cor.ab(comp_df[,-1],comp_df$year,plot_info=plot_info,sp_info=sp_info)
  
}

all_SES_pair_long <- do.call(rbind,lapply(plot_result,function(x) x$SES.pair.long))
species_df <- do.call(rbind,lapply(plot_result,function(x) x$species.df))
#####

all_sp2 <- all_sp2[all_sp2$species %in% gsub("_"," ",species_df$sp),]
tree1 <- V.PhyloMaker2::phylo.maker(all_sp2,GBOTB.extended.WP,nodes.info.1.WP)

phylo_pair <- cov2cor(vcv(tree1$scenario.3))

phylo_pair_df <- expand.grid(rownames(phylo_pair),colnames(phylo_pair))
phylo_pair_df$phy.cor <- c(phylo_pair)
phylo_pair_df$label <- paste0(phylo_pair_df$Var1,"_",phylo_pair_df$Var2)
all_SES_pair_long <- do.call(rbind,lapply(plot_result,function(x) x$SES.pair.long))

all_SES_pair_long <- all_SES_pair_long %>%
  mutate(label=paste0(sp1,"_",sp2)) %>%
  left_join(phylo_pair_df[,c("label","phy.cor")],by=join_by(label))

##### Distribution
unique_sp <- unique(species_df$sp)

distr <- list()
for (k in 1:length(unique_sp)) {
  distr[[k]] <- powo_distribution(unique_sp[[k]])
}

combined_distr <- do.call(rbind,distr)
missing_sp <- unique(combined_distr$sp[is.na(combined_distr$tdwg3_code)])

missing_sp <- data.frame(original=missing_sp,fixed=missing_sp)

missing_sp$fixed<-gsub("Anisantha","Bromus",missing_sp$fixed)
missing_sp$fixed[missing_sp$original == "Anisantha diandra"] <- "Bromus diandrus"
missing_sp$fixed[missing_sp$original == "Anisantha rigida"] <- "Bromus rigidus"

missing_sp$fixed[missing_sp$original == "Agrostis alpina"] <- "Agrostis canina"
missing_sp$fixed[missing_sp$original == "Aphanes arvensis"] <- "Alchemilla arvensis"
missing_sp$fixed[missing_sp$original == "Atrema americanum"] <- "Bifora americana"
missing_sp$fixed[missing_sp$original == "Baccharis pingraea"] <- "Baccharis glutinosa"
missing_sp$fixed[missing_sp$original == "Bellardia viscosa"] <- "Parentucellia viscosa"
missing_sp$fixed[missing_sp$original == "Bombycilaena californica"] <- "Micropus californicus"
missing_sp$fixed[missing_sp$original == "Brassica nigra"] <- "Mutarda nigra"
missing_sp$fixed<-gsub("Bromopsis","Bromus",missing_sp$fixed)
missing_sp$fixed[missing_sp$original == "Ceratocephala orthoceras"] <- "Ranunculus testiculatus"
missing_sp$fixed<-gsub("Ceratochloa","Bromus",missing_sp$fixed)
missing_sp$fixed[missing_sp$original == "Ceratochloa carinata"] <- "Bromus carinatus"
missing_sp$fixed[missing_sp$original == "Ceratochloa cathartica"] <- "Bromus catharticus"

missing_sp$fixed[missing_sp$original == "Cheilanthes sieberi"] <- "Hemionitis sieberi"
missing_sp$fixed[missing_sp$original == "Cocculus carolinus"] <- "Nephroia carolina"
missing_sp$fixed[missing_sp$original == "Coreopsis tripteris"] <- "Gyrophyllum tripteris"
missing_sp$fixed[missing_sp$original == "Crabbea hirsuta"] <- "Crabbea cirsioides"
missing_sp$fixed[missing_sp$original == "Craspedia coolaminica"] <- "Craspedia gracilis"
missing_sp$fixed[missing_sp$original == "Cynapium apiifolium"] <- "Lomatium dissectum"
missing_sp$fixed[missing_sp$original == "Dactylorhiza praetermissa"] <- "Dactylorhiza majalis"
missing_sp$fixed[missing_sp$original == "Deyeuxia montanensis"] <- "Calamagrostis montanensis"
missing_sp$fixed[missing_sp$original == "Dichelostemma capitatum"] <- "Dipterostemon capitatus"
missing_sp$fixed[missing_sp$original == "Escobaria vivipara"] <- "Pelecyphora vivipara"
missing_sp$fixed[missing_sp$original == "Festuca danthonii"] <- "Festuca ambigua"
missing_sp$fixed[missing_sp$original == "Fimbristylis ovata"] <- "Abildgaardia ovata"
missing_sp$fixed[missing_sp$original == "Gamochaeta coarctata"] <- "Gamochaeta americana"
missing_sp$fixed<-gsub("Glandularia","Verbena",missing_sp$fixed)
missing_sp$fixed<-gsub("Kali","Salsola",missing_sp$fixed)
missing_sp$fixed[missing_sp$original == "Linaria vulgaris"] <- "Linaria saxatilis"
missing_sp$fixed[missing_sp$original == "Logfia californica"] <- "Filago californica"
missing_sp$fixed[missing_sp$original == "Nigritella nigra"] <- "Gymnadenia nigra"
missing_sp$fixed[missing_sp$original == "Packera tridenticulata"] <- "Packera thurberi"
missing_sp$fixed[missing_sp$original == "Pascopyrum smithii"] <- "Elymus smithii"
missing_sp$fixed[missing_sp$original == "Paspalidium geminatum"] <- "Setaria geminata"
missing_sp$fixed[missing_sp$original == "Phlox gracilis"] <- "Microsteris gracilis"
missing_sp$fixed[missing_sp$original == "Polygala polygama"] <- "Senega polygama"
missing_sp$fixed[missing_sp$original == "Polygala alba"] <- "Senega alba"

missing_sp$fixed[missing_sp$original == "Podolepis lessonii"] <- "Panaetia lessonii"
missing_sp$fixed[missing_sp$original == "Pseudosclerochloa rupestris"] <- "Puccinellia rupestris"
missing_sp$fixed[missing_sp$original == "Rochelia zeylanica"] <- "Cynoglossum zeylanicum"
missing_sp$fixed[missing_sp$original == "Rubus trivialis"] <- "Rubus flagellaris"
missing_sp$fixed[missing_sp$original == "Ruellia nudiflora"] <- "Ruellia ciliatiflora"
missing_sp$fixed[missing_sp$original == "Sida dregei"] <- "Sida lancifolia"
missing_sp$fixed[missing_sp$original == "Sporobolus densiflorus"] <- "Sporobolus montevidensis"
missing_sp$fixed[missing_sp$original == "Stipa nitida"] <- "Austrostipa nitida"
missing_sp$fixed[missing_sp$original == "Symphyotrichum eatonii"] <- "Symphyotrichum bracteolatum"
missing_sp$fixed[missing_sp$original == "Thrincia saxatilis"] <- "Leontodon saxatilis"
missing_sp$fixed[missing_sp$original == "Viola bicolor"] <- "Viola rafinesquei"
missing_sp$fixed<-gsub("Valerianella","Valeriana",missing_sp$fixed)

new_distr <- list()
for (l in 1:length(missing_sp$fixed)) {
  new_distr[[l]] <- data.frame(sp=missing_sp$original[[l]],powo_distribution(missing_sp$fixed[[l]]))
}


combined_distr <- rbind(combined_distr,do.call(rbind,new_distr)[,c("sp","tdwg3_code")])
combined_distr <- na.omit(combined_distr)
combined_distr$value <- 1

combined_distr <- unique(combined_distr)

combined_distr_wide <- combined_distr %>% 
  pivot_wider(names_from=tdwg3_code,
              values_from=value,
              values_fill=0)

#####
library(adespatial)
distri_sim <- beta.div.comp(combined_distr_wide[,-1])
distri_repl <- distri_sim$repl
distri_rich <- distri_sim$rich

unique_distri_sp <- unique(combined_distr$sp)

unique_distri_sp_df <- expand.grid(unique_distri_sp,unique_distri_sp) 
unique_distri_sp_df$repl <- c(as.matrix(distri_sim$repl))
unique_distri_sp_df$rich <- c(as.matrix(distri_sim$rich))
unique_distri_sp_df$label <- paste0(gsub(" ","_",unique_distri_sp_df$Var1),"_",gsub(" ","_",unique_distri_sp_df$Var2))

all_SES_pair_long <- all_SES_pair_long %>%
  left_join(unique_distri_sp_df[,c("repl","rich","label")],by=join_by(label))

######

all_SES_pair_long$plotID_full <- paste0(all_SES_pair_long$site_code,"_",all_SES_pair_long$block,"_",all_SES_pair_long$plot,"_",all_SES_pair_long$subplot)
all_SES_pair_long$tl_dist <- all_SES_pair_long$repl+all_SES_pair_long$rich
all_SES_pair_long$pairs <- paste0(all_SES_pair_long$local_provenance.x,"-",all_SES_pair_long$local_provenance.y)
all_SES_pair_long$pairs[all_SES_pair_long$pairs == "INT-NAT"] <- "NAT-INT"

all_SES_pair_long <- all_SES_pair_long %>%
  left_join(sr_df,join_by(plotID_full))
library(glmmTMB)
library(lme4)
library(lmerTest)

all_SES_pair_long$prop.n <- (all_SES_pair_long$n.x/all_SES_pair_long$year+all_SES_pair_long$n.y/all_SES_pair_long$year)/2
all_SES_pair_long$pair_cover <- (all_SES_pair_long$avg_cover.x+all_SES_pair_long$avg_cover.y)/2/100

all_SES_pair_long_subset <- all_SES_pair_long[,c("obs.cor","pair_cover","year","repl","tl_dist","rich","phy.cor","plotID_full","label","n.x","n.y")]
all_SES_pair_long_subset <- na.omit(all_SES_pair_long_subset)

threshold <- quantile(all_SES_pair_long_subset$tl_dist,0.5,na.rm=T)
all_SES_pair_long_subset_l <- subset(all_SES_pair_long_subset,pair_cover >= 0 & tl_dist <= threshold)
all_SES_pair_long_subset_h <- subset(all_SES_pair_long_subset,pair_cover >= 0 & tl_dist > threshold) 

# m1_cover <- lmer(obs.cor~pair_cover*phy.cor*tl_dist+
#                      (pair_cover*phy.cor*tl_dist||plotID_full)+(1|label),data=all_SES_pair_long_subset,REML=T)
# summary(m1_cover)

m1_cover_l <- glmmTMB(obs.cor~pair_cover*phy.cor+
                   (pair_cover*phy.cor||plotID_full)+(1|label),data=all_SES_pair_long_subset_l,REML=T)
summary(m1_cover_l)
performance::r2(m1_cover_l,tol=1e-100)

m1_cover_h <- glmmTMB(obs.cor~pair_cover*phy.cor+
                   (pair_cover*phy.cor||plotID_full)+(1|label),data=all_SES_pair_long_subset_h,REML=T)
summary(m1_cover_h)
performance::r2(m1_cover_h)

m1_cover_h <- glmmTMB(obs.cor~pair_cover+phy.cor+
                        (pair_cover+phy.cor||plotID_full)+(1|label),data=all_SES_pair_long_subset_h,REML=T)
summary(m1_cover_h)
min(c(max(all_SES_pair_long_subset_h$pair_cover,na.rm=T),max(all_SES_pair_long_subset_l$pair_cover,na.rm=T)))

library(ggeffects)
predict_df_l <-ggemmeans(m1_cover_l,terms=c("phy.cor[0:1,by=0.01]","pair_cover[0:0.6,by=0.01]"))
predict_df_l$pairs <- "Low"
predict_df_l$group <- as.numeric(as.character(predict_df_l$group))

p_distribution <- ggplot(data=all_SES_pair_long_subset,aes(x=phy.cor,y=pair_cover))+
  geom_point(data=all_SES_pair_long_subset_l,aes(x=phy.cor,y=pair_cover))+
  labs(y="Average cover",
       x = "Phylogenetic correlation")

plot(p_distribution)

ggsave("Figure/p_distribution.tiff", dpi=600)
p_l <- ggplot(data=predict_df_l,aes(x=x,y=group))+
  geom_raster(aes(fill=predicted))+
  geom_point(data=all_SES_pair_long_subset_l,aes(x=phy.cor,y=pair_cover))+
  theme_classic()+
  labs(y="Average cover",
       x = "Phylogenetic correlation",
       fill = "Temporal correlation")+
  theme(legend.position="bottom")+
  scale_fill_continuous(low="#009e73",high="#cc79a7",limits=c(-0.25,0.25))


plot(p_l)
ggsave("Figure/p_l.tiff", dpi=600,width=15,height=15,compression="lzw",unit="cm")

predict_df_h <-ggemmeans(m1_cover_h,terms=c("phy.cor[0:1,by=0.01]","pair_cover[0:0.6,by=0.01]"))
predict_df_h$pairs <- "High"
predict_df_h$group <- as.numeric(as.character(predict_df_h$group))

p_h <- ggplot(data=predict_df_h,aes(x=x,y=group))+
  geom_raster(aes(fill=predicted))+
  geom_point(data=all_SES_pair_long_subset_h,aes(x=phy.cor,y=pair_cover))+
  theme_classic()+
  labs(y="Average cover",
       x = "Phylogenetic correlation",
       fill = "Temporal correlation")+
  theme(legend.position="bottom")+
  scale_fill_continuous(low="#009e73",high="#cc79a7",limits=c(-0.25,0.25))

plot(p_h)
ggsave("Figure/p_h.tiff", dpi=600, width=15,height=15,compression="lzw",unit="cm")

predict_df <- rbind(predict_df_l,predict_df_h)
predict_df$facet_name = paste0("Average abundance= ",predict_df$group)
predict_df <- subset(predict_df, group == 0.01 | group == 0.25 | group == 0.5)
p <- ggplot(data=predict_df,aes(x=x,y=predicted,colour=pairs,fill=pairs))+
  geom_line()+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),colour="transparent",alpha=0.2)+
  facet_wrap(~facet_name)+
  theme_bw()+
  scale_color_manual(values=c("#cc79a7","#009e73"))+
  scale_fill_manual(values=c("#cc79a7","#009e73"))+
  labs(y="Temporal correlation",
       x = "Phylogenetic correlation",
       colour = "Native distribution dissimilarity",
       fill = "Native distribution dissimilarity")+
  theme(legend.position="bottom")
plot(p)

ggsave("Figure/p_line.tiff", dpi=600,width=15,height=15,compression="lzw",unit="cm")

####
library(phytools)

lambda <- signal_sig <- K <- NULL
for (i in 1:length(unique(all_SES_pair_long$plotID_full))){
subset_data <- subset(all_SES_pair_long,plotID_full==unique(all_SES_pair_long$plotID_full)[[i]])
subset_data <- unique(subset_data[,c("sp1","avg_cover.x")])
avg_cover <- subset_data$avg_cover.x
names(avg_cover) <- subset_data$sp1

phylosig_result <- phylosig(tree1$scenario.3,avg_cover,method="lambda",test=TRUE)
phylosig_result_K <- phylosig(tree1$scenario.3,avg_cover,method="K",test=F)

lambda <- c(lambda,phylosig_result$lambda)
K <- c(K,phylosig_result_K)
signal_sig <- c(signal_sig,ifelse(phylosig_result$P<0.05,1,0))
}

t.test(lambda,mu=0.05)
t.test(K)

null_lambda <- NULL
for (i in 1:length(unique(all_SES_pair_long$plotID_full))){
  for (sim in 1:1000) {
  subset_data <- subset(all_SES_pair_long,plotID_full==unique(all_SES_pair_long$plotID_full)[[i]])
  subset_data <- unique(subset_data[,c("sp1","avg_cover.x")])
  avg_cover <- subset_data$avg_cover.x
  names(avg_cover) <- subset_data$sp1[sample(1:length(subset_data$sp1),length(subset_data$sp1))]
  
  phylosig_result <- phylosig(tree1$scenario.3,avg_cover,method="lambda",test=TRUE)
  null_lambda <- rbind(null_lambda,data.frame(null_lambda=phylosig_result$lambda,sim=sim,i=i))
  }
}

###########
m2 <- lmer(pair_cover~tl_dist*log(sr)*phy.cor+(tl_dist*phy.cor||plotID_full),data=all_SES_pair_long_subset,REML=T)
summary(m2)

library(piecewiseSEM)

sem <- psem(m1_cover,m2)
summary(sem)
######
r <- NULL
result_df <- NULL
p_repl <- p_rich <- p_dist<- p5 <- p6 <- NULL
eff_repl <- NULL
p_repl_gam <- p_rich_gam <- p_dist_gam <- NULL
table_low <- table_median <- table_high <- NULL

library(glmmTMB)
all_SES_pair_long$scaled_obs_cor <- (all_SES_pair_long$obs.cor +1)/2

for (l in 1:length(unique(all_SES_pair_long$plotID_full))) {
  subset_all_SES_pair <- subset(all_SES_pair_long,plotID_full == unique(all_SES_pair_long$plotID_full)[[l]])
  
  if(nrow(subset_all_SES_pair) <= 30) {
    next
  }
  
  m <- lm(pair_cover~tl_dist*phy.cor,data=subset_all_SES_pair)
  summary(m)
  
  anova_m <- car::Anova(m,white.adjust=T,type="III")
  result_df <- rbind(result_df,data.frame(var=names(coef(m)),
                                          summary(m)$coefficients,
                                          l=l,
                                          anova_p= na.omit(anova_m$`Pr(>F)`),
                                          min_phy_cor = min(subset_all_SES_pair$phy.cor)))
  
  #m <- glmmTMB(SES.pair~scale(phy.cor)*scale(repl),,data=subset_all_SES_pair)
  #summary(m)
  #result_df <- rbind(result_df,data.frame(var=names(fixef(m)$cond),summary(m)$coefficients$cond,l=l))
  
  # m1 <- rq(SES.pair~scale(phy.cor)*scale(repl),data=subset_all_SES_pair,tau=c(.25,.5,.75))
  # summary(m1)
  # 
  # table_low <- rbind(table_low,data.frame(vars=rownames(summary(m1)[[1]]$coefficients),summary(m1)[[1]]$coefficients,l=l))
  # table_median <- rbind(table_median,data.frame(vars=rownames(summary(m1)[[2]]$coefficients),summary(m1)[[2]]$coefficients,l=l))
  # table_high <- rbind(table_high,data.frame(vars=rownames(summary(m1)[[3]]$coefficients),summary(m1)[[3]]$coefficients,l=l))

  # m <- glmmTMB(scaled_obs_cor ~ phy.cor*repl, data=subset_all_SES_pair, family=beta_family())
  # summary(m)
  # result_df <- rbind(result_df,data.frame(var=names(fixef(m)$cond),summary(m)$coefficients$cond,l=l))
  
  #result_df <- rbind(result_df,var=data.frame(var=names(coef(m)),coef=coefficients(m),l=l))
  r <- c(r,cor.test(subset_all_SES_pair$phy.cor,subset_all_SES_pair$SES.pair)$estimate)
  # p_repl <- c(p_repl,summary(m)$coefficients[3,4])
  # eff_repl <- c(eff_repl,summary(m)$coefficients[3,1])
  # p_rich <- c(p_rich,summary(m)$coefficients[4,4])
  # p_dist <- c(p_dist,summary(m)$coefficients[2,4])
  # p5 <- c(p5,summary(m)$coefficients[5,4])
  # p6<- c(p6,summary(m)$coefficients[6,4])
  # 
  # m <- gam(SES.pair~s(dist,repl,k=5),data=subset_all_SES_pair)
  # summary(m)
  # # 
  # vis.gam(m,view=c("dist","repl"),theta=180,n.grid=50,lwd=0.4)
  # p_dist_gam <- c(p_dist_gam,summary(m)$s.table[1,4])
  # p_repl_gam <- c(p_repl_gam,summary(m)$s.table[2,4])
  #p_rich_gam <- c(p_rich_gam,summary(m)$s.table[1,4])
}

result_df_overall <- result_df %>%
  mutate(sig=ifelse(anova_p < 0.05,1,0)) %>%
  group_by(var) %>%
  dplyr::summarise(sum_sig=sum(sig))

coef_df <- subset(result_df,var=="phy.cor:repl" & anova_p < 0.05)

summary(lm(Estimate~1,weights=1/(Std..Error^2),data=coef_df))

table_low_overall <- table_low %>%
  mutate(sig=ifelse(sign(lower.bd) == sign(upper.bd),1,0)) %>%
  group_by(vars) %>%
  dplyr::summarise(sum_sig=sum(sig))

table_median_overall <- table_median %>%
  mutate(sig=ifelse(sign(lower.bd) == sign(upper.bd),1,0)) %>%
  group_by(vars) %>%
  dplyr::summarise(sum_sig=sum(sig))

table_high_overall <- table_high %>%
  mutate(sig=ifelse(sign(lower.bd) == sign(upper.bd),1,0)) %>%
  group_by(vars) %>%
  dplyr::summarise(sum_sig=sum(sig))

see <- subset(table_low, vars == "scale(phy.cor):scale(repl)")
see <- see[sign(see$lower.bd) == sign(see$upper.bd),]
length(which(p_repl < 0.05))/length(p_repl)

length(which(p_rich<0.05))/length(p_rich)

length(which(p_dist<0.05))/length(p_dist)

length(which(p5<0.05))
length(which(p6<0.05))

sum(eff_repl[which(p_repl < 0.05)] > 0)

###
r <- NULL
result_df <- NULL
p_repl <- p_rich <- p_dist<- p5 <- p6 <- NULL
eff_repl <- NULL
p_repl_gam <- p_rich_gam <- p_dist_gam <- NULL
table_low <- table_median <- table_high <- NULL
library(glmmTMB)
all_SES_pair_long$scaled_obs_cor <- (all_SES_pair_long$obs.cor +1)/2

for (l in 1:length(unique(all_SES_pair_long$plotID_full))) {
  subset_all_SES_pair <- subset(all_SES_pair_long,plotID_full == unique(all_SES_pair_long$plotID_full)[[l]] & n.x > 2 & n.y> 2)
  subset_all_SES_pair <- subset(subset_all_SES_pair, pairs == "NAT-NAT" | pairs == "NAT-INT")
  
  pair_n <- table(subset_all_SES_pair$pairs)
  
  if(length(pair_n) < 2) {
    next
  }
  
  if(pair_n["NAT-NAT"] < 15 | pair_n["NAT-INT"] < 15) {
    next
  }
  
  m <- lm(SES.pair~scale(phy.cor)*pairs,data=subset_all_SES_pair)
  summary(m)
  
  anova_m <- car::Anova(m,white.adjust=T,type="III")
  result_df <- rbind(result_df,data.frame(var=names(coef(m)),
                                          summary(m)$coefficients,
                                          l=l,
                                          anova_p= na.omit(anova_m$`Pr(>F)`),
                                          min_phy_cor = min(subset_all_SES_pair$phy.cor)))
  
  #m <- glmmTMB(SES.pair~scale(phy.cor)*scale(repl),,data=subset_all_SES_pair)
  #summary(m)
  #result_df <- rbind(result_df,data.frame(var=names(fixef(m)$cond),summary(m)$coefficients$cond,l=l))
  
  # m1 <- rq(SES.pair~scale(phy.cor)*scale(repl),data=subset_all_SES_pair,tau=c(.25,.5,.75))
  # summary(m1)
  # 
  # table_low <- rbind(table_low,data.frame(vars=rownames(summary(m1)[[1]]$coefficients),summary(m1)[[1]]$coefficients,l=l))
  # table_median <- rbind(table_median,data.frame(vars=rownames(summary(m1)[[2]]$coefficients),summary(m1)[[2]]$coefficients,l=l))
  # table_high <- rbind(table_high,data.frame(vars=rownames(summary(m1)[[3]]$coefficients),summary(m1)[[3]]$coefficients,l=l))
  
  # m <- glmmTMB(scaled_obs_cor ~ phy.cor*repl, data=subset_all_SES_pair, family=beta_family())
  # summary(m)
  # result_df <- rbind(result_df,data.frame(var=names(fixef(m)$cond),summary(m)$coefficients$cond,l=l))
  
  #result_df <- rbind(result_df,var=data.frame(var=names(coef(m)),coef=coefficients(m),l=l))
  r <- c(r,cor.test(subset_all_SES_pair$phy.cor,subset_all_SES_pair$SES.pair)$estimate)
  # p_repl <- c(p_repl,summary(m)$coefficients[3,4])
  # eff_repl <- c(eff_repl,summary(m)$coefficients[3,1])
  # p_rich <- c(p_rich,summary(m)$coefficients[4,4])
  # p_dist <- c(p_dist,summary(m)$coefficients[2,4])
  # p5 <- c(p5,summary(m)$coefficients[5,4])
  # p6<- c(p6,summary(m)$coefficients[6,4])
  # 
  # m <- gam(SES.pair~s(dist,repl,k=5),data=subset_all_SES_pair)
  # summary(m)
  # # 
  # vis.gam(m,view=c("dist","repl"),theta=180,n.grid=50,lwd=0.4)
  # p_dist_gam <- c(p_dist_gam,summary(m)$s.table[1,4])
  # p_repl_gam <- c(p_repl_gam,summary(m)$s.table[2,4])
  #p_rich_gam <- c(p_rich_gam,summary(m)$s.table[1,4])
}

result_df_overall <- result_df %>%
  mutate(sig=ifelse(anova_p < 0.05,1,0)) %>%
  group_by(var) %>%
  dplyr::summarise(sum_sig=sum(sig))

coef_df <- subset(result_df,var=="scale(phy.cor):pairsNAT-NAT")
coef_df <- subset(result_df,var=="scale(phy.cor)")
coef_df <- subset(result_df,var=="pairsNAT-NAT")
t <- lm(Estimate~1,weights=1/Std..Error,data=coef_df)
summary(t)
###

species_df$plotID_full <- paste0(species_df$site_code,"_",species_df$block,"_",species_df$plot,"_",species_df$subplot)
phylo_sig_result <- NULL
for (l in 1:length(unique(species_df$plotID_full))) {
  message(l)
  subset_species <- subset(species_df,plotID_full == unique(species_df$plotID_full)[[l]])
  tree_subset <- keep.tip(tree1$scenario.3,tip=unique(subset_species$sp))
  
  if (nrow(subset_species) >= 15) {
  x <- subset_species$mean_abundance
  names(x) <- subset_species$sp
  ab_lambda <- phylosig(tree_subset,x,method="K",test=T)
  
  x <- subset_species$var
  names(x) <- subset_species$sp
  var_lambda <- phylosig(tree_subset,x,method="K",test=T)
  
  x <- subset_species$coef.SES
  names(x) <- subset_species$sp
  coef_lambda <- phylosig(tree_subset,x,method="K",test=T)
  
  x <- sqrt(subset_species$var)/subset_species$mean_abundance
  names(x) <- subset_species$sp
  CV_lambda <- phylosig(tree_subset,x,method="K",test=T)
  
  phylo_sig_result <- rbind(phylo_sig_result,
                            data.frame(ab_lambda=ab_lambda$K,
                                       ab_lambda_p = ab_lambda$P,
                                       var_lambda=var_lambda$K,
                                       var_lambda_p = var_lambda$P,
                                       coef_lambda=coef_lambda$K,
                                       coef_lambda_p = coef_lambda$P,
                                       CV_lambda=CV_lambda$K,
                                       CV_lambda_p = CV_lambda$P,
                                       l=l,
                                       n=nrow(subset_species))
  )
  } else {
    next
  }
}

t.test(phylo_sig_result$ab_lambda)
t.test(phylo_sig_result$var_lambda)
t.test(phylo_sig_result$coef_lambda)
t.test(phylo_sig_result$CV_lambda)

