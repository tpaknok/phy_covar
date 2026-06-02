setwd("C:/Users/pakno/OneDrive - University of Toronto/phy_covar")

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
load("C:/Users/pakno/OneDrive - University of Toronto/phy_covar/Data/Plants_WP.rdata")

clean_taxon <- nameMatch(spList=gsub("_"," ",all_sp$species),spSource=database,author=F)
clean_taxon[clean_taxon$Submitted_Name == clean_taxon$Submitted_Genus,"Submitted_Name"] <- paste0(clean_taxon[clean_taxon$Submitted_Name == clean_taxon$Submitted_Genus,"Submitted_Name"]," sp.")
clean_taxon[clean_taxon$Accepted_SPNAME == clean_taxon$Submitted_Genus,"Accepted_SPNAME"] <- paste0(clean_taxon[clean_taxon$Accepted_SPNAME == clean_taxon$Submitted_Genus,"Accepted_SPNAME"]," sp.")
length(unique(clean_taxon$Accepted_SPNAME)) - table(clean_taxon$NOTE)[[1]] #number of species in the analyses

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
  sp_info$var_cover <- apply(comp_df[,-1],2,var)
  sp_info$pop_stab <- sp_info$avg_cover/sqrt(sp_info$var_cover)
  
  plot_result[[i]] <- SES.cor.ab(comp_df[,-1],comp_df$year,plot_info=plot_info,sp_info=sp_info)
  
}

all_SES_pair_long <- do.call(rbind,lapply(plot_result,function(x) x$SES.pair.long))
species_df <- do.call(rbind,lapply(plot_result,function(x) x$species.df))

#####
fam_df <- data.frame(sp=c(all_SES_pair_long$sp1,all_SES_pair_long$sp2),
                     family=c(all_SES_pair_long$Family.x,all_SES_pair_long$Family.y))
fam_df <- unique(fam_df)
fam_df <- fam_df %>%
  left_join(unique(Cover_subset[,c("Accepted_SPNAME","functional_group")]),by=join_by(sp==Accepted_SPNAME))

table(fam_df$functional_group)
#####
library(V.PhyloMaker2)
all_sp2 <- all_sp2[all_sp2$species %in% gsub("_"," ",species_df$sp),]
tree1 <- phylo.maker(all_sp2,GBOTB.extended.WP,nodes.info.1.WP,scenario="S1")
table(tree1$species.list$status)
phylo_pair <- cophenetic.phylo(tree1$scenario.1)

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

all_SES_pair_long$plotID_full <- paste0(all_SES_pair_long$site_code,"_",all_SES_pair_long$block,"_",all_SES_pair_long$plot,"_",all_SES_pair_long$subplot)
all_SES_pair_long$tl_dist <- all_SES_pair_long$repl+all_SES_pair_long$rich
all_SES_pair_long$pairs <- paste0(all_SES_pair_long$local_provenance.x,"-",all_SES_pair_long$local_provenance.y)

######
no_distri_sp <- fam_df[!gsub("_"," ",fam_df$sp) %in% unique_distri_sp,]

no_distri_pair <- data.frame(name=c(all_SES_pair_long[,"sp1"],all_SES_pair_long[,"sp2"]),
                             site_code=rep(all_SES_pair_long[,"site_code"],2),
                             plotID_full=rep(all_SES_pair_long[,"plotID_full"],2))

no_distri_pair <- no_distri_pair[no_distri_pair$name %in% no_distri_sp$sp,]

summary_no_distri_sp <- no_distri_pair %>%   
  group_by(name) %>%
  summarize(n_site = n_distinct(site_code),
            n_plot = n_distinct(plotID_full))
#### synchrony analysis ----
#all_SES_pair_long$pairs[all_SES_pair_long$pairs == "INT-NAT"] <- "NAT-INT"

# all_SES_pair_long <- all_SES_pair_long %>%
#   left_join(sr_df,join_by(plotID_full))
library(glmmTMB)
library(lme4)
library(lmerTest)

all_SES_pair_long$prop.n <- (all_SES_pair_long$n.x/all_SES_pair_long$year+all_SES_pair_long$n.y/all_SES_pair_long$year)/2
all_SES_pair_long$pair_cover <- (all_SES_pair_long$avg_cover.x+all_SES_pair_long$avg_cover.y)/2/100
all_SES_pair_long$min_cover <- apply(all_SES_pair_long[,c("avg_cover.x","avg_cover.y")],1,min)
all_SES_pair_long$cover_diff <- abs(all_SES_pair_long$avg_cover.x/100-all_SES_pair_long$avg_cover.y/100)

all_SES_pair_long_subset <- all_SES_pair_long[,c("sp1","sp2","gamma","obs.cor","log_Vratio","phi","pair_cover","year","repl","tl_dist","rich","phy.cor","plotID_full","label","n.x","n.y","site_code","cover_diff","min_cover","avg_cover.x","var_cover.x","pop_stab.x","pop_stab.y")]
all_SES_pair_long_subset <- na.omit(all_SES_pair_long_subset)

site_df <- Env %>%
  filter(site_code %in% all_SES_pair_long$site_code & trt == "Control") %>%
  group_by(site_code) %>%
  summarize(lat=mean(latitude),
            long=mean(longitude))

# threshold <- quantile(all_SES_pair_long_subset$tl_dist,0.5,na.rm=T)
# all_SES_pair_long_subset$group <- ifelse(all_SES_pair_long_subset$tl_dist == 1 ,"historical allopatry","historical sympatry")
# all_SES_pair_long_subset_l <- subset(all_SES_pair_long_subset,pair_cover >= 0 & tl_dist <= threshold)
# all_SES_pair_long_subset_h <- subset(all_SES_pair_long_subset,pair_cover >= 0 & tl_dist > threshold) 

all_SES_pair_long_clim <- all_SES_pair_long_subset %>%
  left_join(Env[,c("site_code","latitude","MAT_v2","TEMP_VAR_v2","MAP_v2","MAP_VAR_v2","latitude")],by=join_by("site_code"),multiple="first")

all_SES_pair_long_clim$obs.cor.beta <- (all_SES_pair_long_clim$obs.cor+1)/2
all_SES_pair_long_clim$obs.cor.beta <- (all_SES_pair_long_clim$obs.cor.beta*(nrow(all_SES_pair_long_clim)-1)+0.5)/nrow(all_SES_pair_long_clim)
all_SES_pair_long_clim$MAP_v2 <- log(all_SES_pair_long_clim$MAP_v2)
all_SES_pair_long_clim$native_sim <- 1- all_SES_pair_long_clim$tl_dist

all_SES_pair_long_clim <- subset(all_SES_pair_long_clim,phy.cor <= 400)

n <- (nrow(all_SES_pair_long_clim))
all_SES_pair_long_clim$cover_diff_adj <- (all_SES_pair_long_clim$cover_diff * (n - 1) + 0.5) / n

m1_cover_diff <- glmmTMB(cover_diff~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                      (scale(sqrt(phy.cor))*scale(native_sim)||plotID_full),
                    data=all_SES_pair_long_clim,
                    dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                    control = glmmTMBControl(parallel = 6),
                    family= beta_family(link = "logit"),
                    REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_cover_diff)
plot(simulateResiduals(m1_cover_diff))
library(lmtest)
library(sandwich)

m1_sync <- glmmTMB(obs.cor~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                      (scale(sqrt(phy.cor))*scale(native_sim)||plotID_full),
                    data=all_SES_pair_long_clim,
                    dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                    control = glmmTMBControl(parallel = 6),
                    REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_sync)
plot(simulateResiduals(m1_sync))

m1_sync2 <- glmmTMB(obs.cor~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+scale(sqrt(cover_diff))+
                     (scale(sqrt(phy.cor))*scale(native_sim)||plotID_full),
                   data=all_SES_pair_long_clim,
                   dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+scale(sqrt(cover_diff))+(1|plotID_full),
                   control = glmmTMBControl(parallel = 6),
                   REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_sync2)
plot(simulateResiduals(m1_sync2))

all_SES_pair_long_clim$pop_stab_diff <- abs(log(all_SES_pair_long_clim$pop_stab.x)-log(all_SES_pair_long_clim$pop_stab.y))

m1_pop <- glmmTMB(pop_stab_diff~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                     (scale(sqrt(phy.cor))*scale(native_sim)||plotID_full),
                   data=all_SES_pair_long_clim,
                   dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                   control = glmmTMBControl(parallel = 6),
                   family=gaussian,
                   REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_pop)
plot(simulateResiduals(m1_pop))

m1_pop2 <- glmmTMB(pop_stab_diff~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+scale(cover_diff)+
                      (scale(sqrt(phy.cor))*scale(native_sim)||plotID_full),
                    data=all_SES_pair_long_clim,
                    dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+scale(cover_diff)+(1|plotID_full),
                    control = glmmTMBControl(parallel = 6),
                  family=gaussian,
                  REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_pop2)

plot(simulateResiduals(m1_pop2))
# m1_pop <- glmmTMB(pop_stab_diff~scale(phy.cor)*scale(tl_dist)*scale(TEMP_VAR_v2)+
#                        (scale(phy.cor)*scale(tl_dist)||plotID_full)+
#                        (scale(TEMP_VAR_v2)||label),
#                   data=subset(all_SES_pair_long_clim),
#                   family=tweedie,
#                   control = glmmTMBControl(parallel = 6))
# summary(m1_pop)

# m1_pop <- glmmTMB(pop_stab_diff~scale(phy.cor)*scale(tl_dist)*scale(MAP_VAR_v2)+scale(phy.cor)*scale(tl_dist)*scale(MAP_v2)+
#                     (scale(phy.cor)*scale(tl_dist)||plotID_full)+
#                     (scale(MAP_VAR_v2)+scale(MAP_v2)||label),
#                   data=all_SES_pair_long_clim,
#                   family=tweedie,
#                   control = glmmTMBControl(parallel = 6))
# summary(m1_pop)

###
example_df <- Cover_subset %>%
  ungroup() %>%
  filter(ID == "bogong.au_3_21_D") %>%
  filter(Taxon == "Pimelea alpina" | Taxon == "Prasophyllum alpinum")

###
library(ggeffects)

predict_df <-ggemmeans(m1_cover,terms=c("phy.cor[0:781,by=1]","tl_dist[0:1,by=0.5]"),terms_to_colnames=TRUE)
predict_df$predicted <- ((predict_df$predicted*nrow(all_SES_pair_long_clim)-0.5)/(nrow(all_SES_pair_long_clim)-1))*2-1

predict_df$facet <- as.character(predict_df$facet)
predict_df$facet[predict_df$facet == "9"] <- "Precipitation seasonality = 9"
predict_df$facet[predict_df$facet == "52.5"] <- "Precipitation seasonality = 52.5"
predict_df$facet[predict_df$facet == "96"] <- "Precipitation seasonality = 96"
predict_df$facet <- factor(predict_df$facet,levels=c("Precipitation seasonality = 9","Precipitation seasonality = 52.5","Precipitation seasonality = 96"))

p <- ggplot(predict_df,aes(x=x,y=predicted))+
  geom_line(aes(colour=group))+
  facet_wrap(~facet)+
  labs(y="Synchrony",x="Phylogenetic distance",colour="Native distribution dissimilarity")+
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+
  theme_bw()+
  theme(legend.position="bottom")

plot(p)

###

library(ggeffects)

predict_df <-ggemmeans(m1_pop,terms=c("phy.cor[0:781,by=1]","MAP_VAR_v2[9:96,by=43.5]"),terms_to_colnames=TRUE)

p2 <- ggplot(predict_df,aes(x=x,y=predicted))+
  geom_line(aes(colour=group))+
  labs(y="Pairwise population stability difference",x="Phylogenetic distance",colour="Precipitation seasonality")+
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3"))+
  theme_bw()+
  theme(legend.position="bottom")

plot(p2)

###
library(ggeffects)
predict_df_l <-ggemmeans(m1_cover_l,terms=c("phy.cor[0:1,by=0.01]","pair_cover[0:0.6,by=0.01]"))
predict_df_l$pairs <- "Historically sympatric"
predict_df_l$group <- as.numeric(as.character(predict_df_l$group))
predict_df_l$metric <- "Pearson's Correlation Coefficient"
predict_df_Vratio_l <-ggemmeans(m1_cover_l_Vratio,terms=c("phy.cor[0:1,by=0.01]","pair_cover[0:0.6,by=0.01]"))
predict_df_Vratio_l$pairs <- "Historically sympatric"
predict_df_Vratio_l$group <- as.numeric(as.character(predict_df_l$group))
predict_df_Vratio_l$metric <- "log(Variance ratio)"

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
predict_df_h$pairs <- "Historically allopatric"
predict_df_h$group <- as.numeric(as.character(predict_df_h$group))
predict_df_h$metric <- "Pearson's Correlation Coefficient"
predict_df_Vratio_h <-ggemmeans(m1_cover_h_Vratio,terms=c("phy.cor[0:1,by=0.01]","pair_cover[0:0.6,by=0.01]"))
predict_df_Vratio_h$pairs <- "Historically allopatric"
predict_df_Vratio_h$group <- as.numeric(as.character(predict_df_h$group))
predict_df_Vratio_h$metric <- "log(Variance ratio)"

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

predict_df <- rbind(predict_df_l,predict_df_h,predict_df_Vratio_l ,predict_df_Vratio_h )
predict_df$group <- predict_df$group*100
predict_df$facet_name = paste0("Average cover = ",predict_df$group,"%")
predict_df <- subset(predict_df, group == 1 | group == 25 | group == 50)
p <- ggplot(data=predict_df,aes(x=x,y=predicted,colour=pairs,fill=pairs))+
  geom_line()+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),colour="transparent",alpha=0.2)+
  facet_grid(metric~facet_name)+
  theme_bw()+
  scale_color_manual(values=c("#cc79a7","#009e73"))+
  scale_fill_manual(values=c("#cc79a7","#009e73"))+
  labs(y="Synchrony",
       x = "Phylogenetic relatedness",
       colour = "Pairs",
       fill = "Pairs")+
  theme(legend.position="bottom")
plot(p)

ggsave("Figure/p_line.tiff", dpi=600,width=18,height=18,compression="lzw",unit="cm")

#### pop stability analysis ----
library(phytools)
pop_stab_df <- NULL
lambda_avg <- lambda_var <- K_avg <- K_var <- lambda_CV <- K_CV <- NULL

##### observed ----
for (i in 1:length(unique(all_SES_pair_long_clim$plotID_full))){
  lambda_avg_null <- lambda_var_null <- K_avg_null <- K_var_null <- lambda_CV_null <- K_CV_null <- NULL
  
  plotID_full <- unique(all_SES_pair_long_clim$plotID_full)[[i]]
  subset_data <- subset(all_SES_pair_long_clim,plotID_full==unique(all_SES_pair_long_clim$plotID_full)[[i]])

  site_code <- unique(subset_data$site_code)
  subset_env <- subset_data[,c("MAT_v2","TEMP_VAR_v2","MAP_v2","MAP_VAR_v2","tl_dist","year")]
  max.phy.cor <- max(subset_data$phy.cor)
  mean.phy.cor <- mean(subset_data$phy.cor)
  
  tl_dist_sd <- sd(subset_env$tl_dist)
  subset_data <- unique(subset_data[,c("sp1","avg_cover.x","var_cover.x")])

if (nrow(subset_data) < 20) {
  next()
}
  avg_cover <- subset_data$avg_cover.x
  var_cover <- subset_data$var_cover.x
  CV_cover <- sqrt(var_cover)/avg_cover
  
  names(avg_cover) <- subset_data$sp1
  names(CV_cover) <- subset_data$sp1
  names(var_cover) <- subset_data$sp1
  gamma <- length(unique(subset_data$sp1))
  
  phylosig_result_avg <- phylosig(tree1$scenario.1,avg_cover,method="lambda",test=F)
  phylosig_result_K_avg <- phylosig(tree1$scenario.1,avg_cover,method="K",test=F)
  phylosig_result_var <- phylosig(tree1$scenario.1,var_cover,method="lambda",test=F)
  phylosig_result_K_var <- phylosig(tree1$scenario.1,var_cover,method="K",test=F)
  phylosig_result_CV <- phylosig(tree1$scenario.1,log(1/CV_cover),method="lambda",test=F)
  phylosig_result_K_CV <- phylosig(tree1$scenario.1,log(1/CV_cover),method="K",test=F)
  
  lambda_avg <- c(lambda_avg,phylosig_result_avg$lambda)
  lambda_var <- c(lambda_var,phylosig_result_var$lambda)
  K_avg <- c(K_avg,phylosig_result_K_avg)
  K_var <- c(K_var,phylosig_result_K_var)
  lambda_CV <- c(lambda_CV,phylosig_result_CV$lambda)
  K_CV <- c(K_CV,phylosig_result_K_CV)
  
  ##### null lambda ----
  for (iter in 1:100) {
    pos <- sample(1:nrow(subset_data),replace=F)
    names(avg_cover) <- subset_data$sp1[pos]
    names(CV_cover) <- subset_data$sp1[pos]
    names(var_cover) <- subset_data$sp1[pos]
    
    phylosig_result_avg_null <- phylosig(tree1$scenario.1,avg_cover,method="lambda",test=F)
    phylosig_result_K_avg_null <- phylosig(tree1$scenario.1,avg_cover,method="K",test=F)
    phylosig_result_var_null <- phylosig(tree1$scenario.1,var_cover,method="lambda",test=F)
    phylosig_result_K_var_null <- phylosig(tree1$scenario.1,var_cover,method="K",test=F)
    phylosig_result_CV_null <- phylosig(tree1$scenario.1,log(1/CV_cover),method="lambda",test=F)
    phylosig_result_K_CV_null <- phylosig(tree1$scenario.1,log(1/CV_cover),method="K",test=F)
    
    lambda_avg_null <- c(lambda_avg_null,phylosig_result_avg_null$lambda)
    lambda_var_null <- c(lambda_var_null,phylosig_result_var_null$lambda)
    K_avg_null <- c(K_avg_null,phylosig_result_K_avg_null)
    K_var_null <- c(K_var_null,phylosig_result_K_var_null)
    lambda_CV_null <- c(lambda_CV_null,phylosig_result_CV_null$lambda)
    K_CV_null <- c(K_CV_null,phylosig_result_K_CV_null)
  }
  
  lambda_avg_null_mean <- mean(lambda_avg_null)
  lambda_var_null_mean <- mean(lambda_var_null)
  K_avg_null_mean <- mean(K_avg_null)
  K_var_null_mean <- mean(K_var_null)
  lambda_CV_null_mean <- mean(lambda_CV_null)
  K_CV_null_mean <- mean(K_CV_null)
  
  lambda_avg_null_sd <- sd(lambda_avg_null)
  lambda_var_null_sd <- sd(lambda_var_null)
  K_avg_null_sd <- sd(K_avg_null)
  K_var_null_sd <- sd(K_var_null)
  lambda_CV_null_sd <- sd(lambda_CV_null)
  K_CV_null_sd <- sd(K_CV_null)
  
  pop_stab_df <- rbind(pop_stab_df,data.frame(plotID_full = plotID_full,
                                              site_code = site_code,
                                              lambda=phylosig_result_CV$lambda,
                                              K=phylosig_result_K_CV,
                                              SES_lambda= (phylosig_result_CV$lambda-lambda_CV_null_mean)/lambda_CV_null_sd,
                                              SES_K= (phylosig_result_K_CV-K_CV_null_mean)/K_CV_null_sd,
                                              gamma=gamma,
                                              max.phy.cor=max.phy.cor,
                                              mpd = mean.phy.cor,
                                              tl_dist_sd = tl_dist_sd,
                                              as.list(colMeans(subset_env))))
  
}

#### null ----
lambda_avg_null <- lambda_var_null <- K_avg_null <- K_var_null <- lambda_CV_null <- K_CV_null <- NULL
all_lambda_avg_null <- all_lambda_var_null <- all_K_avg_null <- all_K_var_null <- all_lambda_CV_null <- all_K_CV_null <- NULL

for (iter in 1:1000) {
  message(iter)
  for (i in 1:length(unique(all_SES_pair_long_clim$plotID_full))){
    plotID_full <- unique(all_SES_pair_long_clim$plotID_full)[[i]]
    subset_data <- subset(all_SES_pair_long_clim,plotID_full==unique(all_SES_pair_long_clim$plotID_full)[[i]])
    
    site_code <- unique(subset_data$site_code)
    subset_env <- subset_data[,c("MAT_v2","TEMP_VAR_v2","MAP_v2","MAP_VAR_v2","tl_dist","year")]
    max.phy.cor <- max(subset_data$phy.cor)
    
    subset_data <- unique(subset_data[,c("sp1","avg_cover.x","var_cover.x")])
    
    if (nrow(subset_data) < 20) {
      next()
    }
    avg_cover <- subset_data$avg_cover.x
    var_cover <- subset_data$var_cover.x
    CV_cover <- sqrt(var_cover)/avg_cover
    
    pos <- sample(1:nrow(subset_data),replace=F)
    names(avg_cover) <- subset_data$sp1[pos]
    names(CV_cover) <- subset_data$sp1[pos]
    names(var_cover) <- subset_data$sp1[pos]
    
    phylosig_result_avg <- phylosig(tree1$scenario.1,avg_cover,method="lambda",test=TRUE)
    phylosig_result_K_avg <- phylosig(tree1$scenario.1,avg_cover,method="K",test=F)
    phylosig_result_var <- phylosig(tree1$scenario.1,var_cover,method="lambda",test=TRUE)
    phylosig_result_K_var <- phylosig(tree1$scenario.1,var_cover,method="K",test=F)
    phylosig_result_CV <- phylosig(tree1$scenario.1,log(1/CV_cover),method="lambda",test=TRUE)
    phylosig_result_K_CV <- phylosig(tree1$scenario.1,log(1/CV_cover),method="K",test=TRUE)
    
    lambda_avg_null <- c(lambda_avg_null,phylosig_result_avg$lambda)
    lambda_var_null <- c(lambda_var_null,phylosig_result_var$lambda)
    K_avg_null <- c(K_avg_null,phylosig_result_K_avg)
    K_var_null <- c(K_var_null,phylosig_result_K_var)
    lambda_CV_null <- c(lambda_CV_null,phylosig_result_CV$lambda)
    K_CV_null <- c(K_CV_null,phylosig_result_K_CV$K)
  }
  
  all_lambda_avg_null <- c(all_lambda_avg_null,mean(lambda_avg_null))
  all_lambda_var_null <- c(all_lambda_var_null,mean(lambda_var_null))
  all_K_avg_null <- c(all_K_avg_null,mean(K_avg_null))
  all_K_var_null <- c(all_K_var_null,mean(K_var_null))
  all_lambda_CV_null <- c(all_lambda_CV_null,mean(lambda_CV_null))
  all_K_CV_null <- c(all_K_CV_null,mean(K_CV_null))
}

t.test(lambda_avg)
t.test(lambda_var)

t.test(K_avg)
t.test(K_var)

t.test(lambda_CV)
t.test(K_CV)

#### new model 26122025 ----

t.test(pop_stab_df$SES_lambda)
t.test(pop_stab_df$SES_K)

CV_lambda_m <- glmmTMB(lambda~scale(tl_dist)+scale(tl_dist_sd)+scale(MAP_VAR_v2)+scale(max.phy.cor)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(gamma)+(1|site_code),
                       dispformula = ~scale(tl_dist)+scale(tl_dist_sd)+scale(max.phy.cor)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(gamma),
                       data=pop_stab_df,REML=T,family=beta_family())
summary(CV_lambda_m)
plot(simulateResiduals(CV_lambda_m))

CV_lambda_SES_m <- glmmTMB(SES_lambda~scale(tl_dist)+scale(tl_dist_sd)+scale(mpd)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(gamma)+(1|site_code),
                       dispformula = ~scale(tl_dist)+scale(tl_dist_sd)+scale(mpd)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(gamma),
                       data=pop_stab_df,REML=T)
summary(CV_lambda_SES_m)

plot(simulateResiduals(CV_lambda_SES_m))
#### old model ----
CV_lambda_m <- glmmTMB(lambda~scale(tl_dist)*scale(MAP_VAR_v2)+scale(tl_dist)*scale(TEMP_VAR_v2)+scale(tl_dist)*scale(MAT_v2)+scale(tl_dist)*scale(MAP_v2)+(1|site_code),weights=gamma,data=pop_stab_df,family=tweedie)
summary(CV_lambda_m)

library(DHARMa)
plot(simulateResiduals(CV_lambda_m))

CV_lambda_m <- glmmTMB(K~scale(tl_dist)*scale(MAP_v2)+scale(tl_dist)*scale(MAP_VAR_v2)+
                         scale(tl_dist)*scale(MAT_v2)+scale(tl_dist)*scale(TEMP_VAR_v2)+scale(phy.cor)+scale(gamma)+
                         (1|site_code),data=pop_stab_df,REML=T)
summary(CV_lambda_m)


########
physig_df <- data.frame(value=c(lambda_avg,lambda_var,K_avg,K_var))
df_size <-nrow(physig_df)/4
physig_df$metric <- rep(c(rep("Average Cover",df_size),rep("Variance",df_size)),2)
physig_df$signal <- c(rep("Pagel's λ",df_size*2),rep("Bloomberg's K",df_size*2))

library(ggplot2)
physig_p <- ggplot(physig_df)+
  geom_violin(aes(y=value,x=signal))+
  facet_wrap(~metric)+
  labs(x="Phylogenetic signal",
       y="Value")+
  theme_bw()
plot(physig_p)

ggsave("Figure/p_signal.tiff", dpi=600,width=18,height=12,compression="lzw",unit="cm")

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

###
library(tidyverse)

all_SES_pair_long$plotID_full2 <- paste0(all_SES_pair_long$site_code,"_",all_SES_pair_long$block,"_",all_SES_pair_long$plot) #without subplot ID. doesn't make a difference
Map_df <- Env %>%
  mutate(plotID_full = paste0(site_code,"_",block,"_",plot)) %>%
  filter(plotID_full %in% unique(all_SES_pair_long$plotID_full2) & trt == "Control") %>%
  group_by(plotID_full) %>%
  summarize(lat=mean(latitude),
            long=mean(longitude),
            range = max(year)-min(year),
            earliest = min(year),
            latest = max(year)) %>%
  left_join(unique(all_SES_pair_long[,c("plotID_full2","year")]),by=join_by(plotID_full==plotID_full2)) %>%
  mutate(indicator=range+1-year)

library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

map <- ggplot(data = world) +
  geom_sf(colour="transparent",fill="grey75")+
  geom_point(data=Map_df,aes(x=long,y=lat,colour=year),size=2)+
  scale_color_gradient2(low="#f8765c",mid="#982d80",high="#221150",midpoint=13.5)+
  labs(colour="Sampled years")+
  theme_void()+
  theme(legend.position="bottom")
plot(map)

ggsave("Figure/p_map.tiff", dpi=600,compression="lzw",width=8,height=5,bg="white")
