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
#library(lefse)
library(stringr)

library(V.PhyloMaker2)
library(U.Taxonstand)
library(kewr)

source("Script/powo_distribution.R")
source("Script/calc_temporal_metrics.R")

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
plot_result <- null_result <- list()

n_null <- 100

for (i in 1:length(unique(Cover_subset$ID))) {
  message(i)
  
  Cover_df_subset <- Cover_subset %>%
    ungroup() %>%
    filter(ID == unique(Cover_subset$ID)[[i]]) %>%
    mutate(Epithet = word(Accepted_SPNAME, 2, 2, "_")) %>%
    filter(Epithet != "sp.")
  
  Cover_df_subset$Genus <- word(Cover_df_subset$Accepted_SPNAME, 1, sep = "_")
  
  foo_df <- Cover_df_subset %>%
    group_by(Accepted_SPNAME) %>%
    dplyr::count() %>%
    as.data.frame()
  
  Cover_df_subset <- Cover_df_subset %>%
    left_join(foo_df, join_by(Accepted_SPNAME))
  
  comp_df <- Cover_df_subset %>%
    select(year, Accepted_SPNAME, max_cover) %>%
    group_by(year, Accepted_SPNAME) %>%
    dplyr::summarize(max_cover = sum(max_cover), .groups = "drop") %>%
    pivot_wider(id_cols = year,
                names_from  = Accepted_SPNAME,
                values_from = max_cover,
                values_fill = 0)
  
  no_var_sp <- names(which(apply(comp_df[,-1],2,FUN=sd) == 0))
  comp_df   <- comp_df %>% select(-all_of(no_var_sp))
  
  plot_info <- unique(Cover_df_subset[, c("site_code", "block", "plot", "subplot", "trt")])
  plot_info$year  <- length(unique(comp_df$year))
  plot_info$gamma <- ncol(comp_df[, -1])+length(no_var_sp)
  plot_info$alpha <- mean(vegan::specnumber(comp_df[, -1]))
  plot_info$no_var_sp <- length(no_var_sp)
  
  sp_info_base <- unique(Cover_df_subset[, c("Accepted_SPNAME", "Family", "Genus", "n")])
  sp_info_base <- sp_info_base %>% filter(!Accepted_SPNAME %in% no_var_sp)
  
  # observed sp_info with cover stats from real data
  sp_info_obs <- sp_info_base
  sp_info_obs$avg_cover <- colMeans(comp_df[, -1])
  sp_info_obs$var_cover <- apply(comp_df[, -1], 2, var)
  sp_info_obs$pop_stab  <- sp_info_obs$avg_cover / sqrt(sp_info_obs$var_cover)
  
  # --- observed ---
  plot_result[[i]] <- calc_temporal_metrics(
    comp_df[, -1], comp_df$year,
    plot_info = plot_info, sp_info = sp_info_obs
  )
  
  # --- null expectations ---
  null_result[[i]] <- lapply(seq_len(n_null), function(j) {
    
    sp_cols      <- comp_df[, -1]
    shuf_idx     <- sample(ncol(sp_cols))
    sp_cols_null <- sp_cols[, shuf_idx]
    colnames(sp_cols_null) <- colnames(sp_cols)
    
    # null sp_info: cover stats recomputed from shuffled columns
    sp_info_null <- sp_info_base
    sp_info_null$avg_cover <- colMeans(sp_cols_null)
    sp_info_null$var_cover <- apply(sp_cols_null, 2, var)
    sp_info_null$pop_stab  <- sp_info_null$avg_cover / sqrt(sp_info_null$var_cover)
    
    calc_temporal_metrics(
      sp_cols_null, comp_df$year,
      plot_info = plot_info, sp_info = sp_info_null
    )
  })
}

all_SES_pair_long <- do.call(rbind,lapply(plot_result,function(x) x$pair.df))
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
missing_sp$fixed[missing_sp$original == "Acmispon americanum"] <- "Acmispon_americanus"
missing_sp$fixed[missing_sp$original == "Afroaster hispida"] <- "Afroaster_hispidus"
missing_sp$fixed[missing_sp$original == "Andropogon ternarius"] <- "Chrysopogon_argutus"
missing_sp$fixed[missing_sp$original == "Anisantha fasciculata"] <- "Bromus_fasciculatus"
missing_sp$fixed[missing_sp$original == "Calamagrostis epigeios"] <- "Calamagrostis_epigejos"
missing_sp$fixed[missing_sp$original == "Coreopsis tripteris"] <- "Anacis_tripteris"
missing_sp$fixed[missing_sp$original == "Cotula alpina"] <- "Leptinella_alpina"
missing_sp$fixed[missing_sp$original == "Crabbea hirsuta"] <- "Crabbea_nana"
missing_sp$fixed[missing_sp$original == "Croton setigerus"] <- "Croton_setiger"
missing_sp$fixed[missing_sp$original == "Cynosurus echinatus"] <- "Falona_echinata"
missing_sp$fixed[missing_sp$original == "Doellingeria ledophylla"] <- "Eucephalus_ledophyllus"
missing_sp$fixed[missing_sp$original == "Erechtites quadridentata"] <- "Senecio_quadridentatus"
missing_sp$fixed[missing_sp$original == "Conyza ramosissima"] <- "Erigeron_divaricatus"
missing_sp$fixed[missing_sp$original == "Isoetopsis graminifolia"] <- "Minuria_graminifolia"
missing_sp$fixed[missing_sp$original == "Kali tragus"] <- "Salsola tragus"
missing_sp$fixed[missing_sp$original == "Kali turgida"] <- "Salsola kali"
missing_sp$fixed[missing_sp$original == "Microseris troximoides"] <- "Nothocalais_troximoides"
missing_sp$fixed[missing_sp$original == "Nama hispidum"] <- "Nama_hispida"
missing_sp$fixed[missing_sp$original == "Neurachne alopecuroidea"] <- "Neurachne_alopecuroidea"
missing_sp$fixed[missing_sp$original == "Panicum aequinerve"] <- "Trichanthecium_aequinerve"
missing_sp$fixed[missing_sp$original == "Plantago caffra"] <- "Plantago_bellardii"
missing_sp$fixed[missing_sp$original == "Plantago myosurus"] <- "Plantago_myosuros"
missing_sp$fixed[missing_sp$original == "Crabbea hirsuta"] <- "Crabbea_nana"

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
library(tidyverse)

source("Script/calc_overlap_metrics.R")
distri_sim <- beta.div.comp(combined_distr_wide[,-1])
distri_repl <- distri_sim$repl
distri_rich <- distri_sim$rich

unique_distri_sp <- unique(combined_distr$sp)

unique_distri_sp_df <- expand.grid(unique_distri_sp,unique_distri_sp) 
unique_distri_sp_df$repl <- c(as.matrix(distri_sim$repl))
unique_distri_sp_df$rich <- c(as.matrix(distri_sim$rich))
unique_distri_sp_df$label <- paste0(gsub(" ","_",unique_distri_sp_df$Var1),"_",gsub(" ","_",unique_distri_sp_df$Var2))

sym_extent <- calc_overlap_metrics(combined_distr_wide)
sym_extent$label  <- paste0(gsub(" ", "_", sym_extent$sp2), "_", gsub(" ", "_", sym_extent$sp1))
sym_extent$label2 <- paste0(gsub(" ", "_", sym_extent$sp1), "_", gsub(" ", "_", sym_extent$sp2))

# Duplicate rows with swapped label so both orientations are covered
sym_extent_both <- rbind(
  sym_extent %>% mutate(label = label),
  sym_extent %>% mutate(label = label2)
) %>% select(-label2)

all_SES_pair_long <- all_SES_pair_long %>%
  left_join(unique_distri_sp_df[,c("repl","rich","label")],by=join_by(label)) %>%
  left_join(sym_extent_both[,c("sympatry_min","sympatry_union","range_symmetry","label")],by=join_by(label))

all_SES_pair_long$plotID_full <- paste0(all_SES_pair_long$site_code,"_",all_SES_pair_long$block,"_",all_SES_pair_long$plot,"_",all_SES_pair_long$subplot)
all_SES_pair_long$tl_dist <- all_SES_pair_long$repl+all_SES_pair_long$rich

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
library(tidyverse)

source("Script/plot_model_predictions.R")

all_SES_pair_long$prop.n <- (all_SES_pair_long$n.x/all_SES_pair_long$year+all_SES_pair_long$n.y/all_SES_pair_long$year)/2
all_SES_pair_long$pair_cover <- (all_SES_pair_long$avg_cover.x+all_SES_pair_long$avg_cover.y)/2/100
all_SES_pair_long$min_cover <- apply(all_SES_pair_long[,c("avg_cover.x","avg_cover.y")],1,min)
all_SES_pair_long$cover_diff <- abs(all_SES_pair_long$avg_cover.x/100-all_SES_pair_long$avg_cover.y/100)

all_SES_pair_long_subset <- all_SES_pair_long[,c("sp1","sp2","gamma","obs.cor","log_Vratio","phi","pair_cover","year",
                                                 "repl","tl_dist","rich",
                                                 "sympatry_min","sympatry_union","range_symmetry",
                                                 "phy.cor","plotID_full","label",
                                                 "n.x","n.y",
                                                 "site_code",
                                                 "cover_diff","min_cover",
                                                 "avg_cover.x","avg_cover.y",
                                                 "var_cover.x","var_cover.y",
                                                 "pop_stab.x","pop_stab.y")]
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
all_SES_pair_long_clim$async <- all_SES_pair_long_clim$obs.cor*-1
start <- Sys.time()
m1_cover_diff <- glmmTMB(cover_diff_adj~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                           (scale(sqrt(phy.cor))+scale(native_sim)||plotID_full),
                         data=all_SES_pair_long_clim,
                         dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                         control = glmmTMBControl(parallel = 6),
                         family= beta_family(link = "logit"),
                         REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_cover_diff)
end <- Sys.time()
plot(simulateResiduals(m1_cover_diff))

m1_cover_diff <- glmmTMB(cover_diff_adj~scale(sqrt(phy.cor))+scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                           (scale(sqrt(phy.cor))+scale(native_sim)||plotID_full),
                         data=all_SES_pair_long_clim,
                         dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                         control = glmmTMBControl(parallel = 6),
                         family= beta_family(link = "logit"),
                         REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_cover_diff)

p_pair_cover <- plot_model_predictions(m1_cover_diff,data=all_SES_pair_long_clim,family="beta",
                                       y_label="Temporal mean cover difference",linetype_mean = "solid",facet_var = "native_sim",facet_label = "Sympatry extent",nrow=2,ncol=1,dispersion_metric = c("disp"))
dir.create("Result", showWarnings = FALSE)
ggsave("Figure/Result/cover.tiff", p_pair_cover$combined,
       dpi = 600, width = 6, height =6, compression = "lzw")

m1_cover_diff2 <- glmmTMB(cover_diff_adj~scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                           (1|plotID_full),
                         data=all_SES_pair_long_clim,
                         dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                         control = glmmTMBControl(parallel = 6),
                         family= beta_family(link = "logit"),
                         REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_cover_diff2)

cor(predict(m1_cover_diff,re.form=NA),all_SES_pair_long_clim$cover_diff_adj)
cor(predict(m1_cover_diff2,re.form=NA),all_SES_pair_long_clim$cover_diff_adj)

### Synchrony ----
m1_sync <- glmmTMB(async~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                     (scale(sqrt(phy.cor))*scale(native_sim)||plotID_full),
                   data=all_SES_pair_long_clim,
                   dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                   control = glmmTMBControl(parallel = 6),
                   REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_sync)
plot(simulateResiduals(m1_sync))
start<-Sys.time()

m1_sync <- glmmTMB(async~scale(sqrt(phy.cor))+scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                     (scale(sqrt(phy.cor))+scale(native_sim)||plotID_full),
                   data=all_SES_pair_long_clim,
                   dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                   control = glmmTMBControl(parallel = 6),
                   REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_sync)
end<-Sys.time()

p_pair_sync <- plot_model_predictions(m1_sync,data=all_SES_pair_long_clim,family="gaussian",y_label = "Asynchrony",linetype_mean = "dashed",facet_var = "native_sim",facet_label = "Sympatry extent",nrow=2,ncol=1,dispersion_metric = c("disp"))
ggsave("Figure/Result/sync.tiff", p_pair_sync$combined,
       dpi = 600, width = 6, height =6, compression = "lzw")

m1_sync2 <- glmmTMB(async~scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                            (1|plotID_full),
                          data=all_SES_pair_long_clim,
                          dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                          control = glmmTMBControl(parallel = 6),
                          REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_sync2)

cor(predict(m1_sync,re.form=NA),all_SES_pair_long_clim$async)
cor(predict(m1_sync2,re.form=NA),all_SES_pair_long_clim$async)

### population stability ----
all_SES_pair_long_clim$pop_stab_diff <- abs(log(all_SES_pair_long_clim$pop_stab.x)-log(all_SES_pair_long_clim$pop_stab.y))

start<-Sys.time()
m1_pop <- glmmTMB(pop_stab_diff~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                    (scale(sqrt(phy.cor))*scale(native_sim)||plotID_full),
                  data=all_SES_pair_long_clim,
                  dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                  control = glmmTMBControl(parallel = 6),
                  family=gaussian,
                  REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_pop)
end<-Sys.time()

plot(simulateResiduals(m1_pop))

m1_pop <- glmmTMB(pop_stab_diff~scale(sqrt(phy.cor))+scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                    (scale(sqrt(phy.cor))+scale(native_sim)||plotID_full),
                  data=all_SES_pair_long_clim,
                  dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                  control = glmmTMBControl(parallel = 6),
                  family=gaussian,
                  REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_pop)
p_pair_pop <- plot_model_predictions(m1_pop,data=all_SES_pair_long_clim,family="gaussian",y_label = "Population stability difference",facet_var = "native_sim",facet_label = "Sympatry extent", nrow=2,ncol=1,dispersion_metric = c("disp"))
ggsave("Figure/Result/pop.tiff", p_pair_pop$combined,
       dpi = 600, width = 6, height =6, compression = "lzw")

m1_pop2 <- glmmTMB(pop_stab_diff~scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+
                      (1|plotID_full),
                    data=all_SES_pair_long_clim,
                    dispformula = ~scale(sqrt(phy.cor))*scale(native_sim)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|plotID_full),
                    control = glmmTMBControl(parallel = 6),
                    REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_pop2)

cor(predict(m1_pop,re.form=NA),all_SES_pair_long_clim$pop_stab_diff)
cor(predict(m1_pop2,re.form=NA),all_SES_pair_long_clim$pop_stab_diff)

plot(simulateResiduals(m1_pop2))

### community-level
library(phytools)
pop_stab_df <- NULL
lambda_avg <- lambda_sync <- lambda_CV <- K_avg <- K_sync <- K_CV <- I_avg <- I_sync <- I_CV <- NULL
niter = 10

for (i in 1:length(unique(all_SES_pair_long_clim$plotID_full))){
  message(i)
  lambda_avg_null <- lambda_sync_null <- lambda_CV_null <- K_avg_null <- K_sync_null <- K_CV_null <- I_avg_null <- I_sync_null <- I_CV_null <- NULL
  
  plotID_full <- unique(all_SES_pair_long_clim$plotID_full)[[i]]
  subset_data <- subset(all_SES_pair_long_clim,plotID_full==unique(all_SES_pair_long_clim$plotID_full)[[i]])
  
  site_code <- unique(subset_data$site_code)
  subset_env <- subset_data[,c("MAT_v2","TEMP_VAR_v2","MAP_v2","MAP_VAR_v2","native_sim","year")]
  max.phy.cor <- max(subset_data$phy.cor)
  mean.phy.cor <- mean(subset_data$phy.cor)
  VPD <- var(subset_data$phy.cor)
  native_sim_sd <- sd(subset_env$native_sim)
  
  # Mirror the pairs
  subset_mirror <- bind_rows(subset_data, dplyr::rename(subset_data, sp1 = sp2, sp2 = sp1))
  
  matrix_df <-  subset_mirror %>%
    select(sp2,sp1,async) %>%
    pivot_wider(names_from = sp2, values_from = async) %>%
    column_to_rownames("sp1")
  
  sync_sp <- as.data.frame(colMeans(matrix_df,na.rm=T))
  sync_sp$sp1 <- rownames(sync_sp) 
  colnames(sync_sp) <- c("sync","sp1")
  
  subset_data_sp1 <- unique(subset_data[,c("sp1","avg_cover.x","var_cover.x")])
  subset_data_sp2 <- unique(subset_data[,c("sp2","avg_cover.y","var_cover.y")])
  
  subset_data <- rbind(
    subset_data_sp1,
    setNames(
      subset_data_sp2[!subset_data_sp1$sp1 %in% subset_data_sp2$sp2, ],
      names(subset_data_sp1)  # force same column names
    )
  )
  
  subset_data <- subset_data %>%
    left_join(sync_sp,by=join_by(sp1))
  
  subset_data <- unique(subset_data)

    if (nrow(subset_data) < 20) {
    next()
  }
  
  avg_cover <- subset_data$avg_cover.x
  var_cover <- subset_data$var_cover.x
  CV_cover <- sqrt(var_cover)/avg_cover
  sync_cover <- subset_data$sync
  
  names(avg_cover) <- subset_data$sp1
  names(CV_cover) <- subset_data$sp1
  names(sync_cover) <- subset_data$sp1
  gamma <- length(unique(subset_data$sp1))
  
  dist_mat <- cophenetic(tree1$scenario.1)
  
  species_keep <- intersect(rownames(dist_mat), names(avg_cover))
  dist_mat <- dist_mat[species_keep, species_keep]
  # 2. Convert distances to weights (invert so close relatives have high weight)
  weight_mat <- 1 / dist_mat
  diag(weight_mat) <- 0  # diagonal must be zero
  
  weight_mat <- weight_mat[names(avg_cover),names(avg_cover)]
  phylosig_result_avg <- phylosig(tree1$scenario.1,avg_cover,method="lambda",test=F, niter=niter)
  phylosig_result_K_avg <- phylosig(tree1$scenario.1,avg_cover,method="K",test=F, niter=niter)
  phylosig_result_I_avg <- Moran.I(avg_cover, weight_mat,scaled = TRUE)
  
  #phylosig_result_sync <- phylosig(tree1$scenario.1,sync_cover,method="lambda",test=F, niter=niter)
  #phylosig_result_K_sync <- phylosig(tree1$scenario.1,sync_cover,method="K",test=F, niter=niter)
  phylosig_result_CV <- phylosig(tree1$scenario.1,log(1/CV_cover),method="lambda",test=F, niter=niter)
  phylosig_result_K_CV <- phylosig(tree1$scenario.1,log(1/CV_cover),method="K",test=F, niter=niter)
  phylosig_result_I_CV <- Moran.I(log(1/CV_cover), weight_mat,scaled = TRUE)
  
  lambda_avg <- c(lambda_avg,phylosig_result_avg$lambda)
  #lambda_sync <- c(lambda_sync,phylosig_result_sync$lambda)
  lambda_CV <- c(lambda_CV,phylosig_result_CV$lambda)
  
  K_avg <- c(K_avg,phylosig_result_K_avg)
  #K_sync <- c(K_sync,phylosig_result_K_sync)
  K_CV <- c(K_CV,phylosig_result_K_CV)
  
  I_avg <- c(I_avg, phylosig_result_I_avg$observed)
  I_CV <- c(I_CV, phylosig_result_I_CV$observed)
  
  ##### null lambda ----
  for (iter in 1:100) {
    idx <- sample(length(CV_cover))
    
    avg_cover_null  <- setNames(avg_cover[idx],  names(avg_cover))
    CV_cover_null <- setNames(CV_cover[idx], names(CV_cover))
    #names(sync_cover) <- subset_data$sp1[pos]
    
    phylosig_result_avg_null <- phylosig(tree1$scenario.1,avg_cover_null,method="lambda",test=F, niter=niter)
    phylosig_result_K_avg_null <- phylosig(tree1$scenario.1,avg_cover_null,method="K",test=F, niter=niter)
    phylosig_result_I_avg_null <- Moran.I(avg_cover_null, weight_mat,scaled = TRUE)
    
    #phylosig_result_sync_null <- phylosig(tree1$scenario.1,sync_cover,method="lambda",test=F, niter=niter)
    #phylosig_result_K_sync_null <- phylosig(tree1$scenario.1,sync_cover,method="K",test=F, niter=niter)
    phylosig_result_CV_null <- phylosig(tree1$scenario.1,log(1/CV_cover_null),method="lambda",test=F, niter=niter)
    phylosig_result_K_CV_null <- phylosig(tree1$scenario.1,log(1/CV_cover_null),method="K",test=F, niter=niter)
    phylosig_result_I_CV_null <- Moran.I(log(1/CV_cover_null), weight_mat,scaled = TRUE)
    
    lambda_avg_null <- c(lambda_avg_null,phylosig_result_avg_null$lambda)
    #lambda_sync_null <- c(lambda_sync_null,phylosig_result_sync_null$lambda)
    lambda_CV_null <- c(lambda_CV_null,phylosig_result_CV_null$lambda)
    
    K_avg_null <- c(K_avg_null,phylosig_result_K_avg_null)
    #K_sync_null <- c(K_sync_null,phylosig_result_K_sync_null)
    K_CV_null <- c(K_CV_null,phylosig_result_K_CV_null)
    
    I_avg_null <- c(I_avg_null, phylosig_result_I_avg_null$observed)
    I_CV_null <- c(I_CV_null, phylosig_result_I_CV_null$observed)
    
  }
  
  lambda_avg_null_mean <- mean(lambda_avg_null)
  #lambda_sync_null_mean <- mean(lambda_sync_null)
  lambda_CV_null_mean <- mean(lambda_CV_null)

  K_avg_null_mean <- mean(K_avg_null)
  #K_sync_null_mean <- mean(K_sync_null)
  K_CV_null_mean <- mean(K_CV_null)
  
  I_avg_null_mean <- mean(I_avg_null)
  #lambda_sync_null_mean <- mean(lambda_sync_null)
  I_CV_null_mean <- mean(I_CV_null)
  
  lambda_avg_null_sd <- sd(lambda_avg_null)
  #lambda_sync_null_sd <- sd(lambda_sync_null)
  lambda_CV_null_sd <- sd(lambda_CV_null)
  
  K_avg_null_sd <- sd(K_avg_null)
  #K_sync_null_sd <- sd(K_sync_null)
  K_CV_null_sd <- sd(K_CV_null)
  
  I_avg_null_sd <- sd(I_avg_null)
  #lambda_sync_null_mean <- mean(lambda_sync_null)
  I_CV_null_sd <- sd(I_CV_null)
  

  pop_stab_df <- rbind(pop_stab_df,data.frame(plotID_full = plotID_full,
                                              site_code = site_code,
                                              lambda_stability=phylosig_result_CV$lambda,
                                              K_stability=phylosig_result_K_CV,
                                              I_stability = phylosig_result_I_CV$observed,
                                              SES_lambda_stability= (phylosig_result_CV$lambda-lambda_CV_null_mean)/lambda_CV_null_sd,
                                              SES_K_stability= (phylosig_result_K_CV-K_CV_null_mean)/K_CV_null_sd,
                                              SES_I_stability = (phylosig_result_I_CV$observed-I_CV_null_mean)/I_CV_null_sd,
                                              lambda_avg=phylosig_result_avg$lambda,
                                              K_avg=phylosig_result_K_avg,
                                              I_avg = phylosig_result_I_avg$observed,
                                              SES_lambda_avg= (phylosig_result_avg$lambda-lambda_avg_null_mean)/lambda_avg_null_sd,
                                              SES_K_avg= (phylosig_result_K_avg-K_avg_null_mean)/K_avg_null_sd,
                                              SES_I_avg = (phylosig_result_I_avg$observed-I_avg_null_mean)/I_avg_null_sd,
                                              #lambda_sync=phylosig_result_sync$lambda,
                                              #K_sync=phylosig_result_K_sync,
                                              #SES_lambda_sync= (phylosig_result_sync$lambda-lambda_sync_null_mean)/lambda_sync_null_sd,
                                              #SES_K_sync= (phylosig_result_K_sync-K_sync_null_mean)/K_sync_null_sd,
                                              gamma=gamma,
                                              max.phy.cor=max.phy.cor,
                                              mpd = mean.phy.cor,
                                              VPD = VPD,
                                              native_sim_sd = native_sim_sd,
                                              avg_abundance =  mean(subset_data$avg_cover.x),
                                              as.list(colMeans(subset_env))))
  
}

###
pop_stab_df$MAP_v2 <- log(pop_stab_df$MAP_v2)

m1_cover_I <- glmmTMB(SES_I_avg~scale(native_sim)*scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|site_code),
                           data=  pop_stab_df,
                           dispformula = ~scale(native_sim)*scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma)),
                           control = glmmTMBControl(parallel = 6),
                           family=gaussian,
                      REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_cover_I)
plot(simulateResiduals(m1_cover_I))

m1_cover_I <- glmmTMB(SES_I_avg~scale(native_sim)+scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|site_code),
                      data=  pop_stab_df,
                      dispformula = ~scale(native_sim)+scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma)),
                      control = glmmTMBControl(parallel = 6),
                      family=gaussian,
                      REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_cover_I)
p_comm_I_avg <- plot_model_predictions(m1_cover_I,data=pop_stab_df,family="gaussian",
                                       y_label = "Moran's I SES (Average cover)", 
                                       x_var = "VPD",
                                       x_label = "Phylogenetic distance variance",
                                       facet_var = "native_sim",
                                       facet_label = "Average sympatry extent",
                                       ncol= 1,
                                       nrow = 2,
                                       linetype_mean = "dashed",
                                       linetype_var  = "dashed",
                                       log_transform=F)
ggsave("Figure/Result/avg_comm.tiff", p_comm_I_avg$combined,
       dpi = 600, width = 8, height = 8, compression = "lzw")


m1_cover_I2 <- glmmTMB(SES_I_avg~scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|site_code),
                      data=  pop_stab_df,
                      dispformula = ~scale(native_sim)+scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma)),
                      control = glmmTMBControl(parallel = 6),
                      family=gaussian,
                      REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_cover_I2)

cor(predict(m1_cover_I,re.form=NA),pop_stab_df$SES_I_avg)
cor(predict(m1_cover_I2 ,re.form=NA),pop_stab_df$SES_I_avg)
# m1_cover_lambda <- glmmTMB(lambda_avg~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|site_code),
#                          data=  pop_stab_df,
#                          dispformula = ~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma)),
#                          control = glmmTMBControl(parallel = 6),
#                          family=gaussian,REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
# summary(m1_cover_lambda)
# plot(simulateResiduals(m1_cover_lambda))
# 
# p_comm_lambda_avg <- plot_model_predictions(m1_cover_lambda,data=pop_stab_df,family="gaussian",
#                                         y_label = "Pagel's λ (Average cover)", 
#                                         x_var = "native_sim",
#                                         x_label = "Historical sympatry",
#                                         ncol= 2,
#                                         nrow = 1,
#                                         linetype_mean = "dashed",
#                                         linetype_var  = "dashed",
#                                         log_transform=F)
# 
# ggsave("Figure/Result/avg_comm.tiff", p_comm_lambda_avg$combined,
#        dpi = 600, width = 8, height = 4, compression = "lzw")
# 
# m1_cover_K <- glmmTMB(K_avg~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|site_code),
#                            data=  pop_stab_df,
#                            dispformula = ~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma)),
#                            control = glmmTMBControl(parallel = 6),
#                            family=gaussian,REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
# summary(m1_cover_K)
# plot(simulateResiduals(m1_cover_K))
# performance::r2(m1_cover_K)
# 
# p_comm_K_cover <- plot_model_predictions(m1_cover_K,
#                                              data=pop_stab_df,
#                                              family="gaussian",
#                                              y_label = "Bloomberg's K (Mean cover)",
#                                              x_var = "native_sim",
#                                              x_label = "Historical sympatry",
#                                              linetype_mean = "dashed",
#                                              linetype_var  = "dashed",
#                                              ncol= 2,
#                                              nrow = 1,
#                                              dispersion_metric = "disp",
#                                              log_transform=T)

# ggsave("Figure/Result/stab_comm.tiff", p_comm_K_cover$combined,
#        dpi = 600, width = 8, height =4, compression = "lzw")


###
m1_stab_I <- glmmTMB(SES_I_stability~scale(native_sim)*scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|site_code),
                      data=  pop_stab_df,
                      dispformula = ~scale(native_sim)*scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma)),
                      control = glmmTMBControl(parallel = 6),
                      family=gaussian,REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_stab_I)
plot(simulateResiduals(m1_stab_I))

m1_stab_I <- glmmTMB(SES_I_stability~scale(native_sim)+scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|site_code),
                     data=  pop_stab_df,
                     dispformula = ~scale(native_sim)+scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma)),
                     control = glmmTMBControl(parallel = 6),
                     family=gaussian,REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_stab_I)

p_comm_I_stab <- plot_model_predictions(m1_stab_I,data=pop_stab_df,family="gaussian",
                                        y_label = "Moran's I SES (Population stability)", 
                                        x_var = "VPD",
                                        x_label = "Phylogenetic distance variance",
                                        facet_var = "native_sim",
                                        facet_label = "Average sympatry extent",
                                        ncol= 1,
                                        nrow = 2,
                                        linetype_mean = "dashed",
                                        linetype_var  = "dashed",
                                        log_transform=F)
ggsave("Figure/Result/stab_comm.tiff", p_comm_I_stab$combined,
       dpi = 600, width = 8, height = 8, compression = "lzw")

m1_stab_I2 <- glmmTMB(SES_I_stability~scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|site_code),
                     data=  pop_stab_df,
                     dispformula = ~scale(native_sim)+scale(VPD)+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma)),
                     control = glmmTMBControl(parallel = 6),
                     family=gaussian,REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
summary(m1_stab_I2)

cor(predict(m1_stab_I,re.form=NA),pop_stab_df$SES_I_stability)
cor(predict(m1_stab_I2 ,re.form=NA),pop_stab_df$SES_I_stability)
    # m1_cover_lambda <- glmmTMB(lambda_avg~scale(native_sim)+sca
# m1_stability_lambda <- glmmTMB(lambda_stability~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma))+(1|site_code),
#                            data=  pop_stab_df,
#                            dispformula = ~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(log(gamma)),
#                            control = glmmTMBControl(parallel = 6),
#                            family=gaussian,
#                            REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
# summary(m1_stability_lambda)
# plot(simulateResiduals(m1_stability_lambda))
# 
# p_comm_lambda_stab <- plot_model_predictions(m1_stability_lambda,
#                                              data=pop_stab_df,
#                                              family="gaussian",
#                                              y_label = "Pagel's λ (Population stability)",
#                                              x_var = "native_sim",
#                                              x_label = "Historical sympatry",
#                                         linetype_mean = "dashed",
#                                         linetype_var  = "dashed",
#                                         ncol= 2,
#                                         nrow = 1,
#                                         dispersion_metric = "disp",
#                                         log_transform=F)
# 
# # ggsave("Figure/Result/stab_comm.tiff", p_comm_lambda_stab$combined,
# #        dpi = 600, width = 8, height =4, compression = "lzw")
# 
# m1_stability_K <- glmmTMB(log(K_stability)~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(gamma)+(1|site_code),
#                                data=  pop_stab_df,
#                                dispformula = ~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(gamma),
#                                control = glmmTMBControl(parallel = 6),
#                                family=gaussian, REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
# summary(m1_stability_K)
# plot(simulateResiduals(m1_stability_K))
# 
# 
# p_comm_K_stab <- plot_model_predictions(m1_stability_K,
#                                              data=pop_stab_df,
#                                              family="gaussian",
#                                              y_label = "Bloomberg's K (Population stability)",
#                                              x_var = "native_sim",
#                                              x_label = "Historical sympatry",
#                                              linetype_mean = "dashed",
#                                              linetype_var  = "dashed",
#                                              ncol= 2,
#                                              nrow = 1,
#                                              dispersion_metric = "disp",
#                                              log_transform=T)
# 
# ggsave("Figure/Result/stab_comm.tiff", p_comm_K_stab$combined,
#        dpi = 600, width = 8, height =4, compression = "lzw")

# ###
# m1_sync_lambda <- glmmTMB(lambda_sync~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(gamma)+(1|site_code),
#                           data=  pop_stab_df,
#                           dispformula = ~scale(native_sim)+scale(sqrt(mpd))+scale(MAP_VAR_v2)+scale(MAT_v2)+scale(TEMP_VAR_v2)+scale(MAP_v2)+scale(year)+scale(gamma),
#                           control = glmmTMBControl(parallel = 6),
#                           family=Gamma(link="log"),
#                           REML=T) #don't do phi because a pair between an abundant and a rare species will have high value, even though they are asynchronized.
# summary(m1_sync_lambda)
# plot(simulateResiduals(m1_sync_lambda))
