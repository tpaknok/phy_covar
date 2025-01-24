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

Nutnet_year <- Nutnet %>%
  group_by(site_name,site_code,block,plot) %>%
  summarize(year_range = max(year)-min(year))

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
  filter(year_trt >= 1) %>%
  summarize(year_unique = length(unique(year))) %>%
  filter(year_unique >= 10) %>%
  mutate(ID = paste0(site_code,"_",block,"_",plot,"_",subplot)) %>%
  ungroup() %>%
  select(ID)

NAT_Status <- Cover %>%
  group_by(functional_group,local_provenance) %>%
  summarize(n())

Cover_subset <- Cover %>%
  mutate(ID = paste0(site_code,"_",block,"_",plot,"_",subplot)) %>%
  filter(ID %in% Cover_year$ID) %>%
  filter(live == 1) %>%
  filter(Family != "NULL") %>%
  mutate(Genus = word(Taxon,sep=" ",1)) %>%
  filter(Genus != "UNKNOWN") %>%
  #filter(trt == "Control") %>%
  filter(trt != "Fence" & trt != "NPK+Fence" & trt != "Other") %>%
  filter(functional_group != "BRYOPHYTE") %>%
  filter(functional_group != "NULL") %>%
  filter(year_trt >= 1) #%>%
  #filter(local_provenance == "NAT" | local_provenance == "INT")

Biomass_subset <- Biomass %>%
  mutate(ID = paste0(site_code,"_",block,"_",plot,"_",subplot)) %>%
  filter(ID %in% Cover_year$ID) %>%
  filter(live == 1) %>%
  filter(trt != "Fence" & trt != "NPK+Fence" & trt != "Other") %>%
  filter(category != "BRYOPHYTE") %>%
  filter(category != "NULL") %>%
  filter(year_trt >= 1) %>%
  group_by(ID,year) %>%
  summarize(sum_biomass = sum(mass))

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
tree1 <- V.PhyloMaker2::phylo.maker(all_sp2,GBOTB.extended.WP,nodes.info.1.WP)

phylo_pair <- cophenetic(tree1$scenario.3)

fixed_name <- clean_taxon[clean_taxon$Submitted_Name != clean_taxon$Accepted_SPNAME,]
Cover_subset <- Cover_subset %>%
  mutate(Taxon = gsub("_"," ",Cover_subset$Taxon)) %>%
  left_join(clean_taxon[,c("Submitted_Name","Accepted_SPNAME")],by=join_by(Taxon == Submitted_Name))

Cover_subset$Accepted_SPNAME <- gsub(" ","_",Cover_subset$Accepted_SPNAME)

all_cor_df <- all_sp_df <- NULL

mpd_zero <- function(x,phylo_pair_subset,abundance.weighted = T) {
  mpd <- mpd(x,phylo_pair_subset,abundance.weighted)
  mpd <- mpd %>% replace(is.na(.),0)
  mean(mpd)
}

mntd_zero <- function(x,phylo_pair_subset,abundance.weighted = T) {
  mntd <- mntd(x,phylo_pair_subset,abundance.weighted)
  mntd <- mntd %>% replace(is.na(.),0)
  mean(mntd)
}

for (i in 1:length(unique(Cover_subset$ID))) {
  message(i)
  
  Cover_df_subset <- Cover_subset %>%
    ungroup() %>%
    filter(ID == unique(Cover_subset$ID)[[i]])
  
  Biomass_df_subset <- Biomass_subset %>%
    ungroup() %>%
    filter(ID == unique(Cover_subset$ID)[[i]])
  
  comp_df <- Cover_df_subset %>%
    select(year,Accepted_SPNAME,max_cover) %>%
    group_by(year,Accepted_SPNAME) %>%
    summarize(max_cover = sum(max_cover)) %>%
    pivot_wider(id_cols=year,
                names_from= Accepted_SPNAME,
                values_from = max_cover,
                values_fill = 0)
  
  comp_df_biomass <- comp_df
  comp_df_biomass[,-1] <- decostand(comp_df[,-1],"total")
  comp_df_biomass[,-1] <- comp_df_biomass[,-1] * Biomass_df_subset$sum_biomass
  
  sp_sd <- apply(comp_df,2,sd)
  
  #comp_df <- comp_df[,sp_sd > 0]
  
  if (ncol(comp_df) <= 2) next()
  
  phylo_pair_subset <- phylo_pair[rownames(phylo_pair) %in% colnames(comp_df[,-1]),colnames(phylo_pair) %in% colnames(comp_df[,-1])]
  
  mpd_gamma <- mpd(t(as.matrix(colSums(comp_df[,-1]))),phylo_pair_subset,abundance.weighted = T)
  mpd_uw_gamma <- mpd(t(as.matrix(colSums(comp_df[,-1]))),phylo_pair_subset,abundance.weighted = F)
  mntd_gamma <- mntd(t(as.matrix(colSums(comp_df[,-1]))),phylo_pair_subset,abundance.weighted = T)
  mntd_uw_gamma <- mntd(t(as.matrix(colSums(comp_df[,-1]))),phylo_pair_subset,abundance.weighted = F)
  
  mpd_alpha <- mean(mpd_zero(as.matrix(comp_df[,-1]),phylo_pair_subset,abundance.weighted = T))
  mpd_uw_alpha <- mean(mpd_zero(as.matrix(comp_df[,-1]),phylo_pair_subset,abundance.weighted = F))
  mntd_alpha <- mean(mntd_zero(as.matrix(comp_df[,-1]),phylo_pair_subset,abundance.weighted = T))
  mntd_uw_alpha <- mean(mntd_zero(as.matrix(comp_df[,-1]),phylo_pair_subset,abundance.weighted = F))
  
  #phylo_redun <- treeUniqueness(tree_Nutnet,t(as.matrix(colSums(comp_df[,-1]))),index="Shannon")
  gamma <- diversity(t(as.matrix(colSums(comp_df[,-1]))))
  mean_alpha  <- mean(diversity(as.matrix(comp_df[,-1])))
  mean_alpha_uw<- mean(specnumber(as.matrix(comp_df[,-1])))
  pielou <- mean(diversity(as.matrix(comp_df[,-1])))/log(mean(specnumber(as.matrix(comp_df[,-1]))))
  
  sd_test <- ifelse(all(apply(comp_df[,-1],2,sd) > 0),1,0)
  Faith_wPD <- weighted.faith(tree1$scenario.3,as.data.frame(t(colSums(as.matrix(comp_df[,-1])))))

  COV.AB.df <- SES.COV.AB(as.matrix(comp_df[,-1]),comp_df[,1])$SES.analysis
  COV.AB.df2 <- SES.COV.AB(as.matrix(comp_df_biomass[,-1]),comp_df_biomass[,1])$SES.analysis
  
  bioM_variability <- sd(Biomass_df_subset$sum_biomass)/mean(Biomass_df_subset$sum_biomass)
    
  sp_df_subset <- Cover_df_subset %>%
    distinct(Family,Accepted_SPNAME,local_provenance) %>%
    mutate(Genus = word(Accepted_SPNAME,sep="_",1)) %>%
    add_count(Accepted_SPNAME) %>%
    mutate(i = 1)
  
  nat_test <- ifelse(max(sp_df_subset$n) == 1, 1,0)
  
  nat_df_subset <- sp_df_subset %>%
    mutate(local_provenance = ifelse(n==1,local_provenance,"UNK")) %>%
    ungroup() %>%
    distinct(Accepted_SPNAME,local_provenance)
  
  native_sp <- nat_df_subset[nat_df_subset$local_provenance == "NAT","Accepted_SPNAME"]
  introduced_sp <- nat_df_subset[nat_df_subset$local_provenance == "INT","Accepted_SPNAME"]
  
  native_comp <- comp_df[,colnames(comp_df) %in% native_sp]
  int_comp <- comp_df[,colnames(comp_df) %in% introduced_sp]
  
  nat_ratio_alpha <- mean(specnumber(native_comp)/specnumber(comp_df[,-1]))
  nat_ratio_abundance_alpha <- mean(rowSums(native_comp)/rowSums(comp_df[,-1]))
  nat_ratio_gamma <- length(native_sp)/length(colnames(comp_df[,-1]))
  nat_ratio_abundance_gamma <- sum(native_comp)/sum(comp_df[,-1])
  unclassified_ratio <- 1-(length(native_sp)+length(introduced_sp))/length(colnames(comp_df[,-1]))
  unclassified_ratio_abundance <- 1-(sum(native_comp)+sum(int_comp))/sum(comp_df[,-1])
  
  temporal_df_site_subplot <- data.frame(ID=unique(Cover_subset$ID)[[i]],
                                         site.code = unique(Cover_df_subset$site_code),
                                         block = unique(Cover_df_subset$block),
                                         plot = unique(Cover_df_subset$plot),
                                         subplot = unique(Cover_df_subset$subplot),
                                         trt = unique(Cover_df_subset$trt),
                                         year = max(comp_df[,1])-min(comp_df[,1]),
                                         year_sampled = nrow(comp_df),
                                         sd_test = sd_test,
                                         nat_test = nat_test,
                                         i=i,
                                         spp = ncol(comp_df[,-1]),
                                         abundance = mean(rowSums(comp_df[,-1])),
                                         gamma=gamma,
                                         mpd_gamma=mpd_gamma,
                                         mpd_uw_gamma = mpd_uw_gamma,
                                         mntd_gamma=mntd_gamma,
                                         mntd_uw_gamma = mntd_uw_gamma,
                                         mpd_alpha=mpd_alpha,
                                         mpd_uw_alpha = mpd_uw_alpha,
                                         mntd_alpha=mntd_alpha,
                                         mntd_uw_alpha = mntd_uw_alpha,
                                         #Accepted_SPNAME = sort(colnames(comp_df[,-1])),
                                         #phylo_redun = phylo_redun$`Phylogenetic Redundancy`,
                                         Faith_wPD = Faith_wPD,
                                         mean_alpha = mean_alpha,
                                         mean_alpha_uw = mean_alpha_uw,
                                         unclassified_ratio,
                                         unclassified_ratio_abundance,
                                         nat_ratio_alpha,
                                         nat_ratio_abundance_alpha,
                                         nat_ratio_gamma,
                                         nat_ratio_abundance_gamma,
                                         pielou = pielou,
                                         bioM_variability
  )
  
  #temporal_df_site_subplot <- temporal_df_site_subplot %>%
    #left_join(nat_df_subset,by="Accepted_SPNAME")
  
  #SES.coef.ab <- COV.AB.df$SES.coef.ab
  #SES.cov <- COV.AB.df$SES.cov
  #SES.cor <- COV.AB.df$SES.cor
  
  #names(SES.coef.ab) <- names(SES.cov) <- rownames(COV.AB.df)

  # coef.phylo.sig <- phylosig(tree_Nutnet,SES.coef.ab,method="lambda")
  # cov.phylo.sig <- phylosig(tree_Nutnet,SES.cov,method="lambda")
  
  # lambda_result_cov <- geiger::fitContinuous(tree_Nutnet,dat=SES.cov,model=c("lambda"))
  # OU_result_cov <-geiger::fitContinuous(tree_Nutnet,dat=SES.cov,model=c("OU"))
  # lambda_result_ab <- geiger::fitContinuous(tree_Nutnet,dat=SES.coef.ab,model=c("lambda"))
  # OU_result_ab <-geiger::fitContinuous(tree_Nutnet,dat=SES.coef.ab,model=c("OU"))
  # white_result_ab <- geiger::fitContinuous(tree_Nutnet,dat=SES.coef.ab,model=c("white"))
  # white_result_cov <-geiger::fitContinuous(tree_Nutnet,dat=SES.cov,model=c("white"))
  
  temporal_df_site_subplot <- cbind(temporal_df_site_subplot,unique(COV.AB.df[,c("comm.var","w.sp.var","sp.var","sync_w_r","sync_uw_r")]))
  # temporal_df_site_subplot <- cbind(temporal_df_site_subplot,COV.AB.df,
  #                                   lambda_ab=round(lambda_result_ab$opt$lambda,5),
  #                                   OU_ab=round(OU_result_ab$opt$alpha,5),
  #                                   lambda_ab_AICC=round(lambda_result_ab$opt$aicc,5),
  #                                   OU_ab_AICC=round(OU_result_ab$opt$aicc,5),
  #                                   white_ab_AICC = round(white_result_ab$opt$aicc,5),
  #                                   lambda_cov=round(lambda_result_cov$opt$lambda,5),
  #                                   OU_cov=round(OU_result_cov$opt$alpha,5),
  #                                   lambda_cov_AICC=round(lambda_result_cov$opt$aicc,5),
  #                                   OU_cov_AICC=round(OU_result_cov$opt$aicc,5),
  #                                   white_cov_AICC = round(white_result_cov$opt$aicc,5))
  # 
  
  biomass_metric <- data.frame(unique(COV.AB.df2[,c("comm.var","w.sp.var","sp.var","sync_w_r","sync_uw_r")]))
  colnames(biomass_metric) <- c("b.comm.var","b.w.sp.var","b.sp.var","b_sync_w_r","b_sync_uw_r")
  temporal_df_site_subplot <- cbind(temporal_df_site_subplot,biomass_metric)
  all_cor_df <- rbind(all_cor_df,temporal_df_site_subplot)
  
  all_sp_df <- rbind(all_sp_df,sp_df_subset)
}

### global

pooled_df <- all_cor_df %>%
  group_by(Taxon) %>%
  summarize(mean.SES.cov = mean(SES.cov,na.rm=T),
            mean.SES.coef.ab = mean(SES.coef.ab,na.rm=T))

global.SES.cov <- pooled_df$mean.SES.cov
global.SES.coef.ab <- pooled_df$mean.SES.coef.ab
names(global.SES.coef.ab) <- names(global.SES.cov) <- pooled_df$Taxon

coef.phylo.sig <- phylosig(tree_Nutnet,global.SES.coef.ab,method="lambda")
cov.phylo.sig <- phylosig(tree_Nutnet,global.SES.cov,method="lambda")

########
lambda_df <- all_cor_df %>%
  group_by(ID) %>%
  mutate(NAT.SR = sum(local_provenance == "NAT"),
         INT.SR = sum(local_provenance == "INT"),
         ratio = NAT.SR/spp,
         unclassified.ratio = 1-(NAT.SR+INT.SR)/spp) %>%
  group_by(ID,site.code,block,plot,subplot,year,spp,ratio,unclassified.ratio,mpd,gamma,phylo_redun,NAT.SR,INT.SR) %>%
  summarize(AB.phylo = mean(`coef.phylo.sig$lambda`),
            cov.phylo = mean(`cov.phylo.sig$lambda`)) %>%
  filter(spp >=20  & unclassified.ratio <= 0.1) 

plot(lambda_df$cov.phylo,lambda_df$AB.phylo)

lambda_df <- lambda_df %>%
  left_join(distinct(Env[,c("site_code","AI","MAP_v2","MAT_v2","TEMP_VAR_v2","MAP_VAR_v2")]),by=join_by(site.code == site_code))
###
library(DHARMa)
lambda_df$cov.phylo1 <- ifelse(lambda_df$cov.phylo < 0.01,0,lambda_df$cov.phylo)
lambda_df$cov.phylo1 <- ifelse(lambda_df$cov.phylo1 > 1,1,lambda_df$cov.phylo1)

m_phylo_sig <- glmmTMB(cov.phylo1~scale(mpd)*scale(MAT_v2)*scale(ratio)+(1|site.code),
                       family=tweedie,
                       data=lambda_df,
                       control = glmmTMBControl(parallel = 6))
plot(simulateResiduals(m_phylo_sig))
summary(m_phylo_sig)
performance::r2(m_phylo_sig)

library(ggplot2)
library(ggeffects)
predict_df <- as.data.frame(ggemmeans(m_phylo_sig,terms=c("ratio[0.1:1,by=0.05]","MAT_v2[0:19,by=1]")),terms_to_colnames=T)
predict_df$MAT_v2 <- as.numeric(as.character(predict_df$MAT_v2))
p <- ggplot(data=predict_df)+
  geom_point(data=lambda_df,aes(y=cov.phylo1,x=ratio,fill=MAT_v2),shape=21,size=12)+
  geom_line(aes(y=predicted,x=ratio,colour=MAT_v2,group=MAT_v2))+
  scale_fill_continuous(low="blue",high="red")+
  scale_colour_continuous(low="blue",high="red")

plot(p)

p2 <- ggplot(data=predict_df)+
  geom_raster(aes(y=MAT_v2,x=ratio,fill=predicted))+
  geom_point(data=lambda_df,aes(y=MAT_v2,x=ratio,fill=cov.phylo1),shape=21,size=6)+
  scale_fill_continuous(low="blue",high="red",limits=c(0,1))
  
plot(p2)

lambda_df$AB.phylo1 <- ifelse(lambda_df$AB.phylo < 0.01,0,lambda_df$AB.phylo)
m_phylo_sig1 <- glmmTMB(AB.phylo1~scale(spp)*scale(MAT_v2)*scale(ratio)+(1|site.code),
                        data=lambda_df,family=tweedie)
plot(simulateResiduals(m_phylo_sig1))
summary(m_phylo_sig1)
performance::r2(m_phylo_sig1)


p2 <- ggplot()+
  geom_violin(data=lambda_df,aes(y=AB.phylo1,x=1),position=position_dodge())+
  scale_fill_continuous(low="blue",high="red",limits=c(0,1))

plot(p2)
### species-level analyses

all_cor_df <- all_cor_df %>%
  left_join(distinct(Env[,c("site_code","AI","MAP_v2","MAT_v2","TEMP_VAR_v2","MAP_VAR_v2")]),by=join_by(site.code == site_code))

m <- lmer(SES.coef.ab~local_provenance*scale(MAT_v2)*scale(mpd)+(1|Taxon)+(local_provenance||site.code)+(local_provenance||plot),data=subset(all_cor_df,local_provenance != "UNK" & local_provenance != "NULL"),REML=T)
summary(m)
anova(m)

m1 <- lmer(SES.cov~local_provenance*scale(MAT_v2)*scale(mpd)+(1|Taxon)+(local_provenance||site.code)+(1|plot),data=subset(all_cor_df,local_provenance != "UNK" & local_provenance != "NULL"))
summary(m1)
performance::r2(m1)

predict_df <- as.data.frame(ggemmeans(m1,terms=c("MAT_v2","local_provenance")),terms_to_colnames=T)
predict_df$MAT_v2 <- as.numeric(as.character(predict_df$MAT_v2))
p3 <- ggplot(data=predict_df)+
  geom_line(aes(y=predicted,x=local_provenance,colour=MAT_v2,group=MAT_v2))+
  scale_fill_continuous(low="blue",high="red")+
  scale_colour_continuous(low="blue",high="red")

plot(p3)

library(interactions)
johnson_neyman(m1,MAT_2,local_provenance)

### community level analyses
library(lme4)
library(lmerTest)
library(ggeffects)
all_cor_df_pooled <- all_cor_df %>%
  group_by(ID) %>%
  mutate(NAT.SR = sum(local_provenance == "NAT"),
         INT.SR = sum(local_provenance == "INT"),
         ratio = NAT.SR/spp,
         unclassified.ratio = 1-(NAT.SR+INT.SR)/spp) %>%
  group_by(ID,site.code,block,plot,subplot,year,spp,mpd,mpd_uw,gamma,phylo_redun,NAT.SR,INT.SR,unclassified.ratio,ratio,Faith_wPD,mean_alpha,mpd_alpha_uw,mpd_alpha_w) %>%
  summarize(SES.coef.ab = mean(SES.coef.ab,na.rm=T),
            SES.cor = mean(SES.cor,na.rm=T),
            SES.cov = mean(SES.cov, na.rm=T)) %>%
  left_join(distinct(Env[,c("site_code","AI","MAP_v2","MAT_v2","TEMP_VAR_v2","MAP_VAR_v2","PET","elevation","latitude","longitude")]),by=join_by(site.code == site_code)) %>%
  filter(unclassified.ratio < 0.1)

all_cor_df_pooled$ID2 <- paste0(all_cor_df_pooled$site.code,"_",all_cor_df_pooled$block,"_",all_cor_df_pooled$plot)
Env$ID2 <- paste0(Env$site_code,"_",Env$block,"_",Env$plot)

all_cor_df_pooled <- all_cor_df_pooled %>%
  left_join(distinct(Env[,c("ID2","pct_N","ppm_K","ppm_P","trt")]),by="ID2")

all_cor_df_pooled$pct_N <- as.numeric(all_cor_df_pooled$pct_N)
all_cor_df_pooled$ppm_P <- as.numeric(all_cor_df_pooled$ppm_P)
all_cor_df_pooled$ppm_K <- as.numeric(all_cor_df_pooled$ppm_K)

m2 <- lmer(SES.coef.ab~scale(mpd)*scale(ratio)*scale(MAT_v2)+(1|site.code),data=all_cor_df_pooled,REML=T)
summary(m2)
tab_model(m2,vcov.fun = "CR2")
performance::r2(m2)
plot(simulateResiduals(m2))

m2 <- lmer(SES.coef.ab.mean~scale(mpd)*scale(MAT_v2)+(1|site.code),data=all_cor_df_pooled,REML=T)
summary(m2)
performance::r2(m2)

predict_df <- as.data.frame(ggemmeans(m2,terms=c("mpd","MAT_v2[0:20,by=1]")),terms_to_colnames=T)
predict_df$MAT_v2 <- as.numeric(as.character(predict_df$MAT_v2))

p2 <- ggplot(data=predict_df)+
  geom_raster(aes(y=mpd,x=MAT_v2,fill=predicted))+
  labs(fill="SES.coef.ab.mean")+
  geom_point(data=all_cor_df_pooled,aes(y=mpd,x=MAT_v2,fill=SES.coef.ab.mean),shape=21,size=6)+
  scale_fill_gradient2(low="blue",mid="white",high="red")

plot(p2)

library(robustlmm)
library(DHARMa)

m1 <- rlmer(mpd~trt+(1|site.code),
            data=all_cor_df_pooled,REML=T)

summary(m1)

p_mpd <- ggplot(all_cor_df_pooled,aes(x=trt,y=mpd))+
  geom_violin()+
  geom_point()
plot(p_mpd)

m2 <- rlmer(ratio~trt+(1|site.code),
            data=all_cor_df_pooled,REML=T)

summary(m2)

p_ratio <- ggplot(all_cor_df_pooled,aes(x=trt,y=ratio))+
  geom_violin()+
  geom_point()
plot(p_ratio)

m3 <- rlmer(SES.cor~scale(ratio)*scale(log(mpd))*trt+(1|site.code),
            data=all_cor_df_pooled,REML=T)

summary(m3)

all_predict_df <- NULL
for (i in 1:length(unique(all_cor_df_pooled$trt))) {
  trt_data <- subset(all_cor_df_pooled,trt==unique(all_cor_df_pooled$trt)[[i]])
  m3 <- rlmer(SES.cov~ratio*log(mpd)+(1|site.code),
              data=trt_data,REML=T)
  
  summary(m3)
  r2 <- performance::r2(m3)
  
  predict_df <- as.data.frame(ggemmeans(m3,terms=c("ratio[0:1,by=0.5]","mpd[20:220]")),terms_to_colnames=T)
  predict_df$mpd <- as.numeric(as.character(predict_df$mpd))
  predict_df$trt <- unique(all_cor_df_pooled$trt)[[i]]
  predict_df$r2m <- c(r2$R2_marginal)
  all_predict_df <- rbind(all_predict_df,predict_df)
}

all_predict_df$trt <- factor(all_predict_df$trt,levels=c("Control","N","P","K","NP","NK","PK","NPK"))
all_cor_df_pooled$trt <- factor(all_cor_df_pooled$trt,levels=c("Control","N","P","K","NP","NK","PK","NPK"))

all_predict_df$label <- paste0("R2m=",round(all_predict_df$r2m,3))

p1 <- ggplot(data=all_predict_df)+
  geom_line(aes(y=predicted,x=mpd,group=ratio,colour=ratio))+
  #geom_point(data=all_cor_df_pooled,aes(y=SES.cor,x=mpd,colour=ratio))+
  facet_grid(~trt)+
  geom_text(x=Inf,y=Inf,aes(label=all_predict_df$label),vjust=2,hjust=1)+
  labs(y="SES.cor", x = "mpd",colour="native species ratio")+
  scale_colour_continuous(low="blue",high="red")+
  theme_bw()


plot(p1)
ggsave("model_prediction.tiff",width=18,height=6,compression="lzw")

p1 <- ggplot(data=all_predict_df)+
  geom_line(aes(y=predicted,x=mpd,group=ratio,colour=ratio))+
  geom_point(data=all_cor_df_pooled,aes(y=SES.cor,x=mpd,colour=ratio))+
  geom_text(x=Inf,y=Inf,aes(label=all_predict_df$label),vjust=2,hjust=1)+
  facet_grid(~trt)+
  labs(y="SES.cor", x = "mpd",colour="native species ratio")+
  scale_colour_continuous(low="blue",high="red")+
  theme_bw()


plot(p1)

ggsave("with_data.tiff",width=18,height=6,compression="lzw")
predict_df <- as.data.frame(ggemmeans(m3,terms=c("ratio[0:1,by=0.5]","mpd[20:220,by=20]","MAP_VAR_v2[15:95,by=5]")),terms_to_colnames=T)
predict_df$mpd <- as.numeric(as.character(predict_df$mpd))
predict_df$MAP_VAR_v2 <- as.numeric(as.character(predict_df$MAP_VAR_v2))

p1 <- ggplot(data=predict_df)+
  geom_line(aes(y=predicted,x=mpd,group=ratio,colour=ratio))+
  facet_wrap(~MAP_VAR_v2)+
  labs(y="SES.cor.mean", x = "mpd",colour="ratio")+
  scale_colour_continuous(low="blue",high="red")


plot(p1)

p1 <- ggplot(data=predict_df)+
  geom_raster(aes(x=MAP_VAR_v2,y=ratio))+
  geom_point(data=all_cor_df_pooled,aes(x=MAP_VAR_v2,y=ratio))+
  scale_colour_continuous(low="blue",high="red")


plot(p1)

predict_df <- as.data.frame(ggemmeans(m3,terms=c("ratio[0:1,by=0.01]","ppm_K[0:1200,by=50]")),terms_to_colnames=T)
predict_df$ppm_K <- as.numeric(as.character(predict_df$ppm_K))
p1 <- ggplot(data=predict_df)+
  geom_raster(data=predict_df,aes(x=ratio,y=ppm_K,fill=predicted))+
  geom_point(data=na.omit(all_cor_df_pooled),aes(x=ratio,y=ppm_K,colour=SES.cor))+
  scale_colour_continuous(low="blue",high="red")


plot(p1)

p2 <- ggplot(data=predict_df)+
  geom_raster(aes(y=mpd,x=ratio,fill=predicted))+
  labs(fill="SES.cov.mean", x = "Native species ratio")+
  ylim(20,225)+
  geom_point(data=all_cor_df_pooled,aes(y=mpd,x=ratio,fill=SES.cov.mean),shape=21,size=6)+
  scale_fill_gradient2(low="blue",mid="white",high="red",limits=c(-5.5,1.5))

plot(p2)

predict_df <- as.data.frame(ggemmeans(m3,terms=c("MAT_v2[0:20,by=1]","mpd[20:225,by=3]")),terms_to_colnames=T)
predict_df$mpd <- as.numeric(as.character(predict_df$mpd))

p2 <- ggplot(data=predict_df)+
  geom_raster(aes(y=mpd,x=MAT_v2,fill=predicted))+
  labs(fill="SES.cov.mean", x = "MAT")+
  ylim(20,225)+
  geom_point(data=all_cor_df_pooled,aes(y=mpd,x=MAT_v2,fill=SES.cov.mean),shape=21,size=6)+
  scale_fill_gradient2(low="blue",mid="white",high="red",limits=c(-2.5,1.5))

plot(p2)
###

all_cor_df <- all_sp_df <- NULL

phylo_pair <- cophenetic(tree_Nutnet)

  for (i in 1:length(unique(Cover_subset$ID))) {
    message(i)

    Cover_df_subset <- Cover_subset %>%
      ungroup() %>%
      filter(ID == unique(Cover_subset$ID)[[i]])

    comp_df <- Cover_df_subset %>%
      select(Taxon,max_cover,year) %>%
      group_by(year,Taxon) %>%
      summarize(max_cover = sum(max_cover)) %>%
      pivot_wider(id_cols=year,
                  names_from= Taxon,
                  values_from = max_cover,
                  values_fill = 0)
    
    cover_sum <- rowSums(comp_df[,-1])
    TD <- diversity(comp_df[,-1])

    comp_df_pa <- decostand(comp_df,"pa")
    comp_df <- comp_df[,colSums(comp_df_pa) > 2]

    if (ncol(comp_df) <= 2) next()
    sd_test <- ifelse(all(apply(comp_df[,-1],2,sd) > 0),1,0)

    cor_df <- cor(comp_df[,-1])
    inds <- which(lower.tri(cor_df),arr.ind=T)

    sp_df_subset <- Cover_df_subset %>%
      distinct(Family,Taxon,local_provenance) %>%
      mutate(Genus = word(Taxon,sep=" ",1)) %>%
      add_count(Taxon) %>%
      mutate(i = 1)

    nat_test <- ifelse(max(sp_df_subset$n) == 1, 1,0)

    nat_df_subset <- sp_df_subset %>%
      mutate(local_provenance = ifelse(n==1,local_provenance,"UNK")) %>%
      ungroup() %>%
      distinct(Taxon,local_provenance)

    temporal_df_site_subplot <- data.frame(ID=unique(Cover_subset$ID)[[i]],
                                           site.code = unique(Cover_df_subset$site_code),
                                           block = unique(Cover_df_subset$block),
                                           plot = unique(Cover_df_subset$plot),
                                           subplot = unique(Cover_df_subset$subplot),
                                           trt = unique(Cover_df_subset$trt),
                                           year = max(comp_df[,1])-min(comp_df[,1]),
                                           year_sampled = nrow(comp_df),
                                           r=cor_df[inds],
                                           sp1 = rownames(cor_df)[inds[, 1]],
                                           sp2 = colnames(cor_df)[inds[,2]],
                                           sd_test = sd_test,
                                           nat_test = nat_test,
                                           i=i,
                                           spp=ncol(comp_df[,-1])
                                           )

    temporal_df_site_subplot <- temporal_df_site_subplot %>%
      left_join(nat_df_subset,by=join_by(sp1 == Taxon)) %>%
      left_join(nat_df_subset,by=join_by(sp2 == Taxon),suffix=c("_sp1","_sp2"))

    all_cor_df <- rbind(all_cor_df,temporal_df_site_subplot)

    all_sp_df <- rbind(all_sp_df,sp_df_subset)
  }

pooled_cor_df <- all_cor_df %>%
  group_by(ID,site.code,spp) %>%
  summarize(mean_r = mean(r,na.rm=T))

summary(glmmTMB(mean_r~spp,data= pooled_cor_df))
#####
library(rtrees)
library(ape)
library(glmmTMB)

all_sp <- all_sp_df
all_sp$Taxon <- word(all_sp$Taxon,1,2,sep=" ")
all_sp$Taxon<- gsub(" ","_",all_sp$Taxon)
all_sp$Taxon <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", all_sp$Taxon, perl=TRUE)
all_sp$Genus <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", all_sp$Genus, perl=TRUE)

all_sp <- all_sp[!duplicated(all_sp$Taxon),]
colnames(all_sp)[c(1,2,4)] <- c("family","species","genus")

tree_Nutnet <- get_tree(sp_list = all_sp$species,taxon="plant",scenario="at_basal_node")

pair_PD <- cophenetic.phylo(tree_Nutnet)
pair_PD_long <- expand.grid(1:nrow(pair_PD),1:ncol(pair_PD))
pair_PD_long$PD <- c(pair_PD)
pair_PD_long$sp1 <-  rownames(pair_PD)[pair_PD_long [, 1]]
pair_PD_long$sp2 <-  colnames(pair_PD)[pair_PD_long [, 2]]

pair_PD_long$sp.pair <- paste0(pair_PD_long$sp1,"-",pair_PD_long$sp2)

all_cor_df_lmm <- all_cor_df %>%
  mutate(pairs = paste0(local_provenance_sp1,"-",local_provenance_sp2)) %>%
  filter(!is.na(r)) %>%
  mutate(pairs = ifelse(pairs == "INT-NAT","NAT-INT",pairs)) %>%
  filter(pairs == "NAT-NAT" | pairs == "INT-INT" | pairs == "NAT-INT") %>%
  mutate(block = paste0(site.code,"_",block)) %>%
  mutate(plot = paste0(site.code,"_",block,"_",plot)) %>%
  filter(trt == "Control") %>%
  mutate(sp.pair = paste0(sp1,"-",sp2))

all_cor_df_lmm$sp.pair <- gsub(" ","_",all_cor_df_lmm$sp.pair)

Env_df <- Env %>%
  group_by(site_code,latitude,longitude,MAT_v2,MAP_v2,TEMP_VAR_v2,MAP_VAR_v2) %>%
  summarize(n=n())

all_cor_df_lmm <- all_cor_df_lmm %>%
  left_join(pair_PD_long[,c("sp.pair","PD")],by="sp.pair") %>%
  left_join(Env_df,by=join_by(site.code == site_code))

cor(all_cor_df_lmm[,c("MAT_v2","MAP_v2","TEMP_VAR_v2","MAP_VAR_v2")])

all_cor_df_lmm$pairs <- as.factor(all_cor_df_lmm$pairs)

m <- glmmTMB(r~pairs*PD*scale(MAT_v2)+(1|site.code)+(1|plot),data=all_cor_df_lmm)
summary(m)
performance::r2(m)
car::Anova(m)

PD_explore <- all_cor_df_lmm %>% group_by(site.code,pairs) %>% summarize(PD_max = max(PD),PD_min = min(PD))
library(ggeffects)
predict_df <- ggemmeans(m,terms=c("year","pairs","PD"))

r_plot <- ggplot(data=predict_df,aes(x=facet,y=predicted,colour=group,group=group))+
  geom_line()+
  facet_grid(~x)+
  theme_classic()

plot(r_plot)

plot_level_df <- all_cor_df_lmm %>% group_by(site.code,pairs) %>% summarize(PD = mean(PD),
                                                                          year=mean(year),
                                                                          mean_r = mean(r,na.rm=T))

m_plot <- glmmTMB(mean_r~pairs*scale(PD)*scale(year)+(1|site.code),data=plot_level_df)
summary(m_plot)
car::Anova(m_plot)
###
combined_study <- NULL
pairwise_PD <- cophenetic(tree_Nutnet)

for (i in 1:length(unique(Cover_subset$ID))) {
  message(i)

  Cover_df_subset <- Cover_subset %>%
    ungroup() %>%
    filter(ID == unique(Cover_subset$ID)[[i]])

  comp_df <- Cover_df_subset %>%
    select(Taxon,max_cover,year) %>%
    pivot_wider(names_from= Taxon,
                values_from = max_cover,
                values_fill = 0) %>%
    mutate(year = year-min(year))

  comp_df_pa <- decostand(comp_df,"pa")
  comp_df <- comp_df[,colSums(comp_df_pa) > 2]

  mean.cover <- colMeans(comp_df[,-1])

  comp_df[,-1] <- apply(comp_df[,-1],2,function(x) x/max(x))

  subset_pairwise_PD <-pairwise_PD[rownames(pairwise_PD) %in% colnames(comp_df),colnames(pairwise_PD) %in% colnames(comp_df)]
  diag(subset_pairwise_PD) <- NA
  mean.PD <- apply(subset_pairwise_PD,2,mean,na.rm=T)
  min.PD <- apply(subset_pairwise_PD,2,min,na.rm=T)

  sp_year_change_df <- NULL
  for (j in 2:ncol(comp_df)) {
    m <- lm(unlist(comp_df[,j])~year,data=comp_df) #have zero. can't do beta
    summary(m)

    sp_year_change <- data.frame(site.code = unique(Cover_df_subset$site_code),
                                 plot = unique(Cover_df_subset$ID),
                                 Taxon=colnames(comp_df)[[j]],
                                 sp = ncol(comp_df)-1,
                                 year = max(comp_df[,1])-min(comp_df[,1]),
                                 year_sampled = nrow(comp_df),
                                 slope = summary(m)$coefficients[2,1],
                                 se = summary(m)$coefficients[2,2],
                                 z = summary(m)$coefficients[2,3],
                                 p = summary(m)$coefficients[2,4],
                                 r2 = summary(m)$r.squared,
                                 native = Cover_df_subset[Cover_df_subset$Taxon %in% colnames(comp_df)[j],"local_provenance"][[1]]
    )
    sp_year_change_df <- rbind(sp_year_change_df,sp_year_change)
  }
  sp_year_change_df$mean.PD <- mean.PD
  sp_year_change_df$min.PD <- min.PD
  sp_year_change_df$mean.cover <- mean.cover
  combined_study <- rbind(combined_study,sp_year_change_df)
  }


library(spaMM)
combined_study$native_num <- ifelse(combined_study$native == "NAT",1,0)
m2 <- fitme(slope~native_num+mean.PD+mean.cover+(1|Taxon)+(native_num+mean.PD||site.code)+(native_num+mean.PD||plot),data=combined_study,method="REML",init=list(lambda=NaN))
summary(m2)
anova(m2)


m2 <- glmmTMB(slope~native_num+mean.PD+(1|Taxon)+(native_num||site.code)+(native_num||plot),data=combined_study)
summary(m2)
