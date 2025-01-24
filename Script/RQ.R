all_SES_pair_long$plotID_full <- paste0(all_SES_pair_long$site_code,"_",all_SES_pair_long$block,"_",all_SES_pair_long$plot,"_",all_SES_pair_long$subplot)
all_SES_pair_long$tl_dist <- all_SES_pair_long$repl+all_SES_pair_long$rich
all_SES_pair_long$pairs <- paste0(all_SES_pair_long$local_provenance.x,"-",all_SES_pair_long$local_provenance.y)
all_SES_pair_long$pairs[all_SES_pair_long$pairs == "INT-NAT"] <- "NAT-INT"

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
  
  if(nrow(subset_all_SES_pair) <= 30) {
    next
  }
  
  m <- rq(SES.pair~scale(phy.cor)*scale(sqrt(repl)),data=subset_all_SES_pair,tau=c(0.2,0.5,0.8))
  summary(m,se="boot")
  AIC(m)
  m1 <- rq(SES.pair~scale(phy.cor)+scale(sqrt(repl)),data=subset_all_SES_pair,tau=c(0.2,0.5,0.8))
  summary(m1,se="boot")
  AIC(m1)
  
  max_repl <- max(subset_all_SES_pair$repl,na.rm=T)
  min_repl <- min(subset_all_SES_pair$repl,na.rm=T)
  min_phy <- min(subset_all_SES_pair$phy.cor)
  max_phy <- max(subset_all_SES_pair$phy.cor)
  
  AIC_diff <- AIC(m)-AIC(m1)
  rr_cor <- cor(subset_all_SES_pair$rich,subset_all_SES_pair$repl)
  repl_phy_cor <- cor(subset_all_SES_pair$repl,subset_all_SES_pair$phy.cor)
  
  result_df <- rbind(result_df,
                     data.frame(var=rownames(summary(m)[[1]]$coefficients),summary(m,se="boot")[[1]]$coefficients,tau=summary(m)[[1]]$tau,rr_cor,repl_phy_cor,AIC_diff=AIC_diff[[1]],max_repl=max_repl,min_repl=min_repl,min_phy,max_phy,l=l),
                     data.frame(var=rownames(summary(m)[[2]]$coefficients),summary(m,se="boot")[[2]]$coefficients,tau=summary(m)[[2]]$tau,rr_cor,repl_phy_cor,AIC_diff=AIC_diff[[2]],max_repl=max_repl,min_repl=min_repl,min_phy,max_phy,l=l),
                     data.frame(var=rownames(summary(m)[[3]]$coefficients),summary(m,se="boot")[[3]]$coefficients,tau=summary(m)[[3]]$tau,rr_cor,repl_phy_cor,AIC_diff=AIC_diff[[3]],max_repl=max_repl,min_repl=min_repl,min_phy,max_phy,l=l))
  
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

#result_df$sig <- ifelse(sign(result_df$lower.bd)*sign(result_df$upper.bd) == 1, 1,0)
result_df$repl_range <- result_df$max_repl-result_df$min_repl
result_df$phy_range <- result_df$max_phy-result_df$min_phy

library(tidyverse)
result_df_overall <- result_df %>%
  group_by(tau,l) %>%
  dplyr::summarise(AIC=mean(AIC_diff))

t.test(subset(result_df,var=="scale(phy.cor)" & tau == 0.8)$Value)
t.test(subset(result_df,var=="scale(phy.cor):scale(repl)" & tau == 0.8)$Value)

summary(lm(Value~1,weight=1/(Std..Error^2),data=subset(result_df,var=="scale(phy.cor):scale(sqrt(repl))" & tau == 0.8)))

summary(lm(Value~1,weight=1/(Std..Error^2),data=subset(result_df,var=="scale(phy.cor)" & tau == 0.2)))
summary(lm(Value~1,weight=1/(Std..Error^2),data=subset(result_df,var=="scale(sqrt(repl))" & tau == 0.2)))

###
summary(lm(Value~1,weight=1/(Std..Error^2),data=subset(result_df,var=="scale(phy.cor):scale(repl)" & tau == 0.5)))

library(mgcv)
summary(gam(Value~s(phy_range,),weight=1/(Std..Error^2),data=subset(result_df,var=="scale(phy.cor)" & tau == 0.8)))

        