library(lme4)
library(lmerTest)
library(ggeffects)
all_cor_df_pooled <- all_cor_df %>%
  filter(unclassified_ratio < 0.1 & unclassified_ratio_abundance < 0.1)

####

analyzed.df <- subset(all_cor_df_pooled,trt=="Control")
analyzed.df$mean_alpha <- exp(analyzed.df$mean_alpha)

analyzed.df$scaled_sync_w_r <- (analyzed.df$sync_w_r - (-1))/2
analyzed.df$scaled_sync_uw_r <- (analyzed.df$sync_uw_r - (-1))/2

###
library(mgcv)
library(gratia)
library(itsadug)
analyzed.df$site.code <- as.factor(analyzed.df$site.code)

smooth_interact <- gam(mpd_uw_alpha~ s(nat_ratio_alpha)+trt+s(site.code,bs="re"),
                       data = analyzed.df)
summary(smooth_interact)
plot_smooth(smooth_interact,view="nat_ratio_alpha",rm.ranef = TRUE)

smooth_interact <- gam(mntd_uw_alpha~ s(nat_ratio_alpha)+trt+s(site.code,bs="re"),
                       data = analyzed.df)
summary(smooth_interact)
plot_smooth(smooth_interact,view="nat_ratio_alpha",rm.ranef = TRUE)

smooth_interact <- gam(sync_w_r ~ s(mean_alpha_uw)+s(mpd_uw_alpha)+s(mntd_uw_alpha)+s(nat_ratio_alpha)+s(pielou)+s(abundance)+trt+s(site.code,bs="re"),
                        data = analyzed.df)
summary(smooth_interact)

plot_smooth(smooth_interact,view="mpd_uw_alpha",rm.ranef = TRUE)
plot_smooth(smooth_interact,view="nat_ratio_alpha",rm.ranef = TRUE)

smooth_interact <- gam(w.sp.var~ s(mean_alpha_uw)+s(mpd_uw_alpha)+s(mntd_uw_alpha)+s(nat_ratio_alpha)+s(pielou)+s(abundance)+trt+s(site.code,bs="re"),
                        data = analyzed.df)
summary(smooth_interact)

plot_smooth(smooth_interact,view="nat_ratio_alpha",rm.ranef = TRUE)
plot_smooth(smooth_interact,view="mntd_uw_alpha",rm.ranef = TRUE)

predict.gam(smooth_interact$gam)
###
library(DHARMa)
analyzed.df <- subset(all_cor_df_pooled)

analyzed.df <- na.omit(analyzed.df)
analyzed.df$mean_alpha <- exp(analyzed.df$mean_alpha)
cor(analyzed.df[,c("mean_alpha","mntd_alpha","mpd_alpha","nat_ratio_alpha","abundance","mean_alpha_uw","mpd_uw_alpha","pielou")])
cor(analyzed.df[,c("nat_ratio_alpha","mean_alpha_uw","mpd_uw_alpha","pielou")])

analyzed.df$scaled_sync_w_r <- (analyzed.df$sync_w_r - (-1))/2
analyzed.df$scaled_sync_uw_r <- (analyzed.df$sync_uw_r - (-1))/2

plot(analyzed.df$mean_alpha,analyzed.df$mpd_alpha)
###
m1 <- lmer(comm.var~scale(mean_alpha_uw)*scale(mpd_uw_alpha)*trt+(1|site.code),
           data=analyzed.df)
summary(m1)

m1 <- lmer(comm.var~sync_w_r+w.sp.var+(1|site.code),
              data=analyzed.df)
summary(m1)
car::Anova(m1)
performance::r2(m1)
plot(simulateResiduals(m1))

m2 <- lmer(comm.var~scale(mean_alpha)+scale(mntd_alpha)+scale(nat_ratio_alpha)+scale(log(abundance))+(1|site.code),
              data=analyzed.df)
summary(m2)
performance::r2(m2)
plot(simulateResiduals(m2))

m3 <- lmer(sync_w_r~scale(log(mean_alpha))+scale(log(mntd_alpha))+scale(nat_ratio_alpha)+scale(log(abundance))+(1|site.code),
              data=analyzed.df)
summary(m3)
performance::r2(m3)
plot(simulateResiduals(m3))

m4 <- lmer(w.sp.var~scale(log(mean_alpha))+scale(log(mntd_alpha))+scale(nat_ratio_alpha)+scale(log(abundance))+(1|site.code),
              data=analyzed.df)
summary(m4)
performance::r2(m4)
plot(simulateResiduals(m4))

library(piecewiseSEM)

SEM <- psem(m1,m3,m4)
summary(SEM)

plot(SEM)

library(emmeans)
m1 <- glmmTMB(comm.var~scale(log(mean_alpha))*trt+scale(nat_ratio_alpha)*trt+(1|site.code/block),
              data=analyzed.df,family=tweedie)
summary(m1)
car::Anova(m1)
performance::r2(m1)
emt1 <- emtrends(m1, pairwise ~ trt, var = "nat_ratio_alpha",adjust="fdr")

m2 <- glmmTMB(scaled_sync_w_r~scale(log(mean_alpha))*trt+scale(nat_ratio_alpha)*trt+scale(log(mntd_alpha))*trt+(1|site.code/block),
           data=analyzed.df,family=beta_family())
summary(m2)
car::Anova(m2)
performance::r2(m2)
emt2 <- emtrends(m2, pairwise ~ trt, var = "nat_ratio_alpha",adjust="none")
emt2

m3 <- glmmTMB(w.sp.var~scale(log(mean_alpha))*trt+scale(nat_ratio_alpha)*trt+(1|site.code/block),
           data=analyzed.df,family=tweedie)
summary(m3)
car::Anova(m3)
performance::r2(m3)
emt3 <- emtrends(m3, pairwise ~ trt, var = "nat_ratio_alpha",adjust="none")
emt3

###
  
  cor(analyzed.df[,c("nat_ratio_alpha","mean_alpha_uw","mpd_uw_alpha","pielou","mntd_uw_alpha","abundance")])
  all_result <- NULL
  for (i in 1:length(unique(all_cor_df_pooled$trt))) {
    analyzed.df <- all_cor_df_pooled[all_cor_df_pooled$trt==unique(all_cor_df_pooled$trt)[[i]],]
    analyzed.df$mean_alpha <- exp(analyzed.df$mean_alpha)
    
    analyzed.df$scaled_sync_w_r <- (analyzed.df$sync_w_r - (-1))/2
    analyzed.df$scaled_sync_uw_r <- (analyzed.df$sync_uw_r - (-1))/2
    
    library(emmeans)
    m1 <- glmmTMB(comm.var~scale(log(mean_alpha_uw))+
                    scale(nat_ratio_alpha)+
                    scale(log(mpd_uw_alpha))+
                    scale(log(mntd_uw_alpha))+
                    scale(pielou)+
                    (1|site.code),
                  data=analyzed.df,family=tweedie)
    summary(m1)
    car::Anova(m1)
    #performance::r2(m1)
    
    m2 <- glmmTMB(scaled_sync_w_r~scale(log(mean_alpha_uw))+
                    scale(nat_ratio_alpha)+
                    scale(log(mpd_uw_alpha))+
                    scale(log(mntd_uw_alpha))+
                    scale(pielou)+
                    (1|site.code),
                  data=analyzed.df,family=beta_family())
    summary(m2)
    car::Anova(m2)
    #performance::r2(m2)
    
    
    m3 <- glmmTMB(w.sp.var~scale(log(mean_alpha_uw))+
                    scale(nat_ratio_alpha)+
                    scale(log(mpd_uw_alpha))+
                    scale(log(mntd_uw_alpha))+
                    scale(pielou)+
                    (1|site.code),
                  data=analyzed.df,family=tweedie)
    summary(m3)
    car::Anova(m3)
    performance::r2(m3)
    #plot(simulateResiduals(m3))
    
    model_result <- rbind(data.frame(summary(m1)$coefficients$cond,pred=names(fixef(m1)$cond),var="comm.var",trt=unique(all_cor_df_pooled$trt)[[i]]),
                          data.frame(summary(m2)$coefficients$cond,pred=names(fixef(m2)$cond),var="sync.w.r",trt=unique(all_cor_df_pooled$trt)[[i]]),
                          data.frame(summary(m3)$coefficients$cond,pred=names(fixef(m3)$cond),var="w.sp.var",trt=unique(all_cor_df_pooled$trt)[[i]])
    )
    
    all_result <- rbind(all_result,model_result)
  }
  ###
library(ggeffects)

predict_df <- as.data.frame(ggemmeans(m2,terms=c("mpd_uw_alpha[50:335,by=3]")),terms_to_colnames=T)
predict_df$predicted <- predict_df$predicted*2-1
predict_df$conf.low <- predict_df$conf.low*2-1
predict_df$conf.high <- predict_df$conf.high*2-1

p <- ggplot(data=predict_df,aes(x=mpd_uw_alpha,y=predicted))+
  geom_line(aes())+
  geom_point(data=analyzed.df,aes(x=mpd_uw_alpha,y=sync_w_r))+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.1)+
  geom_hline(yintercept=0)+
  labs(y="Sychrony")+
  theme_bw()

plot(p)

###
predict_df <- as.data.frame(ggemmeans(m3,terms=c("nat_ratio_alpha[0:1,by=0.05]")),terms_to_colnames=T)

p2 <- ggplot(data=predict_df,aes(x=nat_ratio_alpha,y=predicted))+
  geom_line()+
  geom_point(data=analyzed.df,aes(x=nat_ratio_alpha,y=w.sp.var))+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.1)+
  labs(y="Species-level variance (Abundance-weighted)")+
  theme_bw()

plot(p2)

predict_df <- as.data.frame(ggemmeans(m3,terms=c("mntd_uw_alpha")),terms_to_colnames=T)


p3 <- ggplot(data=predict_df,aes(x=mntd_uw_alpha,y=predicted))+
  geom_line()+
  geom_point(data=analyzed.df,aes(x=mntd_uw_alpha,y=w.sp.var))+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.1)+
  labs(y="Species-level variance (Abundance-weighted)",x="MNTD")+
  theme_bw()

plot(p3)
