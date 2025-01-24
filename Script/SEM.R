library(piecewiseSEM)

library(lme4)
library(lmerTest)
library(ggeffects)
all_cor_df_pooled <- all_cor_df %>%
  filter(unclassified_ratio < 0.1 & unclassified_ratio_abundance < 0.1)

####

analyzed.df <- subset(all_cor_df_pooled,trt=="NPK")
analyzed.df$mean_alpha <- exp(analyzed.df$mean_alpha)

analyzed.df$scaled_sync_w_r <- (analyzed.df$sync_w_r - (-1))/2
analyzed.df$scaled_sync_uw_r <- (analyzed.df$sync_uw_r - (-1))/2

analyzed.df$site.code <- as.factor(analyzed.df$site.code)
analyzed.df$block.code <- paste0(analyzed.df$site.code,"_",analyzed.df$block)
analyzed.df$block.code <- as.factor(analyzed.df$block.code)
analyzed.df$block <- as.factor(analyzed.df$block)

analyzed.df <- analyzed.df %>%
  group_by(block.code,site.code) %>%
  summarize(nat_ratio_alpha = mean(nat_ratio_alpha),
            mean_alpha_uw = mean(mean_alpha_uw),
            mpd_uw_alpha = mean(mpd_uw_alpha),
            mntd_uw_alpha = mean(mntd_uw_alpha),
            pielou = mean(pielou),
            sync_w_r = mean(sync_w_r),
            sync_uw_r = mean(sync_uw_r),
            w.sp.var = mean(w.sp.var),
            comm.var = mean(comm.var),
            b.comm.var = mean(b.comm.var))

analyzed.df$log_mean_alpha_uw <- log(analyzed.df$mean_alpha_uw)
analyzed.df$log_mpd_uw_alpha<- log(analyzed.df$mpd_uw_alpha)
analyzed.df$log_mntd_uw_alpha<- log(analyzed.df$mntd_uw_alpha)

cor(scale(analyzed.df$nat_ratio_alpha),scale(analyzed.df$nat_ratio_alpha)^2)
cor(analyzed.df[,c("nat_ratio_alpha","log_mean_alpha_uw","log_mpd_uw_alpha","pielou","log_mntd_uw_alpha")])

m1 <- lmer(log_mean_alpha_uw~nat_ratio_alpha+(1|site.code),data=analyzed.df)
summary(m1)

m2 <- lmer(log_mpd_uw_alpha~nat_ratio_alpha+(1|site.code),data=analyzed.df)
summary(m2)

m3 <- lmer(log_mntd_uw_alpha~nat_ratio_alpha+(1|site.code),data=analyzed.df)
summary(m3)

m4 <- lmer(pielou~nat_ratio_alpha+(1|site.code),data=analyzed.df)
summary(m4)

m5 <- lmer(w.sp.var~log_mean_alpha_uw +
                nat_ratio_alpha+
                log_mpd_uw_alpha+
                log_mntd_uw_alpha+
                pielou+
                (1|site.code),
              data=analyzed.df)
summary(m5)
plot(simulateResiduals(m5))


m6 <- lmer(sync_w_r~log_mean_alpha_uw+
                nat_ratio_alpha+
                log_mpd_uw_alpha+
                log_mntd_uw_alpha+
                pielou+
                (1|site.code),
              data=analyzed.df)
summary(m6)
plot(simulateResiduals(m6))

m7 <- lmer(comm.var~sync_w_r+
                w.sp.var+
                log_mean_alpha_uw+
                nat_ratio_alpha+
                log_mpd_uw_alpha+
                log_mntd_uw_alpha+
                pielou+
                (1|site.code),
              data=analyzed.df)

plot(simulateResiduals(m7))

m8 <- lmer(b.comm.var~comm.var+
             (1|site.code),
           data=analyzed.df)

summary(m8)
plot(simulateResiduals(m7))

sem1 <- psem(m5,m6,m7,m8,
             log_mean_alpha_uw %~~% nat_ratio_alpha,
             log_mean_alpha_uw %~~% log_mpd_uw_alpha,
             log_mean_alpha_uw %~~% log_mntd_uw_alpha)
summary(sem1)

####
library(itsadug)
library(mgcv)

m5 <- gam(w.sp.var~ s(mean_alpha_uw)+
            s(nat_ratio_alpha)+
            s(mpd_uw_alpha)+
            s(mntd_uw_alpha)+
            s(pielou)+
            s(site.code,bs="re"),
          data = analyzed.df,
          method="REML")
summary(m5)

###
p <- plot_smooth(m5,view="mntd_uw_alpha",rm.ranef = TRUE)
plot(simulateResiduals(m5))
predict_df <- p$fv

###
library(ggplot2)

p <- ggplot(data=predict_df,aes(y=w.sp.var,x=mntd_uw_alpha))+
  geom_point(data=analyzed.df)+
  geom_ribbon(aes(y=fit,ymin=ll,ymax=ul),alpha=0.25)+
  geom_line(aes(y=fit),linetype=2)+
  labs(x="MNTD",y="Weighted species variability")+
  theme_bw()
plot(p)

ggsave("Figure/sp_var_mntd.tiff",height=4,width=4,dpi=600,compression="lzw")

p <- plot_smooth(m5,view="nat_ratio_alpha",rm.ranef = TRUE)
plot(simulateResiduals(m5))
predict_df <- p$fv

library(ggplot2)

p <- ggplot(data=predict_df,aes(y=w.sp.var,x=nat_ratio_alpha))+
  geom_point(data=analyzed.df)+
  geom_ribbon(aes(y=fit,ymin=ll,ymax=ul),alpha=0.25)+
  geom_line(aes(y=fit))+
  labs(x="Native species richness ratio",y="Weighted species variability")+
  theme_bw()
plot(p)

ggsave("Figure/sp_var_nat.tiff",height=4,width=4,dpi=600,compression="lzw")
###
###

m6 <- gam(sync_w_r~ s(mean_alpha_uw)+
            s(nat_ratio_alpha)+
            s(mpd_uw_alpha)+
            s(mntd_uw_alpha)+
            s(pielou)+
            s(site.code,bs="re"),
          method="REML",
          data = analyzed.df)
summary(m6)
plot(simulateResiduals(m6))
p <- plot_smooth(m6,view="nat_ratio_alpha",rm.ranef = TRUE)
predict_df <- p$fv

p <- ggplot(data=predict_df,aes(y=sync_w_r,x=nat_ratio_alpha))+
  geom_point(data=analyzed.df)+
  geom_ribbon(aes(y=fit,ymin=ll,ymax=ul),alpha=0.25)+
  geom_line(aes(y=fit))+
  labs(x="Native species richness ratio",y="Synchrony")+
  theme_bw()
plot(p)
ggsave("Figure/sync_nat.tiff",height=4,width=4,dpi=600,compression="lzw")

p <- plot_smooth(m6,view="mntd_uw_alpha",rm.ranef = TRUE)
predict_df <- p$fv

p <- ggplot(data=predict_df,aes(y=sync_w_r,x=mntd_uw_alpha))+
  geom_point(data=analyzed.df)+
  geom_ribbon(aes(y=fit,ymin=ll,ymax=ul),alpha=0.25)+
  geom_line(aes(y=fit))+
  labs(x="MNTD",y="Synchrony")+
  theme_bw()
plot(p)
ggsave("Figure/sync_MTND.tiff",height=4,width=4,dpi=600,compression="lzw")

p <- plot_smooth(m6,view="mpd_uw_alpha",rm.ranef = TRUE)
predict_df <- p$fv

p <- ggplot(data=predict_df,aes(y=sync_w_r,x=mpd_uw_alpha))+
  geom_point(data=analyzed.df)+
  geom_ribbon(aes(y=fit,ymin=ll,ymax=ul),alpha=0.25)+
  geom_line(aes(y=fit))+
  labs(x="MPD",y="Synchrony")+
  theme_bw()
plot(p)
ggsave("Figure/sync_MPD_NPK.tiff",height=4,width=4,dpi=600,compression="lzw")

##
m7 <- gam(comm.var~ s(sync_w_r)+
            s(w.sp.var)+
            s(mean_alpha_uw)+
            s(nat_ratio_alpha)+
            s(mpd_uw_alpha)+
            s(mntd_uw_alpha)+
            s(pielou)+
            s(site.code,bs="re"),
          method="REML",
          data = analyzed.df)
summary(m7)

##
m8 <- gam(b.comm.var~ s(comm.var)+
            s(site.code,bs="re"),
          data = analyzed.df,
          method="REML")
summary(m8)

p <- plot_smooth(m8,view="comm.var",rm.ranef = TRUE)
predict_df <- p$fv
p <- ggplot(data=predict_df,aes(y=b.comm.var,x=comm.var))+
  geom_point(data=analyzed.df)+
  geom_ribbon(aes(y=fit,ymin=ll,ymax=ul),alpha=0.25)+
  geom_line(aes(y=fit))+
  labs(x="Cover variabiltiy",y="Biomass variability")+
  theme_bw()
plot(p)

GAM <- psem(m5,m6,m7,m8)
summary(GAM)
