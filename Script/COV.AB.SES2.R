library(picante)

##open git. 

####covariance-abundance SEMs

#year(rows) by species (cols) matrix

#function calculating SES covariance and SES aboundance trend (lm coef), 
#X is time by species matrix

SES.COV.AB<-function(X,year,iter=999, null.model=F){
  set.seed(1000)
  #observed average covariances   
  cov.ob<-cov(X)
  var<-diag(cov.ob)
  cor.ob <- cor(X)
  ave.cov <- mean(cov.ob[lower.tri(cov.ob)],na.rm=T)
  ave.cor <- mean(cor.ob[lower.tri(cor.ob)],na.rm=T)
  comm.var <- sd(rowSums(X))/mean(rowSums(X))
  sp.var <- mean(apply(X,2,function(x) sd(x)/mean(x)))
  w.sp.var <- weighted.mean(apply(X,2,function(x) sd(x)/mean(x)),colSums(X))
  
  #observed abundance trend
  Y <- unlist(year-min(year))
  
  coef.lm<-function(X,Y){
    coefs<-NULL
    for (i in 1:ncol(X)){
      tmp.lm<-lm(X[,i]~Y)
      coefs[i]<-tmp.lm$coefficients[2]
      
    }
    names(coefs)<-colnames(X)
    return(coefs)
  }
  
  coef.obs<-coef.lm(X,Y)
  
  sync_r <- function(X) {
    r_vector <- NULL
    for (j in 1:ncol(X)) {
      sp <- unlist(X[,j])
      other_sp <- rowSums(X[,-j])
      r_vector <- c(r_vector,cor(sp,other_sp))
    }

    RA <- colSums(X)/sum(X)
    sync_w_r <- weighted.mean(r_vector,RA,na.rm=T)
    sync_uw_r <- mean(r_vector,na.rm=T)

    return(c(sync_w_r=sync_w_r,
             sync_uw_r=sync_uw_r))
  }

  sync_r_obs <- sync_r(X)
  
  #randomizer for both nullcov and nullab
  null.cov<-null.cor<-null.comm.var <- null.sp.var <- NULL
  
  null.coef <- sync_r_list <-list()
  
  if (null.model == T) {
      for (i in 1:iter){
      X.null<-randomizeMatrix(X) #i think this is the appropriate null model - if we fix both col and row sums, systematic effect of year on all species would not be detected.
      tmp.null.cov<-cov(X.null)
      tmp.null.cor<-cor(X.null)
  
      null.cov<-c(null.cov,mean(tmp.null.cov[lower.tri(tmp.null.cov)]))
      null.cor<-c(null.cov,mean(tmp.null.cor[lower.tri(tmp.null.cor)]))
      
      # sync_r_list[[i]] <- sync_r(X.null)
      null.coef[[i]] <- data.frame(Taxon=colnames(X),coef.ab=coef.lm(X.null,Y),i=i)
      
      null.comm.var <- c(null.comm.var,sd(rowSums(X.null))/mean(rowSums(X.null)))
      
      null.sp.var.vec <- apply(X.null,2,function(x) mean(x)/sd(x))
      null.sp.var <- c(null.sp.var,mean(null.sp.var.vec))
    }
    
    null.coef.df <- do.call(rbind,null.coef) %>% 
      group_by(Taxon) %>% 
      summarize(null.coef.ab=mean(coef.ab),
                null.coef.ab.sd=sd(coef.ab))
    
    coef.obs <- coef.obs[order(names(coef.obs))]
    null.coef.df$obs.coef <- coef.obs
    
    null.cor.ave<-mean(null.cor)
    null.cor.sd<-sd(null.cor)
    
    null.cov.ave<-mean(null.cov)
    null.cov.sd<-sd(null.cov)
    
    null.sync.r <- as.data.frame(do.call(rbind,sync_r_list))
    null.sync.w.r <- mean(null.sync.r$sync_w_r)
    null.sync.uw.r <-mean(null.sync.r$sync_uw_r)
    null.sync.w.r.sd <-sd(null.sync.r$sync_w_r)
    null.sync.uw.r.sd <- sd(null.sync.r$sync_uw_r)
      
  } else {
    null.cov.ave <- null.cov.sd <- null.cov.ave <- null.cov.sd <- null.cor.ave <- null.cor.sd <- NA
    null.sync.w.r <- null.sync.uw.r <- null.sync.w.r.sd <- null.sync.uw.r.sd <- NA
    null.coef.df <- data.frame(null.coef.ab = NA,null.coef.ab.sd = NA)
  }
  
  return(list(SES.analysis=data.frame(Taxon = order(names(coef.obs)),
                                      ab.coef.obs=coef.obs,
                                      SES.coef.ab=((coef.obs-null.coef.df$null.coef.ab)/null.coef.df$null.coef.ab.sd),
                                      ave.cov,
                                      SES.cov=((ave.cov-null.cov.ave)/null.cov.sd),
                                      ave.cor,
                                      SES.cor=((ave.cor-null.cor.ave)/null.cor.sd),
                                      comm.var = comm.var,
                                      SES.comm.var = (comm.var-mean(null.comm.var))/sd(null.comm.var),
                                      sp.var= sp.var,
                                      w.sp.var= w.sp.var,
                                      sync_w_r = sync_r_obs["sync_w_r"],
                                      sync_uw_r = sync_r_obs["sync_uw_r"],
                                      SES.sync.w.r = (sync_r_obs["sync_w_r"]-null.sync.w.r)/null.sync.w.r.sd,
                                      SES.sync.uw.r = (sync_r_obs["sync_uw_r"]-null.sync.uw.r)/null.sync.uw.r.sd,
                                      var),
              Observed.covariances=cov.ob))
}

# SES.COV.AB<-function(X,iter=999){
#   #observed average covariances   
#   cov.ob<-cov(X)
#   var<-diag(cov.ob)
#   diag(cov.ob)<-NA
#   ave.cov<-apply(cov.ob,1,function(X) mean(X,na.rm=T))
#   
#   #observed abundance trend
#   
#   coef.lm<-function(X){
#     coefs<-NULL
#     Y<-1:nrow(X)
#     for (i in 1:ncol(X)){
#       tmp.lm<-lm(X[,i]~Y)
#       coefs[i]<-tmp.lm$coefficients[2]
#       
#     }
#     names(coefs)<-names(X)
#     return(coefs)
#   }
#   
#   coef.obs<-coef.lm(X)
#   
#   
#   #randomizer for both nullcov and nullab
#   null.cov<-matrix(0,nrow=iter,ncol=length(ave.cov))
#   null.coef<-matrix(0,nrow=iter,ncol=length(ave.cov))
#   
#   for (i in 1:iter){
#     X.null<-randomizeMatrix(X) #i think this is the appropriate null model - if we fix both col and row sums, systematic effect of year on all species would not be detected.
#     tmp.null<-cov(X.null)
#     diag(tmp.null)<-NA
#     null.cov[i,]<-apply(tmp.null,1,function(X) mean(X,na.rm=T))
#     
#     null.coef[i,]<-coef.lm(X.null)
#     
#   }
#   colnames(null.cov)<-names(X)
#   colnames(null.coef)<-names(X)
#   
#   null.cov.ave<-apply(null.cov,2,mean)
#   null.cov.sd<-apply(null.cov,2,sd)
#   
#   null.coef.ave<-apply(null.coef,2,mean)
#   null.coef.sd<-apply(null.coef,2,sd)
#   
#   return(list(SES.analysis=data.frame(ab.coef.obs=coef.obs,
#                                       SES.coef.ab=((coef.obs-null.coef.ave)/null.coef.sd),
#                                       ave.cov,
#                                       SES.cov=((ave.cov-null.cov.ave)/null.cov.sd),
#                                       var),
#               Observed.covariances=cov.ob))
# }
