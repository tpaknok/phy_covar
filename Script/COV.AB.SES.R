library(picante)

##open git. 

####covariance-abundance SEMs

#year(rows) by species (cols) matrix

#function calculating SES covariance and SES aboundance trend (lm coef), 
#X is time by species matrix
SES.cor.ab<-function(X,year,iter=999,seed=999,plot_info,sp_info){
  set.seed(seed)
  X <- as.data.frame(X)
  #observed average covariances
  cor.ob <- cor(X)
  var<-diag(cov(X))
  diag(cor.ob)<-NA

  #observed abundance trend
  
  coef.lm<-function(X,year){
    coefs<-NULL
    Y<-unlist(year-min(year))
    for (i in 1:ncol(X)){
      tmp.lm<-lm(X[,i]~Y)
      coefs[i]<-tmp.lm$coefficients[2]
      
    }
    names(coefs)<-names(X)
    return(coefs)
  }
  
  coef.obs<-coef.lm(X,year)
  
  null.cor.list <- list()
  null.coef <- NULL
  #randomizer for both nullcov and nullab
   for (i in 1:iter){
    X.null<-randomizeMatrix(X) #i think this is the appropriate null model - if we fix both col and row sums, systematic effect of year on all species would not be detected.
    tmp.null<-cor(X.null)
    diag(tmp.null)<-NA
    
    
    null.cor.list[[i]] <- c(tmp.null)
    null.coef <- cbind(null.coef,coef.lm(X.null,year))
   }
  
  null.cor.df <- do.call(cbind,null.cor.list)
  
  coef.SES <- (coef.obs-rowMeans(null.coef))/apply(null.coef,1,sd)
  SES.pair <- (cor.ob - apply(null.cor.df,1,mean))/apply(null.cor.df,1,sd)
  
  sp_name <- expand.grid(colnames(X),colnames(X))
  colnames(sp_name) <- c("sp1","sp2")
  SES.pair.df <- data.frame(plot_info,
                            sp_name,
                            SES.pair=c(SES.pair),
                            expand.grid(1:ncol(X),1:ncol(X))
                            )
  
  
  SES.pair.df <- SES.pair.df %>%
    filter(Var1 > Var2) %>%
    mutate(obs.cor = cor.ob[lower.tri(cor.ob)]) %>%
    mutate(sp1 = as.character(sp1),
           sp2 = as.character(sp2))
  
  SES.pair.df <- SES.pair.df %>%
    left_join(sp_info,by=join_by(sp1==Accepted_SPNAME)) %>%
    left_join(sp_info,by=join_by(sp2==Accepted_SPNAME))
  
  species.df <- data.frame(sp=colnames(X),
                       var=unname(var),
                       mean_abundance=unname(colMeans(X)),
                       SES.pair.avg = unname(colMeans(SES.pair,na.rm=T)),
                       coef.SES = coef.SES,
                       plot_info)
  
  species.df <- species.df %>%
    left_join(sp_info,by=join_by(sp==Accepted_SPNAME))
    
  result <- list(SES.pair.long=SES.pair.df,
                 SEs.pair.matrix=SES.pair,
                 species.df=species.df)
  
  return(result)
}



