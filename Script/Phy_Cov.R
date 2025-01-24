####Phylogeny-covariance
#run gen.cov.R first to gen COVs list
library(EnvStats)
library(picante)
library(V.PhyloMaker)
library(stringr)
library(rlist)
library(RRphylo)
setwd("/Users/marccadotte/Dropbox/2018WorkingFiles/Marc2018/NutNet/Catford meeting/Phy_Cov/")

source("subTree.R")
source("Xunder.R")

#n.tree<-read.tree("Nutnet_tree.txt")

#build tree
int.status.all<-read.csv("~/Dropbox/2018WorkingFiles/Marc2018/NutNet/Catford meeting/COV-ab_SEM/NutNet_species-list-by-origin_Nov2022.csv")
plants<-int.status.all[int.status.all$lifeform!="NULL",]

species<-plants$standard_taxon
family<-plants$Family
genus<-word(species,1)

tax<-data.frame(species, genus, family)

#n.tree<-phylo.maker(tax)
#tr<-n.tree$scenario.3

#tr<-read.tree("Nutnet_tree.txt")
quartz()
plot(tr,show.tip.label = FALSE,type="radial",edge.width=2)

tr$tip.label<-toupper(tr$tip.label)
tr$tip.label<-Xunder(tr$tip.label)

tr.cov<-vcv(tr)

#for-loop or apply through COV list
COV.paired<-list()
  
for (i in 1:length(COVs)) {
sample.time.cov<-COVs[[i]][lower.tri(COVs[[i]])]
sub.phy.cov<-tr.cov[match(row.names(COVs[[i]]),row.names(tr.cov)),
                       match(row.names(COVs[[i]]),row.names(tr.cov))]
sample.phy.cov<-sub.phy.cov[lower.tri(sub.phy.cov)]

COV.paired[[i]]<-data.frame(sample=sample.time.cov,phy=sample.phy.cov)

}

#put below into for loop or into above
#tmp<-COV.paired[[1]]
#outlier<-rosnerTest(tmp$sample,k=1)

#if (outlier$all.stats[,8]==1){
  
#  if (outlier$all.stats[,4]<0){
    
#    tmp[round(tmp$sample,digits=2)!=round(outlier$all.stats[,4],digits=2),]
#  }
  
#  if (outlier$all.stats[,4]>0){
    
 #   tmp[round(tmp$sample,digits=2)!=round(outlier$all.stats[,4],digits=2),]
#  }
#}




#plot(COV.paired[[1]]$phy,COV.paired[[1]]$sample)
#cor(COV.paired[[1]]$phy,COV.paired[[1]]$sample)

cors<-sapply(COV.paired,function(x) {
  x<-na.omit(x)
  cor(x$phy,x$sample)
})

quartz()
hist(cors,main=NULL,
     xlab="Temporal covariance-phylogenetic distance correlation coeficients")
abline(v=mean(cors,na.rm=TRUE),lty="dashed",lwd=3,col="chocolate4")
T<-t.test(cors)
abline(v=T$conf.int,lty="dashed",col="chocolate2")

all.COV.paired<-list.rbind(COV.paired)
all.COV.paired<-na.omit(all.COV.paired)
cor(all.COV.paired)

plot(all.COV.paired$phy,all.COV.paired$sample)

###phylo.signal of variances and ave abundances for sites with more than 15 species

K.var<-NULL
P.var<-NULL
Z.var<-NULL
K.ab<-NULL
P.ab<-NULL
Z.ab<-NULL

for (i in 1:length(SES.out)){
  tmp<-SES.out[[i]]
  tmp<-na.omit(tmp)
  
  if (dim(tmp)[1]<3){
    K.var[i]<-NA
    K.ab[i]<-NA
  }
  
  if (dim(tmp)[1]>2){
    st<-subTree(tr,row.names(tmp))
    tmp<-tmp[match(st$tip.label,row.names(tmp)),]
    vars<-tmp$var
    names(vars)<-row.names(tmp)
    
    if (is.binary(st)==FALSE){
      st<-fix.poly(st,type="resolve")
    }
    PS.var<-phylosignal(vars,st)
    K.var[i]<-PS.var[1,1]
    P.var[i]<-PS.var[1,4]
    Z.var[i]<-PS.var[1,5]
    
    abs<-tmp$ave.ab
    names(abs)<-row.names(tmp)
    
    PS.ab<-phylosignal(abs,st)
    K.ab[i]<-PS.ab[1,1]
    P.ab[i]<-PS.ab[1,4]
    Z.ab[i]<-PS.ab[1,5]
  }
}



quartz()
par(mfrow=c(1,2))
##remove outlier plotting K < 2
hist((K.var[K.var<2]),main=NULL, xlab="Blomberg's K for species temporal variance")
abline(v=mean(K.var,na.rm=TRUE),lty="dashed",lwd=3,col="chocolate4")
Tv<-t.test(K.var)
abline(v=Tv$conf.int,lty="dashed",lwd=1,col="chocolate2")

hist((K.ab[K.ab<2]),main=NULL, xlab="Blomberg's K for species average abundance")
abline(v=mean(K.ab,na.rm=TRUE),lty="dashed",lwd=3,col="chocolate4")
Ta<-t.test(K.ab)
abline(v=Ta$conf.int,lty="dashed",lwd=1,col="chocolate2")

#hist(P.var)
#hist(Z.var)

#hist(P.ab)
#hist(Z.ab)
