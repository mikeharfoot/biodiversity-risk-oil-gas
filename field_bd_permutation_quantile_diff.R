
wants <- c("countrycode","RColorBrewer")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
has <- wants %in% .packages()
for(p in wants[!has]) library(p,character.only = T)

options(stringsAsFactors = F)
use.pa.cats = c(1:7,90,92)
use.pa.cats = list(1,2,3,4,list(c(5,6,7)))
#use.pa.cats = c(0:7,90,92)

#Specify the input directory
indsn = "C:/Users/mikeha/Work/Fossil Fuels and Biodiversity/Fields/"

outdir = "C:/Users/mikeha/Dropbox/Fossil Fuels and biodiversity/Figures/"


# Process the field biodiversity data -------------------------------------

curr.bd.all = read.csv(paste(indsn,"Fields_all_current_moll_bd_no_dup.csv",sep=""))
fut.bd.all = read.csv(paste(indsn,"Fields_all_future_moll_bd_no_dup.csv",sep=""))


curr.bd.all$InvArea[which(curr.bd.all$InvArea < 0)] = 0.0
fut.bd.all$InvArea[which(fut.bd.all$InvArea < 0)] = 0.0

#Convert from per m2 to per 1000 km2
curr.bd.all$InvArea = curr.bd.all$InvArea * 1E9
fut.bd.all$InvArea = fut.bd.all$InvArea * 1E9

bd.colnames = c("SpCount","InvArea")
bd.labels = c("Assessed species richness","Mean range rarity of assessed species")
bd.labels.short = c("SR","RR")
fut.bd.cols = sapply(bd.colnames, function(x) which(colnames(fut.bd.all) == x))
curr.bd.cols = sapply(bd.colnames, function(x) which(colnames(curr.bd.all) == x))

curr.inds = vector(mode = "list",length = 2)
fut.inds = vector(mode = "list",length = 2)
curr.inds[1] = list(which(curr.bd.all$ONSHORE_OF %in% c("Onshore", " ") & curr.bd.all$SpCount > 0))
curr.inds[2] = list(which(curr.bd.all$ONSHORE_OF %in% c("Offshore","Onshore/Offshore") & curr.bd.all$SpCount > 0))
fut.inds[1] = list(which(fut.bd.all$ONSHORE_OF %in% c("Onshore", " ") & fut.bd.all$SpCount > 0))
fut.inds[2] = list(which(fut.bd.all$ONSHORE_OF %in% c("Offshore","Onshore/Offshore") & fut.bd.all$SpCount > 0))

quants = c(0,0.005,0.025,0.25,0.5,0.75,0.975,0.995,1.0)
xlims = c(0,8)

wilcox.p = c(1,1)


# Read in the infrastructure grid -----------------------------------------

indsn="C:/Users/mikeha/Work/Fossil Fuels and biodiversity/"
comb.shp.name = "BD_Realm_OGM_Realm_Continent_no_dup_all_well_count"
master<-read.csv(paste(indsn,"Combined/",comb.shp.name,".csv",sep=""))

#Convert from per m2 to per 1000 km2
master$InvArea = master$InvArea * 1E9

ogm.colnames = c("Count_","sum_length","Count_2")
ogm.labels= c("Wells","Pipeline length","Refineries")
ogm.cols = unlist(sapply(ogm.colnames, function(x) which(colnames(master) == x)))

quantiles = seq(0,1,0.01)

permute_fields_deltabd<-function(f.l,bd,q)
{
  f.l.st = sample(f.l,replace=F)
  q.uf = quantile(bd[which(f.l.st == "unexploited")],probs = q)
  q.ef = quantile(bd[which(f.l.st == "exploited")],probs = q)
  q.uf-q.ef
}


terr.sr.comb = c(curr.bd.all[unlist(curr.inds[1]),curr.bd.cols[1]],fut.bd.all[unlist(fut.inds[1]),fut.bd.cols[1]])
ex.terr.sr.quants = quantile(curr.bd.all[unlist(curr.inds[1]),curr.bd.cols[1]],probs=quantiles)
un.terr.sr.quants = quantile(fut.bd.all[unlist(fut.inds[1]),fut.bd.cols[1]],probs=quantiles)
delta.terr.sr = un.terr.sr.quants-ex.terr.sr.quants
mar.sr.comb = c(curr.bd.all[unlist(curr.inds[2]),curr.bd.cols[1]],fut.bd.all[unlist(fut.inds[2]),fut.bd.cols[1]])
ex.mar.sr.quants = quantile(curr.bd.all[unlist(curr.inds[2]),curr.bd.cols[1]],probs=quantiles)
un.mar.sr.quants = quantile(fut.bd.all[unlist(fut.inds[2]),fut.bd.cols[1]],probs=quantiles)
delta.mar.sr = un.mar.sr.quants-ex.mar.sr.quants


terr.rr.comb = c(curr.bd.all[unlist(curr.inds[1]),curr.bd.cols[2]],fut.bd.all[unlist(fut.inds[1]),fut.bd.cols[2]])
ex.terr.rr.quants = quantile(curr.bd.all[unlist(curr.inds[1]),curr.bd.cols[2]],probs=quantiles)
un.terr.rr.quants = quantile(fut.bd.all[unlist(fut.inds[1]),fut.bd.cols[2]],probs=quantiles)
delta.terr.rr = un.terr.rr.quants-ex.terr.rr.quants
mar.rr.comb = c(curr.bd.all[unlist(curr.inds[2]),curr.bd.cols[2]],fut.bd.all[unlist(fut.inds[2]),fut.bd.cols[2]])
ex.mar.rr.quants = quantile(curr.bd.all[unlist(curr.inds[2]),curr.bd.cols[2]],probs=quantiles)
un.mar.rr.quants = quantile(fut.bd.all[unlist(fut.inds[2]),fut.bd.cols[2]],probs=quantiles)
delta.mar.rr = un.mar.rr.quants-ex.mar.rr.quants


terr.f.labs = c(rep("exploited",length(curr.bd.all[unlist(curr.inds[1]),curr.bd.cols[1]])),
           rep("unexploited",length(fut.bd.all[unlist(fut.inds[1]),fut.bd.cols[1]])))

mar.f.labs = c(rep("exploited",length(curr.bd.all[unlist(curr.inds[2]),curr.bd.cols[1]])),
                rep("unexploited",length(fut.bd.all[unlist(fut.inds[2]),fut.bd.cols[1]])))

n.perms = 10000


n.gt.obs<-function(o,perms)
{
  n = list()
  for(i in 1:length(o))
  {
    n[[i]] = length(which(perms[,i] >= o[i]))
  }
  unlist(n)
}

n.lt.obs<-function(o,perms)
{
  n = list()
  for(i in 1:length(o))
  {
    n[[i]] = length(which(perms[,i] <= o[i]))
  }
  unlist(n)
}


sig.test<-function(x) {
  sig.i = which(x < 0.025)
  return(ifelse(length(sig.i > 0),sig.i,0))
}


terr.sr.perms = t(replicate(n.perms, permute_fields_deltabd(terr.f.labs,terr.sr.comb,quantiles)))
terr.sr.n = array(NA, dim=c(2,length(delta.terr.sr)))
terr.sr.n[1,] = n.gt.obs(delta.terr.sr,terr.sr.perms) 
terr.sr.n[2,] = n.lt.obs(delta.terr.sr,terr.sr.perms) 
terr.sr.sig = terr.sr.n/n.perms
terr.sr.ul = unlist(apply(terr.sr.sig,2,sig.test))

mar.sr.perms = t(replicate(n.perms, permute_fields_deltabd(mar.f.labs,mar.sr.comb,quantiles)))
mar.sr.n = array(NA, dim=c(2,length(delta.mar.sr)))
(mar.sr.n[1,] = n.gt.obs(delta.mar.sr,mar.sr.perms))
(mar.sr.n[2,] = n.lt.obs(delta.mar.sr,mar.sr.perms))
mar.sr.sig = mar.sr.n/n.perms
mar.sr.ul = unlist(apply(mar.sr.sig,2,sig.test))

terr.rr.perms = t(replicate(n.perms, permute_fields_deltabd(terr.f.labs,terr.rr.comb,quantiles)))
terr.rr.n = array(NA, dim=c(2,length(delta.terr.sr)))
terr.rr.n[1,] = n.gt.obs(delta.terr.rr,terr.rr.perms) 
terr.rr.n[2,] = n.lt.obs(delta.terr.rr,terr.rr.perms) 
terr.rr.sig = terr.rr.n/n.perms
terr.rr.ul = unlist(apply(terr.rr.sig,2,sig.test))

mar.rr.perms = t(replicate(n.perms, permute_fields_deltabd(mar.f.labs,mar.rr.comb,quantiles)))
mar.rr.n = array(NA, dim=c(2,length(delta.mar.sr)))
(mar.rr.n[1,] = n.gt.obs(delta.mar.rr,mar.rr.perms))
(mar.rr.n[2,] = n.lt.obs(delta.mar.rr,mar.rr.perms))
mar.rr.sig = mar.rr.n/n.perms
mar.rr.ul = unlist(apply(mar.rr.sig,2,sig.test))
