# Lion fish preference index
#======================
#======================
# 2019 JANUARY 
# DEFINE WORKING DIRECTORY
# 
setwd("~/Desktop/Enclosures/Enclosure statistics/Final")
####################
# LOAD data provided by AMELIA & KIM
# full data
data_all <- read.csv('lionfishdata.csv')
# covariables
data_x <-data.frame(Starve_time = data_all$Starvationtime_hours,
                     cloudcover = data_all$Cloudcover,sex = data_all$Sex,
                     lf_std_length = data_all$Lionfishstandardlength_cm, lf_wetwgt = data_all$LionfishWetWeight_g,
                     mooncycle = data_all$Mooncycle,depthcaptured = data_all$Depthcaptured_m,bodycond = data_all$BodyCondition,Lionfish_totlength = data_all$Lionfishtotallength_cm)
# # fish units consumed
counts <- data.frame(N_chromis = data_all$Numberchromisconsumed, N_goby = data_all$Numbergobyconsumed,N_wrasse = data_all$Numberwrasseconsumed)
# ####################


################# SAVE & LOAD Rdata area ######################
# save.image("lionfishMartin.RData")
# load("lionfishMartin.RData")
###############################################################
#
# # Use zCompositions to replace zeros using bayesian multiplicative methods
 library('zCompositions')
# # total prey consumed by fish
 c_tot <- counts$N_chromis + counts$N_wrasse + counts$N_goby
# # replacement: cmultRepl with SQ prior
 rep_p <- cmultRepl(counts,method = "SQ") # proportions data with zeros replaced
# # pseudo-counts
 rep_count<-rep_p*c_tot
# # names of variables
 colnames(rep_count) =colnames(rep_p) <- colnames(counts)
# ###################################
#
###################################
# T^2 Hotelling test: contrast H_0: mu = TARGET
# MAIL Michele: TARGET proportion: 11chromis, 55 gobies, 11wrasses
# TARGET PROPORTION in %
(c(11,55,11)/sum(c(11,55,11))*100)# 14.28571 71.42857 14.28571
# CENTER of data set: closed (%) geometric mean
g<-exp(apply(log(rep_p),2,mean))
(g<-g/sum(g)*100)# 60.89472  22.93172  16.17356

###################################
# package for T^2 hotelling
 library(rrcov)
#
# function: create log-ratio balances from a SBP
fsbp2ilr<-function(s){
  # given a vector (or a set of vectors) of a SBP formed by 0,+1,-1 
  # returns the ilr-vector (or set of) or balance
  # author: Martin (2015)
  
  bal=s
  
  if(is.null(dim(s))) 
  {
    
    den<-sum(bal==-1)
    num<-sum(bal==1)  
    bal[s==1]<-sqrt(den/((den+num)*num))
    bal[s==-1]<--sqrt(num/((den+num)*den))
    
    
  }
  else
  {
    numsbp=dim(s)[1]
    for (f in 1:numsbp) {
      den<-sum(bal[f,]==-1)
      num<-sum(bal[f,]==1)  
      bal[f,bal[f,]==1]<-sqrt(den/((den+num)*num))
      bal[f,bal[f,]==-1]<--sqrt(num/((den+num)*den))
      
    }
    
  }
  
  return(bal)
}

# AND A FEASIBLE LOG-RATIO BASIS COULD BE
#(1,-1,-1) = "N_chromis"/g("N_goby","N_wrasse")
#(0, 1,-1) = "N_goby"/"N_wrasse" 
tU3=as.matrix(fsbp2ilr(rbind(c(1,-1,-1),c(0,1,-1))))
U3=t(tU3)
# 
ilr_rep_p<-log(as.matrix(rep_p))%*%U3 
ilr_ref<-log(as.vector(c(11,55,11)/sum(c(11,55,11))))%*%U3 
T2.test(as.matrix(ilr_rep_p),mu=as.vector(ilr_ref))
# 
# RESULTS:
# T2 = 80.475, F = 38.120, df1 = 2, df2 = 18, p-value = 3.383e-07
# 
# INTEPRETATION:
# REJECT NULL Hypotesis
# In average, the preference of the data set is different from the NULL HIPOT
#  of NO-PREFERENCE, the (11,55,11) TARGET composition
#
#
##########################################################
# INDEX of PREFERENCE based on Aitchison distance
#########################################################
# Aitchison distance
fdAit<-function(x,y){
  log(x) - mean(log(x))
  return(sqrt(sum((log(x) - mean(log(x)) - log(y) + mean(log(y)))^2)))
}
#
############################################################
# Non-Scaled Index 
fIndexPref<-function (x,xref) 
{
  # xref = target composition
  if (is.vector(x)) {
    return(fdAit(x,xref)^2)
  }
  else {
    return(apply(x,1,fdAit,xref)^2)
  }
}
#
# [0,1]-Scaled Index 
fIndexPrefS<-function (x,xref) 
{
  # xref = target composition
  if (is.vector(x)) {
    d<-fdAit(x,xref)^2
  }
  else {
    d<-apply(x,1,fdAit,xref)^2
  }
  # scaling to [0,1] and return
  return(1-exp(-d))
}
#
###############################################
# example of "no preference"
fIndexPref(c(3,15,3),c(11,55,11))
fIndexPrefS(c(3,15,3),c(11,55,11))
# check for the first fish
fIndexPref(rep_count[1,],c(11,55,11))
fIndexPrefS(rep_count[1,],c(11,55,11))
# Preference index for all 20 fish
(lindex<-fIndexPref(rep_count,c(11,55,11)))
(lindexS<-fIndexPrefS(rep_count,c(11,55,11)))
############################################
# total for each fish
totalN<-apply(rep_count,1,sum)
#

##############################
# BOOSTRAP analysis for the preference
# null Hypo: H_0: multinomial prob = (11,55,11)
# for each fish I generated 10000 samples assuming the total
# consumed and the target probability. Afterwards I calculated
# the preference index for each of the 10000 samples. This
# gives a distribution of the index of preference under the
# null hypo. I calculated the p-value comparing the
# index of preference of the fish with the target distrution
# Therefore, for each fish, we can consider that null Hypo is H_0: "index pref =0"
#
pvalIndex=pvalIndexS<-rep(0,length(totalN))
# # NUMBER of random samples for multinomial
# numboot=10000
# #####
# # loop for all fishes
# for (numfish in (1:length(totalN))){
# # check fish #2  numfish<-2
# lindfish<-lindex[numfish]
# lindSfish<-lindexS[numfish]
# set.seed(10)
# XRand<-t(rmultinom(numboot,totalN[numfish],prob=c(11,55,11)))
# # Zero replacement
# XRand_rep_p <- cmultRepl(XRand,method = "SQ") # proportions data with zeros replaced
# XRand_rep_count<-XRand_rep_p*totalN[numfish]
# # head(XRand_rep_count,20)
# #
# IndexNullHip<-fIndexPref(XRand_rep_count,c(11,55,11))
# IndexSNullHip<-fIndexPrefS(XRand_rep_count,c(11,55,11))
# # P-valor Index
# pvalIndex[numfish]<-(length(IndexNullHip[IndexNullHip>lindfish]))/length(IndexNullHip)
# #(length(IndexNullHip[IndexNullHip>lindfish])+0.5)/length(IndexNullHip)
# 
# # P-valor [0,1]-index
# pvalIndexS[numfish]<-(length(IndexSNullHip[IndexSNullHip>lindSfish]))/length(IndexSNullHip)
# #(length(IndexSNullHip[IndexSNullHip>lindSfish])+0.5)/length(IndexSNullHip)
# ####
# } # END lop for all fishes
# # P-VALUES
# non-scaled index of preference
(pvalIndex)
# [1] 0.0568 0.0000 0.0000 0.0000 0.0000 0.0901 0.0000 0.0000 0.0923 0.0000 0.0000 0.0000 0.0000 0.0000
# [15] 0.4866 1.0000 0.2898 0.0000 0.0923 0.0000
# scaled index of preference
(pvalIndexS)
# [1] 0.0568 0.0000 0.0000 0.0000 0.0000 0.0901 0.0000 0.0000 0.0923 0.0000 0.0000 0.0000 0.0000 0.0000
# [15] 0.4866 1.0000 0.2898 0.0000 0.0923 0.0000
#
# INTERPRETATION:
# First: scaling doesn't affect. 
# Second: for each fish where p-value is below 0.05 we can assume that the fish has an non-zero
# index of preference, it has a SIGNIFICANT index of preference
#############################################
###
####
# CORRELATION BETWEEN bodycondition and index of preference
# not scaled
# cor(data_x$bodycond,lindex) # 0.479191
# plot(data_x$bodycond,lindex)
# # WITH LOGs
# cor(log(data_x$bodycond),log(lindex)) # 0.3974944
# plot(log(data_x$bodycond),log(lindex))
# # scaled
# cor(data_x$bodycond,lindexS) # 0.2113215
# #
#########################################
### Ternary: levels curves for Scaled index
# install.packages("compositions")
library(compositions)
# ##
# # freq repeated
# (rep_count/apply(rep_count,1,sum))
fq<-c(1,1,1,3,3,1,3,3,2,1,3,3,3,3,1,1,1,3,2,1)
# ##
# 
niv=seq(0.1,0.9, by=0.1)
 niv<-c(niv,0.99)
 rad<-sqrt(-log(1-niv))
#############updated from Martin 3.26.19 
 plot(acomp(c(11,55,11)),labels=colnames(counts),col="black",pch="+")
 # plot fish: size proportional to frequency of fishes
 # cluster provided KIM in 2019 January
 #clus<-c(2,1,1,1,1,2,1,1,2,1,1,1,1,1,2,2,2,1,2,1)
 # cluster by email KIM
 # "I just split the data by body condition, so everything 
 # about 0.0138 is red and everything below that number 
 # is blue. I can provide the code as well if that helps.
clus<-ifelse(data_all$BodyCondition<0.0138,1,2)
 # clus is 1 2 2 1 2 1 1 2 2 1 2 1 1 2 1 2 1 2 1 2
 #######
 # cluster 1: "full circles"
plot(acomp(rep_count[clus==1,]),labels=colnames(counts),add=TRUE,pch=20, cex=1*fq[clus==1],col="blue")
 # cluster 2:  full  circles
plot(acomp(rep_count[clus==2,]),labels=colnames(counts),add=TRUE,pch=20, cex=1*fq[clus==2],col="red")
 # plot "*" whole geometric center
plot(acomp(g),labels=colnames(counts),add=TRUE,pch=21,col="black",cex=1.25)
 
 #data_all$BodyCondition
 
 for (i in 1:length(niv)){
   ellipses(acomp(c(11,55,11)),diag(rep(1,3)),r=rad[i],col="black",lty=3,lwd=1)
 }
 legend("topright", inset=c(0,0.1),legend=c("1","2","3"), 
        pch=21,pt.cex=c(0.5,1,1.5),
        xjust=0.5,title="fish number")
 #################################
 #amelia's addition
 #add bc values to counts dataset, visualize overlap
 counts$bc<-data_x$bodycond
 counts$index<-lindex
 counts$group<-clus
 new = formatC(lindex, digits = 3)
 plot(counts$index~jitter(counts$group,0.2))
 plot(counts$index~counts$group)
plot(index~bc, data=counts, col=group)
text(index~bc, data=counts, labels=new)
text(jitter(index,150)~bc, data=counts, labels=new)

 #####################
 
 # plot "+" TARGET
plot(acomp(c(11,55,11)),labels=colnames(counts),col="black",pch="+")
# plot fish: size proportional to frequency of fishes
plot(acomp(rep_count),labels=colnames(counts),add=TRUE,pch=20, cex=1*fq)
# plot "o" geometric center
plot(acomp(g),labels=colnames(counts),add=TRUE,pch="o",col="black",cex=1.25)

#data_all$BodyCondition

for (i in 1:length(niv)){
ellipses(acomp(c(11,55,11)),diag(rep(1,3)),r=rad[i],col="black",lty=3,lwd=2)
}
legend("topright", inset=c(0,0.2),legend=c("1","2","3"), 
       pch=20,pt.cex=c(0.5,1,1.5),xjust=0.5,title="fish number")
#######
# Ternary with clusters
######
niv=seq(0.1,0.9, by=0.1)
niv<-c(niv,0.99)
rad<-sqrt(-log(1-niv))
# plot "+" TARGET
plot(acomp(c(11,55,11)),labels=colnames(counts),col="black",pch="+")
# plot fish: size proportional to frequency of fishes
# cluster provided KIM in 2019 January
clus<-c(2,1,1,1,1,2,1,1,2,1,1,1,1,1,2,2,2,1,2,1)
#######
# cluster 1: "empty circles"
plot(acomp(rep_count[clus==1,]),labels=colnames(counts),add=TRUE,pch=21, cex=1*fq[clus==1])
# cluster 2:  full  circles
plot(acomp(rep_count[clus==2,]),labels=colnames(counts),add=TRUE,pch=20, cex=1*fq[clus==2])
# plot "*" whole geometric center
plot(acomp(g),labels=colnames(counts),add=TRUE,pch="*",col="black",cex=1.25)

#data_all$BodyCondition

for (i in 1:length(niv)){
  ellipses(acomp(c(11,55,11)),diag(rep(1,3)),r=rad[i],col="black",lty=3,lwd=2)
}
legend("topright", inset=c(0,0.2),legend=c("1","2","3"), 
       pch=21,pt.cex=c(0.5,1,1.5),xjust=0.5,title="fish number")
#################################
#################################
# rep_count and other to a CSV file
X=cbind(rep_count,data_all$BodyCondition)
X=cbind(X,X[,4]*100,clus)
colnames(X)=c(colnames(X)[1:3],"BodyCond","BodyCond100","Group")
X
write.table(X,"XLionFish.csv",sep=";")
####################################
# cluster results
####
# DATA: X = XLionFish.csv
#
colnames(X)
# [1] "N_chromis"   "N_goby"      "N_wrasse"    "BodyCond"    "BodyCond100" "Group"  
# geomean in counts and percentage
# group 1
(g1<-exp(apply(log(X[X[,6]==1,1:3]),2,mean)))
# N_chromis    N_goby  N_wrasse 
# 0.9948448 0.1931161 0.2516762
(g1/sum(g1)*100)
# N_chromis    N_goby  N_wrasse 
# 69.10386  13.41422  17.48192 
####
# group 2
(g2<-exp(apply(log(X[X[,6]==2,1:3]),2,mean)))
# N_chromis    N_goby  N_wrasse 
# 1.2688500 1.6358436 0.3688886
(g2/sum(g2)*100)
# N_chromis    N_goby  N_wrasse 
#  38.76029  49.97106  11.26865
#
# RATIO
(g1/sum(g1)*100)/(g2/sum(g2)*100)
# N_chromis    N_goby  N_wrasse 
# 1.7828520 0.2684398 1.5513761 
###
# T^2 : H_0 mu_1=mu_2
x<-filr(X[X[,6]==1,1:3])
y<-filr(X[X[,6]==2,1:3])
library(rrcov)
T2.test(x,y)
# Two-sample Hotelling test
# 
# data:  x and y
# T2 = 69.023, F = 32.594, df1 = 2, df2 = 17, p-value = 1.524e-06
# alternative hypothesis: true difference in mean vectors is not equal to (0,0)
# sample estimates:
#   Var1     Var2
# mean x-vector  1.1591569 0.452989
# mean y-vector -0.1796388 1.112391
# ========================================
# 
# BODY condition
hist(X[,4])
summary(X[,4])
# not skewed, not log is required
(b1<-mean(X[X[,6]==1,4]))
# 0.01447846
(b2<-mean(X[X[,6]==2,4]))
# 0.01245429
t.test(X[X[,6]==1,4],X[X[,6]==2,4])
# Welch Two Sample t-test
# 
# data:  X[X[, 6] == 1, 4] and X[X[, 6] == 2, 4]
# t = 2.9424, df = 13.411, p-value = 0.01112
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.0005426049 0.0035057467
# sample estimates:
#   mean of x  mean of y 
# 0.01447846 0.01245429 
# =========================
# repeat the same for index of preference
hist(lindexS)
summary(lindexS)
# not left skewed, not log is required
(l1<-mean(lindexS[X[,6]==1]))
# 0.9988722
(l2<-mean(lindexS[X[,6]==2]))
# 0.8266204
t.test(lindexS[X[,6]==1],lindexS[X[,6]==2])
# Welch Two Sample t-test
# 
# data:  lindexS[X[, 6] == 1] and lindexS[X[, 6] == 2]
# t = 2.2166, df = 6.0001, p-value = 0.06852
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.01789552  0.36239904
# sample estimates:
#   mean of x mean of y 
# 0.9988722 0.8266204 
t.test(lindexS[X[,6]==1],lindexS[X[,6]==2],var.equal=TRUE)
# Two Sample t-test
# 
# data:  lindexS[X[, 6] == 1] and lindexS[X[, 6] == 2]
# t = 3.0953, df = 18, p-value = 0.006243
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.05533737 0.28916615
# sample estimates:
#   mean of x mean of y 
# 0.9988722 0.8266204 
summary(cbind(lindexS[X[,6]==1],lindexS[X[,6]==2]))
var.test(lindexS[X[,6]==1],lindexS[X[,6]==2])
# F test to compare two variances
# 
# data:  lindexS[X[, 6] == 1] and lindexS[X[, 6] == 2]
# F = 8.2255e-06, num df = 12, denom df = 6, p-value < 2.2e-16
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   1.532815e-06 3.066691e-05
# sample estimates:
#   ratio of variances 
# 8.225457e-06 
wilcox.test(lindexS[X[,6]==1],lindexS[X[,6]==2])
# Wilcoxon rank sum test with continuity correction
# 
# data:  lindexS[X[, 6] == 1] and lindexS[X[, 6] == 2]
# W = 91, p-value = 0.000339
# alternative hypothesis: true location shift is not equal to 0


