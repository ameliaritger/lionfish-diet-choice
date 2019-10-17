# Lion fish preference index
# june 2018
#setwd("~/Documents/1_RECERCA/3_ARTICLES_EN_PREPAR/article_Pierotti_Leon/Kimber/boxduke1314")
#setwd('C:/Users/Kimberly Bourne/Box Sync/Lionfish-Amelia')
setwd('E:/Box Sync/Lionfish-Amelia')
#
###############################################################
# save.image("lionfish.RData")
# load("lionfish.RData")
###############################################################

data_all <- read.csv('lionfishdata.csv')

data_x <-data.frame(Starve_time = data_all$Starvationtime_hours,
                    cloudcover = data_all$Cloudcover,sex = data_all$Sex,
                    lf_std_length = data_all$Lionfishstandardlength_cm, lf_wetwgt = data_all$LionfishWetWeight_g,
                    mooncycle = data_all$Mooncycle,depthcaptured = data_all$Depthcaptured_m,bodycond = data_all$BodyCondition,Lionfish_totlength = data_all$Lionfishtotallength_cm)
counts <- data.frame(N_chromis = data_all$Numberchromisconsumed, N_goby = data_all$Numbergobyconsumed,N_wrasse = data_all$Numberwrasseconsumed)
########################
# Use zCompositions to replace zeros using bayesian multiplicative methods
library('zCompositions')
c_tot <- counts$N_chromis + counts$N_wrasse + counts$N_goby
rep_p <- cmultRepl(counts,method = "SQ") # proportions data with zeros replaced
rep_count<-rep_p*c_tot

colnames(rep_count) <- colnames(rep_p) <- colnames(counts)

rep_biomass <- rep_count*matrix(rep(c(0.641,0.129,0.655),nrow = nrow(rep_count)))
###################################
# MAIL MICHELE:
### 11chromis, 55 gobies, 11wrasses
#############################
# T^2 hotelling
library(rrcov)
#
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

# AND A FEASIBLE LOG-RATIO COULD BE
tU3=as.matrix(fsbp2ilr(rbind(c(1,-1,-1),c(0,1,-1))))
U3=t(tU3)

ilr_rep_p<-log(as.matrix(rep_p))%*%U3 
ilr_ref<-log(as.vector(c(11,55,11)/sum(c(11,55,11))))%*%U3 
#T2.test(as.matrix(ilr_rep_p),mu=as.vector(ilr_ref))
# T2 = 80.475, F = 38.120, df1 = 2, df2 = 18, p-value = 3.383e-07
# 
#  In average, the preference of the data set is different from the NULL HIPOT
#  of NO-PREFERENCE, the (11,55,11) composition
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
  if (is.vector(x)) {
    d<-fdAit(x,xref)^2
  }
  else {
    d<-apply(x,1,fdAit,xref)^2
  }
  return(1-exp(-d))
}
#
###############################################
# no preference
fIndexPref(c(3,15,3),c(11,55,11))
fIndexPrefS(c(3,15,3),c(11,55,11))
# check for the first fish
fIndexPref(rep_count[1,],c(11,55,11))
fIndexPrefS(rep_count[1,],c(11,55,11))
# for all 20 fish with counts
(lindex_count <-fIndexPref(rep_count,c(11,55,11)))
(lindexS_count <-fIndexPrefS(rep_count,c(11,55,11)))
lindex_count <- as.data.frame(lindex_count)
lindexS_count <- as.data.frame(lindexS_count)

# for all 20 fish with biomass
(lindex_biom <-fIndexPref(rep_biomass,c(7.051,7.095,7.205)))
(lindexS_biom <-fIndexPrefS(rep_biomass,c(7.051,7.095,7.205)))
lindex_biom <- as.data.frame(lindex_biom)
lindexS_biom <- as.data.frame(lindexS_biom)
############################################
  # total for each fish
  totalN<-apply(rep_count,1,sum)
  #
  
  ##############################
  
  pvalIndex=pvalIndexS<-rep(0,length(totalN))
  # NUMBER of random samples for multinomial
  numboot=10000
  #####
  # loop for all fishes
  for (numfish in (1:length(totalN))){
  # check fish #2  numfish<-2
  lindfish<-lindex[numfish]
  lindSfish<-lindexS[numfish]
  set.seed(10)
  XRand<-t(rmultinom(numboot,totalN[numfish],prob=c(11,55,11)))
  # Zero replacement
  XRand_rep_p <- cmultRepl(XRand,method = "SQ") # proportions data with zeros replaced
  XRand_rep_count<-XRand_rep_p*totalN[numfish]
  # head(XRand_rep_count,20)
  #
  IndexNullHip<-fIndexPref(XRand_rep_count,c(11,55,11))
  IndexSNullHip<-fIndexPrefS(XRand_rep_count,c(11,55,11))
  # P-valor Index
  pvalIndex[numfish]<-(length(IndexNullHip[IndexNullHip>lindfish]))/length(IndexNullHip)
  #(length(IndexNullHip[IndexNullHip>lindfish])+0.5)/length(IndexNullHip)
  
  # P-valor [0,1]-index
  pvalIndexS[numfish]<-(length(IndexSNullHip[IndexSNullHip>lindSfish]))/length(IndexSNullHip)
  #(length(IndexSNullHip[IndexSNullHip>lindSfish])+0.5)/length(IndexSNullHip)
  ####
  } # END lop for all fishes
  (pvalIndex)
  # [1] 0.0568 0.0000 0.0000 0.0000 0.0000 0.0901 0.0000 0.0000 0.0923 0.0000 0.0000 0.0000 0.0000 0.0000
  # [15] 0.4866 1.0000 0.2898 0.0000 0.0923 0.0000
  (pvalIndexS)
  # [1] 0.0568 0.0000 0.0000 0.0000 0.0000 0.0901 0.0000 0.0000 0.0923 0.0000 0.0000 0.0000 0.0000 0.0000
  # [15] 0.4866 1.0000 0.2898 0.0000 0.0923 0.0000
  #############################################
  ###
  ####
  # CORRELATION BETWEEN bodycondition and index of preference
  # not scaled
  cor(data_x$bodycond,lindex) # 0.479191
  plot(data_x$bodycond,lindex)
  # WITH LOGs
  cor(log(data_x$bodycond),log(lindex)) # 0.3974944
  plot(log(data_x$bodycond),log(lindex))
  # scaled
  cor(data_x$bodycond,lindexS) # 0.2113215
  #
  #########################################
  #########################################
  ait <- numeric()
  for (i in 1:nrow(rep_count)){ait <- c(ait,fdAit(as.numeric(rep_count[i,]),c(11,55,11)))}
  #########################################
  #########################################
  ### levels curves for Scaled index
  library(compositions)
  ##
  # freq repeated
  (rep_count/apply(rep_count,1,sum))
  fq<-c(1,1,1,3,3,1,3,3,2,1,3,3,3,3,1,1,1,3,2,1)
  ##
  
  niv=seq(0.1,0.9, by=0.1)
  niv<-c(niv,0.99)
  rad<-sqrt(-log(1-niv))
  plot(acomp(c(11,55,11)),labels=colnames(counts),col="black",pch="+")
  plot(acomp(rep_count),labels=colnames(counts),add=TRUE,pch=20, cex=1*fq)
  
  for (i in 1:length(niv)){
    ellipses(acomp(c(11,55,11)),diag(rep(1,3)),r=rad[i],col="black",lty=3,lwd=2)
  }
  legend("topright", inset=c(0,0.2),legend=c("1","2","3"), 
         pch=20,pt.cex=c(0.5,1,1.5),xjust=0.5,title="fish number")
  #################################
  
  
  plot(data$bodycond,data$lindex, xlab = "Body Condition", ylab = "Selectivity Index")
  
  clusdata <- data.frame(Bodycondition = data_all$BodyCondition,Index = lindex_count)
  
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra) # clustering algorithms & visualization
  cfit <- kmeans(clusdata, 2)
  
  rownames(clusdata) <- data_all$Sex
  fviz_cluster(cfit, data = clusdata, xlab = "Body Condition",ylab ="Selectivity Index",main = "")
  
  set.seed(123)
  gap_stat <- clusGap(clusdata, FUN = kmeans, nstart = 25,
                      K.max = 10, B = 50)
  # Print the result
  print(gap_stat, method = "firstmax")
  
  fviz_gap_stat(gap_stat)
  ######################################
  
  library(rpart)
  treemodel <- rpart(lindex_count ~ Bodycondition,data = clusdata)
  library(rpart.plot)
  rpart.plot(treemodel,cex = 0.75,digits = 4)
  
  tree.predict <- predict(treemodel)
  plot(clusdata$Bodycondition[which(tree.predict < 4)],clusdata$lindex_count[which(tree.predict < 4)],col = 'blue',
       ylim = c(min(clusdata$lindex_count),max(clusdata$lindex_count)),xlim = c(min(clusdata$Bodycondition),max(clusdata$Bodycondition)),
       ylab = "Selectivity Index",xlab = "Body Condition")
  points(clusdata$Bodycondition[which(tree.predict > 4)],clusdata$lindex_count[which(tree.predict > 4)],col = 'red')
  
######################################
  ## For normalized (0,1) data
  reg_data <-  data.frame(Bodycondition = data_all$BodyCondition,Index = lindexS_count)
  
  breaks <- seq(0.011,0.017,by = 0.0001)
  mse <- numeric()
  for (b in breaks){
    piecewise <- lm(lindexS_count ~ (Bodycondition > b), data = reg_data)
    mse <- rbind(mse,c(b,as.numeric(sqrt(mean(summary(piecewise)$residuals^2))),as.numeric(summary(piecewise)$r.squared)))
  }
  
  b <- mse[which(mse[,2] == min(mse[,2])),1]
  b <- 0.0124
  
  piecewise1 <- lm(lindexS_count ~ Bodycondition, data = reg_data[which(reg_data$Bodycondition < b),])
  piecewise2 <- lm(lindexS_count ~ Bodycondition, data = reg_data[which(reg_data$Bodycondition > b),])
  pred1 <- predict(piecewise1,newdata = data.frame(Bodycondition = breaks[breaks < b]))
  pred2 <- predict(piecewise2,newdata = data.frame(Bodycondition = breaks[breaks > b]))
  
  
  plot(reg_data$Bodycondition[reg_data$Bodycondition < b ],reg_data$lindexS_count[reg_data$Bodycondition < b ],col = "blue",
       xlim=range(reg_data$Bodycondition),ylim = range(reg_data$lindexS_count), xlab = "Body condition", ylab = "Selectivity Index", bty = "n",pch = 16)
  points(reg_data$Bodycondition[reg_data$Bodycondition > b ],reg_data$lindexS_count[reg_data$Bodycondition > b ],col = "red",pch = 17)
  lines(breaks[breaks < b],pred1)
######################################
  ## For non-normalized (1,8) data
  
  normal <- function (d){
    return(1-exp(-d))
  }
  
  reg_data <-  data.frame(Bodycondition = data_all$BodyCondition,Index = lindex_count)
  
  breaks <- seq(0.011,0.017,by = 0.0001)
  mse <- numeric()
  for (b in breaks){
    piecewise <- lm(lindex_count ~ (Bodycondition > b), data = reg_data)
    mse <- rbind(mse,c(b,as.numeric(sqrt(mean(summary(piecewise)$residuals^2))),as.numeric(summary(piecewise)$r.squared)))
  }
  
  b <- mean(mse[which(mse[,2] == min(mse[,2])),1])
  #b <- 0.0138
  piecewise <- lm(lindex_count ~ (Bodycondition > b), data = reg_data)
  
  piecewise1 <- lm(lindex_count ~ Bodycondition, data = reg_data[which(reg_data$Bodycondition < b),])
  piecewise2 <- lm(lindex_count ~ Bodycondition, data = reg_data[which(reg_data$Bodycondition > b),])
  pred1 <- predict(piecewise1,newdata = data.frame(Bodycondition = breaks[breaks < b]))
  pred2 <- predict(piecewise2,newdata = data.frame(Bodycondition = breaks[breaks > b]))
  
  
  plot(reg_data$Bodycondition[reg_data$Bodycondition < b ],normal(reg_data$lindex_count[reg_data$Bodycondition < b ]),col = "blue",
       xlim=range(reg_data$Bodycondition),ylim = range(normal(reg_data$lindex_count)), xlab = "Body condition", ylab = "Selectivity Index", bty = "n",pch = 16)
  points(reg_data$Bodycondition[reg_data$Bodycondition > b ],normal(reg_data$lindex_count[reg_data$Bodycondition > b ]),col = "red",pch = 17)
  lines(breaks[breaks < b],normal(pred1),lwd =  2)
  lines(breaks[breaks > b],normal(pred2),lwd =  2)
  
  plot(reg_data$Bodycondition[reg_data$Bodycondition < b ],reg_data$lindex_count[reg_data$Bodycondition < b ],col = "blue",
       xlim=range(reg_data$Bodycondition),ylim = range(reg_data$lindex_count), xlab = "Body condition", ylab = "Selectivity Index", bty = "n",pch = 16)
  points(reg_data$Bodycondition[reg_data$Bodycondition > b ],reg_data$lindex_count[reg_data$Bodycondition > b ],col = "red",pch = 17)
  lines(breaks[breaks < b],pred1,lwd =  2)
  lines(breaks[breaks > b],pred2,lwd =  2)
  
  plot(reg_data$Bodycondition[reg_data$Bodycondition < b ],reg_data$lindex_count[reg_data$Bodycondition < b ],col = "blue",
       xlim=range(reg_data$Bodycondition),ylim = range(reg_data$lindex_count), xlab = "Body condition", ylab = "Selectivity Index", bty = "n",pch = 16)
  points(reg_data$Bodycondition[reg_data$Bodycondition > b ],reg_data$lindex_count[reg_data$Bodycondition > b ],col = "red",pch = 17)
  lines(breaks,predict(piecewise,newdata = data.frame(Bodycondition = breaks)))
  

###From Amelia
library(ggplot2)
library(tidyverse)
library(grid)
#create new column for each category
norm_data <-  data.frame(Bodycondition = data_all$BodyCondition,Index = lindexS_count)
norm_data$category <- ifelse(norm_data$Bodycondition<b, "A", "B") #or b is 0.0138

#plot
a <- ggplot(norm_data, aes(x=Bodycondition, y=lindexS_count, color=category)) +
  geom_hline(yintercept=0.99, linetype="dashed", color="#F0E442", size=1) + #add prey preference line
  annotate("rect", xmin = -Inf, xmax = 0.0184, ymin = 0.99, ymax = Inf, fill = "#F0E442", alpha = .3, color = NA) + #add prey preference olor
  geom_hline(yintercept=0.03, linetype="dashed", color="#E69F00", size=1) + #add no preference line
  annotate("rect", xmin = -Inf, xmax = 0.0184, ymin = -Inf, ymax = 0.03, fill = "#E69F00", alpha = .2, color = NA) + #add no preference color
  annotate("text", x = 0.0183, y = 0.01, hjust = -0.215, vjust = 0, label="no prey preference", color = "black", size=3.5) +
  annotate("text", x = 0.018, y = 1.0, hjust = -0.215, vjust = 0, label="strong prey preference", color = "black", size=3.5) +
  geom_point(aes(shape=category), size=3, alpha=0.9, show.legend=FALSE) + #plot points
  scale_color_manual(values=c("#0072B2", "#CC79A7")) + #manually change point colors
  labs(x="Body condition", y="Index of Selectivity") + #change axis labels
  geom_segment(x=0.011,xend=b,y=mean(pred1), yend=mean(pred1), colour="#0072B2") + #add group 1 mean line
  geom_segment(x=b,xend=0.018,y=mean(pred2), yend=mean(pred2), colour="#CC79A7") + #add group 2 mean line
  theme_classic() + #remove background crap
  theme(plot.margin=unit(c(1,3.6,0,0), "cm")) + #extend plot area to allow text
  #annotate("text", x=0.02, y=1.0, label="lionfish are not selective", size=2.5) +
  coord_cartesian(xlim = c(0.011, 0.018), clip="off") + #limit plot area
  scale_y_continuous(limits=c(0,1.02), #change min and max values on y axis
                     expand=c(0,0),
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0)) +
  scale_x_continuous(limits=c(0.011, 0.0185),
                     breaks=c(0.010,0.012,0.014,0.016,0.018)) #change min and max values on x axis
  #annotation_custom(grob = linesGrob(), xmin = 0.019, xmax = 0.019, ymin = 0.05, ymax = 0.95) +
  #annotation_custom(grob = linesGrob(), xmin = 0.0195, xmax = 0.019, ymin = 0.05, ymax = 0.15)

# Disable clip-area.
gt <- ggplot_gtable(ggplot_build(a))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
 
# alternative plot
c <- ggplot(norm_data, aes(x=Bodycondition, y=lindexS_count, color=category)) +
  geom_hline(yintercept=0.99, linetype="dashed", color="#F0E442", size=1) + #add prey preference line
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.99, ymax = Inf, fill = "#F0E442", alpha = .3, color = NA) + #add prey preference olor
  geom_hline(yintercept=0.03, linetype="dashed", color="#E69F00", size=1) + #add no preference line
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.03, fill = "#E69F00", alpha = .2, color = NA) + #add no preference color
  annotate("text", x = 0.0175, y = 0.065, label="no prey preference", color = "black", size=3.5) +
  annotate("text", x = 0.0175, y = 0.96, label="strong prey preference", color = "black", size=3.5) +
  geom_point(aes(shape=category), size=3, alpha=0.9, show.legend=FALSE) + #plot points
  scale_color_manual(values=c("#0072B2", "#CC79A7")) + #manually change point colors
  labs(x="Body condition", y="Index of Selectivity") + #change axis labels
  annotate("segment", x=0.0175,
           xend=0.0175,
           y=0.9,
           yend=0.1,
           color="black",
           size=1,
           arrow=arrow(length=unit(0.08,"npc"))) +
#geom_segment(x=0.011,xend=b,y=mean(pred1), yend=mean(pred1), colour="#0072B2") + #add group 1 mean line
  #geom_segment(x=b,xend=0.018,y=mean(pred2), yend=mean(pred2), colour="#CC79A7") + #add group 2 mean line
  theme_classic() + #remove background crap
  theme(plot.margin=unit(c(1,3.6,0,0), "cm")) + #extend plot area to allow text
  #annotate("text", x=0.02, y=1.0, label="lionfish are not selective", size=2.5) +
  scale_y_continuous(limits=c(0,1.02), #change min and max values on y axis
                     expand=c(0,0),
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0)) +
  scale_x_continuous(limits=c(0.011, 0.018),
                     breaks=c(0.010,0.012,0.014,0.016,0.018)) #change min and max values on x axis

# Disable clip-area.
gt <- ggplot_gtable(ggplot_build(c))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

