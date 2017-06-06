#### Clear Workspace, Load packages ####
rm(list = ls())
library(lme4)
library(lmerTest) # will enhance the lme4 summary by df's and p-values
library(pbkrtest) # required by lmerTest for Kenward-Roger's
library(ggplot2)
library(plyr)
library(lsmeans)
library(boot)
library(lsr)
library(pwr)
library(effsize)
library(perm)
library(lmPerm)
library(car)
library(BayesFactor)
#### Define functions ####
# Bootstrapping functions for graphs
# (arguments: data(x) and cases (d))
resamples=10000
# Bootstrapped mean and confidence intervals for plots
bmean <- function(y, d) {
  return(mean(y[d],na.rm=TRUE))
}
booty_lo <- function(y,d) {
  if(length(unique(y)) > 1) {
    b <- boot(data= y[d], statistic= bmean, R= resamples)
    bci <- boot.ci(b,conf = 0.95,type="bca")
    ci<-bci$bca[4]
  }
  else 
  {ci<- bmean(y)
  }
  return(ci)
}
booty_hi <- function(y,d) {
  if(length(unique(y)) > 1) {
    b <- boot(data= y[d], statistic= bmean, R= resamples)
    bci <- boot.ci(b,conf = 0.95,type="bca")
    ci<-bci$bca[5]
  }
  else 
  {ci<- bmean(y)
  }
  return(ci)
}


stats_cont <- function(y,group) {
  #overall mean
  o_mean=signif(bmean(y), digits = 3)
  o_cilo=signif(booty_lo(y), digits = 3)
  o_cihi=signif(booty_hi(y), digits = 3)
  #group stats
  gr_mean=signif(by(y,group,bmean), digits = 3)
  gr_cilo=signif(by(y,group,booty_lo), digits = 3)
  gr_cihi=signif(by(y,group,booty_hi), digits = 3)
  #diff stats
  stat=permTS(y~group,
               alternative = "two.sided",
               method="exact.mc",
               control=permControl(nmc=resamples))
  g=cohen.d(y~group,hedges.correction=TRUE)
  return(list(cbind(o_mean,o_cilo,o_cihi),
              cbind(gr_mean,gr_cilo,gr_cihi),
              stat,
              g))
}


stats_cont_paired <- function(y1,y2) {
  #mean diff
  o_mean_diff=signif(bmean(y1-y2), digits = 3)
  o_cilo_diff=signif(booty_lo(y1-y2), digits = 3)
  o_cihi_diff=signif(booty_hi(y1-y2), digits = 3)

  #diff stats
  stat=permTS(y1,y2,
              paired=TRUE,
              alternative = "two.sided",
              method="exact.mc",
              control=permControl(nmc=resamples))
  g=cohen.d(y1,y2,paired=TRUE,hedges.correction=TRUE)
  return(list(cbind(o_mean_diff,o_cilo_diff,o_cihi_diff),
              stat,
              g))
}
#### #Basic Graphing settings (colours and stufff) ###

####LOAD AND DEFINE DATA#####
# Load data
setwd("/Users/matthiaszunhammer/Dropbox/LDOPA/Paper/L_DOPA_Placebo_Analysis/")
df <- read.table("Data_LDOPA_Long.dat",header=TRUE,sep="\t")

df <- within(df, {
  ids <- factor(ids)
  day <- factor(day, levels = 1:2, labels = c("1","2"))
  testseq <- factor(testseq, levels = 1:2, labels = c("1", "2"))
  placeboFirst  <- factor(placeboFirst, levels = 0:1, labels = c("controlf", "placebof"))
  pla        <- factor(wirk, levels = 0:1, labels = c("con", "pla"))
  male      <- factor(male, levels = 0:1, labels = c("female", "male"))
  ldopa <- factor(ldopa, levels = 0:1, labels = c("Placebo", "Levodopa"))
  #edu <- factor(edu, levels = 0:4, labels = c("none", "basic","middle","high","uni"))
})

# Replace NaN's by NA... otherwise problems.
df[is.na(df)]=NA 

# Create correct medication guesses (guess_med)
df$guess_med=0
df$guess_med[df$subj_LDOPA=='0'&df$ldopa=='Placebo']=1
df$guess_med[df$subj_LDOPA=='1'&df$ldopa=='Levodopa']=1

# Create overall side effects score
df$side_effects_all<-df$side_effects_d1+df$side_effects_d2

# Create overall real HPT
df$HPT_d1<-rowMeans(df[,grep("^HPT._...._pre_d1", colnames(df))],na.rm=1)
df$HPT_d2<-rowMeans(df[,grep("^HPT._...._pre_d2", colnames(df))],na.rm=1)

# Create overall temp
i_temp1=which(colnames(df)=="temps_.1")
df$temp=rowMeans(df[,i_temp1:(i_temp1+14)])

# # To exclude participants with low data quality.
# #df=df[!(df$unclear==1&df$male==1),]
# # z-Transformed versions of continuous data
# df$z_rating=as.numeric(scale(df$rating, center = TRUE, scale = TRUE))
# df$z_temp=as.numeric(scale(df$temp, center = TRUE, scale = TRUE))
# df$z_erwa=as.numeric(scale(df$erwa, center = TRUE, scale = TRUE))
# df$z_dev_from_target_rating<-scale(df$dev_from_target_rating, center = TRUE, scale = TRUE)
# df$z_oxymethyldopa_conz<-scale(df$oxymethyldopa_conz, center = TRUE, scale = TRUE)
# 

# Create contrast data-set contol-placebo for all days
dfd<-df[df$pla=="pla",]
dfd$ratingdiff<-df$rating[df$pla=="con"]-df$rating[df$pla=="pla"]
dfd$tempdiff<-df$temp[df$pla=="con"]-df$temp[df$pla=="pla"]

dfd$expect<-df$erwa[df$pla=="con"] #Erwa/Effi ratings were done before/after con OR pla, depending on which was first
a=df$erwa[df$pla=="pla"]
dfd$expect[is.na(df$erwa[df$pla=="con"])]=a[is.na(df$erwa[df$pla=="con"])]

dfd$actual<-df$effi[df$pla=="con"] #Erwa/Effi ratings were done before/after con OR pla, depending on which was first
a=df$effi[df$pla=="pla"]
dfd$actual[is.na(df$effi[df$pla=="con"])]=a[is.na(df$effi[df$pla=="con"])]

# Create contrast data-set contol-placebo for days 2
df2<-df[df$day=="2",]
dfd2<-dfd[dfd$day=="2",]
dfd2$expect_d1<-dfd$expect[dfd$day=="1"]
dfd2$actual_d1<-dfd$actual[dfd$day=="1"]
dfd2$expect_d2<-dfd$expect[dfd$day=="2"]
dfd2$actual_d2<-dfd$actual[dfd$day=="2"]

dfd2$rating_d1_con<-df$rating[df$day=="1"&df$pla=="con"]
dfd2$rating_d1_pla<-df$rating[df$day=="1"&df$pla=="pla"]
dfd2$rating_d2_con<-df$rating[df$day=="2"&df$pla=="con"]
dfd2$rating_d2_pla<-df$rating[df$day=="2"&df$pla=="pla"]
dfd2$ratingdiff_d1<-dfd$ratingdiff[dfd$day=="1"]
dfd2$ratingdiff_d2<-dfd$ratingdiff[dfd$day=="2"]

dfd2$temp_d1_con<-df$temp[df$day=="1"&df$pla=="con"]
dfd2$temp_d1_pla<-df$temp[df$day=="1"&df$pla=="pla"]
dfd2$temp_d2_con<-df$temp[df$day=="2"&df$pla=="con"]
dfd2$temp_d2_pla<-df$temp[df$day=="2"&df$pla=="pla"]
dfd2$tempdiff_d1<-dfd$tempdiff[dfd$day=="1"]
dfd2$tempdiff_d2<-dfd$tempdiff[dfd$day=="2"]

dfd2$z_ratingdiff_d2<-scale(dfd2$ratingdiff_d2, center = TRUE, scale = TRUE)

# Create LONG Dataset
#allratings<-df[,grep("^ratings_.\\d", colnames(df))]
dfl=reshape(df, 
            varying = list(
              grep("^ratings_.\\d", colnames(df)),
              grep("^duration_.\\d", colnames(df)),
              grep("^temps_.\\d", colnames(df))),
            v.names = c("ratingsall","durationsall","tempsall"),
            timevar = "stimrep", 
            times = 1:15, 
            #new.row.names = 1:1000,
            direction = "long")

dfl2=dfl[dfl$day=="2",]
dfld2=dfl[dfl$day=="2"&dfl$pla=="pla",]
dfld2$rating_d2_pla=dfl$ratingsall[dfl$day=="2"&dfl$pla=="pla"]
dfld2$rating_d2_con=dfl$ratingsall[dfl$day=="2"&dfl$pla=="con"]
dfld2$ratingsalldiff=dfl$ratingsall[dfl$day=="2"&dfl$pla=="con"]-dfl$ratingsall[dfl$day=="2"&dfl$pla=="pla"]
dfld2$durationsalldiff=dfl$durationsall[dfl$day=="2"&dfl$pla=="con"]-dfl$durationsall[dfl$day=="2"&dfl$pla=="pla"]
dfld2$tempssalldiff=dfl$tempsall[dfl$day=="2"&dfl$pla=="con"]-dfl$tempsall[dfl$day=="2"&dfl$pla=="pla"]


dfd2$t_med_to_bloodsample=(dfd2$time_d1_end-dfd2$time_med)*24*60


#### DESCRIPTIVE results, Paragraph1, TABLE S1: ####
# N
by(dfd2$ids,dfd2$ldopa,length)
#by(dfd2$unclear,dfd2$ldopa,sum)

# DEMOGRAFIC
# age
#stats_cont(dfd2$age,dfd2$ldopa)
signif(bmean(dfd2$age),digits=3)# for age min,max is preferred to 95%CIs
signif(min(dfd2$age),digits=3)
signif(max(dfd2$age),digits=3)
signif(by(dfd2$age,dfd2$ldopa,bmean),digits=3) 
signif(by(dfd2$age,dfd2$ldopa,min),digits=3) 
signif(by(dfd2$age,dfd2$ldopa,max),digits=3)

# gender
sex=table(dfd2$male,dfd2$ldopa)
signif(sex[3]/(sex[1]+sex[3]),digits=3) #prop male placebo
signif(sex[4]/(sex[2]+sex[4]),digits=3) #prop male placebo
#chisq.test(sex)

# edu
count(dfd2$edu)
edu=table(dfd2$edu,dfd2$ldopa)
edu
#chisq.test(edu)

# CES-D (german ads-k)
stats_cont(dfd2$CES_D,dfd2$ldopa)

# Placebo First
count(dfd2$placeboFirst)
pla1st=table(dfd2$placeboFirst,dfd2$ldopa)
prop.table(pla1st,2)
chisq.test(pla1st)

# Thermode Location
count(dfd2$loc_wirk_d2)
loc_D2=table(dfd2$loc_wirk_d2,dfd2$ldopa)
loc_D2
prop.table(loc_D2,2)
chisq.test(loc_D2)

# Temperatures
# T40 Day1
stats_cont(dfd2$T40_d1,dfd2$ldopa)
# T60 Day1
#stats_cont(dfd2$T60_d1,dfd2$ldopa)
# T80 Day1
stats_cont(dfd2$T80_d1,dfd2$ldopa)
# T60 Day2
stats_cont(dfd2$T60_d2,dfd2$ldopa)

# SIDE EFFECTS & QUESTIONNAIRES
# side effects score
stats_cont(dfd2$side_effects_all,dfd2$ldopa)
bf0=ttestBF(formula=side_effects_all~ldopa,data=dfd2)
bf0
1/bf0


# medication guess % LDOPA
sum(dfd2$subj_LDOPA)/length(dfd2$subj_LDOPA)
guess_med_l=table(dfd2$subj_LDOPA,dfd2$ldopa)
prop.table(guess_med_l,2)
chisq.test(guess_med_l)
permTS(dfd2$subj_LDOPA~dfd2$ldopa,
            alternative = "two.sided",
            method="exact.mc",
            control=permControl(nmc=resamples))

# medication guess % correct
sum(dfd2$guess_med)/length(dfd2$guess_med)
guess_med_corr=table(dfd2$guess_med,dfd2$ldopa)
prop.table(guess_med_corr,1)
chisq.test(guess_med_corr)

#POMS
# stats_cont(dfd2$POMS_Tension_D1,dfd2$ldopa)
# stats_cont(dfd2$POMS_Depression_D1,dfd2$ldopa)
# stats_cont(dfd2$POMS_Anger_D1,dfd2$ldopa)
# stats_cont(dfd2$POMS_Vigor_D1,dfd2$ldopa)
# stats_cont(dfd2$POMS_Fatigue_D1,dfd2$ldopa)
# stats_cont(dfd2$POMS_Confusion_D1,dfd2$ldopa)
# stats_cont(dfd2$POMS_Total_D1,dfd2$ldopa)
# 
# stats_cont(dfd2$POMS_Tension_D2,dfd2$ldopa)
# stats_cont(dfd2$POMS_Depression_D2,dfd2$ldopa)
# stats_cont(dfd2$POMS_Anger_D2,dfd2$ldopa)
# stats_cont(dfd2$POMS_Vigor_D2,dfd2$ldopa)
# stats_cont(dfd2$POMS_Fatigue_D2,dfd2$ldopa)
# stats_cont(dfd2$POMS_Confusion_D2,dfd2$ldopa)
# stats_cont(dfd2$POMS_Total_D2,dfd2$ldopa)

#### DAY 1 conditioning, Results Paragraph 2, TABLE S2: ####
# DAY 1 RATING RESULTS
# Pain Ratings Day 1
stats_cont(dfd2$ratingdiff_d1,dfd2$ldopa) # Conditioning strength
stats_cont(dfd2$rating_d1_con,dfd2$ldopa)
stats_cont(dfd2$rating_d1_pla,dfd2$ldopa)

stats_cont(dfd2$tempdiff_d1,dfd2$ldopa) # Conditioning strength
stats_cont(dfd2$temp_d1_con,dfd2$ldopa)
stats_cont(dfd2$temp_d1_pla,dfd2$ldopa)

sd(dfd2$ratingdiff_d1)
sd(dfd2$rating_d1_con)
sd(dfd2$rating_d1_pla)

# Expectation and actuation ratings Day1
stats_cont(dfd2$expect_d1,dfd2$ldopa)

stats_cont_paired(dfd2$expect_d1,dfd2$expect_d2)

# Actual Day1
stats_cont(dfd2$actual_d1,dfd2$ldopa)


#### DAY 2 Testing, Results Paragraph 3, TABLE S4: ####
# Pain Ratings Day 2 (within subject)
stats_cont_paired(dfd2$rating_d2_con,dfd2$rating_d2_pla)
# Pain ratings Day 2
stats_cont(dfd2$ratingdiff_d2,dfd2$ldopa) # Conditioning strength
stats_cont(dfd2$rating_d2_con,dfd2$ldopa)
stats_cont(dfd2$rating_d2_pla,dfd2$ldopa)
# Expectations and actuation Day2
stats_cont(dfd2$expect,dfd2$ldopa)
stats_cont(dfd2$actual,dfd2$ldopa)
stats_cont_paired(dfd2$actual_d1,dfd2$actual_d2) # paired change between days

#### MAIN RESULTS LEVODOPA Paragraph 4, Figures 2 and 3 ####

#medication blood levels
stats_cont(dfd2$t_med_to_bloodsample,dfd2$ldopa)
stats_cont(dfd2$ldopa_conz,dfd2$ldopa)
stats_cont(dfd2$oxymethyldopa_conz,dfd2$ldopa)
sd(dfd2$ldopa_conz)

prop_therapeutic=dfd2$ldopa_conz[dfd2$ldopa=="Levodopa"]>.2
prop_therapeutic=sum(prop_therapeutic)/length(prop_therapeutic)

#medication blood levels vs placebo
dfd2$z_ratingdiff_d2=as.numeric(scale(dfd2$ratingdiff_d2, center = TRUE, scale = TRUE))
dfd2$z_ldopa_conz=as.numeric(scale(dfd2$ldopa_conz, center = TRUE, scale = TRUE))

lm1<-lm(z_ratingdiff_d2~z_ldopa_conz,
        data=dfd2[dfd2$ldopa=="Levodopa",],na.action =na.exclude)
summary(lm1)

# MAIN RESULTS: Rating stats
stats_cont_paired(dfd2$rating_d2_con,dfd2$rating_d2_pla) # paired change between days
bf1=ttestBF(dfd2$rating_d2_con,dfd2$rating_d2_pla,paired=TRUE)
bf1
1/bf1

stats_cont(dfd2$ratingdiff_d2,dfd2$ldopa)
bf2=ttestBF(formula=ratingdiff_d2~ldopa,data=dfd2)
bf2
1/bf2

stats_cont(dfd2$rating_d2_con,dfd2$ldopa)
stats_cont(dfd2$rating_d2_pla,dfd2$ldopa)

#### Exploratory Paragraph ####
# Results sans suspected outliers
stats_cont_paired(dfd2$rating_d2_con[!dfd2$unclear],dfd2$rating_d2_pla[!dfd2$unclear]) # paired change between days
bf1_out=ttestBF(dfd2$rating_d2_con[!dfd2$unclear],dfd2$rating_d2_pla[!dfd2$unclear],paired=TRUE)
bf1_out
1/bf1_out

stats_cont(dfd2$ratingdiff_d2[!dfd2$unclear],dfd2$ldopa[!dfd2$unclear])
bf2_out=ttestBF(formula=ratingdiff_d2~ldopa,data=dfd2[!dfd2$unclear,])
bf2_out
1/bf2_out

# Results ANOVA gender
lm1<-lmp(ratingdiff_d2~(ldopa*male),
         perm="Exact",
         data=dfd2,na.action =na.exclude)
Anova(lm1,Type="II")
AIC(lm1)
summary(lm1)
stats_cont(dfd2$ratingdiff_d2[dfd2$male=='male'],dfd2$ldopa[dfd2$male=='male'])
stats_cont(dfd2$ratingdiff_d2[dfd2$male=='female'],dfd2$ldopa[dfd2$male=='female'])


# REPETITIONS
dfld2$z_stimrep=scale(dfld2$stimrep, center = TRUE, scale = TRUE)
dfld2$z_ratingsalldiff=scale(dfld2$ratingsalldiff, center = TRUE, scale = TRUE)

lmm1<-lmer(ratingsalldiff~dfld2$ldopa*z_stimrep+(1+z_stimrep| ids),
           data=dfld2,verbose=TRUE,
           control=lmerControl(optCtrl=list(maxfun=10000)),
           na.action =na.exclude)
Anova(lmm1,type='III')
summary(lmm1)
AIC(lmm1)
plot(dfld2$ratingsalldiff,predict(lmm1))

####################################
#### Placebo stats linear model ####
####################################
#lm1<-lmp(ratingdiff_d2~ldopa_conz+male,
#        data=dfd2,na.action =na.exclude)
lm1<-lmp(ratingdiff_d2~ldopa_conz+male+expect_d2+temp,
         data=dfd2,na.action =na.exclude)
AIC(lm1)
summary(lm1)
anova(lm1)
confint(lm1)
alias(lm1)
Anova(lm1,type='III')
etaSquared(lm1,type=2,anova=TRUE)

cor.test(dfd2$oxymethyldopa_conz[dfd2$ldopa=='Levodopa'],
         dfd2$ratingdiff_d2[dfd2$ldopa=='Levodopa'],
         method='pearson')

etaSquared(lm1,type=2,anova=TRUE)
# MISC
# HPT Day1
#stats_cont(dfd2$HPT_d1,dfd2$ldopa)
# HPT Day2
#stats_cont(dfd2$HPT_d2,dfd2$ldopa)
# alk
#stats_cont(dfd2$alk_p_w,dfd2$ldopa)
# food
#food=table(dfd2$food,dfd2$ldopa)
#prop.table(food)

# Graph L-DOPA VS PLACEBO
#histogram(~dfd2$ratingdiff_d2|dfd2$male:dfd2$ldopa)


####################################
#### Placebo stats linear mixed model ####
####################################

lmm1<-lmer(rating~pla+ldopa_conz+male+stimrep+(1+pla| ids),
        data=dfl2,
        na.action =na.exclude)
summary(lmm1)

#### Plots ####

# For sub-sample selection
#dfd2=dfd2[!(dfd2$ratingdiff_d2<0|dfd2$male=="male"),] #Excluded: Temperature data for day 3 were entered ambiguously/not to be found.

### MAIN RESULTS PLA VS PLA Belief
groupcolor=c(rgb(.5,.5,.5),rgb(.627,.706,1))
grouplabels=c("Placebo",
              "Levodopa")
dodge <- position_dodge(width=0.5)
jdodge <-position_jitterdodge(jitter.width=0.5,dodge.width=0.5)
g=ggplot(dfd2, aes(x=factor(subj_LDOPA),y=ratingdiff_d2, fill=ldopa))+ #,color=male,fill=male
  geom_hline(aes(yintercept=0),color="grey")+
  #geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
  geom_dotplot(position=dodge, binaxis='y',binwidth=1.8, stackdir='center', color=NA)+ #, binwidth=1 ,alpha=0.15
  geom_point(stat = "summary", fun.y=bmean,position=dodge)+
  geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
  labs(y ="Placebo Effect (points VAS) ± 95 % CI",x="Group")+
  scale_y_continuous(breaks=c(-30, -15, 0, 15, 30),limits = c(-30, 30))+
  scale_fill_manual("Group:\n",values =  groupcolor,
                      labels= grouplabels)+
  theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
  theme(axis.text=element_text(size=12), #fix axis font size
        axis.title=element_text(size=12))+ #fix axis title font size
  theme(strip.text.x = element_text(size = 12,face="bold"))+
  theme(axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_blank(), #element_line(colour = "#555555"), #fix grid display
        panel.grid.minor.y = element_blank(), #fix grid display
        panel.grid.major.x = element_blank(),
        axis.title.x= element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        panel.background = element_blank(),
        plot.title = element_text(size = 14,face="bold"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
  theme(legend.position=c(1.75,0.5),                 #place legend
        legend.direction="vertical",                    #stack format
        legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
        legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
        legend.title = element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12))+
  #ggtitle("Placebo vs Group")+# Add inital of investigators first name as graph title
  coord_fixed(ratio=0.20)
g
ggsave(file="/Users/matthiaszunhammer/Dropbox/LDOPA/Paper/Figure_1.eps",g,width = 6,height = 6,units = c("cm"),dpi = 600,scale=2.25)

### MAIN RESULTS PLA VS Levodopa Concentrations
grouplabels=c("Placebo",
              "L-Dopa")
g=ggplot(dfd2, aes(x=ldopa_conz,y=ratingdiff_d2, fill=ldopa))+ #,color=male,fill=male
  geom_hline(aes(yintercept=0),color="grey")+
  #geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
  geom_point()+ #, binwidth=1 ,alpha=0.15
  geom_smooth(method = "lm", se=TRUE, formula = y ~ x,alpha=0.2)+
  
  labs(y ="Placebo Effect (points VAS) ± 95 % CI",x="Oxymethyl-dopa Conzentration")+
  scale_y_continuous(breaks=c(-30, -15, 0, 15, 30),limits = c(-30, 30))+
  scale_fill_manual("Group:\n",values =  groupcolor,
                    labels= grouplabels)+
  theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
  theme(axis.text=element_text(size=12), #fix axis font size
        axis.title=element_text(size=12))+ #fix axis title font size
  theme(strip.text.x = element_text(size = 12,face="bold"))+
  theme(axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_blank(), #element_line(colour = "#555555"), #fix grid display
        panel.grid.minor.y = element_blank(), #fix grid display
        panel.grid.major.x = element_blank(),
        axis.title.x= element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        panel.background = element_blank(),
        plot.title = element_text(size = 14,face="bold"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
  theme(legend.position="bottom",#c(1.75,0.5),                 #place legend
        legend.direction="horizontal",                    #stack format
        legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
        #legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
        legend.title = element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12))+
  #ggtitle("Placebo vs Group")+# Add inital of investigators first name as graph title
  coord_fixed(ratio=0.20)
g

### MAIN RESULTS PLA VS CON

dotcolor=c(rgb(.8,.8,.8),rgb(.627,.706,1))
collabels=c("No Levodopa",
            "Levodopa")
statcolor=c(rgb(.4,.4,.4),rgb(.392,.443,.910))

xlabels=c("Control Site",
              "Placebo Site")
dodge <- position_dodge(width=0.5)
jdodge <-position_jitterdodge(jitter.width=0.5,dodge.width=0.5)
g2=ggplot(df2, aes(x=ldopa,y=rating, color= factor(wirk),fill= factor(wirk)))+ #,color=male,fill=male
  geom_hline(aes(yintercept=0),color="grey")+
  #geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
  geom_dotplot(position=dodge, binaxis='y',binwidth=1.8, stackdir='center',color=NA)+ #, binwidth=1 ,alpha=0.15
  #geom_line(aes(group=ids),color=rgb(.8,.8,.8))+ #, binwidth=1 ,alpha=0.15
  geom_point(stat = "summary", fun.y=bmean,position=dodge)+
  geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
  labs(y ="Pain Ratings (points VAS) ± 95 % CI",x="Group")+
  scale_y_continuous(breaks=c(0, 25, 50, 75, 100),limits = c(0, 100))+
  scale_fill_manual("",values =  dotcolor,
                    labels= xlabels)+#"Group:\n",
  scale_color_manual("",values =  statcolor,
                     labels= xlabels)+
  scale_x_discrete(labels= collabels)+
  theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
  theme(axis.text=element_text(size=12), #fix axis font size
        axis.title=element_text(size=12))+ #fix axis title font size
  theme(strip.text.x = element_text(size = 12,face="bold"))+
  theme(axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_blank(), #element_line(colour = "#555555"), #fix grid display
        panel.grid.minor.y = element_blank(), #fix grid display
        panel.grid.major.x = element_blank(),
        axis.title.x= element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        panel.background = element_blank(),
        plot.title = element_text(size = 14,face="bold"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
  theme(legend.position="bottom",#c(1.75,0.5),                 #place legend
        legend.direction="horizontal",                    #stack format
        legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
        legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
        legend.title = element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12))+
  #ggtitle("Placebo vs Group")+# Add inital of investigators first name as graph title
  coord_fixed(ratio=0.20)
g2
ggsave(file="/Users/matthiaszunhammer/Dropbox/LDOPA/Paper/Figure_1_ersatz.eps",g2,width = 6,height = 6,units = c("cm"),dpi = 600,scale=2.5)

### MAIN RESULTS SIDE EFFECTS
groupcolor=c(rgb(.5,.5,.5),rgb(.627,.706,1))
grouplabels=c("Placebo",
              "Levodopa")
dodge <- position_dodge(width=0.5)
jdodge <-position_jitterdodge(jitter.width=0.5,dodge.width=0.5)
g=ggplot(dfd2, aes(x=ldopa,y=POMS_Fatigue_D1, fill=ldopa))+ #,color=male,fill=male
  geom_hline(aes(yintercept=0),color="grey")+
  #geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
  geom_dotplot(position=dodge, binaxis='y', stackdir='center', color=NA)+ #, binwidth=1 ,alpha=0.15
  geom_point(stat = "summary", fun.y=bmean,position=dodge)+
  geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
  labs(y ="POMS Fatigue (points) ± 95 % CI",x="Group")+
  #scale_y_continuous(breaks=c(0:5),limits = c(0, 5))+
  scale_fill_manual("Group:\n",values =  groupcolor,
                    labels= grouplabels)+
  theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
  theme(axis.text=element_text(size=12), #fix axis font size
        axis.title=element_text(size=12))+ #fix axis title font size
  theme(strip.text.x = element_text(size = 12,face="bold"))+
  theme(axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_blank(), #element_line(colour = "#555555"), #fix grid display
        panel.grid.minor.y = element_blank(), #fix grid display
        panel.grid.major.x = element_blank(),
        axis.title.x= element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        panel.background = element_blank(),
        plot.title = element_text(size = 14,face="bold"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
  theme(legend.position=c(1.75,0.5),                 #place legend
        legend.direction="vertical",                    #stack format
        legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
        legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
        legend.title = element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12))+
  #ggtitle("Placebo vs Group")+# Add inital of investigators first name as graph title
  coord_fixed(ratio=0.20)
g
ggsave(file="/Users/matthiaszunhammer/Dropbox/LDOPA/Paper/Figure_1_side_effects.eps",g,width = 6,height = 6,units = c("cm"),dpi = 600,scale=2.25)




##### ANALYSIS OF PLACEBO-EXTINCTION ######

# Create simple regression coefficients
dfd2$extinction=NA
uniids=unique(dfld2$ids)
for (i in 1:length(uniids)) {
  d=dfld2[dfld2$ids==uniids[i],]
  lm<-lm(data=d[3:nrow(d),], ratingsalldiff~stimrep, na.action =na.exclude)
  #dfd2$extinction[id]=as.numeric(lm$coefficients[2])
  dfd2$extinction[i]=as.numeric(lm$coefficients[2])
}

# Test if regression coefficients differ by ldopa
stats_cont(dfd2$extinction,dfd2$ldopa)
t.test (dfd2$extinction, mu=0)


### Plot Group Extinction ####
groupcolor=c(rgb(.5,.5,.5),rgb(.627,.706,1))
grouplabels=c("Placebo",
              "Levodopa")
dodge <- position_dodge(width=0.5)
jdodge <-position_jitterdodge(jitter.width=0.5,dodge.width=0.5)
c=ggplot(dfd2, aes(x=ldopa,y=extinction, fill=ldopa))+ #,color=male,fill=male
  geom_hline(aes(yintercept=0),color="grey")+
  #geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
  geom_dotplot(position=dodge, binaxis='y', stackdir='center', color=NA)+ #, binwidth=1 ,alpha=0.15
  geom_point(stat = "summary", fun.y=bmean,position=dodge)+
  geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
  labs(y ="Placebo Extinction (points VAS/repetition) ± 95 % CI",x="Group")+
  #scale_y_continuous(breaks=c(-30, -15, 0, 15, 30),limits = c(-30, 30))+
  scale_fill_manual("Group:\n",values =  groupcolor,
                    labels= grouplabels)+
  theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
  theme(axis.text=element_text(size=12), #fix axis font size
        axis.title=element_text(size=12))+ #fix axis title font size
  theme(strip.text.x = element_text(size = 12,face="bold"))+
  theme(axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major.y = element_blank(), #element_line(colour = "#555555"), #fix grid display
        panel.grid.minor.y = element_blank(), #fix grid display
        panel.grid.major.x = element_blank(),
        axis.title.x= element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'transparent'),
        panel.background = element_blank(),
        plot.title = element_text(size = 14,face="bold"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
  theme(legend.position=c(1.75,0.5),                 #place legend
        legend.direction="vertical",                    #stack format
        legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
        legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
        legend.title = element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12))+
  #ggtitle("Placebo vs Group")+# Add inital of investigators first name as graph title
  coord_fixed(ratio=0.20)
c



###### Linear mixed model analysis #####
dfld2$z_ratingsalldiff=(scale(dfld2$ratingsalldiff, center = TRUE, scale = TRUE))
dfld2$zstimrep=(scale(dfld2$stimrep, center = TRUE, scale = TRUE))

lmer1<-lmer(z_ratingsalldiff~ldopa+zstimrep+ldopa:zstimrep+(1+zstimrep|ids),
            data=dfld2,
            na.action =na.omit,
            REML=FALSE,
            control=lmerControl(optimizer="bobyqa",check.conv.singular="warning"))

AIC(lmer1)
summary(lmer1)
anova(lmer1)

#a=summary(lmer1)$coefficients # get coefficients
#Get predicted values
dfld2$FITz_ratingsalldiff=predict(lmer1)
#Get back-transform predicted values
dfld2$FITratingsalldiff=predict(lmer1)*attr(dfld2$z_ratingsalldiff,"scaled:scale")
                          +attr(dfld2$z_ratingsalldiff,"scaled:center")


### PLOT PLACEBO VS STIMULUS REPETITIONS
#dodge <- position_dodge(width=0.5)
#jdodge <-position_jitterdodge(jitter.width=0.5,dodge.width=0.5)
groupcolor=c(rgb(.5,.5,.5),rgb(.627,.706,1))
grouplabels=c("Placebo",
              "Levodopa")
pr=ggplot(dfld2, aes(x=stimrep,y=ratingsalldiff,color=factor(ldopa),fill=factor(ldopa)))+ #,color=male,fill=male
  #geom_jitter()+
  #geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
  #geom_hline(aes(yintercept=0),color="grey")+
  #geom_abline(intercept=a[1],slope=a[3])+
  geom_line(stat = "summary", fun.y=bmean)+
  geom_point(stat = "summary", fun.y=bmean)+
  #geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
  #Regression line for each participant
  #geom_smooth(aes(group=factor(ids)), method = "lm", se=FALSE, formula = y ~ x, linetype = "dashed")+
  #Group regression line
  geom_smooth(method = "lm", se=TRUE, formula = y ~ x, alpha=0.1)+
  #Group regression line from mixed model
  # geom_smooth(aes(y=FITratingsalldiff), method = "lm", se=TRUE, formula = y ~ x, alpha=0)+
  #Regression line for each participant from mixed model
  #geom_smooth(aes(y=FITratingsalldiff,group=factor(ids)), method = "lm", se=FALSE, formula = y ~ x)+
  #geom_dotplot(position=dodge, binaxis='y', stackdir='center',alpha=0.15)+ #, binwidth=1
  labs(y ="Placebo Effect (points VAS)",x="Stimulus repetition")+
  scale_y_continuous(breaks=c(-4, 0, 4, 8))+ #,
  scale_x_continuous(limits = c(1,15), breaks=c(1,8,15))+
  scale_fill_manual("Group: ",values =  groupcolor,
                    labels= grouplabels)+
  scale_color_manual("Group: ",values =  groupcolor,
                    labels= grouplabels)+
  theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
  theme(axis.text=element_text(size=12), #fix axis font size
        axis.title=element_text(size=12))+ #fix axis title font size
  theme(strip.text.x = element_text(size = 12,face="bold"))+
  theme(#axis.ticks =element_blank(),
        #axis.line = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        #panel.grid.major.y = element_blank(), #element_line(colour = "#555555"), #fix grid display
        #panel.grid.minor.y = element_blank(), #fix grid display
        #panel.grid.major.x = element_blank(),
        #axis.title.x= element_blank(),
        #panel.grid.minor = element_blank(),
        #strip.background = element_rect(fill = 'transparent'),
        panel.background = element_blank(),
        plot.title = element_text(size = 14,face="bold"),
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
  theme(legend.position="bottom",#c(1.75,0.5),                 #place legend
        legend.direction="horizontal",                   #stack format
        legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
        legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
        legend.title = element_blank(),#element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12))+
  ggtitle("Placebo vs Group")+# Add inital of investigators first name as graph title
  coord_fixed(ratio=0.20)

pr
ggsave(file="/Users/matthiaszunhammer/Dropbox/LDOPA/Paper/Figure_2.eps",pr,width = 6,height = 6,units = c("cm"),dpi = 600,scale=2.5)




### PLOT RATINGS VS STIMULUS REPETITIONS
groupcolor=c(rgb(.8,.8,.8),rgb(.4,.4,.4),rgb(.627,.706,1),rgb(.392,.443,.910))
grouplabels=c("NoDopa+Con",
              "NoDopa+Pla",
              "L-Dopa+Con",
              "L-Dopa+Pla")
tcr=ggplot(dfl2, aes(x=stimrep,y=ratingsall,color=factor(ldopa):factor(wirk),fill=factor(ldopa):factor(wirk)))+ #,color=male,fill=male
  #geom_jitter()+
  #geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
  #geom_hline(aes(yintercept=0),color="grey")+
  geom_point(stat = "summary", fun.y=bmean)+
  geom_line(stat = "summary", fun.y=bmean)+
  #geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
  #geom_smooth(aes(group=factor(ids)), method = "lm", se=FALSE, formula = y ~ x, color="grey")+
  #geom_smooth(method = "lm", se=TRUE, formula = y ~ x, alpha=0)+
  
  #geom_dotplot(position=dodge, binaxis='y', stackdir='center',alpha=0.15)+ #, binwidth=1
  labs(y ="Pain Ratings (points VAS)", x="Stimulus repetition")+
  scale_fill_manual(values =  groupcolor,labels= grouplabels)+
  scale_color_manual(values =  groupcolor,labels= grouplabels)+
  scale_y_continuous(breaks=c(50, 60, 70, 80))+ #,
  scale_x_continuous(limits = c(1,15), breaks=c(1,8,15))+
  theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
  theme(axis.text=element_text(size=12), #fix axis font size
        axis.title=element_text(size=12))+ #fix axis title font size
  theme(strip.text.x = element_text(size = 12,face="bold"))+
  theme(#axis.ticks =element_blank(),
    #axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #panel.grid.major.y = element_blank(), #element_line(colour = "#555555"), #fix grid display
    #panel.grid.minor.y = element_blank(), #fix grid display
    #panel.grid.major.x = element_blank(),
    #axis.title.x= element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_rect(fill = 'transparent'),
    panel.background = element_blank(),
    plot.title = element_text(size = 14,face="bold"),
    plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
  theme(#legend.position=c(1.75,0.5),                 #place legend
    legend.direction="vertical",                    #stack format
    legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
    legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
    legend.title =element_blank(),# element_text(size = 14,face="bold"),
    legend.text = element_text(size = 12))+
  ggtitle("Placebo vs Group")+# Add inital of investigators first name as graph title
  coord_fixed(ratio=0.20)

tcr
ggsave(file="/Users/matthiaszunhammer/Dropbox/LDOPA/Paper/Time_Course_Ratings_D2.eps",tcr,width = 6,height = 6,units = c("cm"),dpi = 600,scale=2.5)


### Plot Group Side-Effects ####
dodge <- position_dodge(width=0.5)
jitter <- position_jitter(height=0.2)
jdodge <-position_jitterdodge(jitter.width=0.5,dodge.width=0.5)
nw=ggplot(dfd2, aes(x=ratingdiff_d2,y=side_effects_all,color=factor(ldopa),fill=factor(ldopa)))+ #,color=male,fill=male
  geom_jitter(position=jitter)+
  #geom_violin(adjust = 0.5,position=dodge,alpha=0.2,trim = FALSE)+ #
  #geom_hline(aes(yintercept=0),color="grey")+
  #geom_point(stat = "summary", fun.y=bmean)+
  #geom_errorbar(stat = "summary", fun.y=bmean,fun.ymin=booty_lo,  fun.ymax=booty_hi,position=dodge, width=0.25)+
  #Regression line for each participant
  #geom_smooth(aes(group=factor(ids)), method = "lm", se=FALSE, formula = y ~ x, linetype = "dashed")+
  #Group regression line
  geom_smooth(method = "lm", se=1, formula = y ~ x)+
  #Group regression line from mixed model
  # geom_smooth(aes(y=FITratingsalldiff), method = "lm", se=TRUE, formula = y ~ x, alpha=0)+
  #Regression line for each participant from mixed model
  #geom_smooth(aes(y=FITratingsalldiff,group=factor(ids)), method = "lm", se=FALSE, formula = y ~ x)+
  #geom_dotplot(position=dodge, binaxis='y', stackdir='center',alpha=0.15)+ #, binwidth=1
  labs(x ="Placebo Effect (points VAS)",y="Side Effect Score")+
  #scale_y_continuous(breaks=c(-4, 0, 4, 8))+ #,
  #scale_x_continuous(limits = c(1,15), breaks=c(1,8,15))+
  scale_fill_manual("Group: ",values =  groupcolor,
                    labels= grouplabels)+
  scale_color_manual("Group: ",values =  groupcolor,
                     labels= grouplabels)+
  theme(aspect.ratio = 1/sqrt(2))+ #fix plot size
  theme(axis.text=element_text(size=12), #fix axis font size
        axis.title=element_text(size=12))+ #fix axis title font size
  theme(strip.text.x = element_text(size = 12,face="bold"))+
  theme(#axis.ticks =element_blank(),
    #axis.line = element_line(colour = "black"),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    #panel.grid.major.y = element_blank(), #element_line(colour = "#555555"), #fix grid display
    #panel.grid.minor.y = element_blank(), #fix grid display
    #panel.grid.major.x = element_blank(),
    #axis.title.x= element_blank(),
    #panel.grid.minor = element_blank(),
    #strip.background = element_rect(fill = 'transparent'),
    panel.background = element_blank(),
    plot.title = element_text(size = 14,face="bold"),
    plot.margin = unit(c(0.1,0.1,0.1,0.1),"null"))+
  theme(legend.position="bottom",#c(1.75,0.5),                 #place legend
        legend.direction="horizontal",                   #stack format
        legend.key=element_rect(fill = "transparent"),  #remove grey legend bg
        legend.key.size=unit(1.5,"cm"),                 #size of legend bobbles
        legend.title = element_blank(),#element_text(size = 14,face="bold"),
        legend.text = element_text(size = 12))+
  ggtitle("Placebo vs Group")+# Add inital of investigators first name as graph title
  coord_fixed(ratio=0.20)
nw

ggsave(file="/Users/matthiaszunhammer/Dropbox/LDOPA/Paper/Side_Effects_Placebo_Effect.pdf",nw,width = 6,height = 6,units = c("cm"),dpi = 600,scale=2.5)



