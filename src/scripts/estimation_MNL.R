# Script to file for estimation

# Install needed libraries
# install.packages("dplyr")
# install.packages("foreign")
# install.packages("readr")
# install.packages("mosaic")
# install.packages("lazyeval")
# install.packages("DataCombine")
# install.packages("mlogit")

# Load the needed libraries
library(dplyr)
library(foreign)
library(readr)
library(mosaic)
library(lazyeval)
library(DataCombine)
library(broom)
library(mlogit)
library(reshape2)

# Set the working directory
setwd("./Downloads/TRB_ML_Files/")

hh <- read_csv("./data/hhtype.csv")
per <- read_csv("./data/pertype.csv")
access <- read_csv("./data/accessibility.csv")
skims <- read_csv("./data//skims.csv")


# Function to construct utility equations for MLOGIT
CalcUtilityExpr <- function(utiltyExprFile){
  # construct utility equation
  utilFile <- read_csv(utiltyExprFile)
  util <- "choiceBoolean ~ "
  
  
  colIkey <- which( colnames(utilFile)=="ikey" )
  colScript <- which( colnames(utilFile)=="script" )
  colProtected <- which( colnames(utilFile)=="Protected" )
  colOffset <- which( colnames(utilFile)=="Offset" )
  colRef <- which( colnames(utilFile)=="Reference" )
  colStart <- which( colnames(utilFile)=="Start" )
  colEnd <- which( colnames(utilFile)=="End" )
  
  expr1 <- ""
  expr2 <- ""
  utilExpr <-""
  if (nrow(utilFile) > 0){
    for (i in 1:nrow(utilFile)){
      part <- utilFile[i,colIkey]  
      term <- utilFile[i,colScript]  
      prot <- utilFile[i,colProtected]  
      prot <- utilFile[i,colProtected] 
      offset <- utilFile[i,colOffset] 
      ref <- utilFile[i,colRef]
      start <- as.numeric(utilFile[i,colStart])
      end <- as.numeric(utilFile[i,colEnd])
      
      
      if(part == 1 ){
        if(prot == 1) {
          if(nchar(expr1)>0) {
            if(!is.na(offset)){
              expr1 <- paste(expr1," + ","offset(I(",offset,"*",term,"))",sep="") 
            }
            else {
              expr1 <- paste(expr1," + ","I(",term,")",sep="") 
            }
          }
          else
          {
            if(!is.na(offset)){
              expr1 <- paste(expr1,"offset(I(",offset,"*",term,"))",sep="") 
            }
            else{
              expr1 <- paste(expr1,"I(",term,")",sep="") 
            }
          }
        }
        else{
          if(nchar(expr1)>0) {
            if(!is.na(offset)){
              expr1 <- paste(expr1," + ","offset(",offset,"*",term,")",sep="") 
            }
            else {
              expr1 <- paste(expr1," + ",term,sep="") 
            }
          }
          else
          {
            if(!is.na(offset)) {
              expr1 <- paste(expr1,"offset(",offset,"*",term,")",sep="") 
            }
            else {
              expr1 <- paste(expr1,term,sep="") 
            }
            
          }
        }
        
      }
      
      if(part == 3 ){
        for (m in seq(start, end)){
          if( m != ref){
            alt_exp <- paste0("*(altid == ", (m+1),")")
            if(prot == 1) {
              if(nchar(expr1)>0) {
                if(!is.na(offset)){
                  expr1 <- paste(expr1," + ","offset(I(",offset,"*",term, alt_exp, "))",sep="") 
                }
                else {
                  expr1 <- paste(expr1," + ","I(",term,alt_exp,")",sep="") 
                }
              }
              else
              {
                if(!is.na(offset)){
                  expr1 <- paste(expr1,"offset(I(",offset,"*",term,alt_exp,"))",sep="") 
                }
                else{
                  expr1 <- paste(expr1,"I(",term,alt_exp,")",sep="") 
                }
              }
            }
            else{
              if(nchar(expr1)>0) {
                if(!is.na(offset)){
                  expr1 <- paste(expr1," + ","offset(",offset,"*",term,alt_exp,")",sep="") 
                }
                else {
                  expr1 <- paste(expr1," + ",term,alt_exp, sep="") 
                }
              }
              else
              {
                if(!is.na(offset)) {
                  expr1 <- paste(expr1,"offset(",offset,"*",term,alt_exp, ")",sep="") 
                }
                else {
                  expr1 <- paste(expr1,term,alt_exp, sep="") 
                }
                
              }
            }
          }
        }
        
      }
      
      if(part == 2 ){
        if(prot == 1) {
          if(nchar(expr2)>0) {
            if(!is.na(offset)) {
              expr2 <- paste(expr2," + ","offset(I(",offset,"*",term,"))",sep="") 
            }
            else {
              expr2 <- paste(expr2," + ","I(",term,")",sep="") 
            }
          }
          else
          {
            if(!is.na(offset)) {
              expr2 <- paste(expr2,"offset(I(",offset,"*",term,"))",sep="") 
            }
            else {
              expr2 <- paste(expr2,"I(",term,")",sep="") 
            }
            
          }
        }
        else{
          if(nchar(expr1)>0) {
            if(!is.na(offset)) {
              expr2 <- paste(expr2," + offset(",offset,"*",term,")",sep="") 
            }
            else {
              expr2 <- paste(expr2," + ",term,sep="") 
            }
            
          }
          else
          {
            if(!is.na(offset)) {
              expr2 <- paste(expr2,"offset(",offset,"*",term,")",sep="") 
            }
            else {
              expr2 <- paste(expr2,term,sep="") 
            }
            
          }
        }
        
        
      }
      
      
      
    }
    
    if(nchar(expr1)>0 & nchar(expr2)>0) {
      utilExpr <- paste0(util,expr1,"|",expr2)
    } else if(nchar(expr1)>0 & nchar(expr2)==0) {
      utilExpr <- paste0(util,expr1)
    } else if(nchar(expr1)==0 & nchar(expr2)>0) {
      utilExpr <- paste0(util,expr2)
    }
  }
  
  if(nchar(expr1)>0 & nchar(expr2)>0) {
    utilExpr <- paste0(util,expr1,"|",expr2)
  } else if(nchar(expr1)>0 & nchar(expr2)==0) {
    utilExpr <- paste0(util,expr1 ,"|")
  } else if(nchar(expr1)==0 & nchar(expr2)>0) {
    utilExpr <- paste0(util,"0|",expr2)
  }
  
  
  #utilExpr <- paste0(utilExpr,"| 1")
  return (utilExpr)
}


# Variable transformation

hh <- hh %>%
  arrange(sampn)

set.seed(11111)

# random number for sampling uniformly
hh$rand_unif <- runif(nrow(hh), min = 0, max = 1)
hh$rand_bias <- runif(nrow(hh), min = 0, max = 1)
hh$rand_bias_inc <- runif(nrow(hh), min = 0, max = 1)



hh <- hh %>%
  dplyr::mutate(weights_random = ifelse(rand_unif <=0.5, 1, 0),
                weights_bias = ifelse(hhveh == 0 & rand_bias <=0.7, 1,
                                      ifelse(hhveh == 1 & rand_bias <=0.45, 1, 
                                             ifelse(hhveh == 2 & rand_bias <= 0.29, 1, 
                                                    ifelse(hhveh == 3 & rand_bias <= 0.2, 1, 0)))),
                weights_bias_inc = ifelse(incom == 1 & rand_bias_inc <=0.8, 1,
                                      ifelse(incom == 2 & rand_bias_inc <=0.65, 1, 
                                             ifelse(incom == 3 & rand_bias_inc <= 0.45, 1, 
                                                    ifelse(incom == 4 & rand_bias_inc <= 0.3, 1,
                                                           ifelse(incom == 5 & rand_bias_inc <= 0.2, 1, 
                                                                  ifelse(incom >5  & rand_bias_inc <= 0.5, 1,0)))))))


hh <- hh %>%
  mutate(hhveh = ifelse(hhveh >=3, 3, hhveh))


hh_long <- hh[rep(seq_len(nrow(hh)), 4),]

hh_long <- hh_long %>% 
  dplyr::group_by(sampn) %>% 
  dplyr::mutate(altid = sequence(n())) %>%
  dplyr::arrange(sampn, altid)

hh_long <- hh_long %>% 
  dplyr::mutate(choiceBoolean = ifelse(hhveh == altid-1, 1, 0) )

#hh_long <- melt(hh, id.vars=c("sampn"))

hh_long$altid <- as.factor(hh_long$altid)

per$wktaz  <- as.numeric(per$wktaz)
per <- left_join(per, skims, by=c("hhtaz"="ORIG", "wktaz"="DEST"))

per <- per %>%
  dplyr::mutate(auto_dep = ifelse(wktaz>0 & BUSIVT>0, (pmin(BUSIVT, ifelse(AUTODIST<=3, AUTODIST/4.82,9999))-AUTOTIME)/45, 1))

per <- per %>%
  dplyr::mutate(auto_dep = ifelse(wktaz <=0, 0, ifelse(auto_dep > 1, 1, ifelse(BUSIVT == 0, 1, auto_dep))))


### variable transformation
num_adults_df <- per %>% 
  dplyr::group_by(sampn ) %>%
  dplyr::summarise(num_adults = sum(ifelse(age >=18, 1, 0)),
                   hh_auto_dep = sum(auto_dep))

hh_long <- left_join(hh_long, num_adults_df, by="sampn")

# add accessibilities
hh_long <- left_join(hh_long, access[,c("TAZ_NUM","ACCESS4","ACCESS5","ACCESS6")],by=c("hhtaz"="TAZ_NUM"))

hh_long[is.na(hh_long)] <- 0

hh_long <- hh_long %>%
  dplyr::mutate(auto_access = log(1+ACCESS4),
                tran_access = log(1+ACCESS5),
                non_mot_access = log(1+ACCESS6))


utiltyExprFile <- "util_Terms_com_R.csv"

utilExpr <- CalcUtilityExpr(utiltyExprFile)



## ALL sample
# convert into mlogit data
mlogitDataAvail <- mlogit.data(hh_long,
                               shape="long",
                               chid.var = c("sampn"),
                               alt.var = "altid",
                               choice = "choiceBoolean",
                               labels=c("0 Car","1 Car","2 Cars","3+ Cars"))


#res <- mlogit(as.formula(choiceBoolean ~ 0|1), data=mlogitDataAvail, nests=list(zero_car=c('1') ,one_plus=c('2','3','4')), unscaled=FALSE, reflevel = '2',print.level=2, un.nest.el = TRUE)

# Write the model data to file
write.csv(x=hh_long,
          file="./data/mnl_training_data_long_all_sample.csv",
          sep=",",
          row.names=FALSE,
          col.names=TRUE)

# View the formula for the model
print(as.formula(utilExpr))

res <- mlogit(as.formula(utilExpr),
              data=mlogitDataAvail,
              panel=FALSE,
              reflevel = '1',
              print.level=2)

Table_res <- data.frame(summary(res)$CoefTable)
write.csv(Table_res[c(1:4)],file="est_results.csv")


b <- data.frame(summary(res)$logLik)
write.csv(b[c(1:1)],file="logLik.csv")

p <- data.frame(fitted(res,  type = c("probabilities")))

pred_prob <- bind_cols(hh, p)

write_csv(pred_prob, "mnl_prob_all.csv")





## 50-50 random
mlogitDataAvail <- mlogit.data(hh_long,
                               shape="long",
                               chid.var = c("sampn"),
                               alt.var = "altid",
                               choice = "choiceBoolean",
                               labels=c("0 Car","1 Car","2 Cars","3+ Cars"))

res <- mlogit(as.formula(utilExpr),
              data=mlogitDataAvail,
              weights=weights_random,
              panel = FALSE,
              reflevel = '1',
              print.level=1)


Table_res <- data.frame(summary(res)$CoefTable)
write.csv(Table_res[c(1:4)],file="est_results_rand.csv")

b <- data.frame(summary(res)$logLik)
write.csv(b[c(1:1)],file="logLik.csv")

p <- data.frame(fitted(res,  type = c("probabilities")))

pred_prob <- bind_cols(hh, p)

write_csv(pred_prob, "mnl_prob_random.csv")



## 50-50 bias
mlogitDataAvail <- mlogit.data(hh_long,shape="long",chid.var = c("sampn"), alt.var = "altid", choice = "choiceBoolean", labels=c("0 Car","1 Car","2 Cars","3+ Cars"))


res <- mlogit(as.formula(utilExpr),data=mlogitDataAvail,weights=weights_bias, panel = FALSE,reflevel = '1',print.level=1)


Table_res <- data.frame(summary(res)$CoefTable)
write.csv(Table_res[c(1:4)],file="est_results_bias.csv")

b <- data.frame(summary(res)$logLik)
write.csv(b[c(1:1)],file="logLik.csv")

p <- data.frame(fitted(res,  type = c("probabilities")))

pred_prob <- bind_cols(hh, p)

write_csv(pred_prob, "mnl_prob_bias_hhveh.csv")


## 50-50 bias income
mlogitDataAvail <- mlogit.data(hh_long,shape="long",chid.var = c("sampn"), alt.var = "altid", choice = "choiceBoolean", labels=c("0 Car","1 Car","2 Cars","3+ Cars"))


res <- mlogit(as.formula(utilExpr),data=mlogitDataAvail,weights=weights_bias_inc, panel = FALSE,reflevel = '1',print.level=1)


Table_res <- data.frame(summary(res)$CoefTable)
write.csv(Table_res[c(1:4)],file="est_results.csv")

b <- data.frame(summary(res)$logLik)
write.csv(b[c(1:1)],file="logLik.csv")

p <- data.frame(fitted(res,  type = c("probabilities")))

pred_prob <- bind_cols(hh, p)

write_csv(pred_prob, "mnl_prob_bias_inc.csv")



