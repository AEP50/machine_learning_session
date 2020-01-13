# Script to file for estimation

#install.packages("dplyr")
#install.packages("foreign")
#install.packages("readr")
#install.packages("mosaic")
#install.packages("lazyeval")
#install.packages("DataCombine")
#install.packages("mlogit")

library(dplyr)
library(foreign)
library(readr)
library(mosaic)
library(lazyeval)
library(DataCombine)
library(broom)
library(mlogit)
library(reshape2)


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

hh <- hh %>%
  arrange(sampn)


hh <- hh %>%
  mutate(hhveh = ifelse(hhveh >=3, 3, hhveh))

hh <- hh %>%
  mutate(cbd = ifelse(hhtaz <=1000, 1, 0))

# add accessibilities
hh <- left_join(hh, access[,c("TAZ_NUM","ACCESS4","ACCESS5","ACCESS6")],by=c("hhtaz"="TAZ_NUM"))

hh[is.na(hh)] <- 0

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

hh <- left_join(hh, num_adults_df, by="sampn")



hh <- hh %>%
  dplyr::mutate(auto_access = log(1+ACCESS4),
                tran_access = log(1+ACCESS5),
                non_mot_access = log(1+ACCESS6),
                useest = 1,
                recapp = 0)


hh_more_acc <- hh %>%
  mutate(tran_access = tran_access + 0.6,
         useest = 0,
         recapp = 1,
         sampn = sampn + 176338)

hh_less_acc <- hh %>%
  mutate(tran_access = tran_access - 0.6,
         useest = 0,
         recapp = 2,
         sampn = sampn + 2*176338)


hh <- bind_rows(hh, hh_more_acc, hh_less_acc)

hh_long <- hh[rep(seq_len(nrow(hh)), 4),]

hh_long <- hh_long %>% 
  dplyr::group_by(sampn) %>% 
  dplyr::mutate(altid = sequence(n())) %>%
  dplyr::arrange(sampn, altid)

hh_long <- hh_long %>% 
  dplyr::mutate(choiceBoolean = ifelse(hhveh == altid-1, 1, 0) )

#hh_long <- melt(hh, id.vars=c("sampn"))

hh_long$altid <- as.factor(hh_long$altid)





hh_long[is.na(hh_long)] <- 0



utiltyExprFile <- "util_Terms_com_R.csv"

utilExpr <- CalcUtilityExpr(utiltyExprFile)



## ALL sample
# convert into mlogit data
mlogitDataAvail <- mlogit.data(filter(hh_long, useest == 1),shape="long",chid.var = c("sampn"), alt.var = "altid", choice = "choiceBoolean", labels=c("0 Car","1 Car","2 Cars","3+ Cars"))

# Write the model data to file
write.csv(x=filter(hh_long, useest == 1),
          file="./data/mnl_training_data_long_increased_tran_access.csv",
          sep=",",
          row.names=FALSE,
          col.names=TRUE)


#res <- mlogit(as.formula(choiceBoolean ~ 0|1),data=mlogitDataAvail, nests=list(zero_car='1',one_plus=c('2','3','4')), unscaled=TRUE, reflevel = '2',print.level=2)
res <- mlogit(as.formula(utilExpr),data=mlogitDataAvail, weights=useest, panel=FALSE, reflevel = '1',print.level=2)

Table_res <- data.frame(summary(res)$CoefTable)
write.csv(Table_res[c(1:4)],file="est_results.csv")


b <- data.frame(summary(res)$logLik)
write.csv(b[c(1:1)],file="logLik.csv")

p <- data.frame(fitted(res,  type = c("probabilities")))

pred_prob <- bind_cols(filter(hh, useest == 1), p)

#write_csv(pred_prob, "mnl_prob_all_changed_access.csv")

m <- data.frame(predict(res, filter(hh_long, recapp == 1)))

pred_prob_more <- bind_cols(filter(hh, recapp == 1), m)

pred_prob_more <- pred_prob_more %>%
  mutate(sampn = sampn -176338 )

n <- data.frame(predict(res, filter(hh_long, recapp == 2)))

pred_prob_less <- bind_cols(filter(hh, recapp == 2), n)

pred_prob_less <- pred_prob_less %>%
  mutate(sampn = sampn - 2*176338 )


write_csv(pred_prob_more, "mnl_prob_all_more_access.csv")
write_csv(pred_prob_less, "mnl_prob_all_less_access.csv")
