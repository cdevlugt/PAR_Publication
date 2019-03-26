memory.limit(size=(128000+196000)) #phyical RAM + a page file

##### Probability tables #####
subunitDecompress <- function(df, sub=1, drops=c(1:8), drop=T, name="Subunit"){
  # decompresses the data by replicating a column based on read count; can also remove columns
  #
  # Args:
  #   df: data frame containing data
  #   sub:  column containing data to be replicated
  #   drops:  columns to be dropped
  #   drop: should the other subunits be dropped
  #
  # Returns:
  #   df: decompress data frame for a given subunit
  if(drop){
    drops <- drops[drops != sub]#retains subunit column if in drops range
    if(length(drops)!=0){
      df <- df[,-drops]
      sub <- sub - sum(drops < sub)#one more opperation versus replicating than dropping, but stores less stuff in ram during replication
    }
  }
  df <- df[(rep(row.names(df), df[,sub])),]
  df[,sub] <- colnames(df)[sub]
  colnames(df)[sub] <- name
  return (df)
}
subunitSeperate <- function(df, subs=c(1:8), decompress=F, name="Subunit"){
  #Takes a data frame with read counts in columns and then clones the rows seperating the subunits so that they are on seperate rows
  #
  #Args:
  # df: data frame with subunit read counts
  # subs: column numbers of subunits
  # decompress: flag indicating if the subunits should be decompressed
  # name: name for collumn if decompress is T
  #
  #Return:
  # dfFinal: data frame with columns cloned based on subunit
  
  for(i in 1:length(subs)){
    dfTemp <- df[df[,subs[i]]>0,]
    if(decompress){
      dfTemp <- subunitDecompress(dfTemp, sub=subs[i], drops=subs, drop=T, name=name)
    } else {
      dfTemp[,subs[-i]] <- 0
    }
    if(i==1){
      dfFinal <- dfTemp
    } else {
      dfFinal <- rbind(dfFinal, dfTemp)
    }
    rm(dfTemp)
  }
  return(dfFinal)
}
seqComp <- function(df, chars=c('A', 'C', 'G', 'T'), Sequences=10, otherCols=NA, ignore.case=F, Occur=T, Freq=T, All=F, purines=F, strong=F){
  #calculates the relative or absolute portion of a string that is made up of a substring; works to determine the portion of a string made up of
  # RNA or DNA bases/motifs
  #
  #Args:
  # df: dataframe or vector containg the strings
  # chars:  substrings to search for
  # Sequences:  columns containing strings to analyze
  # otherCols:  allows for retention of other data frame columns; you may also cbind results and ignore this feature
  # ignore.case:  ignore case when searching for substrings
  # Occur:  count occurances
  # Freq: count frequency
  # All:  return all data rather than just occurances and/or frequency
  # purines:  determine A and G content of string
  # strong: determine G and C content of string
  #
  #return:
  # cbind(df, dfR): data frame containing sequences and other columns
  # dfR:  dataframe containing sequence composition as specified in Occur and Freq
  
  
  if(!Occur && !Freq){
    stop("Specify either Occur of Freq as true; Occur returns occurances, Freq returns frequency; both may be true")
  }
  
  cols <- Sequences
  Sequences <- as.vector(1:length(Sequences))
  
  if(!is.na(otherCols[1])){
    cols <- c(otherCols, cols)
    otherCols <- 1:length(otherCols)
    Sequences <- as.vector((length(cols)-length(Sequences)+1):length(cols))
  }
  df <- data.frame(df[,cols], stringsAsFactors = F)
  cols <- ncol(df)
  if(Occur){
    for(i in 1:length(Sequences)){
      for(j in 1:length(chars)){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub(chars[j],"",df[,Sequences[i]])))/nchar(chars[j])
        colnames(df)[ncol(df)] <- paste("Abs",chars[j],colnames(df)[Sequences[i]],sep="_")
      }
      if(purines){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("A","",df[,Sequences[i]]))))
        colnames(df)[ncol(df)] <- paste("Abs","Purine",colnames(df)[Sequences[i]],sep="_")
      }
      if(strong){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("C","",df[,Sequences[i]]))))
        colnames(df)[ncol(df)] <- paste("Abs","GC",colnames(df)[Sequences[i]],sep="_")
      }
    }
    dfR <- df[,c((cols+1):ncol(df))]
    df <- df[,1:cols]
  }
  if(Freq){
    for(i in 1:length(Sequences)){
      for(j in 1:length(chars)){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub(chars[j],"",df[,Sequences[i]])))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel",chars[j],colnames(df)[Sequences[i]],sep="_")
      }
      if(purines){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("A","",df[,Sequences[i]]))))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel","Purine",colnames(df)[Sequences[i]],sep="_")
      }
      if(strong){
        df[,(ncol(df)+1)] <- (nchar(df[,Sequences[i]]) - nchar(gsub("G","",gsub("C","",df[,Sequences[i]]))))/nchar(df[,Sequences[i]]) * 100
        colnames(df)[ncol(df)] <- paste("Rel","GC",colnames(df)[Sequences[i]],sep="_")
      }
    }
    
    if(Occur){
      dfR <- cbind(dfR, df[,c((cols+1):ncol(df))])
    } else {
      dfR <- df[,c((cols+1):ncol(df))]
    }
    dfR[is.na(dfR)] <-0
    df <- df[,1:cols]
  }
  
  if(All){
    return(cbind(df, dfR))
  } else {
    return(dfR)
  }
}

#g included dfpp
dfPreProcess <- function(df, strain='', U=1:5, C=6:8) {
  
  #series handing
  df$series <- 'No_Realignment'
  df$series[df$rounds>=1] <- 'Realigned'
  
  #subunit handling
  subs <- c(U, C)
  subs <- unique(subs)
  subs <- subs[order(subs)]
  subNames <- colnames(df)[subs]
  C <- colnames(df)[C]
  
  #set T as U because mRNA
  df$Trim_Sequence <- gsub('T', 'U', df$Trim_Sequence)
  
  #delimit columns
  df <- df[,c(subs, grep('Trim_Sequence', colnames(df)), grep('NT_coupure', colnames(df)), grep('Trim_Len', colnames(df)), grep('series', colnames(df)))]
  gc()
  
  #G+1
  df$NT_coupure <- substr(df$NT_coupure, 1, 1)
  df$NT_coupure <- gsub('[^G]', "", df$NT_coupure)
  
  #Sequence End
  df$Seq_End <- substr(df$Trim_Sequence, nchar(df$Trim_Sequence), nchar(df$Trim_Sequence))
  df$Seq_End <- paste(df$Seq_End, df$NT_coupure, sep='')
  
  #GCA content
  df$Trim_Sequence[df$NT_coupure=='G'] <- paste(df$Trim_Sequence[df$NT_coupure=='G'], "G", sep='')
  df$GCA_content <- 100 - seqComp(df, chars='U', grep('Trim_Sequence', colnames(df)), Occur=F, Freq=T)
  
  #Length
  df$Trim_Len[df$NT_coupure=='G'] <- df$Trim_Len[df$NT_coupure=='G'] + 1
  
  #organize df
  df <- df[,-(grep('Trim_Sequence', colnames(df)))]
  colnames(df)[(length(subs)+1):(length(subs)+2)] <- c('G_Comp', 'Length')
  
  #since data is G_Comp included, drop G_Comp
  df <- df[,-(grep('G_Comp', colnames(df)))]
  
  #factor for glm
  df$series <- factor(df$series, levels=c('No_Realignment', 'Realigned'))
  
  #length is factored because data is of type int, not double; decimal levels are impractical therefore it must be fitted as a factor
  df$Length <- factor(df$Length, levels=c(9:17))
  df$Seq_End <- factor(df$Seq_End, levels =c('U', 'C', 'G', 'A', 'UG', 'CG', 'GG', 'AG'))
  
  df <- subunitSeperate(df, subs=subs, decompress = T)
  rm(subs)
  
  #add template
  df$Template <- "3'U4"
  for(i in 1:length(C)){
    df$Template[df$Subunit==C[i]] <- "3'C4"
  }
  rm(C)
  
  #factor for glm
  df$Template <- factor(df$Template, levels=c("3'U4", "3'C4"))
  df$Subunit <- factor(df$Subunit, levels=c(subNames))
  rm(subNames)
  
  df$strain <- strain
  
  #make it look pretty
  #Subunit|Length|series|Seq_End|GCA_content|Template|strain
  df <- df[,c(3,2,4,6,5,1,7)]
  #series|Length|Seq_End|Template|GCA_content|Subunit|strain
  
  return(df)
}
parStats <- function(df){
  #takes a data frame (df) that has been formatted by dfPreProcess
  #returns the percentage of instances in series 1 or 2 (Realigned vs. No_Realignment) on global and for each template
  #
  #taking the percentage of instances in group 1 or 2 relative to the total number of instances is mathematically equivalent to
  #the logistic regression predictions including higher level interactions, this is faster than performing the regression
  #essentially a better version of:
  #x <- glm(series ~ 1 + Length + Seq_End + Template + Length:Seq_End + Length:Template + Seq_End:Template + Length:Seq_End:Template, data=df, family = binomial(link='logit'))
  
  #create a matrix to store data
  parStats <- matrix(data=0, nrow=length(levels(df$Seq_End)), ncol=length(levels(df$Length)))
  colnames(parStats) <- levels(df$Length)
  rownames(parStats) <- levels(df$Seq_End)
  
  #create a list of matrix with the names as listed
  parStats <- list(parStats, parStats, parStats, parStats, parStats, parStats, parStats, parStats, parStats, parStats, parStats, parStats)
  names(parStats) <- c('Realigned_Counts', 'Realigned_Frequency', 'No_Realignment_Counts', 'No_Realignment_Frequency', 
                       "Realigned_Counts_3'U4", "Realigned_Frequency_3'U4", "No_Realignment_Counts_3'U4", "No_Realignment_Frequency_3'U4",
                       "Realigned_Counts_3'C4", "Realigned_Frequency_3'C4", "No_Realignment_Counts_3'C4", "No_Realignment_Frequency_3'C4")
  
  #obtain counts
  for(i in 1: ncol(parStats[[1]])){ #i is columns
    for(j in 1:nrow(parStats[[1]])){ #j is rows
      #get the counts for each series, and then each series for each template
      parStats[[1]][j,i] <- nrow(df[df$Length==colnames(parStats[[1]])[i] & df$Seq_End==rownames(parStats[[1]])[j] & df$series=='Realigned',])
      parStats[[3]][j,i] <- nrow(df[df$Length==colnames(parStats[[3]])[i] & df$Seq_End==rownames(parStats[[3]])[j] & df$series=='No_Realignment',])
      parStats[[5]][j,i] <- nrow(df[df$Length==colnames(parStats[[5]])[i] & df$Seq_End==rownames(parStats[[5]])[j] & df$series=='Realigned' & df$Template=="3'U4",])
      parStats[[7]][j,i] <- nrow(df[df$Length==colnames(parStats[[7]])[i] & df$Seq_End==rownames(parStats[[7]])[j] & df$series=='No_Realignment'& df$Template=="3'U4",])
      parStats[[9]][j,i] <- nrow(df[df$Length==colnames(parStats[[9]])[i] & df$Seq_End==rownames(parStats[[9]])[j] & df$series=='Realigned' & df$Template=="3'C4",])
      parStats[[11]][j,i] <- nrow(df[df$Length==colnames(parStats[[11]])[i] & df$Seq_End==rownames(parStats[[11]])[j] & df$series=='No_Realignment'& df$Template=="3'C4",])
    }
  }
  
  #get the percentages
  parStats$Realigned_Frequency <- parStats$Realigned_Counts/(parStats$Realigned_Counts + parStats$No_Realignment_Counts)*100
  parStats$No_Realignment_Frequency <- parStats$No_Realignment_Counts/(parStats$Realigned_Counts + parStats$No_Realignment_Counts)*100
  parStats$`Realigned_Frequency_3'U4` <- parStats$`Realigned_Counts_3'U4`/(parStats$`Realigned_Counts_3'U4` + parStats$`No_Realignment_Counts_3'U4`)*100
  parStats$`No_Realignment_Frequency_3'U4` <- parStats$`No_Realignment_Counts_3'U4`/(parStats$`Realigned_Counts_3'U4` + parStats$`No_Realignment_Counts_3'U4`)*100
  parStats$`Realigned_Frequency_3'C4` <- parStats$`Realigned_Counts_3'C4`/(parStats$`Realigned_Counts_3'C4` + parStats$`No_Realignment_Counts_3'C4`)*100
  parStats$`No_Realignment_Frequency_3'C4` <- parStats$`No_Realignment_Counts_3'C4`/(parStats$`Realigned_Counts_3'C4` + parStats$`No_Realignment_Counts_3'C4`)*100
  
  return(parStats)
}

#Brisbane
#setwd("") #set dir if not default
bri <- read.csv("BRI_All_Match.csv", stringsAsFactors = F)
bri <- dfPreProcess(bri, 'Brisbane')
#set bri as df
df <- bri
rm(bri)
df <- df[,c(1:4)]
gc()

#get probabilities and counts, save data
setwd("./Effect Modeling") #folder must exist
x<-parStats(df)
for(z in 1:length(x)){
  write.csv(x[[z]], file = paste('Brisbane ', names(x)[z], '.csv', sep=''))
}

rm(x, df)
gc()

#Hong Kong
setwd("..")
hk <- read.csv("HK_All_Match.csv", stringsAsFactors = F)
hk <- dfPreProcess(hk, 'Hong Kong',  U=1:6, C=7:8)
#set hk as df
df <- hk
rm(hk)
df <- df[,c(1:4)]
gc()

#get probabilities and counts, save data
setwd("./Effect Modeling")
x<-parStats(df)
for(z in 1:length(x)){
  write.csv(x[[z]], file = paste('Hong Kong ', names(x)[z], '.csv', sep=''))
}

rm(x, df)
gc()

#Puerto Rico
setwd("..")
pr8 <- read.csv("PR8_All_Match.csv", stringsAsFactors = F)
pr8 <- dfPreProcess(pr8, 'Puerto Rico')
#set pr8 as df
df <- pr8
rm(pr8)
df <- df[,c(1:4)]
gc()

#get probabilities and counts, save data
setwd("./Effect Modeling")
x<-parStats(df)
for(z in 1:length(x)){
  write.csv(x[[z]], file = paste('Puerto Rico ', names(x)[z], '.csv', sep=''))
}

rm(x, df)
gc()

#WSN
setwd("..")
wsn <- read.csv("wsn_All_Match.csv", stringsAsFactors = F)
wsn <- dfPreProcess(wsn, 'WSN')
#set wsn as df
df <- wsn
rm(wsn)
df <- df[,c(1:4)]
gc()

#get probabilities and counts, save data
setwd("./Effect Modeling")
x<-parStats(df)
for(z in 1:length(x)){
  write.csv(x[[z]], file = paste('WSN ', names(x)[z], '.csv', sep=''))
}

rm(x, df)
gc()
setwd('..')

#Luciferase
#setwd("") #set dir if not default
luc <- read.csv("LUC_All_Match.csv", stringsAsFactors = F)
luc <- dfPreProcess(luc, 'Luciferase', U=c(2, 4), C=c(1, 3))
#set luc as df
df <- luc
rm(luc)
df <- df[,c(1:4)]
gc()

#get probabilities and counts, save data
setwd("./Effect Modeling") #folder must exist
x<-parStats(df)
for(z in 1:length(x)){
  write.csv(x[[z]], file = paste('Luciferase ', names(x)[z], '.csv', sep=''))
}

rm(x, df)
gc()
setwd('..')

##### Predictive Model Heatmap #####
require(ggplot2)
require(cowplot)
require(doParallel)

#import strain PAR
#setwd() #if not default directory, set it here
setwd('./Overview')
parav <- read.csv('PAR rates.csv', stringsAsFactors = F, row.names = 1)[2,]
stPAR <- as.numeric(rowMeans(parav[1,c(1,2,3)])) #because WSN length is invalid it is omitted

#splits for graph
splits <- c(0, mean(c(mean(c(stPAR, 0)),0)), mean(c(mean(c(stPAR, 0)),stPAR)), #below
            stPAR, #PAR rate
            mean(c(mean(c(stPAR, 100)),stPAR)), mean(c(mean(c(stPAR, 100)),100)), 100) #above
splits <- splits/100 #fix for heatmap

#colour scheme
#this is set up so that the colours are similar to the distance from the mean PAR rate of all 3 strains 
colours = c('green', '#40E0D0', 'light grey',
            'white',
            'yellow', 'orange', 'red')

setwd('..')
setwd('./Effect Modeling')

#All data heatmaps
#read data in parallel
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
series <- c('Realigned', 'No_Realignment')
cluster <- 8
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar%
  read.csv(paste(strains[floor((j-1)/2)+1], ' ', series[floor(j/2)-ceiling(j/2)+2], '_counts.csv', sep=''), stringsAsFactors = F, row.names = 1)
stopCluster(cl)
gc()

#unlist data
varNames <- c('pre', 'pno', 'pal',
              'hre', 'hno', 'hal',
              'wre', 'wno', 'wal',
              'bre', 'bno', 'bal')
for(j in 1:4){
  assign(varNames[(j-1)*3+1], dfs[[(j-1)*2+1]])
  assign(varNames[(j-1)*3+2], dfs[[(j-1)*2+2]])
  assign(varNames[(j-1)*3+3],(dfs[[(j-1)*2+1]] + dfs[[(j-1)*2+2]]))
}

#threshold
#6400th of reads as threshold
threshold <- 6400
pal[pal<=sum(pal)/threshold] <- 0
hal[hal<=sum(hal)/threshold] <- 0
wal[wal<=sum(wal)/threshold] <- 0
bal[bal<=sum(bal)/threshold] <- 0

#obtain percentages
pr8 <- as.matrix(pre/pal*100)
hk <- as.matrix(hre/hal*100)
wsn <- as.matrix(wre/wal*100)
bri <- as.matrix(bre/bal*100)

#drop inf
for(i in 1:ncol(pr8)){
  pr8[,i][is.infinite(pr8[,i])] <- NaN
  hk[,i][is.infinite(hk[,i])] <- NaN
  wsn[,i][is.infinite(wsn[,i])] <- NaN
  bri[,i][is.infinite(bri[,i])] <- NaN
}

#create a matrix with the average of all matrix in varNames
#sums ignoring NaN
varNames <- c('pr8', 'hk', 'bri')
all <- matrix(NaN, nrow = nrow(pr8), ncol = ncol(pr8))
row.names(all) <- row.names(pr8)
colnames(all) <- colnames(pr8)
for(i in 1:ncol(pr8)){
  for(j in 1:nrow(pr8)){
    val <- 0 #placeholder value
    cnt <- 0 #count of vals
    for(k in 1:length(varNames)){
      if(!is.na(get(varNames[k])[j,i])){
        val <- val + get(varNames[k])[j,i]
        cnt <- cnt + 1
      }
    }
    if(cnt>0){
      all[j,i] <- val/cnt
    }
  }
}

#make data ggplot compatible
Seq_Ends <- rep(row.names(all), ncol(all))
lengths <- rep(colnames(all)[1], nrow(all))
PAR <- all[,1]
for(i in 2:ncol(all)){
  lengths <- c(lengths, rep(colnames(all)[i], nrow(all)))
  PAR <- c(PAR, all[,i])
}

#put ggplot friendly data into a df
hmData <- data.frame(length = lengths, Sequence_End = Seq_Ends, PAR = PAR, stringsAsFactors = F)
hmData$length <- as.numeric(as.character(gsub('X', '', hmData$length)))
hmData$Sequence_End <- factor(hmData$Sequence_End, levels = rev(row.names(all)))
rm(PAR, Seq_Ends, lengths)
gc()

#heatmap
hmA <- ggplot(data=hmData, aes(x=length, y=Sequence_End, fill=PAR)) +
  geom_tile(colour="white",size=0.25) +
  geom_text(aes(label=formatC(PAR, digits = 3)), size=1.5) +
  scale_fill_gradientn(limits=c(0,100),
                       colours = colours,
                       values = splits,
                       na.value='#000000') +
  theme_bw() +
  ggtitle("All Conserved Sequences") +
  theme(plot.title = element_text(size=8,face="bold", hjust=0.5)) +
  scale_x_continuous(breaks = c(9:17), labels = c(9:17)) +
  xlab('Length (nt)') +
  ylab("Sequence End") +
  theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
  theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) +
  theme(legend.title =element_blank(), legend.text=element_text(size=7))
legend <- get_legend(hmA)
hmA <- hmA + theme(legend.position='none')

#3'U4 data heatmaps
#read data in parallel
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
series <- c('Realigned', 'No_Realignment')
cluster <- 8
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar%
  read.csv(paste(strains[floor((j-1)/2)+1], ' ', series[floor(j/2)-ceiling(j/2)+2], "_counts_3'U4.csv", sep=''), stringsAsFactors = F, row.names = 1)
stopCluster(cl)
gc()

#unlist data
varNames <- c('pre', 'pno', 'pal',
              'hre', 'hno', 'hal',
              'wre', 'wno', 'wal',
              'bre', 'bno', 'bal')
for(j in 1:4){
  assign(varNames[(j-1)*3+1], dfs[[(j-1)*2+1]])
  assign(varNames[(j-1)*3+2], dfs[[(j-1)*2+2]])
  assign(varNames[(j-1)*3+3],(dfs[[(j-1)*2+1]] + dfs[[(j-1)*2+2]]))
}

#threshold
#1200th of reads as threshold
threshold <- 1200
pal[pal<=sum(pal)/threshold] <- 0
hal[hal<=sum(hal)/threshold] <- 0
wal[wal<=sum(wal)/threshold] <- 0
bal[bal<=sum(bal)/threshold] <- 0

#calculate PAR rates
pr8 <- as.matrix(pre/pal*100)
hk <- as.matrix(hre/hal*100)
wsn <- as.matrix(wre/wal*100)
bri <- as.matrix(bre/bal*100)

#drop inf
for(i in 1:ncol(pr8)){
  pr8[,i][is.infinite(pr8[,i])] <- NaN
  hk[,i][is.infinite(hk[,i])] <- NaN
  wsn[,i][is.infinite(wsn[,i])] <- NaN
  bri[,i][is.infinite(bri[,i])] <- NaN
}

#create a matrix with the average of all matrix in varNames
#sums ignoring NaN
varNames <- c('pr8', 'hk', 'bri')
all <- matrix(NaN, nrow = nrow(pr8), ncol = ncol(pr8))
row.names(all) <- row.names(pr8)
colnames(all) <- colnames(pr8)
for(i in 1:ncol(pr8)){
  for(j in 1:nrow(pr8)){
    val <- 0 #placeholder value
    cnt <- 0 #count of vals
    for(k in 1:length(varNames)){
      if(!is.na(get(varNames[k])[j,i])){
        val <- val + get(varNames[k])[j,i]
        cnt <- cnt + 1
      }
    }
    if(cnt>0){
      all[j,i] <- val/cnt
    }
  }
}

#make data ggplot compatible
Seq_Ends <- rep(row.names(all), ncol(all))
lengths <- rep(colnames(all)[1], nrow(all))
PAR <- all[,1]
for(i in 2:ncol(all)){
  lengths <- c(lengths, rep(colnames(all)[i], nrow(all)))
  PAR <- c(PAR, all[,i])
}

#tabulate ggplot data
hmData <- data.frame(length = lengths, Sequence_End = Seq_Ends, PAR = PAR, stringsAsFactors = F)
hmData$length <- as.numeric(as.character(gsub('X', '', hmData$length)))
hmData$Sequence_End <- factor(hmData$Sequence_End, levels = rev(row.names(all)))
rm(PAR, Seq_Ends, lengths)
gc()

hmU <- ggplot(data=hmData, aes(x=length, y=Sequence_End, fill=PAR)) +
  geom_tile(colour="white",size=0.25) +
  geom_text(aes(label=formatC(PAR, digits = 3)), size=1.5) +
  scale_fill_gradientn(colours = colours,
                       values = splits,
                       na.value='#000000') +
  theme_bw() +
  ggtitle("3'U4 Conserved Sequences") +
  theme(plot.title = element_text(size=8,face="bold", hjust=0.5)) +
  scale_x_continuous(breaks = c(9:17), labels = c(9:17)) +
  xlab('Length (nt)') +
  ylab("Sequence End") +
  theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
  theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) + 
  theme(legend.title =element_blank(), legend.text=element_text(size=7))
hmU <- hmU + theme(legend.position='none')

#3'C4 data heatmaps
#read data in parallel
strains <- c('Puerto Rico', 'Hong Kong', 'WSN', 'Brisbane')
series <- c('Realigned', 'No_Realignment')
cluster <- 8
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar%
  read.csv(paste(strains[floor((j-1)/2)+1], ' ', series[floor(j/2)-ceiling(j/2)+2], "_counts_3'C4.csv", sep=''), stringsAsFactors = F, row.names = 1)
stopCluster(cl)
gc()

#unlist data
varNames <- c('pre', 'pno', 'pal',
              'hre', 'hno', 'hal',
              'wre', 'wno', 'wal',
              'bre', 'bno', 'bal')
for(j in 1:4){
  assign(varNames[(j-1)*3+1], dfs[[(j-1)*2+1]])
  assign(varNames[(j-1)*3+2], dfs[[(j-1)*2+2]])
  assign(varNames[(j-1)*3+3],(dfs[[(j-1)*2+1]] + dfs[[(j-1)*2+2]]))
}

#threshold
#12800th of reads as threshold
threshold <- 12800
pal[pal<=sum(pal)/threshold] <- 0
hal[hal<=sum(hal)/threshold] <- 0
wal[wal<=sum(wal)/threshold] <- 0
bal[bal<=sum(bal)/threshold] <- 0

#calculate PAR frequencies
pr8 <- as.matrix(pre/pal*100)
hk <- as.matrix(hre/hal*100)
wsn <- as.matrix(wre/wal*100)
bri <- as.matrix(bre/bal*100)

#drop inf
for(i in 1:ncol(pr8)){
  pr8[,i][is.infinite(pr8[,i])] <- NaN
  hk[,i][is.infinite(hk[,i])] <- NaN
  wsn[,i][is.infinite(wsn[,i])] <- NaN
  bri[,i][is.infinite(bri[,i])] <- NaN
}

#create a matrix with the average of all matrix in varNames
#sums ignoring NaN
varNames <- c('pr8', 'hk', 'bri')
all <- matrix(NaN, nrow = nrow(pr8), ncol = ncol(pr8))
row.names(all) <- row.names(pr8)
colnames(all) <- colnames(pr8)
for(i in 1:ncol(pr8)){
  for(j in 1:nrow(pr8)){
    val <- 0 #placeholder value
    cnt <- 0 #count of vals
    for(k in 1:length(varNames)){
      if(!is.na(get(varNames[k])[j,i])){
        val <- val + get(varNames[k])[j,i]
        cnt <- cnt + 1
      }
    }
    if(cnt>0){
      all[j,i] <- val/cnt
    }
  }
}

#make data ggplot compatible
Seq_Ends <- rep(row.names(all), ncol(all))
lengths <- rep(colnames(all)[1], nrow(all))
PAR <- all[,1]
for(i in 2:ncol(all)){
  lengths <- c(lengths, rep(colnames(all)[i], nrow(all)))
  PAR <- c(PAR, all[,i])
}

#tabulate
hmData <- data.frame(length = lengths, Sequence_End = Seq_Ends, PAR = PAR, stringsAsFactors = F)
hmData$length <- as.numeric(as.character(gsub('X', '', hmData$length)))
hmData$Sequence_End <- factor(hmData$Sequence_End, levels = rev(row.names(all)))
rm(PAR, Seq_Ends, lengths)
gc()

hmC <- ggplot(data=hmData, aes(x=length, y=Sequence_End, fill=PAR)) +
  geom_tile(colour="white",size=0.25) +
  geom_text(aes(label=formatC(PAR, digits = 3)), size=1.5) +
  scale_fill_gradientn(colours = colours,
                       values = splits,
                       na.value='#000000') +
  theme_bw() +
  ggtitle("3'C4 Conserved Sequences") +
  theme(plot.title = element_text(size=8,face="bold", hjust=0.5)) +
  scale_x_continuous(breaks = c(9:17), labels = c(9:17)) +
  xlab('Length (nt)') +
  ylab("Sequence End") +
  theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
  theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) + 
  theme(legend.title =element_blank(), legend.text=element_text(size=7))
hmC <- hmC + theme(legend.position='none')

#Save Plots
hm <- plot_grid(hmA, hmU, hmC, legend, ncol=4, nrow=1, rel_widths = c(1,1,1,0.3))
save_plot('Effect Modeling Heatmap.png', hm, base_height = 2.3, base_width = 8, dpi=900)

setwd('..')

##### Strain Probability Heatmaps #####
#requires: 8 - Effect Modeling Probability Tables
require(ggplot2)
require(cowplot)
require(doParallel)

#import strain PAR
#setwd() #if not default directory, set it here
setwd('./Overview')
parav <- read.csv('PAR rates.csv', stringsAsFactors = F, row.names = 1)[2,]
setwd('..')


#need to pass parav and the tables
strainEffectModels <- function(strain, stPAR){
  #splits for graph
  splits <- c(0, mean(c(mean(c(stPAR, 0)),0)), mean(c(mean(c(stPAR, 0)),stPAR)), #below
              stPAR, #PAR rate
              mean(c(mean(c(stPAR, 100)),stPAR)), mean(c(mean(c(stPAR, 100)),100)), 100) #above
  splits <- splits/100 #fix for heatmap
  
  #colour scheme
  #this is set up so that the colours are similar to the distance from the mean PAR rate of all 3 strains 
  colours = c('green', '#40E0D0', 'light grey',
              'white',
              'yellow', 'orange', 'red')
  
  #read data
  rea <- read.csv(paste(strain, 'Realigned_Counts.csv'), stringsAsFactors = F, row.names = 1)
  noR <- read.csv(paste(strain, 'No_Realignment_Counts.csv'), stringsAsFactors = F, row.names = 1)
  All <- rea + noR
  
  #threshold
  #6400th of reads as threshold
  threshold <- 6400
  All[All<=sum(All)/threshold] <- 0
  
  #obtain percentages
  df <- as.matrix(rea/All*100)
  
  #drop inf
  for(i in 1:ncol(df)){
    df[,i][is.infinite(df[,i])] <- NaN
  }
  
  #create a matrix with the average of all matrix in varNames
  #sums ignoring NaN
  varNames <- c('df')
  all <- matrix(NaN, nrow = nrow(df), ncol = ncol(df))
  row.names(all) <- row.names(df)
  colnames(all) <- colnames(df)
  for(i in 1:ncol(df)){
    for(j in 1:nrow(df)){
      val <- 0 #placeholder value
      cnt <- 0 #count of vals
      for(k in 1:length(varNames)){
        if(!is.na(get(varNames[k])[j,i])){
          val <- val + get(varNames[k])[j,i]
          cnt <- cnt + 1
        }
      }
      if(cnt>0){
        all[j,i] <- val/cnt
      }
    }
  }
  
  #make data ggplot compatible
  Seq_Ends <- rep(row.names(all), ncol(all))
  lengths <- rep(colnames(all)[1], nrow(all))
  PAR <- all[,1]
  for(i in 2:ncol(all)){
    lengths <- c(lengths, rep(colnames(all)[i], nrow(all)))
    PAR <- c(PAR, all[,i])
  }
  
  #put ggplot friendly data into a df
  hmData <- data.frame(length = lengths, Sequence_End = Seq_Ends, PAR = PAR, stringsAsFactors = F)
  hmData$length <- as.numeric(as.character(gsub('X', '', hmData$length)))
  hmData$Sequence_End <- factor(hmData$Sequence_End, levels = rev(row.names(all)))
  rm(PAR, Seq_Ends, lengths)
  gc()
  
  #heatmap
  hmA <- ggplot(data=hmData, aes(x=length, y=Sequence_End, fill=PAR)) +
    geom_tile(colour="white",size=0.25) +
    geom_text(aes(label=formatC(PAR, digits = 3)), size=1.5) +
    scale_fill_gradientn(limits=c(0,100),
                         colours = colours,
                         values = splits,
                         na.value='#000000') +
    theme_bw() +
    ggtitle("All Conserved Sequences") +
    theme(plot.title = element_text(size=8,face="bold", hjust=0.5)) +
    scale_x_continuous(breaks = c(9:17), labels = c(9:17)) +
    xlab('Length (nt)') +
    ylab("Sequence End") +
    theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
    theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) +
    theme(legend.title =element_blank(), legend.text=element_text(size=7))
  legend <- get_legend(hmA)
  hmA <- hmA + theme(legend.position='none')
  
  #####3'U4 data heatmaps
  #read data
  rea <- read.csv(paste(strain, "Realigned_Counts_3'U4.csv"), stringsAsFactors = F, row.names = 1)
  noR <- read.csv(paste(strain, "No_Realignment_Counts_3'U4.csv"), stringsAsFactors = F, row.names = 1)
  All <- rea + noR
  
  #threshold
  #1200th of reads as threshold
  threshold <- 1200
  All[All<=sum(All)/threshold] <- 0
  
  #obtain percentages
  df <- as.matrix(rea/All*100)
  
  #drop inf
  for(i in 1:ncol(df)){
    df[,i][is.infinite(df[,i])] <- NaN
  }
  
  #create a matrix with the average of all matrix in varNames
  #sums ignoring NaN
  varNames <- c('df')
  all <- matrix(NaN, nrow = nrow(df), ncol = ncol(df))
  row.names(all) <- row.names(df)
  colnames(all) <- colnames(df)
  for(i in 1:ncol(df)){
    for(j in 1:nrow(df)){
      val <- 0 #placeholder value
      cnt <- 0 #count of vals
      for(k in 1:length(varNames)){
        if(!is.na(get(varNames[k])[j,i])){
          val <- val + get(varNames[k])[j,i]
          cnt <- cnt + 1
        }
      }
      if(cnt>0){
        all[j,i] <- val/cnt
      }
    }
  }
  
  #make data ggplot compatible
  Seq_Ends <- rep(row.names(all), ncol(all))
  lengths <- rep(colnames(all)[1], nrow(all))
  PAR <- all[,1]
  for(i in 2:ncol(all)){
    lengths <- c(lengths, rep(colnames(all)[i], nrow(all)))
    PAR <- c(PAR, all[,i])
  }
  
  #tabulate ggplot data
  hmData <- data.frame(length = lengths, Sequence_End = Seq_Ends, PAR = PAR, stringsAsFactors = F)
  hmData$length <- as.numeric(as.character(gsub('X', '', hmData$length)))
  hmData$Sequence_End <- factor(hmData$Sequence_End, levels = rev(row.names(all)))
  rm(PAR, Seq_Ends, lengths)
  gc()
  
  hmU <- ggplot(data=hmData, aes(x=length, y=Sequence_End, fill=PAR)) +
    geom_tile(colour="white",size=0.25) +
    geom_text(aes(label=formatC(PAR, digits = 3)), size=1.5) +
    scale_fill_gradientn(colours = colours,
                         values = splits,
                         na.value='#000000') +
    theme_bw() +
    ggtitle("3'U4 Conserved Sequences") +
    theme(plot.title = element_text(size=8,face="bold", hjust=0.5)) +
    scale_x_continuous(breaks = c(9:17), labels = c(9:17)) +
    xlab('Length (nt)') +
    ylab("Sequence End") +
    theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
    theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) + 
    theme(legend.title =element_blank(), legend.text=element_text(size=7))
  hmU <- hmU + theme(legend.position='none')
  
  #####3'C4 data heatmaps
  #read data
  rea <- read.csv(paste(strain, "Realigned_Counts_3'C4.csv"), stringsAsFactors = F, row.names = 1)
  noR <- read.csv(paste(strain, "No_Realignment_Counts_3'C4.csv"), stringsAsFactors = F, row.names = 1)
  All <- rea + noR
  
  #threshold
  #12800 of reads as threshold
  threshold <- 12800
  All[All<=sum(All)/threshold] <- 0
  
  #obtain percentages
  df <- as.matrix(rea/All*100)
  
  #drop inf
  for(i in 1:ncol(df)){
    df[,i][is.infinite(df[,i])] <- NaN
  }
  
  #create a matrix with the average of all matrix in varNames
  #sums ignoring NaN
  varNames <- c('df')
  all <- matrix(NaN, nrow = nrow(df), ncol = ncol(df))
  row.names(all) <- row.names(df)
  colnames(all) <- colnames(df)
  for(i in 1:ncol(df)){
    for(j in 1:nrow(df)){
      val <- 0 #placeholder value
      cnt <- 0 #count of vals
      for(k in 1:length(varNames)){
        if(!is.na(get(varNames[k])[j,i])){
          val <- val + get(varNames[k])[j,i]
          cnt <- cnt + 1
        }
      }
      if(cnt>0){
        all[j,i] <- val/cnt
      }
    }
  }
  
  #make data ggplot compatible
  Seq_Ends <- rep(row.names(all), ncol(all))
  lengths <- rep(colnames(all)[1], nrow(all))
  PAR <- all[,1]
  for(i in 2:ncol(all)){
    lengths <- c(lengths, rep(colnames(all)[i], nrow(all)))
    PAR <- c(PAR, all[,i])
  }
  
  #tabulate
  hmData <- data.frame(length = lengths, Sequence_End = Seq_Ends, PAR = PAR, stringsAsFactors = F)
  hmData$length <- as.numeric(as.character(gsub('X', '', hmData$length)))
  hmData$Sequence_End <- factor(hmData$Sequence_End, levels = rev(row.names(all)))
  rm(PAR, Seq_Ends, lengths)
  gc()
  
  hmC <- ggplot(data=hmData, aes(x=length, y=Sequence_End, fill=PAR)) +
    geom_tile(colour="white",size=0.25) +
    geom_text(aes(label=formatC(PAR, digits = 3)), size=1.5) +
    scale_fill_gradientn(colours = colours,
                         values = splits,
                         na.value='#000000') +
    theme_bw() +
    ggtitle("3'C4 Conserved Sequences") +
    theme(plot.title = element_text(size=8,face="bold", hjust=0.5)) +
    scale_x_continuous(breaks = c(9:17), labels = c(9:17)) +
    xlab('Length (nt)') +
    ylab("Sequence End") +
    theme(axis.text.x = element_text(size=7,face="bold"), axis.title.x = element_text(size=7, face='bold')) + 
    theme(axis.text.y = element_text(size=7,face="bold"), axis.title.y = element_text(size=7, face='bold')) + 
    theme(legend.title =element_blank(), legend.text=element_text(size=7))
  hmC <- hmC + theme(legend.position='none')
  
  #####Save Plots
  hm <- plot_grid(hmA, hmU, hmC, legend, ncol=4, nrow=1, rel_widths = c(1,1,1,0.3))
  save_plot(paste(strain, 'Effect Modeling Heatmap.png'), hm, base_height = 2.3, base_width = 8, dpi=900)
}

setwd('./Effect Modeling')
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase')
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot')) %dopar% 
  strainEffectModels(strains[j], as.numeric(parav[,j]))
stopCluster(cl)

setwd('..')
