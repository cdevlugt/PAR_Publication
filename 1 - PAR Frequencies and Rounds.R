##### PAR Frequencies #####
require(doParallel)

parFreq <- function(df, strain=''){
  #obtain subunits
  if(strain=='Luciferase') {
    subs <- 1:4
  } else {
    subs <- 1:8
  }
  
  #obtain the number of reads that are realigned or not realigned
  Realigned <- sum(df[,c(subs)][df$rounds>0,])
  No_Realignment <- sum(df[,c(subs)][df$rounds==0,])
  
  #tabulate as percent
  parFreq <- data.frame(x = c(
    No_Realignment / (Realigned + No_Realignment) *100,
    Realigned / (Realigned + No_Realignment) *100
  ))
  row.names(parFreq) <- c("No Realignment", "Realigned")
  colnames(parFreq) <- strain
  
  #garbage collect
  rm(df, subs, Realigned, No_Realignment)
  gc()
  
  #return
  return(parFreq)
}

#setwd("") #location of the data if not default directory
names <- c('pr8', 'hk', 'bri', 'wsn', 'luc')
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase')
cluster <- length(names)
cl <- makeCluster(cluster)
registerDoParallel(cl)
parRates <- foreach(j = 1:cluster, .combine = cbind) %dopar% 
  parFreq(read.csv(paste(toupper(names[j]), '_All_Match.csv', sep=''), stringsAsFactors = F), strains[j])
stopCluster(cl)

parRates$Average <- rowMeans(parRates[,1:4])
parRates[,(ncol(parRates)+1)] <- rowMeans(parRates[,1:3])
colnames(parRates)[ncol(parRates)] <- "Average Without WSN"

#save
setwd("./Overview")
write.csv(parRates, 'PAR rates.csv')
setwd("..")

##### PAR Rounds #####
require(doParallel)

parRounds <- function(df, strain=''){
  #obtain subunits
  if(strain=='Luciferase') {
    subs <- 1:4
  } else {
    subs <- 1:8
  }
  
  #obtain the maximum number of PAR rounds
  mr <- max(df$rounds)
  
  #obtain the number of reads that are realigned or not realigned
  Reads <- sum(df[,c(subs)])
  
  #obtain the percentage that passes through each round
  parRounds <- matrix(ncol=1, nrow = (mr+1))
  for(i in 1:nrow(parRounds)){
    parRounds[i,1] <- sum(df[,c(subs)][df$rounds>=(i-1),]) / Reads *100
  }
  
  rm(df)
  gc()
  
  return(parRounds)
}

#setwd("") #location of the data if not default directory
names <- c('pr8', 'hk', 'bri', 'wsn', 'luc')
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase')
cluster <- length(names)
cl <- makeCluster(cluster)
registerDoParallel(cl)
parRnd <- foreach(j = 1:cluster) %dopar% 
  parRounds(read.csv(paste(toupper(names[j]), '_All_Match.csv', sep=''), stringsAsFactors = F), strains[j])
stopCluster(cl)

#get number of rows for combined data frame for rounds
rounds <- 1
for(i in 1:length(names)){
  if(nrow(parRnd[[i]])>rounds){
    rounds <- nrow(parRnd[[i]])
  }
}

#tabulate data
parMulti <- matrix(nrow = rounds, ncol=length(names), data=0)
for(i in 1:length(names)){
  for(j in 1:nrow(parRnd[[i]])){
    parMulti[j,i] <- parRnd[[i]][j,1]
  }
}

parMulti <- as.data.frame(parMulti)
colnames(parMulti) <- strains
row.names(parMulti) <- paste('Round', 0:(rounds-1))

#average
parMulti$Average <- rowMeans(parMulti[,1:4])

#average without WSN
parMulti$`Average Without WSN` <- rowMeans(parMulti[,1:3])

#save
setwd("./Overview")
write.csv(parMulti, 'PAR rounds.csv')
setwd("..")

##### Rounds Graphs #####
require(ggplot2)
require(cowplot)
require(doParallel)

roundsGraphsAndModels <- function(df, num = 1, roundLimit=4){
  #generates several graphs
  #Primers Salvaged is the percentage of primers which do not undergo PAR in a given round relative to all primers in that SPECIFIC round
  #Removal is the percentage of the total data set which is used for transcription (does not PAR) in a given round
  #Rounds is the percentage of primers that undergo PAR in that round, including those that undergo PAR in a subsequent round
  
  #get strain name
  strain <- colnames(df)[num]
  strain <- gsub("\\.", ' ', strain)
  
  #calculate salvage
  salvage <- data.frame(round=0:(nrow(df)-1), percent=df[,num], stringsAsFactors = F)
  salvage <- salvage[salvage$percent!=0,]
  
  salvage$salvage <- 0
  
  #calculate salvage
  for(i in 1:(nrow(salvage)-1)){
    salvage$salvage[i] <- ((salvage$percent[i] - sum(salvage$percent[(i+1):nrow(salvage)])) / salvage$percent[i] * 100)
  }
  salvage$salvage[salvage$salvage==0] <- 100
  salvage <- salvage[,-2]
  salvage$lnSalvage <- log(salvage$salvage)
  salvage <- salvage[salvage$round<=roundLimit,]
  
  #add strain and labels
  salvage$strain <- strain
  salvage$label <- as.character(signif(salvage$salvage, 3))
  
  #graph
  salvagePlot <- ggplot(data=salvage, aes(x=round, y=salvage, group=strain, colour=strain)) +
    geom_point(size=1) +
    scale_color_manual(breaks=strain,
                       values='red') +
    geom_smooth(data = salvage[salvage$round>=1,],
                aes(x=round, y=salvage, group=strain, colour=strain),
                method='lm', se=F, formula = y~x, size=0.5, colour='black', fullrange = T) +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(salvage, 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(0, roundLimit, by=1), limits = c(0,roundLimit+.5)) +
    scale_y_continuous(breaks = seq(75, 100, by=5), limits=c(75, 105)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    labs(x='Rounds of PAR', y='Primers Rescued (%)') +
    theme(legend.position = 'none') +
    theme(aspect.ratio = 1:1)
  
  #save data
  save_plot(paste(strain, 'Primers Salvaged.png'), plot=salvagePlot, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  
  
  #r2 values and p-values
  salvageStats <- data.frame(r2 = summary(lm(data = salvage[salvage$round>=1,], round ~ salvage))$r.squared,
                            pval = summary(lm(data = salvage[salvage$round>=1,], round ~ salvage))$coefficients[2,4],
                            average = mean(salvage$salvage[salvage$round!=0]))
  
  #save
  write.csv(file = paste(strain, 'Salvage Rate R2 and pvalue.csv'), salvageStats, row.names = F)
  
  
  
  
  #get primers that undergo transcription in a given round
  removal <- data.frame(round=0:(nrow(df)-1),percent=df[,num])
  for(i in 1:(nrow(df)-1)){
    removal[i,2] <- removal[i,2] - removal[(i+1),2]
  }
  removal <- removal[1:(roundLimit+1),]
  removal <- removal[removal$percent>0,]
  removal$lnPercent <- log(removal$percent) #may warning if a round has 0 or less primers
  
  #add strain
  removal$strain <- strain
  
  #r2 values and p-values
  r2val <- summary(lm(data = removal[removal$round>=1,], round ~ lnPercent))$r.squared
  pval <- summary(lm(data = removal[removal$round>=1,], round ~ lnPercent))$coefficients[2,4]
  
  #caluculate PAR frequency from lm 
  randPAR <- as.numeric(exp(coef(lm(data=removal[removal$round>=1,], formula = lnPercent ~ round)))[2]) *100
  
  #write
  write.csv(file=paste(strain, 'Removal R2 and pvalue.csv', sep=' '), data.frame(r2 = r2val, p = pval), row.names = F)
  write.csv(file=paste(strain, 'Removal Random PAR rate.csv', sep = ' '), data.frame(PAR = randPAR), row.names = F)
  
  #graph
  #generate removal correlation plot
  removalPlot <- ggplot(data=removal, aes(x=round, y=percent, group=strain, colour=strain)) +
    geom_point(size=1) +
    scale_color_manual(breaks=strain,
                       values='red') +
    geom_smooth(data = removal[removal$round>=1,],
                aes(x=round, y=percent, group=strain, colour=strain),
                method='glm', se=F, formula = y~x, size=0.5, method.args = list(family = gaussian(link = 'log')), 
                fullrange=T, colour='black') +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(percent, 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(0, roundLimit, by=1), limits=c(0,roundLimit+.5)) +
    scale_y_continuous(breaks = seq(0, 100, by=20), limits=c(0, 100)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    labs(x='Rounds of PAR', y='Primers (%)') +
    theme(legend.position = 'none') +
    theme(aspect.ratio = 1:1)
  
  #generate removal correlation plot on log scale
  logRemovalPlot <- ggplot(removal[removal$round>=0,], aes(x=round, y=lnPercent, group=strain, colour=strain)) +
    geom_point(size=1) +
    scale_color_manual(breaks=strain,
                       values='red') +
    geom_smooth(data = removal[removal$round>=1,],
                aes(x=round, y=lnPercent, group=strain, colour=strain),
                method='lm', se=F, formula = y~x, size=0.5, colour='black', fullrange = T) +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(lnPercent, 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(0, roundLimit, by=1), limits=c(0, roundLimit+.5)) +
    scale_y_continuous(breaks = seq(-6, 6, by=2), limits=c(-6.5,6.5)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    labs(x='Rounds of PAR', y='Primers (ln(%))') +
    theme(legend.position = 'none') +
    theme(aspect.ratio = 1:1)
  
  #save graphs
  save_plot(paste(strain, 'Removal.png'), plot=removalPlot, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  save_plot(paste(strain, 'Removal Log.png'), plot=logRemovalPlot, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  
  #delimit data
  df <- df[1:(roundLimit+1),]
  df <- df[df[,num]!=0,]
  
  #get strain name
  strain <- colnames(df)[num]
  strain <- gsub("\\.", ' ', strain)
  
  #extract rounds and percentage
  df <- df[,num][df[,num]!=0] #remove useless data
  roundsPlotdf <- data.frame(round = 0:(length(df)-1), percent = df, lnPercent = log(df), stringsAsFactors = F)
  
  #r2 values and p-values
  r2val <- summary(lm(data = roundsPlotdf[roundsPlotdf$round>=1,], round ~ lnPercent))$r.squared
  pval <- summary(lm(data = roundsPlotdf[roundsPlotdf$round>=1,], round ~ lnPercent))$coefficients[2,4]
  
  #caluculate PAR frequency from lm 
  randPAR <- as.numeric(exp(coef(lm(data=roundsPlotdf[roundsPlotdf$round>=1,], formula = lnPercent ~ round)))[2]) *100
  
  #write
  write.csv(file=paste(strain, 'Rounds R2 and pvalue.csv', sep=' '), data.frame(r2 = r2val, p = pval), row.names = F)
  write.csv(file=paste(strain, 'Random PAR rate.csv', sep = ' '), data.frame(PAR = randPAR), row.names = F)
  
  #graph
  #generate rounds correlation plot
  roundsPlot <- ggplot(data=roundsPlotdf, aes(x=round, y=percent, group=strain, colour=strain)) +
    geom_point(size=1) +
    scale_color_manual(breaks=strain,
                       values='red') +
    geom_smooth(data = roundsPlotdf[roundsPlotdf$round>=1,],
                aes(x=round, y=percent, group=strain, colour=strain),
                method='glm', se=F, formula = y~x, size=0.5, method.args = list(family = gaussian(link = 'log')), 
                fullrange=T, colour='black') +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(percent, 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(0, roundLimit, by=1), limits = c(0, roundLimit+.5)) +
    scale_y_continuous(breaks = seq(0, 100, by=20), limits=c(0, 110)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    labs(x='Rounds of PAR', y='Primers (%)') +
    theme(legend.position = 'none') +
    theme(aspect.ratio = 1:1)
  
  #generate rounds correlation plot on log scale
  logRoundsPlot <- ggplot(roundsPlotdf[roundsPlotdf$round>=0,], aes(x=round, y=lnPercent, group=strain, colour=strain)) +
    geom_point(size=1) +
    scale_color_manual(breaks=strain,
                       values='red') +
    geom_smooth(data = roundsPlotdf[roundsPlotdf$round>=1,],
                aes(x=round, y=lnPercent, group=strain, colour=strain),
                method='lm', se=F, formula = y~x, size=0.5, colour='black', fullrange = T) +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(lnPercent, 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(0, roundLimit, by=1), limits = c(0, roundLimit+.5)) +
    scale_y_continuous(breaks = seq(-6, 6, by=2), limits=c(-6.5,6.5)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    labs(x='Rounds of PAR', y='Primers (ln(%))') +
    theme(legend.position = 'none') +
    theme(aspect.ratio = 1:1)
  
  #save graphs
  save_plot(paste(strain, 'Rounds.png'), plot=roundsPlot, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  save_plot(paste(strain, 'Rounds Log.png'), plot=logRoundsPlot, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
}

#setwd("") #location of the data if not default directory
setwd('./Overview')
df <- read.csv('PAR rounds.csv', stringsAsFactors = F, row.names = 1)

cluster <- ncol(df)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:cluster, .inorder = F, .packages = c('ggplot2', 'cowplot')) %dopar%
  roundsGraphsAndModels(df, j, roundLimit = 3)
stopCluster(cl)
setwd('..')


##### Subunit PAR Frequencies #####
require(doParallel)

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
subunitPARFrequencies <- function(df, strain='', U=1:5, C=6:8, round2plusGroup=T){
  
  #subunit handling from strain... lets me doParallel the pre-processing
  if(is.na(U) && is.na(C) && strain=='Luciferase') {
    U = c(2,4)
    C = c(1,3)
  } else if (is.na(U) && is.na(C) && strain=='Hong Kong') {
    U = c(1:6)
    C = c(7:8)
  } else if (is.na(U) && is.na(C)) {
    U = c(1:5)
    C = c(6:8)
  }
  
  subs <- 1:length(c(U, C))
  
  #delimit data
  df <- df[,c(U,C,grep('rounds', colnames(df)))]
  df <- subunitSeperate(df, subs)
  
  #convert U and C to the names of the respective columns
  U <- colnames(df[,1:length(U)])
  C <- colnames(df[,(length(U)+1):length(subs)])
  
  if(round2plusGroup){
    #create matrix for number of reads per subunit per round
    subPAR <- matrix(ncol = (length(subs)+3), nrow = 3) #3 added for Template and Average
    colnames(subPAR) <- c(U, C, "3'U4", "3'C4", "All")
    
    #populate table
    for(i in 1:length(subs)){
      #the number of reads for a given subunit (i) matched in the given rounds
      subPAR[1,i] <- sum(df[,i][df$rounds==(0)]) #no trimming
      subPAR[2,i] <- sum(df[,i][df$rounds==(1)]) #one round
      subPAR[3,i] <- sum(df[,i][df$rounds>=(2)]) #2+ rounds
      subPAR[,i] <- subPAR[,i]/sum(df[,i]) * 100 #convert to percentage
    }
    
    #calculate the averages
    #this may error if the number of U, C, or total subunits is 1 or less
    #0
    subPAR[1,(length(subs)+1)] <- sum(df[,grep(paste(U, collapse = '|'), colnames(df))][df$rounds==(0),]) #3'U4
    subPAR[1,(length(subs)+2)] <- sum(df[,grep(paste(C, collapse = '|'), colnames(df))][df$rounds==(0),]) #3'C4
    subPAR[1,(length(subs)+3)] <- sum(df[,subs][df$rounds==(0),]) #All
    
    #1
    subPAR[2,(length(subs)+1)] <- sum(df[,grep(paste(U, collapse = '|'), colnames(df))][df$rounds==(1),]) #3'U4
    subPAR[2,(length(subs)+2)] <- sum(df[,grep(paste(C, collapse = '|'), colnames(df))][df$rounds==(1),]) #3'C4
    subPAR[2,(length(subs)+3)] <- sum(df[,subs][df$rounds==(1),]) #All
    
    #2
    subPAR[3,(length(subs)+1)] <- sum(df[,grep(paste(U, collapse = '|'), colnames(df))][df$rounds>=(2),]) #3'U4
    subPAR[3,(length(subs)+2)] <- sum(df[,grep(paste(C, collapse = '|'), colnames(df))][df$rounds>=(2),]) #3'C4
    subPAR[3,(length(subs)+3)] <- sum(df[,subs][df$rounds>=(2),]) #All
    
    row.names(subPAR) <- c('Round 0', 'Round 1', 'Round 2+')
    
  } else { #by total rounds
    #obtain the max rounds for data frame construction
    maxRounds <- max(df$rounds) + 1 #one added to account for round 0
    
    #create matrix for number of reads per subunit per round
    subPAR <- matrix(ncol = (length(subs)+3), nrow = maxRounds) #3 added for Template and Average
    colnames(subPAR) <- c(U, C, "3'U4", "3'C4", "All")
    
    #populate table
    for(i in 1:length(subs)){
      for(j in 1:nrow(subPAR)){
        subPAR[j,i] <- sum(df[,i][df$rounds==(j-1)]) #the number of reads for a given subunit (i) matched in the given round (j-1)
      }
      subPAR[,i] <- subPAR[,i]/sum(df[,i]) * 100 #convert to percentage
    }
    
    #calculate the averages
    #this may error if the number of U, C, or total subunits is 1 or less
    for(j in 1:nrow(subPAR)){
      subPAR[j,(length(subs)+1)] <- sum(df[,grep(paste(U, collapse = '|'), colnames(df))][df$rounds==(j-1),]) #3'U4
      subPAR[j,(length(subs)+2)] <- sum(df[,grep(paste(C, collapse = '|'), colnames(df))][df$rounds==(j-1),]) #3'C4
      subPAR[j,(length(subs)+3)] <- sum(df[,subs][df$rounds==(j-1),]) #All
    }
    row.names(subPAR) <- paste('Round', 0:(maxRounds-1))
  }
  
  #calculate percentages
  subPAR[,(length(subs)+1)] <- subPAR[,(length(subs)+1)]/sum(df[,grep(paste(U, collapse = '|'), colnames(df))]) * 100 #3'U4
  subPAR[,(length(subs)+2)] <- subPAR[,(length(subs)+2)]/sum(df[,grep(paste(C, collapse = '|'), colnames(df))]) * 100 #3'C4
  subPAR[,(length(subs)+3)] <- subPAR[,(length(subs)+3)]/sum(df[,subs]) * 100 #All
  
  rm(df)
  gc()
  
  return(subPAR)
}

#setwd("") #location of the data if not default directory
names <- c('pr8', 'hk', 'bri', 'wsn', 'luc')
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase')
cluster <- length(names)
cl <- makeCluster(cluster)
registerDoParallel(cl)
subPAR <- foreach(j = 1:cluster) %dopar% 
  subunitPARFrequencies(read.csv(paste(toupper(names[j]), '_All_Match.csv', sep=''), stringsAsFactors = F), strains[j], NA, NA, T)
stopCluster(cl)

#average
subPAR[[6]] <- subPAR[[1]]
for(i in 2:4){
  subPAR[[6]] <- subPAR[[6]] + subPAR[[i]]
}
subPAR[[6]] <- subPAR[[6]]/4

#average without wsn
subPAR[[7]] <- subPAR[[1]]
for(i in 2:3){
  subPAR[[7]] <- subPAR[[7]] + subPAR[[i]]
}
subPAR[[7]] <- subPAR[[7]]/3

#save
strains <- c(strains, 'Average', 'Average Without WSN')
setwd("./Overview")
for(i in 1:7){
  write.csv(subPAR[[i]], file=paste(strains[i], "Subunit PAR.csv"), row.names = T)
}
setwd('..')

##### Subunit PAR Graphs #####
require(ggplot2)
require(cowplot)
require(doParallel)

subunitPARGraphs <- function(df, parav, strain){
  
  #graph the data with the templates 
  colnames(df) <- gsub('\\.', '', colnames(df))
  colnames(df)[grep("U4", colnames(df))] <- "3'U4"
  colnames(df)[grep("C4", colnames(df))] <- "3'C4"
  
  graphData <- data.frame(round=rep(row.names(df),ncol(df)),
                          subunit=rep(colnames(df),each=3),
                          frequency=0)
  
  for(i in 1:ncol(df)){
    graphData[(3*(i-1)+1):(3*(i-1)+3),3] <- df[,i]
  }
  graphData$subunit <- factor(graphData$subunit, levels = c(colnames(df)))
  
  graphWithTemplate <- ggplot() +
    geom_bar(data = graphData,aes(group=subunit, y=frequency, x=subunit, fill=round), stat='identity') +
    scale_fill_manual(values = c(hcl(h=seq(15,375, length=(3+1))[3], c=100, l=65),
                        hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65),
                        hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65))) +
    geom_hline(yintercept = parav, colour = 'black', linetype = 'dashed', size=0.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), panel.border = element_blank()) +
    scale_y_continuous(limits = c(0, 100.001), breaks=seq(0, 100, by=20), expand = c(0,0)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text.y = element_text(size=8), 
          axis.text.x = element_text(size=8, angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_text(size=8))+
    labs(x='viral mRNAs', y='Percentage of mRNA (%)')
  
  #graph without templates
  colnames(df) <- gsub('\\.', '', colnames(df))
  
  df <- df[,-c(grep(".4", colnames(df)))] #remove the templates
  graphData <- data.frame(round=rep(row.names(df),ncol(df)),
                          subunit=rep(colnames(df),each=3),
                          frequency=0)
  
  for(i in 1:ncol(df)){
    graphData[(3*(i-1)+1):(3*(i-1)+3),3] <- df[,i]
  }
  graphData$subunit <- factor(graphData$subunit, levels = c(colnames(df)))
  
  graph <- ggplot() +
    geom_bar(data = graphData,aes(group=subunit, y=frequency, x=subunit, fill=round), stat='identity') +
    scale_fill_manual(values = c(hcl(h=seq(15,375, length=(3+1))[3], c=100, l=65),
                                 hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65),
                                 hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65))) +
    geom_hline(yintercept = parav, colour = 'black', linetype = 'dashed', size=0.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), panel.border = element_blank()) +
    scale_y_continuous(limits = c(0, 100.001), breaks=seq(0, 100, by=20), expand = c(0,0)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text.y = element_text(size=8), 
          axis.text.x = element_text(size=8, angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_text(size=8))+
    labs(x='viral mRNAs', y='Percentage of mRNA (%)')
    
  
  #save and make a legend
  if(strain=='Average'){
    legend <- get_legend(graph)
    save_plot("Subunit PAR Legend.png", plot=legend, dpi=900)
  }
  graph <- graph + 
    theme(legend.position = 'none')
  graphWithTemplate <- graphWithTemplate +
    theme(legend.position = 'none')
  
  save_plot(paste(strain, "Subunit PAR.png"), plot=graph, base_height = 3.5, base_width = 2, dpi = 900)
  save_plot(paste(strain, "Subunit and Template PAR.png"), plot=graphWithTemplate, base_height = 3.5, base_width = 2, dpi = 900)
}

setwd("./Overview")
#obtain strain PAR frequencies
parav <- read.csv("PAR rates.csv", stringsAsFactors = F, row.names=1)[1,]

strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase', 'Average', 'Average Without WSN')
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot')) %dopar% 
  subunitPARGraphs(read.csv(paste(strains[j], 'Subunit PAR.csv', sep=' '), stringsAsFactors = F, row.names = 1), as.numeric(parav[1,j]), strains[j])
stopCluster(cl)

setwd('..')

