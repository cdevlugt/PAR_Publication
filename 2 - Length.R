##### Obtain Length Data #####
require(doParallel)
memory.limit(size=(128000+100000))
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
dfPreProcess <- function(df, strain='', U=1:5, C=6:8) {
  
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
  
  #delimit the data
  df <- df[,c(U, C, grep('rounds', colnames(df)), grep('Length', colnames(df)), grep('Trim_Len', colnames(df)), grep('NT_coupure', colnames (df)))]
  U <- 1:length(U)
  C <- length(U) + 1:length(C)
  
  #obtain G at +1
  df$NT_coupure <- substr(df$NT_coupure, 1, 1)
  
  #modify length
  df$Trim_Len[grep('G', df$NT_coupure)] <- df$Trim_Len[grep('G', df$NT_coupure)] + 1
  df <- df[,-(grep('NT_coupure', colnames(df)))]
  
  #series handing
  df$series <- 'No_Realignment'
  df$series[df$rounds>=1] <- 'Realigned'
  
  #handle post realignment
  df2 <- df[df$series=='Realigned',]
  
  #fix the length of both data frames
  df <- df[,-(grep('Length', colnames(df)))]
  colnames(df)[grep('Trim_Len', colnames(df))] <- 'Length'
  df2 <- df2[,-(grep('Trim_Len', colnames(df2)))]
  
  #modify series name in df2
  df2$series <- 'Post-Realignment'
  
  #Merge and drop data
  df <- rbind(df, df2)
  rm(df2)
  gc()
  
  #subunit handling
  subs <- c(U, C)
  subs <- unique(subs)
  subs <- subs[order(subs)]
  subNames <- colnames(df)[subs]
  C <- colnames(df)[C]
  
  df <- subunitSeperate(df, subs=subs, decompress = T)
  rm(subs)
  
  #add template
  df$Template <- "3'U4"
  for(i in 1:length(C)){
    df$Template[df$Subunit==C[i]] <- "3'C4"
  }
  rm(C)
  
  #factor data
  df$series <- factor(df$series, levels=c('No_Realignment', 'Realigned', 'Post-Realignment'))
  df$Template <- factor(df$Template, levels=c("3'U4", "3'C4"))
  
  df$strain <- strain
  gc()
  return(df)
}

#setwd("") #location of the data if not default directory
names <- c('pr8', 'hk', 'bri', 'wsn', 'luc')
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase')
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
dfs <- foreach(j = 1:cluster) %dopar% 
  dfPreProcess(read.csv(paste(toupper(names[j]), '_All_Match.csv', sep=''), stringsAsFactors = F), strains[j], NA, NA)
stopCluster(cl)

#length populations tables
setwd("./Length")
lengthPopulations <- function(df, lower = 9, upper = 17){
  #requires use of dfPreProcess in this script
  
  #limit to lengths between upper and lower args
  df <- df[df$Length>=lower & df$Length<=upper,]
  
  #make table
  x <- data.frame(table(df$Length, df$series))
  x[,1] <- as.numeric(as.character(x[,1])) #fix to numeric
  x[,3] <- as.numeric(as.character(x[,3])) #fix to numeric
  x <- x[,c(2,1,3)] #reorder
  colnames(x) <- c('series', 'Length', 'Frequency')
  
  #total reads, ignoring the duplicated realigned series
  totalReads <- sum(x[,3][x$series!='Post-Realignment'])
  
  x[,3] <- x[,3] / totalReads *100 #calculate percentage
  
  #label strain
  x$strain <- df$strain[1]
  
  rm(df)
  gc()
  return(x)
}

cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
pops <- foreach(j = 1:cluster) %dopar% 
  lengthPopulations(dfs[[j]], 9, 17)
stopCluster(cl)

#calculate average
pops[[6]] <- pops[[1]]
pops[[6]]$Frequency <- 0
pops[[6]]$strain <- 'Average'

#calculate average
for(i in 1:4){
  pops[[6]]$Frequency <- pops[[6]]$Frequency + pops[[i]]$Frequency
}
pops[[6]]$Frequency <- pops[[6]]$Frequency/4

#calculate average without wsn
pops[[7]] <- pops[[1]]
pops[[7]]$Frequency <- 0
pops[[7]]$strain <- 'Average Without WSN'

#calculate average
for(i in 1:3){
  pops[[7]]$Frequency <- pops[[7]]$Frequency + pops[[i]]$Frequency
}
pops[[7]]$Frequency <- pops[[7]]$Frequency/3

#save
write.csv(pops[[6]], 'Average Length Populations.csv', row.names = F)
write.csv(pops[[7]], 'Average Without WSN Length Populations.csv', row.names = F)
for(i in 1:(length(pops)-1)){
  write.csv(pops[[i]], paste(pops[[i]]$strain[1], 'Length Populations.csv'), row.names = F)
}

rm(pops)
gc()

#length PAR frequencies
lengthPAR <- function(df, lower = 9, upper = 17){
  #requires use of dfPreProcess in this script
  
  #limit to lengths between upper and lower args
  df <- df[df$Length>=lower & df$Length<=upper,]
  df <- df[df$series!='Post-Realignment',]
  df$series <- factor(df$series, levels=c('No_Realignment', 'Realigned'))
  strain <- df$strain[1]
  
  #make table
  x <- data.frame(table(df$Length, df$series))
  x[,1] <- as.numeric(as.character(x[,1])) #fix to numeric
  x[,3] <- as.numeric(as.character(x[,3])) #fix to numeric
  x <- x[,c(2,1,3)]
  colnames(x) <- c('series', 'Length', 'Frequency')
  
  rm(df)
  gc()
  
  #par frequencies
  par <- x[c(1:length(lower:upper)),]
  par$series <- 'PAR'
  
  for(i in 1:length(lower:upper)){
    par[i,3] <- x[length(lower:upper)+i,3] / (x[i,3] + x[length(lower:upper)+i,3]) * 100
  }
  
  #label strain
  par$strain <- strain
  
  return(par)
}

cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
pars <- foreach(j = 1:cluster) %dopar% 
  lengthPAR(dfs[[j]], 9, 17)
stopCluster(cl)

#get average
pars[[6]] <- pars[[1]]
pars[[6]]$Frequency <- 0
pars[[6]]$strain <- 'Average'

#calculate average
for(i in 1:4){
  pars[[6]]$Frequency <- pars[[6]]$Frequency + pars[[i]]$Frequency
}
pars[[6]]$Frequency <- pars[[6]]$Frequency/4

#get average without WSN
pars[[7]] <- pars[[1]]
pars[[7]]$Frequency <- 0
pars[[7]]$strain <- 'Average Without WSN'

#calculate average
for(i in 1:3){
  pars[[7]]$Frequency <- pars[[7]]$Frequency + pars[[i]]$Frequency
}
pars[[7]]$Frequency <- pars[[7]]$Frequency/3

#save
write.csv(pars[[6]], 'Average Length PAR.csv', row.names = F)
write.csv(pars[[7]], 'Average Without WSN Length PAR.csv', row.names = F)
for(i in 1:(length(pars)-1)){
  write.csv(pars[[i]], paste(pars[[i]]$strain[1], 'Length PAR.csv'), row.names = F)
}

rm(dfs, cl, cluster, pars, pops, names, strains)
gc()
setwd('..')

##### Populations #####
require(ggplot2)
require(cowplot)
require(doParallel)

#Import data to graph
setwd("./Length")
lengthPopulationsGraph <- function(df, min=8, max=17){
  #make data on opposite side for realigned data
  df$Frequency[df$series!='No_Realignment'] <- df$Frequency[df$series!='No_Realignment'] *-1
  df$series <- gsub('-', '_', df$series)
  
  #isoalte strain name
  strain <- df$strain[1]
  
  #check for missing data between min and max, and add these as 0
  for(i in 1:length(min:max)){
    if(nrow(df[df$Length==(min:max)[i],])==0){
      df <- rbind(df,
                  data.frame(
                    series=c('No_Realignment', 'Realigned', 'Post_Realignment'),
                    Length = (min:max)[i],
                    Frequency = 0,
                    strain = df$strain[1],
                    stringsAsFactors = F
                  ))
    }
  }
  df$series <- factor(df$series, levels = c('No_Realignment', 'Realigned', 'Post_Realignment'))
  df <- df[order(df$series, df$Length),]
  
  lengthPopGraph <- ggplot(data=df, aes(group=series, x=Length, y=Frequency, fill=series)) + 
    geom_density(alpha=.9, stat='identity', position='identity') +
    ylab("Frequency (%)") + 
    scale_y_continuous(breaks=seq(-10, 30, by=10),
                       limits=c(-10, 35),
                       labels=c("10", "0", "10", "20", "30")) +
    scale_x_continuous(breaks=seq(min,max, by=1),
                       limits=c(min,max)) +
    xlab("Length (nt)") +
    coord_flip() +
    labs(fill="Sequence Type") +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_fill_manual(name='legend', breaks=c('No_Realignment', 'Realigned', 'Post_Realignment'), 
                      values=c(No_Realignment=hcl(h=seq(15,375, length=(3+1))[2], c=100, l=65),
                               Realigned=hcl(h=seq(15,375, length=(3+1))[1], c=100, l=65),
                               Post_Realignment=hcl(h=seq(15,375, length=(3+1))[3], c=100, l=65)),
                      labels=c(No_Realignment='NoPAR',
                               Realigned='Pre-PAR',
                               Post_Realignment='Post-PAR')) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    theme(aspect.ratio = 1:1)
  
  #steal legend
  lengthPopLegend <- get_legend(lengthPopGraph)
  
  #remove legend
  lengthPopGraph <- lengthPopGraph + theme(legend.position = 'none')
  
  #save
  save_plot(paste(strain,'Length Populations.png'), plot=lengthPopGraph, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  
  #save legend
  if(strain=='Average'){
    save_plot('Length Populations Legend.png', plot=lengthPopLegend, base_height = 3.5, base_width = 2, dpi = 900)
  }
  rm(df, lengthPopGraph, lengthPopLegend, min, max)
  gc()
}

#import and graph data
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase', 'Average', 'Average Without WSN')
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot')) %dopar% 
  lengthPopulationsGraph(read.csv(paste(strains[j], 'Length Populations.csv', sep=' '), stringsAsFactors = F))
stopCluster(cl)

setwd('..')

##### Correlations Plot #####
require(ggplot2)
require(cowplot)

correlationsPlotAndModels <- function(df, parav, min=9, max=17){
  
  #get strain name
  strain <- df$strain[1]
  
  #get R^2 and pvalue 9 to 13
  r2val <- summary(lm(data = df[df$Length<=13,], log(Frequency) ~ Length))$r.squared
  pval <- summary(lm(data = df[df$Length<=13,], log(Frequency) ~ Length))$coefficients[2,4]
  
  #generate length correlation plot
  lengthCorrelations <- ggplot(data=df[df$Frequency!=0,], aes(x=Length, y=Frequency, group=strain, colour=strain)) +
    geom_hline(yintercept = parav, colour = 'blue', linetype = 'dashed', size=0.5) +
    geom_point(size=1) +
    scale_color_manual(values='red') +
    geom_smooth(data = df,
                aes(x=Length, y=Frequency, group=strain, colour=strain),
                method='glm', se=F, formula = y~x, size=0.5, method.args = list(family = poisson(link = 'log')), 
                fullrange=T, colour='black') +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(Frequency, 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(min, max, by=1), limits = c(min, max+1)) +
    scale_y_continuous(breaks = seq(0, 100, by=20), limits=c(0, 100)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    labs(x='Length (nt)', y='Prime-and-Realign Frequency (%)') +
    theme(aspect.ratio = 1:1) + 
    theme(legend.position = 'none')
  
  #generate length correlation plot on log scale
  logLengthCorrelations <- ggplot(data=df[!is.infinite(log(df$Frequency)),], aes(x=Length, y=log(Frequency), group=strain, colour=strain)) +
    geom_hline(yintercept = log(parav), colour = 'blue', linetype = 'dashed', size=0.5) +
    geom_point(size=1) +
    scale_color_manual(values='red') +
    geom_smooth(data = df[df$Length<=13,],
                aes(x=Length, y=log(Frequency), group=strain, colour=strain),
                method='lm', se=F, formula = y~x, size=0.5, 
                fullrange=F, colour='black') +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(log(Frequency), 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(min, max, by=1), limits = c(min, max+1)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    labs(x='Length (nt)', y='Prime-and-Realign Frequency ln(%)') +
    theme(aspect.ratio = 1:1) + 
    theme(legend.position = 'none')
  
  save_plot(paste(strain, 'Length and PAR.png'), plot=lengthCorrelations, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  save_plot(paste(strain, 'Length and ln PAR.png'), plot=logLengthCorrelations, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  write.csv(data.frame(stat = c('R^2', 'p'),
                       val = c(r2val, pval),
                       stringsAsFactors = F),
            file=paste(strain, 'Length PAR p and R2 9 to 13.csv'), row.names = F)
  
  #full range
  r2val <- summary(lm(data = df[!is.infinite(log(df$Frequency)),], log(Frequency) ~ Length))$r.squared
  pval <- summary(lm(data = df[!is.infinite(log(df$Frequency)),], log(Frequency) ~ Length))$coefficients[2,4]
  
  write.csv(data.frame(stat = c('R^2', 'p'),
                       val = c(r2val, pval),
                       stringsAsFactors = F),
            file=paste(strain, 'Length PAR p and R2 9 to 17.csv'), row.names = F)
  
}

#get PAR rate for strains
setwd("./Overview")
parav <- read.csv("PAR rates.csv", stringsAsFactors = F, row.names=1)[2,]
setwd('..')

#Import data to graph
setwd("./Length")
df <- read.csv('Puerto Rico Length PAR.csv')

#import and graph data
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase', 'Average', 'Average Without WSN')
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot')) %dopar% 
  correlationsPlotAndModels(read.csv(paste(strains[j], 'Length PAR.csv', sep=' '), stringsAsFactors = F), as.numeric(parav[1,j]))
stopCluster(cl)

setwd('..')
