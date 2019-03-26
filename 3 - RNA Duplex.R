##### Obtain Sequence End and Template Data #####
require(doParallel)
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
  df <- df[,c(U, C, grep('rounds', colnames(df)), grep('Trim_Sequence', colnames(df)), grep('NT_coupure', colnames (df)))]
  U <- 1:length(U)
  C <- length(U) + 1:length(C)
  
  #3'U1 directed nucleotide
  df$Trim_Sequence <- substr(df$Trim_Sequence, nchar(df$Trim_Sequence), nchar(df$Trim_Sequence))
  df$Trim_Sequence <- gsub('T', 'U', df$Trim_Sequence) #DNA to RNA
  
  #obtain G at +1
  df$NT_coupure <- substr(df$NT_coupure, 1, 1)
  df$NT_coupure <- gsub('[^G]', '', df$NT_coupure) #get rid of all non-guanine
  
  #sequence end, the sequence used to prime transcription
  df$Seq_End <- paste(df$Trim_Sequence, df$NT_coupure, sep='')
  
  #series handing
  df$series <- 'No_Realignment'
  df$series[df$rounds>=1] <- 'Realigned'
  df$series <- factor(df$series, levels=c('No_Realignment', 'Realigned'))
  
  #drop useless columns
  df <- df[,-c(grep('rounds', colnames(df)), grep('Trim_Sequence', colnames(df)), grep('NT_coupure', colnames (df)))]
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
  
  #factor data
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

#Sequence Ends data
setwd("./RNA Duplex")
sequenceEndsPAR <- function(df){
  #requires use of dfPreProcess in this script
  strain <- df$strain[1] #obtain strain name
  
  #make table for 3'U1 only
  x <- data.frame(table(substr(df$Seq_End, 1, 1), df$series))
  x[,1] <- as.character(x[,1]) #fix to character
  x[,3] <- as.numeric(as.character(x[,3])) #fix to numeric
  x <- x[,c(2,1,3)] #reorder
  colnames(x) <- c('series', 'Seq_End', 'Frequency')
  
  #make table for 3'U1 and 3'C2
  y <- data.frame(table(df$Seq_End, df$series))
  y[,1] <- as.character(y[,1]) #fiy to character
  y[,3] <- as.numeric(as.character(y[,3])) #fiy to numeric
  y <- y[,c(2,1,3)] #reorder
  colnames(y) <- c('series', 'Seq_End', 'Frequency')
  
  rm(df)
  gc()
  
  #chars
  cx <- unique(x$Seq_End) #3'U1
  cy <- unique(y$Seq_End) #factoring in internal
  
  #par frequencies 3'U1 only
  parx <- x[c(1:length(cx)),]
  parx$series <- 'PAR'
  
  for(i in 1:length(cx)){
    parx[i,3] <- x[length(cx)+i,3] / (x[i,3] + x[length(cx)+i,3]) * 100
  }
  parx$Internal <- F #no distinction for 3'C2 directed G
  
  #par frequencies factoring in the G
  pary <- y[c(1:length(cy)),]
  pary$series <- 'PAR'
  
  for(i in 1:length(cy)){
    pary[i,3] <- y[length(cy)+i,3] / (y[i,3] + y[length(cy)+i,3]) * 100
  }
  pary$Internal <- T #distinction for 3'C2 directed G
  pary <- pary[order(nchar(pary$Seq_End), pary$Seq_End),]
  
  #merge tables
  par <- rbind(parx, pary)
  
  #label strain
  par$strain <- strain
  
  #fix series
  par$series <- 'Terminal'
  par$series[par$Internal==F] <- "3'U1 Only"
  par$series[nchar(par$Seq_End)==2] <- "Internal"
  
  par <- par[,-c(grep('Internal', colnames(par)))]
  
  #label maximum energy
  #A = -1.73, G = -1.57, C = -1.40, U = -1.2
  par$Energy <- c(-1.73, -1.40, -1.57, -1.2)#kcal/mol
  
  return(par)
}

#obtain PAR frequencies for sequence ends
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
sePAR <- foreach(j = 1:cluster) %dopar% 
  sequenceEndsPAR(dfs[[j]])
stopCluster(cl)

#obtain average
sePAR[[6]] <- sePAR[[1]]
sePAR[[6]]$Frequency <- 0
sePAR[[6]]$strain <- 'Average'

#calculate average
for(i in 1:4){ #4 is up to WSN
  sePAR[[6]]$Frequency <- sePAR[[6]]$Frequency + sePAR[[i]]$Frequency
}
sePAR[[6]]$Frequency <- sePAR[[6]]$Frequency/4

#obtain average without wsn
sePAR[[7]] <- sePAR[[1]]
sePAR[[7]]$Frequency <- 0
sePAR[[7]]$strain <- 'Average Without WSN'

#calculate average
for(i in 1:4){ #4 is up to WSN
  sePAR[[7]]$Frequency <- sePAR[[7]]$Frequency + sePAR[[i]]$Frequency
}
sePAR[[7]]$Frequency <- sePAR[[7]]$Frequency/3

#save
for(i in 1:length(sePAR)){
  write.csv(sePAR[[i]], paste(sePAR[[i]]$strain[1], 'Priming PAR.csv'), row.names = F)
}
rm(sePAR)
gc()

#vRNA template data
templatePAR <- function(df){
  #requires use of dfPreProcess in this script
  strain <- df$strain[1] #obtain strain name
  
  #make table for genome segments 
  x <- data.frame(table(df$Subunit, df$series, df$Template))
  x[,1] <- as.character(x[,1]) #fix to character
  x[,4] <- as.numeric(as.character(x[,4])) #fix to numeric
  x <- x[x[,4]!=0,] #drop zero values
  x <- x[,c(2,1,3,4)] #reorder
  colnames(x) <- c('series', 'Gene', 'Template', 'Frequency')
  
  #make table for 3'U1 and 3'C2
  y <- data.frame(table(df$Template, df$series))
  y[,1] <- as.character(y[,1]) #fiy to character
  y[,3] <- as.numeric(as.character(y[,3])) #fiy to numeric
  y <- y[,c(2,1, 1, 3)] #reorder
  colnames(y) <- c('series', 'Gene', 'Template', 'Frequency')
  
  rm(df)
  gc()
  
  #join tables
  z <- rbind(x, y)
  z <- z[order(z$series),]
  
  #genes
  cz <- unique(z$Gene)
  
  #par frequencies 3'U1 only
  par <- z[c(1:length(cz)),]
  par$series <- 'PAR'
  
  for(i in 1:length(cz)){
    par[i,4] <- z[length(cz)+i,4] / (z[i,4] + z[length(cz)+i,4]) * 100
  }
  
  #label strain
  par$strain <- strain
  
  return(par)
}

#obtain PAR frequencies for vRNA templates
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
tPAR <- foreach(j = 1:cluster) %dopar% 
  templatePAR(dfs[[j]])
stopCluster(cl)

#obtain average
tPAR[[6]] <- tPAR[[1]]
tPAR[[6]]$Frequency <- 0
tPAR[[6]]$strain <- 'Average'

#calculate average
for(i in 1:4){ #first 4 strains
  tPAR[[6]]$Frequency <- tPAR[[6]]$Frequency + tPAR[[i]]$Frequency
}
tPAR[[6]]$Frequency <- tPAR[[6]]$Frequency/4

#obtain average without WSN
tPAR[[7]] <- tPAR[[1]]
tPAR[[7]]$Frequency <- 0
tPAR[[7]]$strain <- 'Average Without WSN'

#calculate average
for(i in 1:4){ #first 4 strains
  tPAR[[7]]$Frequency <- tPAR[[7]]$Frequency + tPAR[[i]]$Frequency
}
tPAR[[7]]$Frequency <- tPAR[[7]]$Frequency/3

#save
for(i in 1:length(tPAR)){
  write.csv(tPAR[[i]], paste(tPAR[[i]]$strain[1], 'Template PAR.csv'), row.names = F)
}

#All genes
templatePopulations <- tPAR[[1]]
if(length(tPAR)>1){
  for(i in 2:length(tPAR)){
    templatePopulations <- rbind(templatePopulations, tPAR[[i]])
  }
}
write.csv(templatePopulations, "Template Populations PAR.csv", row.names = F)

rm(tPAR, dfs, templatePopulations)
gc()

setwd('..')

##### Sequence Ends Correlations #####
require(ggplot2)
require(cowplot)
require(doParallel)

sequenceEndCorrelationPlotsAndModels <- function(df, parav){
  
  #strain
  strain <- df$strain[1]
  
  #get r2 and p values; ordered as No Distinction, Terminal, Internal
  r2val <- c(summary(lm(data = df[df$series=="3'U1 Only",], Energy ~ Frequency))$r.squared,
             summary(lm(data = df[df$series=="Terminal",], Energy ~ Frequency))$r.squared,
             summary(lm(data = df[df$series=="Internal",], Energy ~ Frequency))$r.squared,
             NA)
  pval <- c(summary(lm(data = df[df$series=="3'U1 Only",], Energy ~ Frequency))$coefficients[2,4],
            summary(lm(data = df[df$series=="Terminal",], Energy ~ Frequency))$coefficients[2,4],
            summary(lm(data = df[df$series=="Internal",], Energy ~ Frequency))$coefficients[2,4],
            t.test(df[df$series=='Terminal',3], df[df$series=='Internal',3], paired=T)$p.value)
  cors <- c("Average", 'Terminal', 'Internal', 'Terminal vs. Internal')
  
  write.csv(data.frame(type=cors, r2 = r2val, p = pval), file = paste(strain, "Priming Nucleotide r2 and p.csv"), row.names = F)
  
  #use the Terminal and Internal lm() to obtain the free energy difference from priming internally
  tef <- lm(data = df[df$series=="Terminal",], Energy ~ Frequency) 
  tfe <- lm(data = df[df$series=="Terminal",], Frequency ~ Energy)
  ief <- lm(data = df[df$series=="Internal",], Energy ~ Frequency) 
  ife <- lm(data = df[df$series=="Internal",], Frequency ~ Energy)
  
  #use the predicted PAR frequencies from the Internal Energy ~ Frequency model
  #multiply by the slope of the Terminal Energy ~ Frequency model
  #add the intercept of the Terminal Energy ~ Frequency model
  #subtract the predicted Energies of the Internal Energy ~ Frequency model
  #this gives the difference between the two linear models
  #Simply: PARInternal * EnergyTerminal/PARTerminal + InterceptTerminal - EnergyInternal #distance in kcal/mol between 2 models
  gEnergy <- mean((predict(ife)*tef$coefficients[2] + tef$coefficients[1])-predict(ief))
  write.csv(gEnergy, file=paste(strain, "Energy bonus of priming with a 3'C2 directed G.csv"), row.names = F)
  
  #3'U1 only plot
  U1only <- ggplot(data=df[df$series=="3'U1 Only",], aes(x=Energy, y=Frequency, group=strain, colour=strain)) +
    geom_hline(yintercept = parav, colour = 'blue', linetype = 'dashed', size=0.5) +
    geom_point(size=1) +
    scale_color_manual(values='red') +
    geom_smooth(method='lm', se=F, formula = y~x, size=0.5, colour='black') +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(Frequency, 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(limits = c(min(df$Energy)-.05, max(df$Energy)+.05),
                       breaks=c(-1.73, -1.57, -1.4, -1.20),
                       labels =paste(c('A', 'G', 'C', 'U'), ' (', c('-1.73', '-1.57', '-1.40', '-1.20'), ')', sep='')) +
    scale_y_continuous(limits = c(0, 45), breaks=seq(0, 40, by=10)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text.y = element_text(size=8), 
          axis.text.x = element_text(size=8),
          axis.title = element_text(size=10))+
    labs(x=bquote(Delta~G[37]^{o}*' (kcal/mol)'), y='Prime-and-Realign Frequency (%)') +
    theme(aspect.ratio = 1:1) +
    theme(legend.position = 'none')
  
  #G distinction
  gDist <- ggplot(data=df[df$series!="3'U1 Only",], aes(x=Energy, y=Frequency, group=series, colour=series)) +
    geom_hline(yintercept = parav, colour = 'blue', linetype = 'dashed', size=0.5) +
    geom_point(size=1) +
    geom_line(stat = 'smooth', method='lm', se=F, formula = y~x, size=0.5, alpha=0.4) +
    scale_color_manual(values=c("#009292", '#d98000')) +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(Frequency, 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(limits = c(min(df$Energy)-.05, max(df$Energy)+.05),
                       breaks=c(-1.73, -1.57, -1.4, -1.20),
                       labels =paste(c('A', 'G', 'C', 'U'), ' (', c('-1.73', '-1.57', '-1.40', '-1.20'), ')', sep='')) +
    scale_y_continuous(limits = c(0, 45), breaks=seq(0, 40, by=10)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text.y = element_text(size=8), 
          axis.text.x = element_text(size=8),
          axis.title = element_text(size=10))+
    labs(x=bquote(Delta~G[37]^{o}*' (kcal/mol)'), y='Prime-and-Realign Frequency (%)') +
    theme(legend.position = 'none') +
    theme(aspect.ratio = 1:1) 
  
  #all 3
  types <- unique(df$series)
  all3 <- ggplot(data=df, aes(x=Energy, y=Frequency, group=series, colour=series)) +
    geom_hline(yintercept = parav, colour = 'blue', linetype = 'dashed', size=0.5) +
    geom_point(size=1) +
    geom_line(stat = 'smooth', method='lm', se=F, formula = y~x, size=0.5, alpha=0.4) +
    geom_text(data=df[df$series=="Terminal",], aes(x=Energy, y=Frequency, group=series, colour=series, label = gsub("\\.$", "",as.character(formatC(Frequency, 3, flag="###")))),
              vjust=-.5, hjust=-0.2, size =2.5, show.legend = F) +
    geom_text(data=df[df$series=="3'U1 Only",], aes(x=Energy, y=Frequency, group=series, colour=series, label = gsub("\\.$", "",as.character(formatC(Frequency, 3, flag="###")))),
              vjust=0.5, hjust=-0.2, nudge_x = -.07, size =2.5, show.legend = F) +
    geom_text(data=df[df$series=="Internal",], aes(x=Energy, y=Frequency, group=series, colour=series, label = gsub("\\.$", "",as.character(formatC(Frequency, 3, flag="###")))),
              vjust=1, hjust=-0.2, size =2.5, show.legend = F) +
    scale_color_manual(breaks = c("3'U1 Only", "Terminal", "Internal"),
                       values=c('#ff6db6', "#009292", '#d98000'), 
                       labels = c("Average", "Terminal", "Internal")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(limits = c(min(df$Energy)-.07, max(df$Energy)+.05),
                       breaks=c(-1.73, -1.57, -1.4, -1.20),
                       labels =paste(c('A', 'G', 'C', 'U'), ' (', c('-1.73', '-1.57', '-1.40', '-1.20'), ')', sep='')) +
    scale_y_continuous(limits = c(0, 45), breaks=seq(0, 40, by=10)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text.y = element_text(size=8), 
          axis.text.x = element_text(size=8),
          axis.title = element_text(size=10))+
    labs(x=bquote(Delta~G[37]^{o}*' (kcal/mol)'), y='Prime-and-Realign Frequency (%)') +
    theme(aspect.ratio = 1:1)
  
  #steal legend
  primingLegend <- get_legend(all3)
  
  #remove legend
  all3 <- all3 + theme(legend.position = 'none')
  
  #save
  save_plot(paste(strain, "3'U1 Directed Nucleotide.png"), plot=U1only, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  save_plot(paste(strain, 'Priming Sequence.png'), plot=gDist, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  save_plot(paste(strain, "3'U1 directed Nucleotide and Priming Sequence.png"), plot=all3, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  if(strain=='Average'){
    save_plot("Priming Legend.png", plot=primingLegend, base_height = 3.5, base_width = 2, base_aspect_ratio = 1:1, dpi = 900)
  }
}

#setwd("") #if not default directory

#get PAR rate for strains
setwd("./Overview")
parav <- read.csv("PAR rates.csv", stringsAsFactors = F, row.names=1)[2,]
setwd('..')

#analyze data
setwd("./RNA Duplex")
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase', 'Average', 'Average Without WSN')
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot')) %dopar% 
  sequenceEndCorrelationPlotsAndModels(read.csv(paste(strains[j], 'Priming PAR.csv', sep=' '), stringsAsFactors = F), as.numeric(parav[1,j]))
stopCluster(cl)

setwd('..')

##### Template Populations #####
require(ggplot2)
require(cowplot)
require(doParallel)

templatePopulationsAndStats <- function(df, parav, strain){
  
  #label subunits
  df$shape <- 49
  df$shape[df$Gene=="PB1"] <- 50
  df$shape[df$Gene=="PA"] <- 51
  df$shape[df$Gene=="HA"] <- 52
  df$shape[df$Gene=="NP"] <- 53
  df$shape[df$Gene=="NA."] <- 54
  df$shape[df$Gene=="M"] <- 55
  df$shape[df$Gene=="NS"] <- 56
  df$shape[df$Gene=="NPhetero"] <- 65
  df$shape[df$Gene=="PB1hetero"] <- 66
  df$shape[df$Gene=="NPhomo"] <- 67
  df$shape[df$Gene=="PB1homo"] <- 68
  
  #set colours based on strain
  df$colour <- "#ffff6d" #luc
  df$colour[df$strain=='Puerto Rico'] <- "#009292" #pr8
  df$colour[df$strain=='Hong Kong'] <- "#24ff24" #hk
  df$colour[df$strain=='WSN'] <- "#b66dff" #wsn
  df$colour[df$strain=='Brisbane'] <- "#ff6db6" #bri
  
  if(strain=='Average'){
    df <- df[df$strain != 'Average' & df$strain != 'Luciferase' & df$strain != 'Average Without WSN',]
  } else if(strain=='Average Without WSN'){
    df <- df[df$strain != 'Average' & df$strain != 'Luciferase' & df$strain != 'WSN'& df$strain != 'Average Without WSN',]
  } else {
    df <- df[df$strain == strain,]
  }
  
  #population p-values
  pval <- t.test(x=df$Frequency[df$Template=="3'U4"], y=df$Frequency[df$Template=="3'C4"])$p.value
  write.csv(pval, file = paste(strain, "Template Groups p-value.csv"), row.names = F)
  
  #factor Template to use as x-units
  df$Template <- factor(df$Template, levels = c("3'U4", "3'C4"))
  
  #add x-units
  df$x <- 1
  df$x[df$Template=="3'C4"] <- 2
  
  stats <- data.frame(series = c("3'U4", "3'C4"),
                      mean = c(mean(df$Frequency[df$Template=="3'U4"]), mean(df$Frequency[df$Template=="3'C4"])),
                      sd = c(sd(df$Frequency[df$Template=="3'U4"]), sd(df$Frequency[df$Template=="3'C4"])),
                      x = c(1, 2))
  
  #get colours and strains
  cols <- unique(df$colour)
  strains <- unique(df$strain)
  subunits <- unique(df$Gene)
  subunits <- gsub("\\.", "", subunits)
  shape <- as.character(unique(df$shape))
  subunits <- subunits[order(shape)]
  
  genePlot <- ggplot() +
    geom_hline(yintercept = parav, colour = 'blue', linetype = 'dashed', size=0.5, show.legend = F) +
    geom_errorbar(data=stats, aes(ymin=mean-2*sd, ymax=mean+2*sd, x=x), size=0.5, show.legend = F) +
    stat_summary(fun.y='mean', data=stats, colour="#000000", geom="errorbar", aes(group=series,x=x, y=mean,ymax=..y.., ymin=..y..), width=0.5, linetype="dashed", show.legend=F) +
    geom_point(data=df, aes(x=x, y=Frequency, group=Template, colour=colour, shape=shape), size=4, position = position_jitterdodge(jitter.width = 0.8)) +
    scale_continuous_identity(aesthetics = 'shape', guide = 'legend',
                              name = "Gene",
                              breaks=min(shape):max(shape),
                              labels=subunits) +
    scale_colour_identity(name = "Strain", guide = "legend", breaks = cols, labels = strains)+
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks=c(1,2), limit=c(0.5,2.5), labels = c("3'U4 (-5.09 kcal/mol)", "3'C4 (-5.71 kcal/mol)")) +
    scale_y_continuous(limits = c(0, 40), breaks=seq(0, 40, by=10)) +
    theme(legend.text = element_text(size=8),
          axis.text.y = element_text(size=8), 
          axis.text.x = element_text(size=8),
          axis.title = element_text(size=10))+
    labs(x="vRNA Template", y='Prime-and-Realign Frequency (%)') +
    theme(aspect.ratio = 1:1)
  genePlot
  
  geneLegend <- get_legend(genePlot)
  
  genePlot <- genePlot +
    theme(legend.position = 'none')
  
  #save
  save_plot(paste(strain, "Template PAR.png"), plot=genePlot, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  if(strain=="Average" | strain=="Luciferase"){
    save_plot(paste(strain, "Template PAR Legend.png"), plot=geneLegend, base_height = 4, base_width = 2, base_aspect_ratio = 1:1, dpi = 900)
  }
  
  rm(df, stats, genePlot, geneLegend)
  gc()
}

#setwd("") #if not default directory

#get PAR rate for strains
setwd("./Overview")
parav <- read.csv("PAR rates.csv", stringsAsFactors = F, row.names=1)[2,]
setwd('..')

#read data
setwd("./RNA Duplex")
df <- read.csv('Template Populations PAR.csv', stringsAsFactors = F)

#drop templates as genes
df <- df[df$Gene!="3'U4" & df$Gene!="3'C4",]

strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase', 'Average')
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot')) %dopar% 
  templatePopulationsAndStats(df, as.numeric(parav[1,j]), strains[j])
stopCluster(cl)

gc()
setwd('..')
