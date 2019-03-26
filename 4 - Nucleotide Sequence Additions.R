##### Obtain Nucleotide Sequence Additions Data #####
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
dfPreProcess <- function(df, strain='', U=1:5, C=6:8, guanineLength = T) {
  
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
  df <- df[,c(U, C, grep('rounds', colnames(df)), grep('Trim_Len', colnames(df)), grep('NT_coupure', colnames (df)), grep('Seq_Trim_R', colnames(df)))]
  U <- 1:length(U)
  C <- length(U) + 1:length(C)
  df <- df[df$rounds>0,]
  gc() #free up RAM from parallel
  
  #obtain G at +1
  df$NT_coupure <- substr(df$NT_coupure, 1, 1)
  
  #subunit handling
  subs <- c(U, C)
  subs <- unique(subs)
  subs <- subs[order(subs)]
  
  #seperate reads by subunit
  df <- subunitSeperate(df, subs)
  
  #label templates
  df$Template <- "3'U4"
  df$Template[rowMeans(df[,c(C)])>0] <- "3'C4"
  
  #add strain
  df$strain <- strain
  
  #obtain rounds
  rounds <- max(df$rounds)
  
  #seperate by rounds
  x <- df[1,] #clone df
  x <- x[-1,] #remove all rows
  x <- x[, - grep('Seq_Trim_R[^1]', colnames(x))] #remove addition rounds columns
  colnames(x)[grep('Seq_Trim_R1', colnames(x))] <- 'Addition' #lable that column as addition
  
  #Construct df
  for(i in 1:rounds){
    #create temporay df "y" with data to add for the given round
    y <- df[df$rounds>=i & df[,grep(paste('Seq_Trim_R', i, sep=''), colnames(df))]!="",]
    y$rounds <- i #lower rounds to the ith round
    
    if(i == 2){ #length changes
      y$Trim_Len <- y$Trim_Len + nchar(y[,c(grep('Seq_Trim_R1', colnames(y)))])
    } else if(i > 2){
      y$Trim_Len <- y$Trim_Len + nchar(apply(y[,c(grep('Seq_Trim_R1', colnames(y)):grep(paste('Seq_Trim_R', (i - 1), sep=''), colnames(y)))], 1, paste, collapse=""))
    } else if(guanineLength && i == 1) {
      y$Trim_Len[y$NT_coupure=="G"] <- y$Trim_Len[y$NT_coupure=="G"] +1
    }
    
    #remove columns and rename the sequence added in that round to match the colnames of x
    y <- y[,-(grep(paste("Seq_Trim_R[^", i, "]", sep=''), colnames(y)))]
    colnames(y)[grep(paste("Seq_Trim_R", i, sep=''), colnames(y))] <- 'Addition'
    
    #merge and delete data, gc() to restore RAM
    x <- rbind(x, y)
    rm(y)
    gc()
  }
  
  colnames(x)[grep('Trim_Len', colnames(x))] <- 'Length'
  
  return(x)
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

sequenceAdditions <- function(df, target, range = NA, addAll = F){
  #returns a table of nucleotide sequence additions based on the target column
  #requires the data be processed by dfPreProcess
  
  #subunits
  if(df$strain[1]=='Luciferase') {
    subs = 1:4
  } else {
    subs = 1:8
  }
  
  #target column
  if(is.character(target)){
    target <- grep(target, colnames(df))
  }
  varName <- colnames(df)[target]
  
  #range
  if(is.na(range)[1]){
    range <- unique(df[,target])
    range <- range[order(range)]
  }
  
  #construct data.frame divisions based on target
  sa <- data.frame(target=df[1,target], Addition='G', reads = 0, stringsAsFactors = F)
  sa <- sa[-1,]
  
  #construct data.frame for all reads
  all <- data.frame(target='All',
                    Addition = c('G', 'GC', 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA'),
                    reads= 0, stringsAsFactors = F)
  
  #add data do data.frame
  for(i in 1:length(range)){
    y <- data.frame(target=range[i],
                    Addition = c('G', 'GC', 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA'),
                    reads= 0, stringsAsFactors = F)
    for(j in 1:nrow(y)){
      if(nrow(df[df[,target]==range[i] & df$Addition==y$Addition[j],])>=1){
        y$reads[j] <- sum(df[,c(subs)][df[,target]==range[i] & df$Addition==y$Addition[j],])
      }
    }
    
    #merge with sa
    sa <- rbind(sa, y)
    
    #add reads to all
    all$reads <- all$reads + y$reads
    
    #manage the RAM
    rm(y)
    gc()
  }
  
  #add columns for byTarget, byAddition, and byReads
  sa$byTarget <- 0
  sa$byAddition <- 0
  sa$byReads <- sa$reads / sum(sa$reads)
  
  all$byTarget <- all$reads / sum(all$reads)
  all$byAddition <- 1 #(dividing all additions by addition will give 1)
  all$byReads <- all$byTarget #(these are the same when applied to all data)
  
  for(i in 1:length(range)){
    sa$byTarget[sa$target==range[i]] <- sa$reads[sa$target==range[i]] / sum(sa$reads[sa$target==range[i]])
    sa$byAddition[sa$target==range[i]] <- sa$reads[sa$target==range[i]] / all$reads
  }
  
  if(addAll){
    sa <- rbind(sa, all)
  }
  
  #rename target column
  colnames(sa)[1] <- varName
  
  #keep strain name
  sa$strain <- df$strain[1]
  
  rm(df)
  gc()
  
  return(sa)
}

#obtain sequence additions from length, template, and rounds
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
saLengths <- foreach(j = 1:cluster) %dopar% 
  sequenceAdditions(dfs[[j]], 'Length', 9:17, F)
saTemplate <- foreach(j = 1:cluster) %dopar% 
  sequenceAdditions(dfs[[j]], 'Template', NA, T)
saRounds <- foreach(j = 1:cluster) %dopar% 
  sequenceAdditions(dfs[[j]], 'rounds', 1:4, F)
stopCluster(cl)

#get average
saLengths[[6]] <- saLengths[[1]]
saTemplate[[6]] <- saTemplate[[1]]
saRounds[[6]] <- saRounds[[1]]

saLengths[[6]]$strain <- 'Average'
saTemplate[[6]]$strain <- 'Average'
saRounds[[6]]$strain <- 'Average'

#calculate average
if(length(strains)>1){
  for(i in 2:4){
    saLengths[[6]][,3:6] <- saLengths[[6]][,3:6] + saLengths[[i]][,3:6]
    saTemplate[[6]][,3:6] <- saTemplate[[6]][,3:6] + saTemplate[[i]][,3:6]
    saRounds[[6]][,3:6] <- saRounds[[6]][,3:6] + saRounds[[i]][,3:6]
  }
}
#average all except reads (column 3)
saLengths[[6]][,4:6] <- saLengths[[6]][,4:6] / 4
saTemplate[[6]][,4:6] <- saTemplate[[6]][,4:6] / 4
saRounds[[6]][,4:6] <- saRounds[[6]][,4:6] / 4

#get average without WSN
saLengths[[7]] <- saLengths[[1]]
saTemplate[[7]] <- saTemplate[[1]]
saRounds[[7]] <- saRounds[[1]]

saLengths[[7]]$strain <- 'Average Without WSN'
saTemplate[[7]]$strain <- 'Average Without WSN'
saRounds[[7]]$strain <- 'Average Without WSN'

#calculate average
if(length(strains)>1){
  for(i in 2:3){
    saLengths[[7]][,3:6] <- saLengths[[7]][,3:6] + saLengths[[i]][,3:6]
    saTemplate[[7]][,3:6] <- saTemplate[[7]][,3:6] + saTemplate[[i]][,3:6]
    saRounds[[7]][,3:6] <- saRounds[[7]][,3:6] + saRounds[[i]][,3:6]
  }
}
#average all except reads (column 3)
saLengths[[7]][,4:6] <- saLengths[[7]][,4:6] / 3
saTemplate[[7]][,4:6] <- saTemplate[[7]][,4:6] / 3
saRounds[[7]][,4:6] <- saRounds[[7]][,4:6] / 3

#save
setwd('./Nucleotide Sequence Additions')
strains <- c(strains, "Average", "Average Without WSN")
for(i in 1:length(strains)){
  write.csv(file = paste(strains[i], 'Additions by Length.csv'), x = saLengths[[i]], row.names = F)
  write.csv(file = paste(strains[i], 'Additions by Template.csv'), x = saTemplate[[i]], row.names = F)
  write.csv(file = paste(strains[i], 'Additions by Rounds.csv'), x = saRounds[[i]], row.names = F)
}

setwd('..')

##### Length and Nucleotide Sequence Addition #####
require(doParallel)
require(ggplot2)
require(ggpubr)
require(cowplot)

sequenceAdditionLengthPlotCorrelation <- function(df){
  #extract strain
  strain <- df$strain[1]
  
  #factor addition for ordering
  df$Addition <- factor(df$Addition, levels = c('G', 'GC', 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA'))
  
  #graph margins to bind graphs
  t <- -0.995
  b <- -0.995
  l <- -1.88
  r <- -1.88
  
  hml <- ggplot(data=df, aes(x=Length, y=Addition, fill=byReads)) +
    geom_tile(colour="white",size=0.25) +
    geom_text(data=df[df$byReads>=0.001 & or(df$byReads <=0.1, df$byReads >=0.8),],
              aes(x=Length, y=Addition, label=as.character(formatC(signif(byReads*100, 3), 3, flag='###'))),
              size=2, colour = "#000000", show.legend = F) +
    geom_text(data=df[df$byReads>=0.001 & and(df$byReads >=0.1, df$byReads <=0.8),],
              aes(x=Length, y=Addition, label=as.character(formatC(signif(byReads*100, 3), 3, flag='###'))),
              size=2, colour = "#FFFFFF", show.legend = F) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme_bw() +
    scale_x_continuous(labels = seq(min(df$Length), max(df$Length), by=1),
                       breaks = seq(min(df$Length), max(df$Length), by=1)) +
    xlab("Length (nt)") +
    ylab("Sequence Added") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=8))+
    guides(fill = guide_colorbar(title.vjust = -20)) +
    theme(plot.margin=unit(c(t, r,0,0), 'cm')) #top, right, bottom, left
  
  #get legend
  hmlegend <- get_legend(hml)
  hml <- hml + theme(legend.position = 'none')
  
  #x (length) bar plot
  xdf <- data.frame(x=unique(df$Length), y=0)
  for(i in 1:nrow(xdf)){
    xdf$y[i] <- sum(df$byReads[df$Length==xdf$x[i]])
  }
  
  xbar <- ggplot(data=xdf, aes(x=x, y=y, fill=y)) + 
    geom_bar(stat='identity') +
    scale_y_continuous(limits=c(0,1)) +
    ylab("") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    clean_theme() +
    theme(plot.margin = unit(c(t,r,b,l), 'cm')) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme(legend.position='none')
  
  #y (additions) bar plot
  ydf <- data.frame(x=levels(df$Addition), y=0)
  for(i in 1:nrow(ydf)){
    ydf$y[i] <- sum(df$byReads[df$Addition==ydf$x[i]])
  }
  ydf$x <- factor(ydf$x, levels = levels(df$Addition))
  
  ybar <- ggplot(data=ydf, aes(x=x, y=y, fill=y)) + 
    geom_bar(stat='identity') +
    scale_y_continuous(limits=c(0,1)) +
    ylab("") +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    clean_theme() +
    theme(plot.margin = unit(c(t,r,b,l), 'cm')) + 
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme(legend.position='none') +
    coord_flip()
  
  graph <- plot_grid(xbar, NULL, NULL, hml, ybar, NULL, ncol=3, nrow=2, rel_widths=c(2.5, 1, 0.35, 2.5, 1, 0.35), rel_heights=c(1,1,1,2,2,2), align = 'hv')
  save_plot(paste(strain, 'Length and Nucleotide Sequence Additions.png'), graph, base_height = 2, base_width = 4)
  if(strain=='Average'){
    save_plot(paste(strain, "Sequence Additions Heatmaps Legend.png"), plot=hmlegend, base_height = 3.5, base_width = 1, base_aspect_ratio = 1:1, dpi = 900)
  }
  
  #length and addition length correlations
  cdf <- data.frame(Length = rep(unique(df$Length), 1, each=6),
                    Addition = rep(1:6, length(unique(df$Length))),
                    reads=0)
  
  #unfactor df$Addition so nchar will work
  df$Addition <- as.character(df$Addition)
  
  for(i in 1:nrow(cdf)){
    cdf$reads[i] <- sum(df$reads[df$Length==cdf$Length[i] & nchar(df$Addition)==cdf$Addition[i]])
  }
  
  #remove zero values
  cdf <- cdf[cdf$reads!=0,]
  cdf <- cdf[rep(row.names(cdf), cdf$reads),]
  
  lenAddCor <- cor.test(cdf$Length, cdf$Addition, method = 'pearson')
  
  write.csv(file = paste(strain, 'Length and Addition Length Correlation.csv'), x = data.frame(r2 = lenAddCor$estimate, pval = lenAddCor$p.value))
  
  rm(cdf, df, ydf, xdf, xbar, ybar, hml, hmlegend)
  gc()
}

#analyze data
setwd('./Nucleotide Sequence Additions')
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase', 'Average', "Average Without WSN")
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot', 'ggpubr')) %dopar% 
  sequenceAdditionLengthPlotCorrelation(read.csv(paste(strains[j], 'Additions by Length.csv', sep=' '), stringsAsFactors = F))
stopCluster(cl)

setwd('..')

##### Template and Nucleotide Sequence Additions #####
require(ggplot2)
require(ggpubr)
require(cowplot)
require(doParallel)

sequenceAdditionTemplatePlot <- function(df){
  #get strain
  strain <- df$strain[1]
  
  #graph margins to bind graphs
  t <- -0.095
  b <- t
  
  #factor data for ggplot ordering
  df$Addition <- factor(df$Addition, levels = c('G', 'GC', 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA'))
  df$Template <- factor(df$Template, levels = c("All", "3'U4", "3'C4"))
  df$x <- 1
  df$x[df$Template=="3'U4"] <- 2
  df$x[df$Template=="3'C4"] <- 3
  
  hmt <- ggplot(data=df, aes(x=x, y=Addition, fill=byTarget)) +
    geom_tile(colour="white",size=0.25) +
    geom_text(data=df[df$byTarget>=0.001 & or(df$byTarget <=0.1, df$byTarget >=0.8),],
              aes(x=x, y=Addition, label=as.character(formatC(signif(byTarget*100, 3), 3, flag='###'))),
              size=2, colour = "#000000", show.legend = F) +
    geom_text(data=df[df$byTarget>=0.001 & and(df$byTarget >=0.1, df$byTarget <=0.8),],
              aes(x=x, y=Addition, label=as.character(formatC(signif(byTarget*100, 3), 3, flag='###'))),
              size=2, colour = "#FFFFFF", show.legend = F) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme_bw() +
    scale_x_continuous(labels = levels(df$Template),
                       breaks = 1:3) +
    xlab("vRNA Template") +
    ylab("Sequence Added") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=8))+
    guides(fill = guide_colorbar(title.vjust = -20)) +
    theme(plot.margin = unit(c(t,0,0,0), 'cm'))
  
  xdf <- data.frame(x = c(1, 2, 3),
                    y = c(0, 
                          sum(df$reads[df$Template=="3'U4"])/ sum(df$reads[df$Template=="All"]),
                          sum(df$reads[df$Template=="3'C4"])/ sum(df$reads[df$Template=="All"])))
  
  xbar <- ggplot(data=xdf, aes(x=x, y=y, fill=y)) + 
    geom_bar(stat='identity') +
    scale_y_continuous(limits=c(0,1), breaks = c(0,1)) +
    ylab("") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    clean_theme() +
    theme(plot.margin = unit(c(0,0,b,0), 'cm')) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme(legend.position='none')
  
  #get legend
  hmlegend <- get_legend(hmt)
  hmt <- hmt + theme(legend.position = 'none')
  
  graph <- plot_grid(xbar,hmt, ncol=1, nrow=2, rel_widths=c(2.5, 2.5), rel_heights=c(1,2), align = 'v')
  save_plot(paste(strain, 'Template and Nucleotide Sequence Additions.png'), graph, base_height = 2, base_width = 1.6)
  
  rm(df, hmt, hmlegend, graph)
  gc()
}

#analyze data
setwd('./Nucleotide Sequence Additions')
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase', 'Average', "Average Without WSN")
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot', 'ggpubr')) %dopar% 
  sequenceAdditionTemplatePlot(read.csv(paste(strains[j], 'Additions by Template.csv', sep=' '), stringsAsFactors = F))
stopCluster(cl)

setwd('..')

##### Rounds and Nucleotide Sequence Additions #####
require(ggplot2)
require(ggpubr)
require(cowplot)
require(doParallel)

sequenceAdditionRoundsPlot <- function(df){
  #get strain
  strain <- df$strain[1]
  
  #graph margins to bind graphs
  t <- -0.095
  b <- t
  
  #factor data for ggplot ordering
  df$Addition <- factor(df$Addition, levels = c('G', 'GC', 'GCA', 'GCG', 'GCAA', 'GCGA', 'GCAAA', 'GCGAA', 'GCAAAA', 'GCGAAA'))
  
  hmr <- ggplot(data=df, aes(x=rounds, y=Addition, fill=byTarget)) +
    geom_tile(colour="white",size=0.25) +
    geom_text(data=df[df$byTarget>=0.001 & or(df$byTarget <=0.1, df$byTarget >=0.8),],
              aes(x=rounds, y=Addition, label=gsub("\\.$", "", as.character(formatC(signif(byTarget*100, 3), 3, flag='###')))),
              size=2, colour = "#000000", show.legend = F) +
    geom_text(data=df[df$byTarget>=0.001 & and(df$byTarget >=0.1, df$byTarget <=0.8),],
              aes(x=rounds, y=Addition, label=gsub("\\.$", "", as.character(formatC(signif(byTarget*100, 3), 3, flag='###')))),
              size=2, colour = "#FFFFFF", show.legend = F) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme_bw() +
    xlab("Rounds of PAR") +
    ylab("Sequence Added") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=8))+
    guides(fill = guide_colorbar(title.vjust = -20)) +
    theme(plot.margin = unit(c(t,0,0,0), 'cm'))
  
  xdf <- data.frame(x = c(1:max(df$rounds)),
                    y = 0)
  
  for(i in 1:max(df$rounds)){
    xdf[i, 2] <- sum(df$reads[df$rounds==i])/ sum(df$reads)
  }
  
  xbar <- ggplot(data=xdf, aes(x=x, y=y, fill=y)) + 
    geom_bar(stat='identity') +
    scale_y_continuous(limits=c(0,1), breaks = c(0,1)) +
    ylab("") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    clean_theme() +
    theme(plot.margin = unit(c(0,0,b,0), 'cm')) +
    scale_fill_gradientn(limits=c(0,1), values=c(0,0.03125,0.0625,0.125,0.5,1), colours = c('#FFFFFF', 'lightgrey','lightblue', 'steelblue', 'purple', 'red'), na.value='#FFFFFF') +
    theme(legend.position='none')
  
  #get legend
  hmlegend <- get_legend(hmr)
  hmr <- hmr + theme(legend.position = 'none')
  
  graph <- plot_grid(xbar, hmr, ncol=1, nrow=2, rel_widths=c(2.5, 2.5), rel_heights=c(1,2), align = 'v')
  save_plot(paste(strain, 'Rounds and Nucleotide Sequence Additions.png'), graph, base_height = 2, base_width = 1.8)
  
  rm(df, hmr, hmlegend)
  gc()
}

setwd('./Nucleotide Sequence Additions')
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase', 'Average', "Average Without WSN")
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot', 'ggpubr')) %dopar% 
  sequenceAdditionRoundsPlot(read.csv(paste(strains[j], 'Additions by Rounds.csv', sep=' '), stringsAsFactors = F))
stopCluster(cl)

setwd('..')

##### Addition Length and Frequency #####
require(doParallel)
require(ggplot2)
require(cowplot)
additionLengthFrequency <- function(df){
  #get strain
  strain <- df$strain[1]
  
  #limit to all
  df <- df[df$Template=='All',]
  
  #addition lengths
  al <- data.frame(Length = c(1:6), frequency = 0, stringsAsFactors = F)
  
  for(i in 1:6){
    al$frequency[i] <- sum(df$byReads[nchar(df$Addition)==i]) * 100
  }
  
  #generate addition length correlation plot
  addLengthCorrelations <- ggplot(data=al[al$frequency!=0,], aes(x=Length, y=frequency, group=strain, colour=strain)) +
    geom_point(size=1) +
    scale_color_manual(values='red') +
    geom_smooth(data = al[3:6,],
                aes(x=Length, y=frequency, group=strain, colour=strain),
                method='glm', se=F, formula = y~x, size=0.5, method.args = list(family = poisson(link = 'log')), 
                fullrange=F, colour='black') +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(frequency, 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(1, 6, by=1), limits = c(1, 7)) +
    scale_y_continuous(breaks = seq(0, 80, by=20), limits=c(0, 80)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    labs(x='Addition Length (nt)', y='Frequency (%)') +
    theme(aspect.ratio = 1:1) + 
    theme(legend.position = 'none')
  
  #generate addition length correlation plot on log scale
  logAddLengthCorrelations <- ggplot(data=al[!is.infinite(log(al$frequency)),], aes(x=Length, y=log(frequency), group=strain, colour=strain)) +
    geom_point(size=1) +
    scale_color_manual(values='red') +
    geom_smooth(data = al[3:6,],
                aes(x=Length, y=log(frequency), group=strain, colour=strain),
                method='lm', se=F, formula = y~x, size=0.5, 
                fullrange=F, colour='black') +
    geom_text(aes(label = gsub("\\.$", "",as.character(formatC(log(frequency), 3, flag="###")))),
              vjust=0, hjust=-0.2, size =2.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks = seq(1, 6, by=1), limits = c(1, 7)) +
    scale_y_continuous(breaks = seq(-1, 5, by=1), limits = c(-1, 5)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=8),
          axis.text = element_text(size=8), 
          axis.title = element_text(size=10))+
    labs(x='Addition Length (nt)', y='Frequency ln(%)') +
    theme(aspect.ratio = 1:1) + 
    theme(legend.position = 'none')
  
  save_plot(paste(strain, 'Addition Length and Frequency.png'), plot=addLengthCorrelations, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  save_plot(paste(strain, 'Addition Length and ln Frequency.png.png'), plot=logAddLengthCorrelations, base_height = 3.5, base_width = 3.5, base_aspect_ratio = 1:1, dpi = 900)
  
  al$lf <- log(al$frequency)
  alLog <- lm(data = al[al$Length>=3,], lf ~ Length)
  
  write.csv(data.frame(r2=summary(alLog)$r.squared, pval = summary(alLog)$coefficients[2,4]),
            paste(strain, "Addition Length Frequency.csv"), row.names = F)
}

#analyze data
setwd('./Nucleotide Sequence Additions')
strains <- c('Puerto Rico', 'Hong Kong', 'Brisbane', 'WSN', 'Luciferase', 'Average', "Average Without WSN")
cluster <- length(strains)
cl <- makeCluster(cluster)
registerDoParallel(cl)
foreach(j = 1:length(strains), .inorder = F, .packages = c('ggplot2', 'cowplot', 'ggpubr')) %dopar% 
  additionLengthFrequency(read.csv(paste(strains[j], 'Additions by Template.csv', sep=' '), stringsAsFactors = F))
stopCluster(cl)

setwd('..')
