##### Prime and Realign Determinations #####
#one of the first scripts I programmed for my MSc

##### Strains (PR8, HK, WSN, Brisbane) #####
require(compiler)
enableJIT(3)

tslClean <- function(table){ 
  # Cleans up the transcript_support_level (tsl) of a given table to numeric values only
  #
  # Args:
  #	table:	30 column table containing the data
  #
  # Returns:
  #	table:	the table as specified in args with the tsl adjusted to numeric values
  
  #assign all NA or "" values to 6
  table$transcript_support_level[is.na(table$transcript_support_level)] <- 6
  table$transcript_support_level[table$transcript_support_level==""] <- 6
  table$transcript_support_level[table$transcript_support_level==""] <- 6
  table$transcript_support_level[startsWith(table$transcript_support_level,"N")] <- 6
  
  #cuts down annotation to a single number
  table$transcript_support_level[startsWith(table$transcript_support_level,"5")] <- 5	
  table$transcript_support_level[startsWith(table$transcript_support_level,"4")] <- 4
  table$transcript_support_level[startsWith(table$transcript_support_level,"3")] <- 3
  table$transcript_support_level[startsWith(table$transcript_support_level,"2")] <- 2
  table$transcript_support_level[startsWith(table$transcript_support_level,"1")] <- 1
  return(table)
}
keepR <- function(table, col=1, val=0, val2=NA, t="==", chr=FALSE){
  # Keeps Rows (keepR) within a specified range
  #
  # Args:
  #	table:	a data.frame to perform operations on, this data.frame will have only rows
  #			except that meet the specified parameters
  #	col: 	column number which the range of values kept will be applied to
  #	val: 	value to keep, or rnage limiter number 1
  #	val2:	range limiter number 2; will be set as val by default if needed
  #	t:  type of operation to be performed to keep rows where this value is TRUE
  # chr: spefies if the values kept are to be relted to string length
  #
  # Returns:
  #	table:	the table passed following removal of rows where "t" argued opperation is false
  if(!chr) {
    if(t=="=="){
      table <- table[table[,col]==val,]
    } else if(t=="!=") {
      table <- table[table[,col]!=val,]
    } else if(t=="<") {
      table <- table[table[,col]<val,]
    } else if(t=="<=") {
      table <- table[table[,col]<=val,]
    } else if(t==">") {
      table <- table[table[,col]>val,]
    } else if(t==">=") {
      table <- table[table[,col]>=val,]
    } else if(t=="na") {
      table <- table[is.na(table[,col]),]
    } else if(t=="!na") {
      table <- table[!(is.na(table[,col])),]
    } else if(t=="<>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]<val) & (table[,col]>val2),]
    } else if(t=="><") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]>val) | (table[,col]<val2),]
    } else if(t=="<=>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]<=val) & (table[,col]>=val2),]
    } else if(t==">=<") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]>=val) | (table[,col]<=val2),]
    } else {
      print("invalid type (t)")
    }
  } else if(chr) {
    if(t=="=="){
      table <- table[nchar(table[,col])==val,]
    } else if(t=="!=") {
      table <- table[nchar(table[,col])!=val,]
    } else if(t=="<") {
      table <- table[nchar(table[,col])<val,]
    } else if(t=="<=") {
      table <- table[nchar(table[,col])<=val,]
    } else if(t==">") {
      table <- table[nchar(table[,col])>val,]
    } else if(t==">=") {
      table <- table[nchar(table[,col])>=val,]
    } else if(t=="<>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])<val) & (nchar(table[,col])>val2),]
    } else if(t=="><") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])>val) | (nchar(table[,col])<val2),]
    } else if(t=="<=>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])<=val) & (nchar(table[,col])>=val2),]
    } else if(t==">=<") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])>=val) | (nchar(table[,col])<=val2),]
    } else {
      print("invalid type (t)")
    }
  }
  return(table)
}
trimData <- function(df, prel=10, tsl=4, art=10){
  # trims data tables to a given limitations for further analysis
  #
  # Args:
  #	df:	a data.frame to perform operations on, removing data above the cut-offs
  #	prel: 	cut-off value based on distance to the transcription start site
  #	tsl: 	transcription support level cut-off (requires cleaned up TSL data)
  #	art:	number of reads for a given subunit as to not be considered background
  #
  # Returns:
  #	df:	the table passed following removal of rows outside of cut-offs
  
  df <- tslClean(df)	#converts transcript_support_level data tonumeric
  df <- keepR(df, col=30, val= tsl, t="<=")	#keeps only rows with a transcript_support_level less than the tsl arg.
  df <- keepR(df, col=6, val=prel, t="<=>")	#keeps only rows with within +/- prel arg. of a transcription start site
  df <- df[(df[,10]>=art | df[,11]>=art | df[,12]>=art | df[,13]>=art | df[,14]>=art | df[,15]>=art | df[,16]>=art | df[,17]>=art),]	#removes rows where no single read count is greater than the art arg.
  df$Round <- df$Round - 1  #decrease rounds to reflect number of trims
  df[,10:17][df[,10:17] < art] <- 0	#sets any values below the art arg. to 0 in rows that have not been removed
  colnames(df)[9] <- "rounds"
  colnames(df)[3] <- "Sequence"
  df$Trim_Sequence <- df$Sequence
  df$Added_Sequence <- ""
  df$Length <- nchar(df$Sequence)
  df$Trim_Len <- nchar(df$Trim_Sequence)
  df$Len_Adjust <- 0
  subs <- 10:17
  others <- (1:ncol(df))[!((1:ncol(df)) %in% subs)]
  df <- df[,c(subs, others)]
  
  return(df)
}
numVectorReverse <- function(df){
  #Takes a vector or data frame of numbers and reverses the order forcing all zero values to the end
  # ex., takes 3, 2, 1, 0 and returns 1, 2, 3, 0
  #
  #Args:
  # df: a data frame or vector
  #
  #Returns:
  # df: dataframe containing rows with reversed numbers
  # as.vector(c(t(tb[1,]))):  vector containing reversed numbers
  
  if(is.vector(df)){  #determins arg type and return type
    retVect = T
    df <- data.frame(t(df), stringsAsFactors = F)
    names <- NA
  } else {  #data frame return info
    retVect = F
    names <- c(colnames(df))
  }
  
  row.names(df) <- c(1:nrow(df))  #keeps row names
  tb <- data.frame(df[,ncol(df)]) #creates data frame
  tb <- tb[,-1] #empties data frame
  z <- ncol(df) #ncol of the passed df, for 0 handling
  allData <- F  #boolean for deterimining if all data has been examined
  
  while(!(allData)){  #until all columns, or all data
    if(length(data.frame(df[,colSums(df)>0], stringsAsFactors = F)[1,])==0){  #if there are no non-zero values left
      for(i in ncol(tb):(z-1)){ #adds columns containing 0 to the left side of the data frame
        tb[,(ncol(tb)+1)] <- 0
      }
      allData = T #all data has been added to the data frame
    } else {  #non-zero values remain
      df <- data.frame(df[,colSums(df)>0], stringsAsFactors = F)  #drops all the 0-only containing columns
      tb[,(ncol(tb)+1)] <- df[,ncol(df)]  #takes all data in last column
      df[,ncol(df)] <- 0  #assigns all values in that column to 0 as they have been passed
      if(ncol(df)==1){  #checks to see if this is the last column, if it is, no further action is needed on passed zeros
        allData = T #all data has been added to the data frame
      } else {
        for(i in 1:(ncol(df)-1)){ #gets all the next non-zero values
          x <- row.names(tb)[tb[,ncol(tb)]==0 & rowSums(df)!=0] #the rows that have another non-zero value
          tb[x,ncol(tb)] <- df[x,(ncol(df)-i)]  #assigns those values to the data frame
          df[x,(ncol(df)-i)] <- 0 #assigns those values to 0 as they have been transfered
        } #repeats the checking process for all columns (it is a bit slower this way, but it is checking for zero things)
      }
    }
  }
  if(ncol(tb)<z){
    for(i in ncol(tb):(z-1)){
      tb[,(ncol(tb)+1)] <- 0
    }
  }
  if(!(is.na(names[1]))){
    colnames(tb) <- names
  }
  if(!(retVect)){
    return(tb)
  } else {
    return(as.vector(c(t(tb[1,]))))
  }
}
charVectorReverse <- function(df){
  #Takes a vector or data frame of characters and reverses the order forcing all "" values to the left
  # ex., takes c, b, a, "" and returns a, b, c, ""
  #
  #Args:
  # df: a data frame or vector
  #
  #Returns:
  # df: dataframe containing rows with reversed numbers
  # as.vector(c(t(tb[1,]))):  vector containing reversed numbers
  
  if(is.vector(df)){  #determins arg type and return type
    retVect = T
    df <- data.frame(t(df), stringsAsFactors = F)
    names <- NA
  } else {  #data frame return info
    retVect = F
    names <- c(colnames(df))
  }
  
  row.names(df) <- c(1:nrow(df))  #keeps row names
  tb <- data.frame(df[,ncol(df)]) #creates data frame
  tb <- tb[,-1] #empties data frame
  z <- ncol(df) #ncol of the passed df, for "" handling
  allData <- F  #boolean for deterimining if all data has been examined
  
  while(!(allData)){  #until all columns, or all data
    if((all(colSums(df!="")==0))){  #if there are no non-"" values left
      for(i in ncol(tb):(z-1)){ #adds columns containing "" to the left side of the data frame
        tb[,(ncol(tb)+1)] <- ""
      }
      allData = T #all data has been added to the data frame
    } else {  #non-"" values remain 
      df <- data.frame(df[,colSums(df!="")>0], stringsAsFactors = F)  #drops all the ""-only containing columns
      tb[,(ncol(tb)+1)] <- df[,ncol(df)]  #takes all data in last column
      df[,ncol(df)] <- ""  #assigns all values in that column to "" as they have been passed
      if(ncol(df)==1){  #checks to see if this is the last column, if it is, no further action is needed on passed ""s
        allData = T #all data has been added to the data frame
      } else {
        for(i in 1:(ncol(df)-1)){ #gets all the next non-"" values
          x <- row.names(tb)[tb[,ncol(tb)]=="" & rowSums(df!="")!=0] #the rows that have another non-"" value
          tb[x,ncol(tb)] <- df[x,(ncol(df)-i)]  #assigns those values to the data frame
          df[x,(ncol(df)-i)] <- "" #assigns those values to "" as they have been transfered
        } #repeats the checking process for all columns (it is a bit slower this way, but it is checking for "" things)
      }
    }
  }
  if(ncol(tb)<z){
    for(i in ncol(tb):(z-1)){
      tb[,(ncol(tb)+1)] <- ""
    }
  }
  if(!(is.na(names[1]))){
    colnames(tb) <- names
  }
  if(!(retVect)){
    return(tb)
  } else {
    return(as.vector(c(t(tb[1,]))))
  }
}
trimPrimeAndRealignSeq <- function(df, cols=NA, Length=30, Trim_Len=31, Len_Adjust=32, lmin=1, trimRounds=1, reorder=T, fromTrim=T, setRounds=F) {
  # Removes and tracks possible sequences added via the prime and realign mechanism
  #   "G","GC","GCA","GCAA", "GCAAA", "GCAAAA", "GCG", "GCGA", "GCGAA", "GCGAAA"
  #
  # Args:
  #	df: a data frame containg columns named "Trim_Sequence", "rounds", as well as the starting sequence
  # cols: leave as NA, number of columns in data frame
  # Length: column containing sequence length
  # Trim_Len: column containing Trim_Length
  # Len_Adjust: column containing length adjustments
  # lmin: minimum length where a sequence can still be trimmed further
  # trimRounds: variable to keep track of recursion instances; do not change via user input else it will trim only once
  # reorder:  change output so that the sequences and numbers are grouped together
  # fromTrim: if true, returns the additions starting from the maximally trimmed part; if false, returns the additions 
  #   as they are removed; ex., GCGGCA: true= GCG, GCA; false= GCA, GCG
  #
  # Returns:
  # trimPrimeAndRealignSeq(df, lmin=lmin, trimRounds=(trimRounds + 1)):
  #  Recursive trimming of added sequences (potentially) from prime and realign
  #	df:	data frame containing all passed data in addition to the sequences trimmed
  #   Sequence
  #   Trim_Sequence
  #   Added_Sequence
  #   variable region of Seq_Trim_R#
  #   variable region + Length
  #   variable region + Trim_Len
  #   variable region + Len_Adjust
  #   variable region + Velocity
  #   variable region + variable region of Num_Trim_R#
  
  if(is.na(cols)){
    cols = ncol(df)
  }
  
  #creates columns to determine if the seqence ends in a possible prime and realign sequence
  df[,(ncol(df)+1):(ncol(df)+10)] <- FALSE
  colnames(df)[(ncol(df)-9):ncol(df)] <- c("G","GC","GCA","GCAA", "GCAAA", "GCAAAA", "GCG", "GCGA", "GCGAA", "GCGAAA")#adds columns with possible sequence outcomes
  
  #prevents over-trimming
  df$Max_cut <- nchar(df$Trim_Sequence) - lmin#creates a column forcing max trim length
  
  #determines if a sequence ends in "G","GC","GCA","GCAA", "GCAAA", "GCAAAA", "GCG", "GCGA", "GCGAA", or "GCGAAA"
  #designed to trim the maximum possible length each round
  df$GC[((df$Max_cut>=2) & (endsWith(df$Trim_Sequence,"GC")))] <- TRUE
  df$GCA[((df$Max_cut>=3) & (endsWith(df$Trim_Sequence,"GCA")))] <- TRUE
  df$GCAA[((df$Max_cut>=4) & (endsWith(df$Trim_Sequence,"GCAA")))] <- TRUE
  df$GCAAA[((df$Max_cut>=5) & (endsWith(df$Trim_Sequence,"GCAAA")))] <- TRUE
  df$GCAAAA[((df$Max_cut>=6) & (endsWith(df$Trim_Sequence,"GCAAAA")))] <- TRUE
  df$GCG[((df$Max_cut>=3) & (endsWith(df$Trim_Sequence,"GCG")))] <- TRUE
  df$G[((df$Max_cut>=1) & (endsWith(df$Trim_Sequence,"G")) & !(df$GCG))] <- TRUE  #done like this to prevent 2 "TRUE" values (if a sequence ends in GCG forces GCG to be trimmed rather than G followed by GC the next round)
  df$GCGA[((df$Max_cut>=4) & (endsWith(df$Trim_Sequence,"GCGA")))] <- TRUE
  df$GCGAA[((df$Max_cut>=5) & (endsWith(df$Trim_Sequence,"GCGAA")))] <- TRUE
  df$GCGAAA[((df$Max_cut>=6) & (endsWith(df$Trim_Sequence,"GCGAAA")))] <- TRUE
  
  if(!setRounds){
    df$rounds<-  df$rounds + rowSums(df[(ncol(df)- 10):(ncol(df)-1)])#increases number of rounds trimmed if any of the nt sequences can be removed
  } else {
    df[,c((ncol(df)-10):(ncol(df)-1))][df$rounds < trimRounds,] <- FALSE
  }
  
  df[,(ncol(df)+1)] <- "" #creates a column to store the removed sequence 
  colnames(df)[ncol(df)] <- paste("Seq_Trim_R",trimRounds, sep="") #renames the column based on the round "Seq_Trim_R1"
  for(i in (ncol(df)-11):(ncol(df)-2)) {
    df[,ncol(df)][df[,i]==TRUE] <- colnames(df)[i] #writes the sequence removed in the Seq_Trim_R# column
  }
  df <- df[,-((ncol(df)- 11):(ncol(df) - 1))] #removes the prime and realign columns containing true/false, also removes Max_cut
  df[,(ncol(df)+1)] <- nchar(df[,ncol(df)]) #adds a column for the number of nt removed
  colnames(df)[ncol(df)] <- paste("Num_Trim_R",trimRounds, sep="")#renames the column based on the round
  df$Trim_Sequence <- substr(df$Trim_Sequence, 1, (nchar(df$Trim_Sequence)-df[,ncol(df)])) #removes the trimmed sequence to generate the new Trim_Sequence
  
  if(max(df$rounds)>=trimRounds) { #recursion if more rounds possible
    return(trimPrimeAndRealignSeq(df, cols=cols, Length=Length, Trim_Len=Trim_Len, Len_Adjust=Len_Adjust, lmin=lmin, trimRounds=(trimRounds + 1),reorder=reorder, fromTrim=fromTrim, setRounds=setRounds))
  } else { #ends if no nt were trimmed from any sequence in that round
    df <- df[,-((ncol(df)-1):(ncol(df)))]#dropping columns when nothing is trimmed that are created to store the sequence and number for that round
    df$Trim_Len <- nchar(df$Trim_Sequence)
    df$Len_Adjust <- df$Length - df$Trim_Len
    df$Added_Sequence <- substr(df$Sequence, start=(nchar(df$Trim_Sequence) + 1), stop=nchar(df$Sequence))
    df$Velocity <- df$Len_Adjust/df$rounds
    df$Velocity[is.na(df$Velocity)] <- 0
    if(reorder){
      df <- df[,c(1:(Length-1),seq((cols+1),ncol(df)-1,by=2), Length, Trim_Len, Len_Adjust, ncol(df), seq((cols+2),ncol(df)-1,by=2))]
      if(fromTrim){
        df[,((ncol(df)-trimRounds+2):(ncol(df)))] <- numVectorReverse(df[,((ncol(df)-trimRounds+2):(ncol(df)))])
        df[,((ncol(df)-2*(trimRounds-1)+1-4):(ncol(df)-(trimRounds-1)-4))] <- charVectorReverse(df[,((ncol(df)-2*(trimRounds-1)+1-4):(ncol(df)-(trimRounds-1)-4))])
      }
    }
    return(df)
  }
}

#you can set your working directory here if not using default
#setwd()
HK <- read.csv(file="C:/Users/HardCor/Desktop/IAV Data/seqNh-HK-NhL9_bigTable_allRounds-norc.csv", header=TRUE, sep="\t", stringsAsFactors=FALSE)	#HK
BR <- read.csv(file="C:/Users/HardCor/Desktop/IAV Data/seqNh_Brisbane-NhL9_bigTable_allRounds-norc.csv", header=TRUE, sep="\t", stringsAsFactors=FALSE)	#Brisbane
PR <- read.csv(file="C:/Users/HardCor/Desktop/IAV Data/seqNh-PR8-NhL9_bigTable_allRounds-norc.csv", header=TRUE, sep="\t", stringsAsFactors=FALSE)	#PR8
WSN <- read.csv(file="C:/Users/HardCor/Desktop/IAV Data/seqNh-wisconsin-NhL9_bigTable_allRounds-norc.csv", header=TRUE, sep="\t", stringsAsFactors=FALSE)	#WSN
WSN <- WSN[,-18]	#trims duplicate

cutOff_pos_relative <- 101		#removes all values above this number and below the negative value of this number; 101 means no cut-off
cutOff_transcript_support_level <- 8	#removes values above this match type; 8 means no cut-off
cutOff_count_artifact <- 10		#removes all values below this number; 0 means no cut-off

HK <- trimData(HK, prel=cutOff_pos_relative, tsl=cutOff_transcript_support_level, art=cutOff_count_artifact)
BR <- trimData(BR, prel=cutOff_pos_relative, tsl=cutOff_transcript_support_level, art=cutOff_count_artifact)
PR <- trimData(PR, prel=cutOff_pos_relative, tsl=cutOff_transcript_support_level, art=cutOff_count_artifact)
WSN <- trimData(WSN, prel=cutOff_pos_relative, tsl=cutOff_transcript_support_level, art=cutOff_count_artifact)

#length limits
lmin <- 9 #cuts all mRNA below this length
lmax <- 16 #cuts all mRNA above this length
ltrim <- function(df, lmin=8, lmax=17){
  #cuts length around the min and max lengths
  
  df <- df[nchar(df$Seq_Bowtie) >=lmin,]
  df <- df[nchar(df$Seq_Bowtie) <=lmax,]
  return(df)
}
HK <- ltrim(HK, lmin=lmin, lmax=lmax)
PR <- ltrim(PR, lmin=lmin, lmax=lmax)
BR <- ltrim(BR, lmin=lmin, lmax=lmax)
WSN <- ltrim(WSN, lmin=lmin, lmax=lmax)

#extract PAR additions
lmin=1
HK <- trimPrimeAndRealignSeq(HK, cols=NA, Length=grep("Length",colnames(HK)), Trim_Len=grep("Trim_Len",colnames(HK)), Len_Adjust=grep("Len_Adjust",colnames(HK)), lmin=lmin, reorder=T, fromTrim=T, setRounds=T)
BR <- trimPrimeAndRealignSeq(BR, cols=NA, Length=grep("Length",colnames(BR)), Trim_Len=grep("Trim_Len",colnames(BR)), Len_Adjust=grep("Len_Adjust",colnames(BR)), lmin=lmin, reorder=T, fromTrim=T, setRounds=T)
PR <- trimPrimeAndRealignSeq(PR, cols=NA, Length=grep("Length",colnames(PR)), Trim_Len=grep("Trim_Len",colnames(PR)), Len_Adjust=grep("Len_Adjust",colnames(PR)), lmin=lmin, reorder=T, fromTrim=T, setRounds=T)
WSN <- trimPrimeAndRealignSeq(WSN, cols=NA, Length=grep("Length",colnames(WSN)), Trim_Len=grep("Trim_Len",colnames(WSN)), Len_Adjust=grep("Len_Adjust",colnames(WSN)), lmin=lmin, reorder=T, fromTrim=T, setRounds=T)

#change saving directory here if you want to
#setwd("")
write.csv(file="HK_All_Match.csv", HK, row.names=FALSE)
write.csv(file="BRI_All_Match.csv", BR, row.names=FALSE)
write.csv(file="PR8_All_Match.csv", PR, row.names=FALSE)
write.csv(file="WSN_All_Match.csv", WSN, row.names=FALSE)

##### Luciferase #####
tslClean <- function(table){ 
  # Cleans up the transcript_support_level (tsl) of a given table to numeric values only
  #
  # Args:
  #	table:	30 column table containing the data
  #
  # Returns:
  #	table:	the table as specified in args with the tsl adjusted to numeric values
  
  #assign all NA or "" values to 6
  table$transcript_support_level[is.na(table$transcript_support_level)] <- 6
  table$transcript_support_level[table$transcript_support_level==""] <- 6
  table$transcript_support_level[table$transcript_support_level==""] <- 6
  table$transcript_support_level[startsWith(table$transcript_support_level,"N")] <- 6
  
  #cuts down annotation to a single number
  table$transcript_support_level[startsWith(table$transcript_support_level,"5")] <- 5	
  table$transcript_support_level[startsWith(table$transcript_support_level,"4")] <- 4
  table$transcript_support_level[startsWith(table$transcript_support_level,"3")] <- 3
  table$transcript_support_level[startsWith(table$transcript_support_level,"2")] <- 2
  table$transcript_support_level[startsWith(table$transcript_support_level,"1")] <- 1
  return(table)
}
keepR <- function(table, col=1, val=0, val2=NA, t="==", chr=FALSE){
  # Keeps Rows (keepR) within a specified range
  #
  # Args:
  #	table:	a data.frame to perform operations on, this data.frame will have only rows
  #			except that meet the specified parameters
  #	col: 	column number which the range of values kept will be applied to
  #	val: 	value to keep, or rnage limiter number 1
  #	val2:	range limiter number 2; will be set as val by default if needed
  #	t:  type of operation to be performed to keep rows where this value is TRUE
  # chr: spefies if the values kept are to be relted to string length
  #
  # Returns:
  #	table:	the table passed following removal of rows where "t" argued opperation is false
  if(!chr) {
    if(t=="=="){
      table <- table[table[,col]==val,]
    } else if(t=="!=") {
      table <- table[table[,col]!=val,]
    } else if(t=="<") {
      table <- table[table[,col]<val,]
    } else if(t=="<=") {
      table <- table[table[,col]<=val,]
    } else if(t==">") {
      table <- table[table[,col]>val,]
    } else if(t==">=") {
      table <- table[table[,col]>=val,]
    } else if(t=="na") {
      table <- table[is.na(table[,col]),]
    } else if(t=="!na") {
      table <- table[!(is.na(table[,col])),]
    } else if(t=="<>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]<val) & (table[,col]>val2),]
    } else if(t=="><") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]>val) | (table[,col]<val2),]
    } else if(t=="<=>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]<=val) & (table[,col]>=val2),]
    } else if(t==">=<") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(table[,col]>=val) | (table[,col]<=val2),]
    } else {
      print("invalid type (t)")
    }
  } else if(chr) {
    if(t=="=="){
      table <- table[nchar(table[,col])==val,]
    } else if(t=="!=") {
      table <- table[nchar(table[,col])!=val,]
    } else if(t=="<") {
      table <- table[nchar(table[,col])<val,]
    } else if(t=="<=") {
      table <- table[nchar(table[,col])<=val,]
    } else if(t==">") {
      table <- table[nchar(table[,col])>val,]
    } else if(t==">=") {
      table <- table[nchar(table[,col])>=val,]
    } else if(t=="<>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])<val) & (nchar(table[,col])>val2),]
    } else if(t=="><") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])>val) | (nchar(table[,col])<val2),]
    } else if(t=="<=>") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])<=val) & (nchar(table[,col])>=val2),]
    } else if(t==">=<") {
      if(is.na(val2)){val2 <- val * -1}
      table <- table[(nchar(table[,col])>=val) | (nchar(table[,col])<=val2),]
    } else {
      print("invalid type (t)")
    }
  }
  return(table)
}
trimDataLuc <- function(df, prel=10, tsl=4, art=10){
  # trims data tables to a given limitations for further analysis, specifically for the data set used for prediction
  #
  # Args:
  #	df:	a data.frame to perform operations on, removing data above the cut-offs
  #	prel: 	cut-off value based on distance to the transcription start site
  #	tsl: 	transcription support level cut-off (requires cleaned up TSL data)
  #	art:	number of reads for a given subunit as to not be considered background
  #
  # Returns:
  #	df:	the table passed following removal of rows outside of cut-offs
  
  df <- tslClean(df)	#converts transcript_support_level data tonumeric
  df <- keepR(df, col=grep('transcript_support_level', colnames(df)), val= tsl, t="<=")	#keeps only rows with a transcript_support_level less than the tsl arg.
  df <- keepR(df, col=6, val=prel, t="<=>")	#keeps only rows with within +/- prel arg. of a transcription start site
  df <- df[(df[,10]>=art | df[,11]>=art | df[,12]>=art | df[,13]>=art),]	#removes rows where no single read count is greater than the art arg.
  df$Round <- df$Round - 1  #decrease rounds to reflect number of trims
  df[,10:13][df[,10:13] < art] <- 0	#sets any values below the art arg. to 0 in rows that have not been removed
  colnames(df)[9] <- "rounds"
  colnames(df)[3] <- "Sequence"
  df$Trim_Sequence <- df$Sequence
  df$Added_Sequence <- ""
  df$Length <- nchar(df$Sequence)
  df$Trim_Len <- nchar(df$Trim_Sequence)
  df$Len_Adjust <- 0
  subs <- 10:13
  others <- (1:ncol(df))[!((1:ncol(df)) %in% subs)]
  df <- df[,c(subs, others)]
  
  return(df)
}
numVectorReverse <- function(df){
  #Takes a vector or data frame of numbers and reverses the order forcing all zero values to the left
  # ex., takes 3, 2, 1, 0 and returns 1, 2, 3, 0
  #
  #Args:
  # df: a data frame or vector
  #
  #Returns:
  # df: dataframe containing rows with reversed numbers
  # as.vector(c(t(tb[1,]))):  vector containing reversed numbers
  
  if(is.vector(df)){  #determins arg type and return type
    retVect = T
    df <- data.frame(t(df), stringsAsFactors = F)
    names <- NA
  } else {  #data frame return info
    retVect = F
    names <- c(colnames(df))
  }
  
  row.names(df) <- c(1:nrow(df))  #keeps row names
  tb <- data.frame(df[,ncol(df)]) #creates data frame
  tb <- tb[,-1] #empties data frame
  z <- ncol(df) #ncol of the passed df, for 0 handling
  allData <- F  #boolean for deterimining if all data has been examined
  
  while(!(allData)){  #until all columns, or all data
    if(length(data.frame(df[,colSums(df)>0], stringsAsFactors = F)[1,])==0){  #if there are no non-zero values left
      for(i in ncol(tb):(z-1)){ #adds columns containing 0 to the left side of the data frame
        tb[,(ncol(tb)+1)] <- 0
      }
      allData = T #all data has been added to the data frame
    } else {  #non-zero values remain
      df <- data.frame(df[,colSums(df)>0], stringsAsFactors = F)  #drops all the 0-only containing columns
      tb[,(ncol(tb)+1)] <- df[,ncol(df)]  #takes all data in last column
      df[,ncol(df)] <- 0  #assigns all values in that column to 0 as they have been passed
      if(ncol(df)==1){  #checks to see if this is the last column, if it is, no further action is needed on passed zeros
        allData = T #all data has been added to the data frame
      } else {
        for(i in 1:(ncol(df)-1)){ #gets all the next non-zero values
          x <- row.names(tb)[tb[,ncol(tb)]==0 & rowSums(df)!=0] #the rows that have another non-zero value
          tb[x,ncol(tb)] <- df[x,(ncol(df)-i)]  #assigns those values to the data frame
          df[x,(ncol(df)-i)] <- 0 #assigns those values to 0 as they have been transfered
        } #repeats the checking process for all columns (it is a bit slower this way, but it is checking for zero things)
      }
    }
  }
  if(ncol(tb)<z){
    for(i in ncol(tb):(z-1)){
      tb[,(ncol(tb)+1)] <- 0
    }
  }
  if(!(is.na(names[1]))){
    colnames(tb) <- names
  }
  if(!(retVect)){
    return(tb)
  } else {
    return(as.vector(c(t(tb[1,]))))
  }
}
charVectorReverse <- function(df){
  #Takes a vector or data frame of characters and reverses the order forcing all "" values to the left
  # ex., takes c, b, a, "" and returns a, b, c, ""
  #
  #Args:
  # df: a data frame or vector
  #
  #Returns:
  # df: dataframe containing rows with reversed numbers
  # as.vector(c(t(tb[1,]))):  vector containing reversed numbers
  
  if(is.vector(df)){  #determins arg type and return type
    retVect = T
    df <- data.frame(t(df), stringsAsFactors = F)
    names <- NA
  } else {  #data frame return info
    retVect = F
    names <- c(colnames(df))
  }
  
  row.names(df) <- c(1:nrow(df))  #keeps row names
  tb <- data.frame(df[,ncol(df)]) #creates data frame
  tb <- tb[,-1] #empties data frame
  z <- ncol(df) #ncol of the passed df, for "" handling
  allData <- F  #boolean for deterimining if all data has been examined
  
  while(!(allData)){  #until all columns, or all data
    if((all(colSums(df!="")==0))){  #if there are no non-"" values left
      for(i in ncol(tb):(z-1)){ #adds columns containing "" to the left side of the data frame
        tb[,(ncol(tb)+1)] <- ""
      }
      allData = T #all data has been added to the data frame
    } else {  #non-"" values remain 
      df <- data.frame(df[,colSums(df!="")>0], stringsAsFactors = F)  #drops all the ""-only containing columns
      tb[,(ncol(tb)+1)] <- df[,ncol(df)]  #takes all data in last column
      df[,ncol(df)] <- ""  #assigns all values in that column to "" as they have been passed
      if(ncol(df)==1){  #checks to see if this is the last column, if it is, no further action is needed on passed ""s
        allData = T #all data has been added to the data frame
      } else {
        for(i in 1:(ncol(df)-1)){ #gets all the next non-"" values
          x <- row.names(tb)[tb[,ncol(tb)]=="" & rowSums(df!="")!=0] #the rows that have another non-"" value
          tb[x,ncol(tb)] <- df[x,(ncol(df)-i)]  #assigns those values to the data frame
          df[x,(ncol(df)-i)] <- "" #assigns those values to "" as they have been transfered
        } #repeats the checking process for all columns (it is a bit slower this way, but it is checking for "" things)
      }
    }
  }
  if(ncol(tb)<z){
    for(i in ncol(tb):(z-1)){
      tb[,(ncol(tb)+1)] <- ""
    }
  }
  if(!(is.na(names[1]))){
    colnames(tb) <- names
  }
  if(!(retVect)){
    return(tb)
  } else {
    return(as.vector(c(t(tb[1,]))))
  }
}
trimPrimeAndRealignSeq <- function(df, cols=NA, Length=30, Trim_Len=31, Len_Adjust=32, lmin=1, trimRounds=1, reorder=T, fromTrim=T, setRounds=F) {
  # Removes and tracks possible sequences added via the prime and realign mechanism
  #   "G","GC","GCA","GCAA", "GCAAA", "GCAAAA", "GCG", "GCGA", "GCGAA", "GCGAAA"
  #
  # Args:
  #	df: a data frame containg columns named "Trim_Sequence", "rounds", as well as the starting sequence
  # cols: leave as NA, number of columns in data frame
  # Length: column containing sequence length
  # Trim_Len: column containing Trim_Length
  # Len_Adjust: column containing length adjustments
  # lmin: minimum length where a sequence can still be trimmed further
  # trimRounds: variable to keep track of recursion instances; do not change via user input else it will trim only once
  # reorder:  change output so that the sequences and numbers are grouped together
  # fromTrim: if true, returns the additions starting from the maximally trimmed part; if false, returns the additions 
  #   as they are removed; ex., GCGGCA: true= GCG, GCA; false= GCA, GCG
  #
  # Returns:
  # trimPrimeAndRealignSeq(df, lmin=lmin, trimRounds=(trimRounds + 1)):
  #  Recursive trimming of added sequences (potentially) from prime and realign
  #	df:	data frame containing all passed data in addition to the sequences trimmed
  #   Sequence
  #   Trim_Sequence
  #   Added_Sequence
  #   variable region of Seq_Trim_R#
  #   variable region + Length
  #   variable region + Trim_Len
  #   variable region + Len_Adjust
  #   variable region + Velocity
  #   variable region + variable region of Num_Trim_R#
  
  if(is.na(cols)){
    cols = ncol(df)
  }
  
  #creates columns to determine if the seqence ends in a possible prime and realign sequence
  df[,(ncol(df)+1):(ncol(df)+10)] <- FALSE
  colnames(df)[(ncol(df)-9):ncol(df)] <- c("G","GC","GCA","GCAA", "GCAAA", "GCAAAA", "GCG", "GCGA", "GCGAA", "GCGAAA")#adds columns with possible sequence outcomes
  
  #prevents over-trimming
  df$Max_cut <- nchar(df$Trim_Sequence) - lmin#creates a column forcing max trim length
  
  #determines if a sequence ends in "G","GC","GCA","GCAA", "GCAAA", "GCAAAA", "GCG", "GCGA", "GCGAA", or "GCGAAA"
  #designed to trim the maximum possible length each round
  df$GC[((df$Max_cut>=2) & (endsWith(df$Trim_Sequence,"GC")))] <- TRUE
  df$GCA[((df$Max_cut>=3) & (endsWith(df$Trim_Sequence,"GCA")))] <- TRUE
  df$GCAA[((df$Max_cut>=4) & (endsWith(df$Trim_Sequence,"GCAA")))] <- TRUE
  df$GCAAA[((df$Max_cut>=5) & (endsWith(df$Trim_Sequence,"GCAAA")))] <- TRUE
  df$GCAAAA[((df$Max_cut>=6) & (endsWith(df$Trim_Sequence,"GCAAAA")))] <- TRUE
  df$GCG[((df$Max_cut>=3) & (endsWith(df$Trim_Sequence,"GCG")))] <- TRUE
  df$G[((df$Max_cut>=1) & (endsWith(df$Trim_Sequence,"G")) & !(df$GCG))] <- TRUE  #done like this to prevent 2 "TRUE" values (if a sequence ends in GCG forces GCG to be trimmed rather than G followed by GC the next round)
  df$GCGA[((df$Max_cut>=4) & (endsWith(df$Trim_Sequence,"GCGA")))] <- TRUE
  df$GCGAA[((df$Max_cut>=5) & (endsWith(df$Trim_Sequence,"GCGAA")))] <- TRUE
  df$GCGAAA[((df$Max_cut>=6) & (endsWith(df$Trim_Sequence,"GCGAAA")))] <- TRUE
  
  if(!setRounds){
    df$rounds<-  df$rounds + rowSums(df[(ncol(df)- 10):(ncol(df)-1)])#increases number of rounds trimmed if any of the nt sequences can be removed
  } else {
    df[,c((ncol(df)-10):(ncol(df)-1))][df$rounds < trimRounds,] <- FALSE
  }
  
  df[,(ncol(df)+1)] <- "" #creates a column to store the removed sequence 
  colnames(df)[ncol(df)] <- paste("Seq_Trim_R",trimRounds, sep="") #renames the column based on the round "Seq_Trim_R1"
  for(i in (ncol(df)-11):(ncol(df)-2)) {
    df[,ncol(df)][df[,i]==TRUE] <- colnames(df)[i] #writes the sequence removed in the Seq_Trim_R# column
  }
  df <- df[,-((ncol(df)- 11):(ncol(df) - 1))] #removes the prime and realign columns containing true/false, also removes Max_cut
  df[,(ncol(df)+1)] <- nchar(df[,ncol(df)]) #adds a column for the number of nt removed
  colnames(df)[ncol(df)] <- paste("Num_Trim_R",trimRounds, sep="")#renames the column based on the round
  df$Trim_Sequence <- substr(df$Trim_Sequence, 1, (nchar(df$Trim_Sequence)-df[,ncol(df)])) #removes the trimmed sequence to generate the new Trim_Sequence
  
  if(max(df$rounds)>=trimRounds) { #recursion if more rounds possible
    return(trimPrimeAndRealignSeq(df, cols=cols, Length=Length, Trim_Len=Trim_Len, Len_Adjust=Len_Adjust, lmin=lmin, trimRounds=(trimRounds + 1),reorder=reorder, fromTrim=fromTrim, setRounds=setRounds))
  } else { #ends if no nt were trimmed from any sequence in that round
    df <- df[,-((ncol(df)-1):(ncol(df)))]#dropping columns when nothing is trimmed that are created to store the sequence and number for that round
    df$Trim_Len <- nchar(df$Trim_Sequence)
    df$Len_Adjust <- df$Length - df$Trim_Len
    df$Added_Sequence <- substr(df$Sequence, start=(nchar(df$Trim_Sequence) + 1), stop=nchar(df$Sequence))
    df$Velocity <- df$Len_Adjust/df$rounds
    df$Velocity[is.na(df$Velocity)] <- 0
    if(reorder){
      df <- df[,c(1:(Length-1),seq((cols+1),ncol(df)-1,by=2), Length, Trim_Len, Len_Adjust, ncol(df), seq((cols+2),ncol(df)-1,by=2))]
      if(fromTrim){
        df[,((ncol(df)-trimRounds+2):(ncol(df)))] <- numVectorReverse(df[,((ncol(df)-trimRounds+2):(ncol(df)))])
        df[,((ncol(df)-2*(trimRounds-1)+1-4):(ncol(df)-(trimRounds-1)-4))] <- charVectorReverse(df[,((ncol(df)-2*(trimRounds-1)+1-4):(ncol(df)-(trimRounds-1)-4))])
      }
    }
    return(df)
  }
}

#setwd("")
luc <- read.csv("seqNh-luc-NhL9_bigTable_allRounds-norc.csv", header=TRUE, sep="\t", stringsAsFactors=FALSE) #prediction set

cutOff_pos_relative <- 101		#removes all values above this number and below the negative value of this number
cutOff_transcript_support_level <- 8	#removes values above this match type
cutOff_count_artifact <- 10		#removes all values below this number

luc <- trimDataLuc(luc, prel=cutOff_pos_relative, tsl=cutOff_transcript_support_level, art=cutOff_count_artifact)

lmin <- 9 #cuts all mRNA below this length
lmax <- 16 #cuts all mRNA above this length

ltrim <- function(df, lmin=8, lmax=17){
  #cuts length around the min and max lengths
  
  df <- df[nchar(df$Seq_Bowtie) >=lmin,]
  df <- df[nchar(df$Seq_Bowtie) <=lmax,]
  return(df)
}

luc <- ltrim(luc, lmin=lmin, lmax=lmax)

lmin=1

luc <- trimPrimeAndRealignSeq(luc, cols=NA, Length=grep("Length",colnames(luc)), Trim_Len=grep("Trim_Len",colnames(luc)), Len_Adjust=grep("Len_Adjust",colnames(luc)), lmin=lmin, reorder=T, fromTrim=T, setRounds=T)

write.csv(file="LUC_All_Match.csv", luc, row.names=FALSE)
