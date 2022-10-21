#!/usr/bin/env Rscript
library('Xmisc')

#This was created by Step 4 and has the data needed about each candidate.
load('Step4.RData')

#Default parameters
maxReadMin <- 10
minStartRead <- 1
maxPerEndRead <- .15
distanceTEMax <- 2500

#User argument parsing code
parser <- ArgumentParser$new()

parser$add_argument('-r',type='character',help='Minimum number of uniquely mapping reads supporting TE')
parser$add_argument('-s',type='character',help='Minimum number of start reads in a single file')
parser$add_argument('-e',type='character',help='Maximum percent of files that cna have a read going from gene to TE which suggests an exonization rather than a new promoter event')
parser$add_argument('-d',type='character',help='The minimum total number of samples that a candidate can be present in to be included')
argsparse <- parser$get_args()

#Argument parse to change defaults
if (!identical(argsparse$r, character(0))){
  maxReadMin <- as.numeric(argsparse$r)
}

if (!identical(argsparse$s, character(0))){
  minStartRead <- as.numeric(argsparse$s)
}

if (!identical(argsparse$e, character(0))){
  maxPerEndRead <- as.numeric(argsparse$e)
}

if (!identical(argsparse$d, character(0))){
  distanceTEMax <- as.numeric(argsparse$d)
}

print("Arguments")
print(paste0("Read Support Minimum in 1 File: ", maxReadMin))
print(paste0("PE Chimeric Read Support: ", minStartRead))
print(paste0("Exonization Percentage Maximum: ", maxPerEndRead))
print(paste0("Distance Upstream TE must be from start of reference transcript: ", distanceTEMax))

readstats <- read.delim('filter_read_stats.txt', sep = "\t", stringsAsFactors = FALSE, header = FALSE)
colnames(readstats) <- c('uniqid_and_file_name', 'read', 'startread', 'endread', 'filetype')

readstats$uniqidfile <- gsub("--","_",readstats$uniqid_and_file_name)
readstats$uniqidfile <- basename(readstats$uniqidfile)
readstats$uniqidfile <- gsub(".stats","",readstats$uniqidfile)

filter_combined_table_final_stat <- merge(filter_combined_table,readstats, by = c('uniqidfile'), all.x = TRUE, all.y = FALSE)
filter_combined_table_final_stat$uniqidfile <- NULL
filter_combined_table_final_stat$uniqid_and_file_name <- NULL

stataggregate <- aggregate(filter_combined_table_final_stat[,c("read",'startread', 'endread')], by=list(filter_combined_table_final_stat[,c("uniqid")]), FUN=sum, na.rm=TRUE)
names(stataggregate)[names(stataggregate) == 'Group.1'] <- 'uniqid'

res2stat<- merge(resexonlength, stataggregate, by=c('uniqid'), all.x = TRUE, all.y = FALSE)

readAggregateStats <- function(testtable){
  totalsamples <- nrow(testtable)
  totalsamplesstart <- sum(testtable$startread > 0)
  totalsamplesend <- sum(testtable$endread > 0)
  sesamples <- sum(testtable$filetype == 'se')
  startreadper <- totalsamplesstart/totalsamples
  endreadper <- totalsamplesend/totalsamples
  seper <- sesamples/totalsamples
  return(c(startreadper, endreadper, seper))
}

resultsinitial <- apply(res2stat[,c('gene2','uniqid')], 1, function(x) readAggregateStats(filter_combined_table_final_stat[filter_combined_table_final_stat$uniqid == x[2], ]))

resultsinitial  <- data.frame(t(resultsinitial))
colnames(resultsinitial) <- c("startreadper", "endreadper", "seper")

res2stat<- cbind(res2stat, resultsinitial)
res2stat$distanceTE <- apply(res2stat[,c('strand','transcriptstart2','transcriptend2','startTE','endTE')],1,function(x) if (x[1] == '+'){return(as.numeric(x[2]) - as.numeric(x[5]))}else{return(as.numeric(x[4])-as.numeric(x[3]))})

res2statfil <- res2stat[res2stat$endreadper <= maxPerEndRead, ]
res2statfil <- res2statfil[res2statfil$startread >= minStartRead,]
res2statfil <- res2statfil[res2statfil$distanceTE > distanceTEMax | res2statfil$distanceTE < 0,]

res2statfil$maxReadsSample <- apply(res2statfil[,c('uniqid','read')],1,function(x) max(filter_combined_table_final_stat[filter_combined_table_final_stat$uniqid == x[1],c('read')]))

res2statfil <- res2statfil[res2statfil$maxReadsSample >= maxReadMin,]

chooseTopIsoform <- function(testtable){
  transcriptids <- c()
  transcoords <- c()
  covtransvec <- c()
  transcriptstrucs <- c()
  intron1start <- c()
  intron1end <- c()
  intronanno <- c()
  for (i in 1:nrow(testtable)){
    row1 <- testtable[i,]
    structurestring <- row1$transcoord
    strandexample <- strsplit(structurestring, ",")[[1]][1]
    elementsvec <- as.numeric(tail(strsplit(structurestring, ",")[[1]],-2))
    elementsvec <- sort(elementsvec)
    
    #Get rid of the arbitrary 5' start and 3' end site of the transcript
    elementsvec <- elementsvec[-1]
    elementsvec <- elementsvec[-length(elementsvec)]
    transcriptstrucs <- c(transcriptstrucs,paste(elementsvec, collapse = ","))
    transcoords <- c(transcoords, row1$transcoord)
    covtransvec <-c(covtransvec, row1$covtrans)
    transcriptids <- c(transcriptids, row1$id2)
    intronanno <- c(intronanno, row1$intronanno)
    if (row1$strand == '+'){
      intron1start <- c(intron1start, elementsvec[1] + 1)
      intron1end <- c(intron1end, elementsvec[2] - 1)
    } else {
      intron1start <- c(intron1start, elementsvec[length(elementsvec)-1] + 1)
      intron1end <- c(intron1end, elementsvec[length(elementsvec)] - 1)
    }
  }
  
  df.temp <- data.frame(transcriptstrucs = transcriptstrucs, transcoords = transcoords, covtransvec = covtransvec, transcriptids = transcriptids, intron1start = intron1start, intron1end = intron1end, intronanno = intronanno, stringsAsFactors = FALSE)
  df.temp <- df.temp[order(-covtransvec),]
  struccount <- as.data.frame(table(transcriptstrucs))
  struccount$transcriptstrucs <- as.character(struccount$transcriptstrucs)
  
  df.temp <- df.temp[df.temp$transcriptstrucs %in% struccount[struccount$Freq == max(struccount$Freq),c('transcriptstrucs')],]
  return(c(df.temp$transcoords[1], df.temp$transcriptids[1], df.temp$intron1start[1], df.temp$intron1end[1],df.temp$intronanno[1]))
}

resultsinitial <- apply(res2statfil[,c('gene2','uniqid')], 1, function(x) chooseTopIsoform (filter_combined_table_final_stat[filter_combined_table_final_stat$uniqid == x[2], ]))

resultsinitial  <- data.frame(t(resultsinitial))
colnames(resultsinitial) <- c("transcoords", "transcriptid", "intron1start", "intron1end", "intronanno")

res2statfil<- cbind(res2statfil, resultsinitial)

i <- sapply(res2statfil, is.factor)
res2statfil[i] <- lapply(res2statfil[i], as.character)

#GFF3 Creation
getGff <- function(uniqid, structurestring){
  longvector <- c()
  chrexample <- strsplit(structurestring, ",")[[1]][2]
  strandexample <- strsplit(structurestring, ",")[[1]][1]
  elementsvec <- as.numeric(tail(strsplit(structurestring, ",")[[1]],-2))
  elementsvec <- sort(elementsvec)
  numexons <- length(elementsvec)/2
  longvector <- c(longvector, chrexample)
  longvector  <- c(longvector, "curated")
  longvector  <- c(longvector, "mRNA")
  longvector  <- c(longvector, min(elementsvec))
  longvector  <- c(longvector, max(elementsvec))
  longvector  <- c(longvector, ".")
  longvector  <- c(longvector, strandexample)
  longvector  <- c(longvector, ".")
  longvector  <- c(longvector, paste0('ID=',uniqid))
  for (j in 1:numexons){
    longvector <- c(longvector, chrexample)
    longvector <- c(longvector, "curated")
    longvector <- c(longvector, "exon")
    longvector <- c(longvector, elementsvec[j*2-1])
    longvector <- c(longvector, elementsvec[j*2])
    longvector <- c(longvector, ".")
    longvector <- c(longvector, strandexample)
    longvector <- c(longvector, ".")
    longvector <- c(longvector, paste0('Parent=',uniqid))
  }
  return(longvector)
}

gffexamples <- apply(res2statfil[,c("uniqid", "transcoords")], 1, function(x) getGff(x[1],x[2]))
gffexamples <- unlist(gffexamples)
gffexamples.m <- matrix(gffexamples, nrow=9)
gffexamples.m <- as.data.frame(t(gffexamples.m))


write.table(gffexamples.m, "candidate_transcripts.gff3", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

write.table(res2statfil, "read_filtered_candidates.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
save.image('Step6.RData')

