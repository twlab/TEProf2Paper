#!/usr/bin/env Rscript

library('Xmisc')

#Get the directory of the rnapipeline scripts
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

#Default parameters
parseTreatment <- ''
exon1LengthMax <- 2588
exonSkipMax <- 2
totalMin <- 1
totalTreatmentMin <- 0
treatmentexclusive <- 'no'
keepNone <- 'no'
filterTE <- 'yes'
argumentFile <- paste0(script.basename,'/arguments.txt')

#User argument parsing code
parser <- ArgumentParser$new()

parser$add_argument('-e',type='character',help='A character of parsing the treatment samples')
parser$add_argument('-l',type='character',help='A number of the exon 1 length max')
parser$add_argument('-s',type='character',help='A number representing the maximum numbers of exons that can be skipped by transcript')
parser$add_argument('-n',type='character',help='The minimum total number of samples that a candidate can be present in to be included')
parser$add_argument('-t',type='character',help='The minimum number of treatment samples the candidate has to be present in')
parser$add_argument('-x',type='character',help='If yes is entered, then only candidates exclusive to the treatment samples will be in the output')
parser$add_argument('-k',type='character',help='If yes is entered, then candidates that do not end up splicing into a gene will be included')
parser$add_argument('-f',type='character',help='If yes is entered, only the following TE classes will be used (LINE, SINE, SINE?, LTR, LTR?, DNA, DNA?, and Retroposon in Repeatmasker)')
parser$add_argument('-a',type='character',help='The arguments.txt file that is being used for this run.')
argsparse <- parser$get_args()

#Argument parse to change defaults
if (!identical(argsparse$e, character(0))){
  parseTreatment <- argsparse$e
}

if (!identical(argsparse$l, character(0))){
  exon1LengthMax <- as.numeric(argsparse$l)
}

if (!identical(argsparse$s, character(0))){
  exonSkipMax <- as.numeric(argsparse$s)
}

if (!identical(argsparse$n, character(0))){
  totalMin <- as.numeric(argsparse$n)
}

if (!identical(argsparse$t, character(0))){
  totalTreatmentMin <- as.numeric(argsparse$t)
}

if (!identical(argsparse$x, character(0))){
  treatmentexclusive <- argsparse$x
}

if (!identical(argsparse$k, character(0))){
  keepNone <- argsparse$k
}

if (!identical(argsparse$f, character(0))){
  filterTE <- argsparse$f
}

if (!identical(argsparse$a, character(0))){
  argumentFile <- argsparse$a
}

print("Arguments")
print(paste0("Treatment Parse: ", parseTreatment))
print(paste0("Exon 1 Length Max: ", exon1LengthMax))
print(paste0("Exon Skip Max: ", exonSkipMax))
print(paste0("Minimum Total Samples: ", totalMin))
print(paste0("Minimum Treatment Samples: ", totalTreatmentMin))
print(paste0("Treatment Exclusive: ", treatmentexclusive))
print(paste0("Keep None Splice Target: ", keepNone))
print(paste0("Filter for only TE classes: ", filterTE))
print(paste0("Argument File Location: ", argumentFile))

#Get the arugments from the arguments.txt file that are needed

argTable <- read.delim(argumentFile, header=FALSE, stringsAsFactors = FALSE)
repeatAnnotationFile <- argTable[argTable[,1] == 'rmskannotationfile',c(2)]


print("")
print("")
print("Files Detected:")
print("")
files_annotated <- list.files(pattern = "_c$")
print(files_annotated)

firsttime = TRUE
i = 0
for (file_name in files_annotated){

  if (file.size(file_name) == 0){
    next
  }
  filtered_table = read.delim(file_name, sep = "\t", stringsAsFactors = FALSE, header = FALSE)

  columnlabels <- c("transcriptname", "type1", "gene1", "chr1", "start1", "end1", "exonintron1", "number1" ,"startcodon1" , "transcriptstart1", "transcriptend1", "elements1", "id1", "splicing", "type2", "gene2", "chr2", "start2", "end2", "exonintron2", "number2" ,"startcodon2", "transcriptstart2", "transcriptend2", "elements2", "id2", "chromtrans", "starttrans", "endtrans", "transcoord", "intronanno", "chrTE", "startTE", "endTE", "subfamTE", "numTE", "strandTE", "strand","covtrans", "covexon", "tpm", "maxtranscov", "maxtpm", "totaltpm")

  colnames(filtered_table) <- columnlabels
  
  filtered_table$file_label <- rep(strsplit(file_name, split = ".gtf")[[1]][1], length(filtered_table$strand))

  if (grepl(parseTreatment,file_name)){
    filtered_table$tumorcount <- rep(1,length(filtered_table$type1))
    filtered_table$normalcount <- rep(0,length(filtered_table$type1))
  } else {
    filtered_table$tumorcount <- rep(0,length(filtered_table$type1))
    filtered_table$normalcount <- rep(1,length(filtered_table$type1))
  }

  if (firsttime == TRUE){
    combined_table = filtered_table
    firsttime = FALSE
  } else {
    combined_table <- rbind(combined_table,filtered_table,stringsAsFactors = FALSE)
  }

  #print(i)
  i = i + 1
}

filter_combined_table = combined_table

filter_combined_table <- filter_combined_table[filter_combined_table$exonintron1 != "exon",]
# If the transcript starts from an exon it is removed


filter_combined_table$permaxcov <- filter_combined_table$covtrans/filter_combined_table$maxtranscov
filter_combined_table$pertottpm <- filter_combined_table$tpm/filter_combined_table$totaltpm


TEreftable <- read.delim(repeatAnnotationFile, skip=1, header=FALSE, stringsAsFactors = FALSE)

indexsubfams <- match(filter_combined_table$subfamTE, TEreftable$V1)
TEclass_v1 <- TEreftable$V2[indexsubfams]
TEfamily_v1 <- TEreftable$V3[indexsubfams]

filter_combined_table$classTE <- TEclass_v1
filter_combined_table$familyTE <- TEfamily_v1

filter_combined_table$uniqidfile <- paste(filter_combined_table$subfamTE,filter_combined_table$startTE,filter_combined_table$gene1,filter_combined_table$exonintron1, filter_combined_table$number1, filter_combined_table$gene2, filter_combined_table$exonintron2, filter_combined_table$number2,filter_combined_table$transcriptstart2, filter_combined_table$file_label, sep = "_")
freqmultipletable <- as.data.frame(table(filter_combined_table$uniqidfile))

#Fix the total score so that it counts based on tumors NOT based on different isoforms. Thus, only one TE to gene isoform is selected based on its tpm.

i <- sapply(freqmultipletable, is.factor)
freqmultipletable[i] <- lapply(freqmultipletable[i], as.character)

freqmultipletable <- freqmultipletable[freqmultipletable$Freq >1,]

getNewStats <- function(transcoords){
  elements <- strsplit(transcoords, ",")[[1]]
  strandtrans <- elements[1]
  if (strandtrans == "+"){
    sortedvec <- sort(as.numeric(elements[3:length(elements)]))
    exon1length <- sortedvec[2] - sortedvec[1] + 1
    intron1start <- sortedvec[2] + 1
    intron1end <- sortedvec[3] - 1
  } else {
    sortedvec <- sort(as.numeric(elements[3:length(elements)]), decreasing = TRUE)
    exon1length <- sortedvec[1] - sortedvec[2] + 1
    intron1start <- sortedvec[3] + 1
    intron1end <-  sortedvec[2] - 1
  }
  return(c(exon1length, intron1start, intron1end))
}

resultsStats <- apply(filter_combined_table[,c('transcoord', 'tpm')],1,function(x) getNewStats(x[1]))
resultsStats  <- data.frame(t(resultsStats))
colnames(resultsStats) <- c("exon1length", "intron1start", "intron1end")

filter_combined_table <- cbind(filter_combined_table, resultsStats)

#There can be multiple trnscirpts for each uniq id combination. Those are thrown out based ont he highest tpm value
filter_combined_table <- filter_combined_table[order(-filter_combined_table$tpm),]
filter_combined_table <- filter_combined_table[!duplicated(filter_combined_table$uniqidfile),]

counts_table = filter_combined_table[,c('subfamTE', 'classTE', 'familyTE', 'chrTE', 'startTE', 'endTE', 'gene1', 'gene2', 'exonintron1', 'number1', 'exonintron2', 'number2', 'transcriptstart2', 'transcriptend2','strand','exon1length', 'tumorcount', 'normalcount', 'permaxcov', 'pertottpm','tpm')]
resexonlength <- aggregate(. ~ subfamTE + classTE + familyTE + chrTE + startTE + endTE + gene1 + gene2 + exonintron1 + number1 + exonintron2 + number2 + transcriptstart2 + transcriptend2 + strand, counts_table, FUN = sum)

i <- sapply(resexonlength, is.factor)
resexonlength[i] <- lapply(resexonlength[i], as.character)

#resexonlength$exonlengthaverage <- as.numeric(resexonlength$exon1length)/(as.numeric(resexonlength$normalcount) + as.numeric(resexonlength$tumorcount))

calcExonSkip <- function(gene1, gene2, number1, number2){
  if (gene1 != gene2 | (gene1 == 'None' & gene2 == 'None')){
    exonskipped <- -1
  } else {
    exonskipped <- as.numeric(number2)-as.numeric(number1)
  }
  return(exonskipped)
}

resexonlength$exonskipped <- apply(resexonlength[,c('gene1','gene2', 'number1', 'number2')],1,function(x) calcExonSkip(x[1],x[2],x[3],x[4]))

resexonlength <- resexonlength[resexonlength$exonskipped <= exonSkipMax,]

if (treatmentexclusive == 'yes'){
  resexonlength <- resexonlength[resexonlength$normalcount == 0,]
}
if (keepNone == 'no'){
  resexonlength <- resexonlength[resexonlength$gene2 != 'None',]
  resexonlength <- resexonlength[resexonlength$exonintron2 != 'intron',]
}

# Treatment minimum filter
resexonlength <- resexonlength[resexonlength$tumorcount >= totalTreatmentMin,]

# Treatment minimum filter
resexonlength <- resexonlength[(resexonlength$tumorcount + resexonlength$normalcount) >= totalMin,]

#TEfilter
if(filterTE == 'yes'){
  resexonlength <- resexonlength[resexonlength$classTE %in% c('LINE', 'SINE', 'SINE?', 'LTR', 'LTR?', 'DNA', 'DNA?', 'Retroposon'),]
}

resexonlength$ExonlengthAverage <- resexonlength$exon1length/(resexonlength$tumorcount + resexonlength$normalcount)
resexonlength<- resexonlength[resexonlength$ExonlengthAverage < exon1LengthMax,]

resexonlength$PerMaxCov <- resexonlength$permaxcov/(resexonlength$tumorcount + resexonlength$normalcount)
resexonlength$PerTotTpm <- resexonlength$pertottpm/(resexonlength$tumorcount + resexonlength$normalcount)
resexonlength$TpmAverage <- resexonlength$tpm/(resexonlength$tumorcount + resexonlength$normalcount)

resexonlength$uniqid <- paste(resexonlength$subfamTE,resexonlength$startTE,resexonlength$gene1,resexonlength$exonintron1, resexonlength$number1, resexonlength$gene2, resexonlength$exonintron2, resexonlength$number2,resexonlength$transcriptstart2,sep = "_")

filter_combined_table$uniqid <- paste(filter_combined_table$subfamTE,filter_combined_table$startTE,filter_combined_table$gene1,filter_combined_table$exonintron1, filter_combined_table$number1, filter_combined_table$gene2, filter_combined_table$exonintron2, filter_combined_table$number2,filter_combined_table$transcriptstart2,sep = "_")
filter_combined_table <- filter_combined_table[filter_combined_table$uniqid %in% resexonlength$uniqid,]

write.table(filter_combined_table, "filter_combined_candidates.tsv", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(resexonlength, "initial_candidate_list.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
save.image('Step4.RData')