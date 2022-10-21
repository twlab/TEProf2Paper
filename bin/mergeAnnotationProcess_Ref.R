#!/usr/bin/env Rscript
library('Xmisc')

#Get the directory of the rnapipeline scripts
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

#Default parameters
exonSkipMax <- 2
keepNone <- 'no'
filterTE <- 'yes'
argumentFile <- paste0(script.basename,'/arguments.txt')
annotationFile <- 'reference_merged_candidates.gff3_annotated_filtered_test_all'

#User argument parsing code
parser <- ArgumentParser$new()

parser$add_argument('-f',type='character',help='The reference annotation file to use for processing')
parser$add_argument('-s',type='character',help='A number representing the maximum numbers of exons that can be skipped by transcript')
parser$add_argument('-k',type='character',help='If yes is entered, then candidates that do not end up splicing into a gene will be included')
parser$add_argument('-t',type='character',help='If yes is entered, only the following TE classes will be used (LINE, SINE, SINE?, LTR, LTR?, DNA, DNA?, and Retroposon in Repeatmasker)')
parser$add_argument('-a',type='character',help='The arguments.txt file that is being used for this run.')
argsparse <- parser$get_args()

#Argument parse to change defaults
if (!identical(argsparse$f, character(0))){
  annotationFile <- argsparse$f
}

if (!identical(argsparse$s, character(0))){
  exonSkipMax <- as.numeric(argsparse$s)
}

if (!identical(argsparse$k, character(0))){
  keepNone <- argsparse$k
}

if (!identical(argsparse$t, character(0))){
  filterTE <- argsparse$t
}

if (!identical(argsparse$a, character(0))){
  argumentFile <- argsparse$a
}

print("Arguments")
print(paste0("Reference annotation name: ", annotationFile ))
print(paste0("Exon Skip Max: ", exonSkipMax))
print(paste0("Keep None Splice Target: ", keepNone))
print(paste0("Filter for only TE classes: ", filterTE))
print(paste0("Argument File Location: ", argumentFile))

argTable <- read.delim(argumentFile, header=FALSE, stringsAsFactors = FALSE)
repeatAnnotationFile <- argTable[argTable[,1] == 'rmskannotationfile',c(2)]

annotatedcufftranscripts <-read.delim(annotationFile, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
columnlabels <- c("transcriptname", "type1", "gene1", "chr1", "start1", "end1", "exonintron1", "number1" ,"startcodon1" , "transcriptstart1", "transcriptend1", "elements1", "id1", "splicing", "type2", "gene2", "chr2", "start2", "end2", "exonintron2", "number2" ,"startcodon2", "transcriptstart2", "transcriptend2", "elements2", "id2", "chromtrans", "starttrans", "endtrans", "transcoord", "intronanno", "chrTE", "startTE", "endTE", "subfamTE", "numTE", "strandTE", "strand")
colnames(annotatedcufftranscripts) <- columnlabels

#This is one exaple where a different name for an isoofrm of a gene is causing issues in tabulation. 

annotatedcufftranscripts[annotatedcufftranscripts == "RP11-545J16.1"] <- 'SLCO1B3'
annotatedcufftranscripts <- annotatedcufftranscripts[annotatedcufftranscripts$exonintron1 != 'exon',]

annotatedcufftranscripts$uniqid <- paste(annotatedcufftranscripts$subfamTE,annotatedcufftranscripts$startTE,annotatedcufftranscripts$gene1,annotatedcufftranscripts$exonintron1, annotatedcufftranscripts$number1, annotatedcufftranscripts$gene2, annotatedcufftranscripts$exonintron2, annotatedcufftranscripts$number2,annotatedcufftranscripts$transcriptstart2,sep = "_")

#Annotate the TEs with class and family information

TEreftable <- read.delim(repeatAnnotationFile, skip=1, header=FALSE, stringsAsFactors = FALSE)

indexsubfams <- match(annotatedcufftranscripts$subfamTE, TEreftable$V1)
TEclass_v1 <- TEreftable$V2[indexsubfams]
TEfamily_v1 <- TEreftable$V3[indexsubfams]

annotatedcufftranscripts$classTE <- TEclass_v1
annotatedcufftranscripts$familyTE <- TEfamily_v1

#Calculate additional statistics on the transcripts and filter based on exon skipping and exon 1 length. 

calcExonSkip <- function(gene1, gene2, number1, number2){
  if (gene1 != gene2 | (gene1 == 'None' & gene2 == 'None')){
    exonskipped <- -1
  } else {
    exonskipped <- as.numeric(number2)-as.numeric(number1)
  }
  return(exonskipped)
}

annotatedcufftranscripts$exonskipped <- apply(annotatedcufftranscripts[,c('gene1','gene2', 'number1', 'number2')],1,function(x) calcExonSkip(x[1],x[2],x[3],x[4]))

annotatedcufftranscripts <- annotatedcufftranscripts[annotatedcufftranscripts$exonskipped <= exonSkipMax,]

if (keepNone == 'no'){
  annotatedcufftranscripts <- annotatedcufftranscripts[annotatedcufftranscripts$gene2 != 'None',]
  annotatedcufftranscripts <- annotatedcufftranscripts[annotatedcufftranscripts$exonintron2 != 'intron',]
}

if(filterTE == 'yes'){
  annotatedcufftranscripts <- annotatedcufftranscripts[annotatedcufftranscripts$classTE %in% c('LINE', 'SINE', 'SINE?', 'LTR', 'LTR?', 'DNA', 'DNA?', 'Retroposon'),]
}

#Obtain the main intron junction that goes from the TE-transcript to the gene. Subsequent steps can use coverage of this intron as an additional filter. 

getIntronJun <- function(structurestring, exonstart, exonend){
  elementsvec <- as.numeric(tail(strsplit(structurestring, ",")[[1]],-2))
  elementsvec <- sort(elementsvec)
  strand <- strsplit(structurestring, ",")[[1]][1]
  
  if (length(elementsvec) == 4){
    intron1start = elementsvec[2] + 1
    intron1end = elementsvec[3] - 1
    return(c(intron1start,intron1end))
  }
  
  if (strand == "+"){
    for (i in seq(4,length(elementsvec), by = 2)){
      if (elementsvec[i] >= exonstart){
        intron1start = elementsvec[i-2] + 1
        intron1end = elementsvec[i-1] - 1
        return(c(intron1start,intron1end))
        break
      }
    }
    
  } else {
    elementsvec <- sort(elementsvec, decreasing = TRUE)
    for (i in seq(4,length(elementsvec), by = 2)){
      if (elementsvec[i] <= exonend){
        intron1start = elementsvec[i-1] + 1
        intron1end = elementsvec[i-2] - 1
        return(c(intron1start,intron1end))
        break
      }
    }
  }
}

results <- apply(annotatedcufftranscripts[,c("transcoord","start2", "end2")],1,function(x) getIntronJun(x[1],as.numeric(x[2]),as.numeric(x[3])))
results <- as.data.frame(t(results))
colnames(results) <- c("intronjunstart","intronjunend")
annotatedcufftranscripts <- cbind(annotatedcufftranscripts, results)

grepintrons <- paste0(annotatedcufftranscripts$chrTE, "\t", annotatedcufftranscripts$strand,"\t",annotatedcufftranscripts$intronjunstart,"\t",annotatedcufftranscripts$intronjunend)
grepintrons <- paste(grepintrons, collapse = '\n')

write(grepintrons, file = "candidate_introns.txt")
write(annotatedcufftranscripts$transcriptname, file = "candidate_names.txt")
save.image('Step10.RData')

