#!/usr/bin/env Rscript
library('Xmisc')

load('Step6.RData')

#Default parameters
annotationFile <- 'reference_merged_candidates.gff3_annotated_filtered_test_all'

#User argument parsing code
parser <- ArgumentParser$new()

parser$add_argument('-f',type='character',help='The reference annotation file to use for processing')
argsparse <- parser$get_args()

#Argument parse to change defaults
if (!identical(argsparse$f, character(0))){
  annotationFile <- argsparse$f
}

#The Rsession from previous steps is reloaded

annotatedcufftranscripts <-read.delim(annotationFile, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
columnlabels <- c("transcriptname", "type1", "gene1", "chr1", "start1", "end1", "exonintron1", "number1" ,"startcodon1" , "transcriptstart1", "transcriptend1", "elements1", "id1", "splicing", "type2", "gene2", "chr2", "start2", "end2", "exonintron2", "number2" ,"startcodon2", "transcriptstart2", "transcriptend2", "elements2", "id2", "chromtrans", "starttrans", "endtrans", "transcoord", "intronanno", "chrTE", "startTE", "endTE", "subfamTE", "numTE", "strandTE", "strand")
colnames(annotatedcufftranscripts) <- columnlabels

annotatedcufftranscripts[annotatedcufftranscripts == "RP11-545J16.1"] <- 'SLCO1B3'
annotatedcufftranscripts <- annotatedcufftranscripts[annotatedcufftranscripts$exonintron1 != 'exon',]

annotatedcufftranscripts$uniqid <- paste(annotatedcufftranscripts$subfamTE,annotatedcufftranscripts$startTE,annotatedcufftranscripts$gene1,annotatedcufftranscripts$exonintron1, annotatedcufftranscripts$number1, annotatedcufftranscripts$gene2, annotatedcufftranscripts$exonintron2, annotatedcufftranscripts$number2,annotatedcufftranscripts$transcriptstart2,sep = "_")

#Only candidates passing the previous thressholds will be considered. 
annotatedcufftranscripts <- annotatedcufftranscripts[annotatedcufftranscripts$uniqid %in% res2statfil$uniqid,]


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

