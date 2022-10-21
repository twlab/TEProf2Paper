#!/usr/bin/env Rscript
library('Xmisc')

#This was created by Step 4 and has the data needed about each candidate.
load("Step12.RData")

#Get the directory of the rnapipeline scripts
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

#Default parameters
argumentFile <- paste0(script.basename,'/arguments.txt')
genomeName <- "BSgenome.Hsapiens.UCSC.hg38"

#User argument parsing code
parser <- ArgumentParser$new()

parser$add_argument('-g',type='character',help='A character string of the genome to use (Defualt: BSgenome.Hsapiens.UCSC.hg38)')
argsparse <- parser$get_args()

#Argument parse to change defaults
if (!identical(argsparse$g, character(0))){
  genomeName <- argsparse$g
}

detach("package:Xmisc", unload=TRUE)

print("Arguments")
print(paste0("Genome Name ", genomeName))

#Load the genome file of the species that is being studied
library('reshape2')
library('biomaRt')
library(genomeName, character.only = TRUE)

codingpotential_results <- read.delim('candidates_cpcout.fa', sep = "\t", stringsAsFactors = FALSE, header = TRUE)

coding_label <- c()
start_location <- c()
for (candidateid in annotatedcufftranscripts.fil$uniqid){
  coding_label <- c(coding_label, codingpotential_results[codingpotential_results$X.ID == candidateid, c('label')])
  start_location <- c(start_location, codingpotential_results[codingpotential_results$X.ID == candidateid, c('start_codon_position')])
}
annotatedcufftranscripts.fil$coding_label <- coding_label
annotatedcufftranscripts.fil$start_location <- start_location

#This is now the part of the code that will analyze the CPC2 longest openin reading frame results as well

calORFCPC2 <- function(newtranscript, oldtranscript, oldstartcodon, strongstart, codingpotential){
  
  strongstart <- as.numeric(strongstart) + 1
  newtlist = strsplit(newtranscript, ",")[[1]]
  strand = newtlist[1]
  chromosome = newtlist[2]
  for (i in 3:length(newtlist)) {
    if (!i %% 2){
      next
    } else{
      if (strand == "+"){
        if (i==3){
          exons = c(strtoi(newtlist[i]),strtoi(newtlist[i + 1]))
        } else {
          exons = rbind(exons, c(strtoi(newtlist[i]),strtoi(newtlist[i + 1])))
        }
      }else {
        if (i==3){
          exons = c(strtoi(newtlist[i+1]),strtoi(newtlist[i]))
        } else {
          exons = rbind(exons, c(strtoi(newtlist[i+1]),strtoi(newtlist[i])))
        }
      }
    }
  }
  if (is.null(dim(exons))){
    exons = t(as.matrix(exons))
  }
  numexons = dim(exons)[1]
  if (numexons > 1){
    exons <- exons[order(exons[,1]),]
  }
  locationexons = c()
  sequenceexons = ""
  for (i in 1:numexons){
    locationexons = append(locationexons,seq.int(exons[i,1], exons[i,2]))
    sequenceexons = paste0(sequenceexons, toString(getSeq(Hsapiens, chromosome, exons[i,1], exons[i,2])))
  }
  
  if (strand == "-"){
    sequenceexons = toString(reverseComplement(DNAString(sequenceexons)))
    locationexons = rev(locationexons)
  }
  strongproteintrans = "None"
  strongstartloc = "None"
  strongendloc = "None"
  
  # Find the start codon
  
  if (codingpotential == "coding") {
    proteintrans = toString(translate(DNAString(substring(sequenceexons,strongstart,nchar(sequenceexons)))))
    endloc = gregexpr(pattern ='\\*',proteintrans)[[1]][1]
    strongproteintrans = proteintrans
    fullproteintrans = proteintrans
    strongendlocfactor = endloc*3
    strongstartloc = locationexons[strongstart]
  } else {
    strongstart = "None"
  }
          
  
  if (strongstart != "None"){
    
    # This will generate a tuple list of all the dimensions of the exons of the old transcript (exons2) this also includes the introns
    oldtlist = strsplit(oldtranscript, ",")[[1]]
    if(oldstartcodon != "None"){
      oldstartbeg <- strtoi(strsplit(oldstartcodon, ";")[[1]][2])
      oldstartend <- strtoi(strsplit(oldstartcodon, ";")[[1]][3])
      
      indexcorrection <- 0
      if (strand == "-"){
        
        if ((oldstartend - oldstartbeg) < 2){
          indexcorrection <- (2-(oldstartend - oldstartbeg))
        }
        
        oldstart = oldstartend
      } else {
        
        if ((oldstartend - oldstartbeg) < 2){
          indexcorrection <- (2-(oldstartend - oldstartbeg))
        }
        
        oldstart = oldstartbeg
        
      }
    } else{
      oldstart = "None"
      indexcorrection <- 0
    }
    if (length(oldtlist) == 2){
      exons2 = t(as.matrix(c(strtoi(oldtlist[1]),strtoi(oldtlist[2]))))
    } else {
      for (i in 1:length(oldtlist)) {
        if (!i %% 2){
          next
        } else {
          if (i==1){
            exons2 = c(strtoi(oldtlist[i]),strtoi(oldtlist[i + 1]))
          } else {
            exons2 = rbind(exons2, c(strtoi(oldtlist[i]),strtoi(oldtlist[i + 1])))
          }
        }
      }
    }
    
    
    # Find the last intron/exon junction for NMD analysis
    if (length(oldtlist) > 2){
      if (strand == "-"){
        exonsmd = exons2[order(-exons2[,1]),]
      } else {
        exonsmd = exons2[order(exons2[,1]),]
      }
    } else {
      exonsmd = exons2
    }
    
    # Find the junction that is last
    lastjunc = "None"
    if (dim(exonsmd)[1] > 1){
      lastjunc = as.integer(exonsmd[dim(exonsmd)[1] - 2,2])
    }
    
    # This generates the sequence
    numexons = dim(exons2)[1]
    if (dim(exonsmd)[1] > 1){
      exons2 <- exons2[order(exons2[,1]),]
    }
    locationexonsorig = c()
    sequenceexonsorig = ""
    for (i in 1:numexons){
      if (!i %% 2){
        next
      } else{
        locationexonsorig = append(locationexonsorig,seq.int(exons2[i,1], exons2[i,2]))
        #print(exons2)
        sequenceexonsorig = paste0(sequenceexonsorig, toString(getSeq(Hsapiens, chromosome, strtoi(exons2[i,1]), strtoi(exons2[i,2]))))
      }
    }
    
    if (strand == "-"){
      sequenceexonsorig = toString(reverseComplement(DNAString(sequenceexonsorig)))
      locationexonsorig = rev(locationexonsorig)
    }
    
    lastjuncloc = "None"
    if (dim(exonsmd)[1] > 1){
      lastjuncloc = match(lastjunc, locationexonsorig)
    }
    
    # Create a vector with frame informtion (1 being start of codon) to determine the frame of the new transcript
    
    indexoldstart = match(oldstart, locationexonsorig) - indexcorrection
    if (is.na(indexoldstart)){
      returnstrongstart = "None"
      returnmatchloc = "None"
      returnproteintrans = "None"
      returnfullproteintrans = "None"
      returnframe = "None"
      returnlabeltype = "None"
      returnnmddistance = "None"
      returnnmdcall = "None"
      returnproteinleft = "None"
      returnframecall = "None"
    } else {
      
      # Find the 'frame' of the transcript using the modulus of the index of the start within the vector.
      modoldstart = indexoldstart %% 3
      numcodon = ceiling(nchar(sequenceexonsorig)/3)
      if (modoldstart == 0){
        codonvec = rep(c(2,3,1), numcodon)
      } else if (modoldstart == 1){
        codonvec = rep(c(1,2,3), numcodon)
      } else if (modoldstart == 2){
        codonvec = rep(c(3,1,2), numcodon)
      }
      
      firstmatch = "None"
      frameinfo = "None"
      i = 0
      for (loc in locationexonsorig){
        i = i + 1
        if (is.na(match(loc, locationexons))){
          next
        } else {
          frameinfo = codonvec[i]
          firstmatch = match(loc, locationexons)
          origlocmatch = match(loc, locationexonsorig)
          break
        }
      }
      
      if (firstmatch != "None"){
        # Calculate the frame
        if (frameinfo == 1){
          inframe = firstmatch
          # the codons are the same in this case so there will not be a difference
          oldcodon = "None"
          newcodon = "None"
        } else if (frameinfo == 2){
          inframe = firstmatch + 2
          # in case it is in-frame chimeric, this will be used to determine if this codon is new
          oldcodon = substring(sequenceexonsorig, origlocmatch-1, origlocmatch+1)
          newcodon = substring(sequenceexons, firstmatch-1, firstmatch+1)
        } else if (frameinfo == 3) {
          inframe = firstmatch + 1
          # in case it is in-frame chimeric, this will be used to determine if this codon is new
          oldcodon = substring(sequenceexonsorig, origlocmatch-2, origlocmatch)
          newcodon = substring(sequenceexons, firstmatch-2, firstmatch)
        }
        
        modnewstart = inframe %% 3
        numcodon2 = ceiling(nchar(sequenceexons)/3)
        
        if (modnewstart == 0){
          codonvec2 = rep(c(2,3,1), numcodon2)
        } else if (modnewstart == 1){
          codonvec2 = rep(c(1,2,3), numcodon2)
        } else if (modnewstart == 2){
          codonvec2 = rep(c(3,1,2), numcodon2)
        }
        
        if (codonvec2[strongstart] == 1){
          framenew = "in-frame"
        } else {
          framenew = "out-of-frame"
        }
        
        
        # Figure out if the protein is the same, the same truncated, chimeric same, chimeric truncated.
        
        labeltype = "None"
        nmdcall = "None"
        
        
        disnew = firstmatch - strongstart
        disold = origlocmatch - indexoldstart
        
        chimericcheck = FALSE
        
        if (framenew == "in-frame"){
          if (disnew > 0){
            if (disold > 0){
              labeltype = "chimeric truncated"
              chimericcheck = TRUE
            } else {
              labeltype = "chimeric normal"
              chimericcheck = TRUE
            }
          }else {
            if (disold > 0){
              labeltype = "truncated"
            } else {
              # The original startcodon should be chosen if both were present in the original transcript
              framenew = "in-frame"
              strongstart2 = match(oldstart, locationexons) 
              if (is.na(strongstart2)){
                strongstart = "None"
                strongproteintrans = "None"
                fullproteintrans = "None"
              } else {
                strongstart = strongstart2
                strongproteintrans = toString(translate(DNAString(substring(sequenceexons,strongstart,nchar(sequenceexons)))))
                fullproteintrans = toString(translate(DNAString(substring(sequenceexons,strongstart,nchar(sequenceexons)))))
              }
              labeltype = "normal"
            }
          }
        } else {
          if (disnew <= 0 & disold <= 0){
            locationstartmatch <- match(oldstart, locationexons)
            if (is.na(locationstartmatch)){
              labeltype = "out-of-frame"
            } else {
              diffstarts <- strongstart - locationstartmatch 
              if ((diffstarts %% 3) == 0){
                if (strongstart >= locationstartmatch){
                  strongproteintrans = toString(translate(DNAString(substring(sequenceexons,strongstart,nchar(sequenceexons)))))
                  labeltype = "normal"
                  framenew = "in-frame"
                } else {
                  strongproteintrans = toString(translate(DNAString(substring(sequenceexons,strongstart,nchar(sequenceexons)))))
                  labeltype = "chimeric normal"
                  chimericcheck = TRUE
                  framenew = "in-frame"
                }
              } else {
                labeltype = "out-of-frame"
              }
            }
          } else {
            labeltype = "out-of-frame"
          }
        }
        
        
        aaCorrection = 0
        if (strongstart != "None"){
          # Truncating protein sequences for the chimeric ones. Adding 1
          numexons2 = ceiling((inframe + 30 - strongstart)/3)
          
          # In the case that the TE- transcript includes exon 1 and it starts before the previous start codon, the following adjustment needs to be made
          if (labeltype == "chimeric normal"){
            numexons2 = ceiling((firstmatch - disold + 30 - strongstart)/3)
            oldcodon = "None"
          }
          if (framenew == "in-frame" && chimericcheck){
            if (oldcodon != "None"){
              aaold = toString(translate(DNAString(oldcodon)))
              aanew = toString(translate(DNAString(newcodon)))
              if (aaold == aanew){
                strongproteintrans = substring(strongproteintrans,1,(numexons2-1))
                aaCorrection = -3
              } else {
                strongproteintrans = substring(strongproteintrans,1,(numexons2))
              }
            } else{
              strongproteintrans = substring(strongproteintrans,1,(numexons2-1))
            }
            endloc = gregexpr(pattern ='\\*',strongproteintrans)[[1]][1]
            # This isn't perfect if it is only chimeric with like 9-10 aa of the original but that would be unusual
            if (endloc != -1 && (nchar(strongproteintrans)-endloc > 9)){
              labeltype = "novel non-chimeric"
            }
          }
        }
        
        endloc = gregexpr(pattern ='\\*',strongproteintrans)[[1]][1]
        if (endloc != -1){
          strongproteintrans2 = substring(strongproteintrans, 1, endloc - 1)
          strongproteintrans = strongproteintrans2
        }
        # See if this peptide will likely go through nonsense mediated decay
        nmdthresh = 55
        distancelastjun = "None"
        if (dim(exons2)[1] > 1){
          if (strongendlocfactor == -3){
            nmdcall = "Inconclusive"
          } else {
            #print(c(lastjuncloc, origlocmatch, disnew, strongendlocfactor))
            distancelastjun = lastjuncloc - origlocmatch + disnew - strongendlocfactor
            if (distancelastjun > nmdthresh){
              nmdcall = "Yes"
            } else {
              nmdcall = "No"
            }
          }
        } else {
          nmdcall = "No"
        }
        
        framecorrection <- codonvec[origlocmatch]
        if(framecorrection == 1){
          correctionFactor = 0
        } else if(framecorrection == 2){
          correctionFactor = 2
        } else if(framecorrection == 3){
          correctionFactor = 1
        }
        
        origProteinTranslate <- toString(translate(DNAString(substring(sequenceexonsorig,origlocmatch + correctionFactor + aaCorrection ,nchar(sequenceexonsorig)))))
        
        returnstrongstart = strongstart
        returnmatchloc = firstmatch
        returnproteintrans = strongproteintrans
        returnfullproteintrans = fullproteintrans
        returnframe = framenew
        returnlabeltype = labeltype
        returnnmddistance = distancelastjun
        returnnmdcall = nmdcall
        returnproteinleft = nchar(strsplit(origProteinTranslate, '\\*')[[1]][1])
        returnframecall = codonvec[1]
      } else{
        returnstrongstart = strongstart
        returnmatchloc = firstmatch
        returnproteintrans = strongproteintrans
        returnfullproteintrans = fullproteintrans
        returnframe = "None"
        returnlabeltype = "None"
        returnnmddistance = "None"
        returnnmdcall = "None"
        returnproteinleft = "None"
        returnframecall = "None"
      }
    }
  } else {
    returnstrongstart = "None"
    returnmatchloc = "None"
    returnproteintrans = "None"
    returnfullproteintrans = "None"
    returnframe = "None"
    returnlabeltype = "None"
    returnnmddistance = "None"
    returnnmdcall = "None"
    returnproteinleft = "None"
    returnframecall = "None"
  }
  
  return(c(returnstrongstart, returnmatchloc, returnproteintrans, returnfullproteintrans, returnframe, returnlabeltype, returnnmddistance, returnnmdcall, returnproteinleft, returnframecall))
}

# Do analysis on if the sequence is in frame or out of frame relative to the main transcript. This is how many threads. 

results <- apply(annotatedcufftranscripts.fil[,c('transcoord','elements2','startcodon2', 'start_location','coding_label')], 1, function(x) calORFCPC2(x[1], x[2], x[3], x[4], x[5]))
#stopCluster(clus)
results = data.frame(t(results))
columnnames <- colnames(annotatedcufftranscripts.fil)
annotatedcufftranscripts.fil = cbind.data.frame(annotatedcufftranscripts.fil, results)
columnlabels <- append(columnnames, c("cpc2_strongstartindex", "cpc2_matchbetweenindex", "cpc2_proteinseq", "cpc2_proteinseqfull", "cpc2_frame", "cpc2_type", "cpc2_nmddistance", "cpc2_nmdcall", "cpc2_origprotleft", "cpc2_frameNum"))
colnames(annotatedcufftranscripts.fil) <- columnlabels

annotatedcufftranscripts.fil$cpc2_proteinseqfull2 <- apply(annotatedcufftranscripts.fil[,c('cpc2_proteinseqfull','uniqid')],1, function(x) strsplit(x[1],"\\*")[[1]][1])

#Adjusting the labels for new definitions. 

newLabelsKozak <- apply(annotatedcufftranscripts.fil[,c("strongstartindex", "matchbetweenindex", "type", "proteinseqfull2")],1, function(x){
  if (x[3] == "chimeric normal" | x[3] == "chimeric truncated" | x[3] == "truncated" | x[3] == "normal" | x[3] == "None"){
    if (x[4] == "None"){
      return("None")
    } else {
      return(x[3])
    }
  } else if(x[3] == "novel non-chimeric" | x[3] == "out-of-frame"){
    return("out-of-frame")
  }
})

newLabelsCPC2 <- apply(annotatedcufftranscripts.fil[,c("cpc2_strongstartindex", "cpc2_matchbetweenindex", "cpc2_type", "cpc2_proteinseqfull2")],1, function(x){
  if (x[3] == "chimeric normal" | x[3] == "chimeric truncated" | x[3] == "truncated" | x[3] == "normal" | x[3] == "None"){
    if (x[4] == "None"){
      return("None")
    } else {
      return(x[3])
    }
  } else if(x[3] == "novel non-chimeric" | x[3] == "out-of-frame"){
    return("out-of-frame")
  }
})

annotatedcufftranscripts.fil$type_final <- newLabelsKozak 
annotatedcufftranscripts.fil$cpc2_type_final <- newLabelsCPC2

#Final code that will output a fasta with teh sequences of teh potentially antigenic candidates. Those are the ones that are chimeric or out of frame

annotatedcufftranscripts.fil$samecall <- apply(annotatedcufftranscripts.fil[,c("type", "cpc2_type")], 
                                                 1,
                                                 function(x){
                                                   if (x[1] == x[2]){
                                                     return("Yes")
                                                   } else {
                                                     return("No")
                                                   }
                                                 })

allowedtypes <- c("chimeric normal","chimeric truncated","out-of-frame","novel non-chimeric")

annotatedcufftranscripts.fil.c <- annotatedcufftranscripts.fil

fastaoutputvector <- c()

for (i in 1:nrow(annotatedcufftranscripts.fil.c)){
  if (annotatedcufftranscripts.fil.c$samecall[i] == "Yes" & annotatedcufftranscripts.fil.c$type[i] %in% allowedtypes){
     fastaoutputvector <- c(fastaoutputvector, paste0(">tr|",annotatedcufftranscripts.fil.c$uniqid[i],"_both|",annotatedcufftranscripts.fil.c$uniqid[i],"_HUMAN_both ",annotatedcufftranscripts.fil.c$uniqid[i],"_HUMAN_both OS=Homo sapiens OX=9606 GN=",annotatedcufftranscripts.fil.c$gene2[i]))
        fastaoutputvector <- c(fastaoutputvector,annotatedcufftranscripts.fil.c$cpc2_proteinseqfull2[i])
  } else {
      if (annotatedcufftranscripts.fil.c$type[i] %in% allowedtypes){
        fastaoutputvector <- c(fastaoutputvector, paste0(">tr|",annotatedcufftranscripts.fil.c$uniqid[i],"_og|",annotatedcufftranscripts.fil.c$uniqid[i],"_HUMAN_og ",annotatedcufftranscripts.fil.c$uniqid[i],"_HUMAN_og OS=Homo sapiens OX=9606 GN=",annotatedcufftranscripts.fil.c$gene2[i]))
        fastaoutputvector <- c(fastaoutputvector,annotatedcufftranscripts.fil.c$proteinseqfull2[i])
      }
    
      if (annotatedcufftranscripts.fil.c$cpc2_type[i] %in% allowedtypes){
        fastaoutputvector <- c(fastaoutputvector, paste0(">tr|",annotatedcufftranscripts.fil.c$uniqid[i],"_cpc2|",annotatedcufftranscripts.fil.c$uniqid[i],"_HUMAN_cpc2 ",annotatedcufftranscripts.fil.c$uniqid[i],"_HUMAN_cpc2 OS=Homo sapiens OX=9606 GN=",annotatedcufftranscripts.fil.c$gene2[i]))
        fastaoutputvector <- c(fastaoutputvector,annotatedcufftranscripts.fil.c$cpc2_proteinseqfull2[i])
      }
  }
}

fileConn<-file("candidates_proteinseq.fa")
writeLines(fastaoutputvector, fileConn)
close(fileConn)

save.image('Step13.RData')