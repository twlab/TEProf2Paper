#!/usr/bin/env Rscript
library('Xmisc')

#This was created by Step 4 and has the data needed about each candidate.
load("Step11_FINAL.RData")

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

#Filtering the transcripts
#Translation can be expensive, and thus there is ususally a step of filtering the transcripts that are exclusive to certain samples. 
#However, for now, no filtering is done to allow maximum flexibility per user

annotatedcufftranscripts.fil <- annotatedcufftranscripts

#Function to call ORF using the Kozak method
calORF <- function(newtranscript, oldtranscript, oldstartcodon){
  
  #Testing Purposes
  #annotatedcufftranscripts.lin28b <- annotatedcufftranscripts[annotatedcufftranscripts$gene2 == 'NLRC4', ]
  #newtranscript <- as.character(annotatedcufftranscripts.lin28b$transcoord[1])
  #oldtranscript <- as.character(annotatedcufftranscripts.lin28b$elements2[1])
  #oldstartcodon <- as.character(annotatedcufftranscripts.lin28b$startcodon2[1])
  
  library('BSgenome.Hsapiens.UCSC.hg38')
  
  #print(oldtranscript)
  #print(oldstartcodon)
  newtlist = strsplit(newtranscript, ",")[[1]]
  strand = newtlist[1]
  chromosome = newtlist[2]
  for (i in 3:length(newtlist)) {
    if (!i %% 2){
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
  startcodons <- gregexpr(pattern ='ATG',sequenceexons)[[1]]
  strongstart = "None"
  strongproteintrans = "None"
  strongstartloc = "None"
  strongendloc = "None"
  
  # Find the start codon
  for (startcodon in startcodons){
    
    substringcodon = substr(sequenceexons,startcodon-6,startcodon+3)
    substringcodonvec = c()
    for (i in 1:nchar(substringcodon)){
      substringcodonvec = append(substringcodonvec, substr(substringcodon,i,i))
    }
    kosak1 = c("G","C","C","G","C","C","A","T","G","G")
    kosak2 = c("G","C","C","A","C","C","A","T","G","G")
    trusum1 = sum(kosak1 == substringcodonvec)
    trusum2 = sum(kosak2 == substringcodonvec)
    score = max(trusum1, trusum2)
    
    # 80% consensus was too much, so I decided the purine position which is important and having 60% consensus was enough
    if (score >= 4){
      if (substr(substringcodon,4,4) == "A" || substr(substringcodon,4,4) == "G" || substr(substringcodon,10,10) == "G"){
        proteintrans = toString(translate(DNAString(substring(sequenceexons,startcodon,nchar(sequenceexons)))))
        endloc = gregexpr(pattern ='\\*',proteintrans)[[1]][1]
        if (endloc > 50 || (endloc == -1 && nchar(proteintrans) > 50)){
          strongstart = startcodon
          strongproteintrans = proteintrans
          fullproteintrans = proteintrans
          strongendlocfactor = endloc*3
          strongstartloc = locationexons[strongstart]
          break
        }
      }
    }
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
        
        # This is going to note for if the 
        
        aamatch = "no"
        
        # Figure out if the protein is the same, the same truncated, chimeric same, chimeric truncated.
        
        labeltype = "None"
        nmdcall = "None"
        
        
        disnew = firstmatch - startcodon
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
                  framenew = "in-frame"
                  chimericcheck = TRUE
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
            # This isn't perfect if it is only chimeric with 9-10 aa of the original but that would be unusual
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


#Left some old code commented out that allowed this to be run in parallel. Initial version will keep to single thread. 

# clus <- makeCluster(1)
# clusterEvalQ(clus, .libPaths( "/bar/nshah/R/x86_64-pc-linux-gnu-library/3.5"))
# clusterExport(clus, "calORF") <This was here for me since versionsing issues made it so parallel parApply only worked with a specific version of R. May not need this so commented out. 
#results <- parApply(clus,annotatedcufftranscripts.fil[,c('transcoord','elements2','startcodon2')], 1, function(x) calORF(x[1], x[2], x[3]))

results <- apply(annotatedcufftranscripts.fil[,c('transcoord','elements2','startcodon2')], 1, function(x) calORF(x[1], x[2], x[3]))

#stopCluster(clus)
results = data.frame(t(results))
columnnames <- colnames(annotatedcufftranscripts.fil)
annotatedcufftranscripts.fil = cbind.data.frame(annotatedcufftranscripts.fil, results)
columnlabels <- append(columnnames, c("strongstartindex", "matchbetweenindex", "proteinseq", "proteinseqfull", "frame", "type", "nmddistance", "nmdcall", "origprotleft","frameNum"))
colnames(annotatedcufftranscripts.fil) <- columnlabels

#Now you have the protein sequence information from Kozak method
annotatedcufftranscripts.fil$proteinseqfull2 <- apply(annotatedcufftranscripts.fil[,c('proteinseqfull','uniqid')],1, function(x) strsplit(x[1],"\\*")[[1]][1])

#This function will allow you to get the RNA sequence from coordinates of the transcript
getSequenceFromCoords <- function(newtranscript){
  
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
  return(sequenceexons)
}

#clus <- makeCluster(8)

#clusterEvalQ(clus, .libPaths( "/bar/nshah/R/x86_64-pc-linux-gnu-library/3.5"))
#clusterExport(clus, "getSequenceFromCoords")
#annotatedcufftranscripts.fil$rnasequence <- parApply(clus, annotatedcufftranscripts.fil[,c('uniqid','transcoord')], 1, function(x) getSequenceFromCoords(x[2]))
annotatedcufftranscripts.fil$rnasequence <- apply(annotatedcufftranscripts.fil[,c('uniqid','transcoord')], 1, function(x) getSequenceFromCoords(x[2]))
#stopCluster(clus)

#Now we will output the RNA sequences into a fast faile. 
fastaoutputvector <- c()

for (i in 1:nrow(annotatedcufftranscripts.fil)){
  fastaoutputvector <- c(fastaoutputvector, paste0(">",annotatedcufftranscripts.fil$uniqid[i]))
  fastaoutputvector <- c(fastaoutputvector,annotatedcufftranscripts.fil$rnasequence[i])
}

fileConn<-file("candidates.fa")
writeLines(fastaoutputvector, fileConn)
close(fileConn)

#Final Image File
save.image('Step12.RData')