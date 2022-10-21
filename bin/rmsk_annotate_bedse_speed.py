#!/usr/bin/python
# programmer : nshah 
# usage: To be used in order to annotate a GTF file for presence of transcripts that begin in transposable elements. 

from __future__ import division #integer division has to be //, now all / are floating point
import sys

allowedchr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

TEloc = sys.argv[1]
TEstart = int(TEloc.split(',')[1])
TEend = int(TEloc.split(',')[2])
geneElements = sys.argv[2]
elementsVec = [int(x) for x in geneElements.split(',')]
elementsVec.sort()
strand = sys.argv[3]
totelements = len(elementsVec)
refIntronNum = sys.argv[4]
refExonNum = sys.argv[5]

# The file name of the peak file has information on the chromosome, start, end, and strand of the peak
totreads = 0
startreads = 0
endreads = 0


for eachline in sys.stdin:
    
    #print elementsVec
    temp = eachline.strip().split("\t")
    
    chr_read = temp[0]
    
    if chr_read not in allowedchr:
        continue
    
    start_read = min([int(temp[i]) for i in [1,2]])
    end_read = max([int(temp[i]) for i in [1,2]])
    
    startinTE = False
    endinTE = False
    
    if start_read >= TEstart and start_read <= TEend:
        startinTE = True
        
    if end_read >= TEstart and end_read <= TEend:
        endinTE = True
    
    #print(startinTE)
    #print(endinTE)
    
    #If none of the start or end of the paired end are in the TE then discard  
    if not (startinTE or endinTE):
        continue
    
    #If both the start and end are in the TE, then that means that the read is within the TE and should only count for the TE read count
    elif startinTE and endinTE:
        totreads += 1
        
    #If only one end is in the TE then it is possible this could be a start or end read. 
    else:
        
        totreads += 1
        
        if strand == '+':
            
            # Only the start is within the TE so there is a change that the end is in the exons of the splice target transcript giving a start read
            if startinTE:
                
                indexOfGap = 0
                c = 1
                for i,k in zip(elementsVec[0::2], elementsVec[1::2]):
                    if end_read >= i and end_read <= k:
                        indexOfGap = c
                        break
                    c = c + 1
                
                if indexOfGap % 2 != 0:
                    exonNum = (indexOfGap + 1)/2
                    if exonNum >= int(refExonNum):
                        startreads += 1
            
            # Only the end of the read is within the TE. 
            else: 
                
                if refIntronNum == "None":
                    # This is intergenic so an end-read is not possible
                    continue
                
                indexOfGap = 0
                c = 1
                for i,k in zip(elementsVec[0::2], elementsVec[1::2]):
                    if start_read >= i and start_read <= k:
                        indexOfGap = c
                        break
                    c = c + 1
                
                if indexOfGap % 2 != 0:
                    exonNum = (indexOfGap + 1)/2
                    if exonNum <= int(refIntronNum):
                        endreads += 1  
        
        if strand == '-':
            
            # Only the start is within the TE so there is a change that the end is in the exons of the splice target transcript giving a start read
            if startinTE:
                
                if refIntronNum == "None":
                    # This is intergenic so an end-read is not possible
                    continue
                
                indexOfGap = 0
                c = 1
                for i,k in zip(elementsVec[0::2], elementsVec[1::2]):
                    if end_read >= i and end_read <= k:
                        indexOfGap = c
                        break
                    c = c + 1
                
                if indexOfGap % 2 != 0:
                    exonNum = (totelements + 2)/4 - (indexOfGap + 1)/2 + 1
                    if exonNum <= int(refIntronNum):
                        endreads += 1
            
            # Only the end of the read is within the TE. 
            else:   
                
                indexOfGap = 0
                c = 1
                for i,k in zip(elementsVec[0::2], elementsVec[1::2]):
                    if start_read >= i and start_read <= k:
                        indexOfGap = c
                        break
                    c = c + 1

                if indexOfGap % 2 != 0:
                    exonNum = (totelements+2)/4 - (indexOfGap + 1)/2 + 1
                    if exonNum >= int(refExonNum):
                        startreads += 1
                        
                         
print('\t'.join([str(totreads), str(startreads), str(endreads), 'se']))