#!/usr/bin/python
# programmer : nshah using code from dli
# usage: To be used in order to create a dictionary of genecode transcripts, exons, and introns. This can be quickly accessed in future programs
import sys
import re
import random
import string
import time
import cPickle as pickle

# These are the dictionaries that will hold the transcripts for the plus strand and the minus strand for Gencode
plus={}
minus={}

# The argument should be a GTF file from gencode
with open(sys.argv[1]) as genecode:
# First generate a dictionary for every transcript and its corresonding exons
    for each in genecode:

        temp=each.strip().split('\t')           
        #How a genecode entry normally looks:
        #chr1   	HAVANA 	gene   	11869  	14412  	.      	+      	.      	gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";
        #I split it by tabs into an array to get the information
        #temp[0] chromosome
        #temp[1] ENSEMBL/HAVANA
        #temp[2] transcript/gene/start_codon/etc
        #temp[3] start
        #temp[4] end
        #temp[6] strand
        #temp[-1] transcript ID
        
        tempstr= temp[-1]
                
                # Only those transcripts tht are identified as the principle ones and are KNOWN (Might change this)
        if (('appris_principal' in temp[-2])):
                
                # Needs to reset the transcript label every time a new transcript is reached
            if temp[2]=='transcript':
                current_temp = temp
                current_tempstr = tempstr
                if temp[6]=='+':
                    
                    # When the chromosome is not in the dictionary, I make a new chromosome key with a empty dictionary as the result 
                    if temp[0] not in plus:
                        plus[temp[0]]={}
                    
                    # Within the dictionary of the current chromosome, I make a new transcript key (start, end, ID) with a empty dictionary as the result 
                    plus[temp[0]][(int(temp[3]),int(temp[4]),tempstr)]={}
                else:
                    if temp[0] not in minus:
                        minus[temp[0]]={}
                        
                    minus[temp[0]][(int(temp[3]),int(temp[4]),tempstr)]={}


            elif temp[2]=='exon':
                    
                # When there is an exon, I make a new key in the current transcript's dictionary with the (start, end) as the key, and then the result is a
                # comma-delimitted list of all the metadata for that exon in gencode
                if temp[6]=='+':
                    plus[temp[0]][(int(current_temp[3]),int(current_temp[4]),current_tempstr)][(int(temp[3]),int(temp[4]))]=",".join(temp).replace("\n","")
                else:
                    minus[temp[0]][(int(current_temp[3]),int(current_temp[4]),current_tempstr)][(int(temp[3]),int(temp[4]))]=",".join(temp).replace("\n","")
                    
            elif temp[2]=='start_codon':
                
                # When there is an start codon, I make a new key in the current transcript's dictionary with "start_codon" as the key, and then the result is a
                # comma-delimitted list of all the metadata for that start codon in gencode
                if temp[6]=='+':
                    plus[temp[0]][(int(current_temp[3]),int(current_temp[4]),current_tempstr)]['start_codon']=",".join(temp).replace("\n","")
                else:
                    minus[temp[0]][(int(current_temp[3]),int(current_temp[4]),current_tempstr)]['start_codon']=",".join(temp).replace("\n","")
                            
            elif temp[2]=='stop_codon':
                
                # Does not actually work. I will fix this later if the stop codon turns out to be important for us to know in downstream analysis
                if temp[6]=='+':
                    plus[temp[0]][(int(current_temp[3]),int(current_temp[4]),current_tempstr)]['stop_codon']=",".join(temp).replace("\n","")
                else:
                    minus[temp[0]][(int(current_temp[3]),int(current_temp[4]),current_tempstr)]['stop_codon']=",".join(temp).replace("\n","")
genecode.close()

# The problem with gencode is that it does not give any entries for the introns. instead of having to have
# downstream transcripts take care of determing the intron locations and their characteristics, I add in entries
# for the regions betweent he exons as introns and number them


# Plus strand intron insertion
for chrom in plus.keys():
    for transcript in plus[chrom].keys():
        exons = plus[chrom][transcript].keys()
            
            # i only want exon entries
        if 'start_codon' in exons:
            exons.remove('start_codon')
        if 'stop_codon' in exons:
            exons.remove('start_codon')
                    
        n_exons = len(exons)
        
        # I add one more entry into the metadata which is the number of exons for each exon
        for  i in range(0,n_exons):
            label = plus[chrom][transcript][exons[i]].split(",")
            label = label[:]
            label.append(str(n_exons))
            labelstr = ",".join(label)

            plus[chrom][transcript][exons[i]] = labelstr.replace("\n","")
        
        
        if n_exons > 1:
                
            # Sort the exons by the start site
            exons.sort(key=lambda tup: int(tup[0]))
            
            for  i in range(0,n_exons-1):
                start = exons[i][1]+1
                end = exons[i+1][0]-1
                intron_tup = (start, end)
                label = plus[chrom][transcript][exons[i]].split(",")
                label = label[:]
                label[2] = 'intron'
                label[3] = str(start)
                label[4] = str(end)

                annotation = label[-3].split("; ")
                annotation[8] = 'intron_number {}'.format(i+1) 
                annotationintron = "; ".join(annotation)
                label[-3] = annotationintron
                labelstr = ",".join(label)
    
                                # Note the exon ID will still be of the exon before the intron. I am using the metadata of the previous exon
                                # but just adding in the intron number
                plus[chrom][transcript][intron_tup] = labelstr.replace("\n","")

for chrom in minus.keys():
    for transcript in minus[chrom].keys():
        exons = minus[chrom][transcript].keys()
        
        if 'start_codon' in exons:
            exons.remove('start_codon')
        if 'stop_codon' in exons:
            exons.remove('start_codon')
        
        n_exons = len(exons)
        
                        # Need to add number of exons
        for  i in range(0,n_exons):
            label = minus[chrom][transcript][exons[i]].split(",")
            label = label[:]
            label.append(str(n_exons))
            labelstr = ",".join(label)
    
            minus[chrom][transcript][exons[i]] = labelstr.replace("\n","")
        
        if n_exons > 1:
                
            # Sort the exons by the start site
            exons.sort(key=lambda tup: int(tup[0]))
            
            for  i in range(0,n_exons-1):
                start = exons[i][1]+1
                end = exons[i+1][0]-1
                intron_tup = (start, end)
                label = minus[chrom][transcript][exons[i]].split(",")
                
                label = label[:]
                label[2] = 'intron'
                label[3] = str(start)
                label[4] = str(end)
    
                annotation = label[-3].split("; ")
                annotation[8] = 'intron_number {}'.format(i+1) 
                annotationintron = "; ".join(annotation)
                label[-3] = annotationintron
                labelstr = ",".join(label)
    
                # Note the exon ID will still be of the exon before the intron. I am using the metadata of the previous exon
                # but just adding in the intron number
                minus[chrom][transcript][intron_tup] = labelstr.replace("\n","")


with open("genecode_plus.dic",'wb') as gp:
    pickle.dump(plus,gp)

with open("genecode_minus.dic",'wb') as gm:
    pickle.dump(minus,gm)

gp.close()
gm.close()




