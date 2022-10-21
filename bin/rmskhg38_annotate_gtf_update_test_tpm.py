#!/usr/bin/python
# programmer : nshah 
# usage: To be used in order to annotate a GTF file for presence of transcripts that begin in transposable elements.
# argument structure: gtfAnnotator.py <gtffile> <argumentfile.txt>

from __future__ import division #integer division has to be //, now all / are floating point
import cPickle as pickle
import os.path
import sys
import subprocess as sp
import re
import math
import copy
import tabix

# If the arguments.txt file is not specified then default the one in the directory as the transcript will be used
if len(sys.argv) < 3:
    argfile = sys.path[0] + '/arguments.txt'
else:
    argfile = sys.argv[2]
    
print argfile

# Import files needed to perform annotation. Check to make sure all the files do exist
argdic = {}
with open(argfile) as f:
    for line in f:
        (key, val) = line.split('\t')
        argdic[key] = val.strip()
  
# Check if all required arguments are present in the dictionary
#
# rmsk: tabix formatted bed6 files with repeatmasker or other file that user wants to check start location for
# gencodeplusdic: Dictionary of all gencode elements including introns and exons for the plus (+) strand
# gencodeminusdic: Dictionary of all gencode elements including introns and exons for the minus (-) strand

requiredkeys = ['rmsk', 'gencodeplusdic', 'gencodeminusdic']
if not (requiredkeys[0] in argdic.keys() and requiredkeys[1] in argdic.keys() and requiredkeys[2] in argdic.keys()):
    sys.exit('All required arguments (rmsk, gencodeplusdic, and gencodeminusdic) are not present')
else:
    for keyfile in requiredkeys:
        if not os.path.isfile(argdic[keyfile]):
            sys.exit('{} does not exist!'.format(argdic[keyfile]))
            
# Check if the 3 optional arguments are present and files exist
#
# focusgenes: The program has two outputs (1) on a focus set of genes (2) with all genes. This file lists the genes that the user wants to filter for originally
# plusintron: Tabix file of all the plus strand introns 
# minusintron: Tabix file of all the minus strand introns

optionalkeys = ['focusgenes', 'plusintron', 'minusintron']

filtergenes = False
annotateintrons = False
if optionalkeys[0] in argdic.keys():
    if os.path.isfile(argdic[optionalkeys[0]]):
        filtergenes = True
if optionalkeys[1] in argdic.keys() and optionalkeys[2] in argdic.keys():
    if os.path.isfile(argdic[optionalkeys[1]]) and os.path.isfile(argdic[optionalkeys[2]]):
        annotateintrons = True

# Add in information for the repeatmasker annotation (user can use any bed6 file)

rmsk = tabix.open(argdic['rmsk'])
input = sys.argv[1]
peaktotalsize = 0 #Will be added up from the bed/peak file

if filtergenes == True:
    oncogenelist = []
    # These are the oncogenes that I want to focus on. This can be excluded in future versions where a global analysis is desired.
    with open(argdic['focusgenes'],"r") as fin:
            for eachline in fin:
                    temp = eachline.strip()
                    oncogenelist = oncogenelist + [temp]
                
# Legacy code allowing for TE filtering within this script, but it was decided that it makes more sense to do this in downstream analysis
#
#TE={}
#TEsub={}
#TEdic={}
# Format: Subfamily     Class   Family
#with open("/bar/nshah/programs/Pyscript/TE.lst","r") as fin:
#       for eachline in fin:
#                temp = eachline.strip().split('\t')
#                TE[temp[0]]=temp[1]
#                TEsub[temp[0]]=temp[2]
#                TEdic[temp[2]]=temp[1]

# Function that annotates a location based on the Gencode dictionary made in genecode_to_dic.py
def annotate(start, DIC):
    startelement = "None"
    start_codon_coor = "None"
    stop_codon_coor = "None"
    transcriptstart = "None"
    transcriptend = "None"
    elementlist = "None"
    overlappercentage = "None"
    previous = "None"
    
    allannotation = []
    for transcript in DIC.keys():
        if start >= transcript[0] and start <= transcript[1]:
            elementlist = []
            for element in DIC[transcript].keys():
                if element != "start_codon" and element != "stop_codon":
                    elementlist.append(str(element[0]))
                    elementlist.append(str(element[1]))
                if element == "start_codon":
                    start_codon = DIC[transcript][element].split(",")
                    start_codon_coor = ";".join([start_codon[0], start_codon[3], start_codon[4]])
                elif element == "stop_codon":
                    stop_codon = DIC[transcript][element].split(",")
                    stop_codon_coor = ";".join([stop_codon[0], stop_codon[3], stop_codon[4]])
                elif start >= element[0] and start <= element[1]:
                    startelement = DIC[transcript][element]
            transcriptstart = str(transcript[0])
            transcriptend = str(transcript[1])
            elementlist = ",".join(elementlist)
            
            coordinateanno = ("^".join([startelement, start_codon_coor, stop_codon_coor, transcriptstart, transcriptend, elementlist, overlappercentage]))
            allannotation.append(coordinateanno)
    
    #If nothing overlaps a placeholder is needed for all the variables which will be "None"
    if elementlist == "None":          
        coordinateanno = ("^".join([startelement, start_codon_coor, stop_codon_coor, transcriptstart, transcriptend, elementlist, overlappercentage]))
        allannotation.append(coordinateanno)
   
    return [allannotation, start, "None"]

def annotateintron(chromosome, start, end, strand):

    if strand == "+":
        pIntrons = tabix.open(argdic['plusintron'])
    else:
        pIntrons = tabix.open(argdic['minusintron'])
    res = []
    tmp = ''
    try:
        tmp=pIntrons.query(chromosome,start,end)
    except:
        tmp = ''
    
    for i in tmp:
        res.append(i)

    if res:
        introns = res
        dels = []
        dist = []
        intronnums = []
        #Remove those not in the correct range
        for i in introns:
            i_start = int(i[1])
            i_end = int(i[2])
            
            if start-1 <= i_start and end+1 >= i_end:
                dist.append(i_end-i_start)
            else:
                dist.append(0)
            
            # This allows you to extract the intron number from the description
            intronnums.append(int(i[3].split("; ")[8].split(" ")[1]))
                
        
        maximumdist = max(dist)
        
        maxindices = [i for i,x in enumerate(dist) if x == maximumdist ]
        
        intronmin = 99999
        for maxindex in maxindices:  
                if intronnums[maxindex] < intronmin:
                        indexi = maxindex
                        intronmin = intronnums[maxindex]
        
        intronmax = introns[indexi]
        i_startmax = int(intronmax[1])
        i_endmax = int(intronmax[2])
        
        #Has to have a maximum distance location that is at least 90% of the intron
        if maximumdist != 0 and abs(float(i_endmax-i_startmax)/float(end-start))>.9:
            indexi = dist.index(maximumdist)
            return ",".join(intronmax)
        else:
            return "None"
    else:
        return "None"

# Simple function to calculate the overlap in Bp of two sequences
def getOverlap(a ,b):
    return max(0, min(a[1] + 1, b[1] + 1) - max(a[0], b[0]))

# Function that annotates a location based on the Gencode dictionary made in genecode_to_dic.py
# previous is a variable that holds the transcript ID of a previous annotation (if a previous annotation has been made)
# If there has been a previous annotation, then that transcript ID will be searched for first before anything else.
def annotateregion(start, end, DIC):
    
    #print 'regionannotate'
    startelement = "None"
    start_codon_coor = "None"
    stop_codon_coor = "None"
    transcriptstart = "None"
    transcriptend = "None"
    elementlist = "None"
    overlappercentage = "None"
    
    allannotation = []
    
    for transcript in DIC.keys():
        overlapbw = getOverlap([start, end], [transcript[0], transcript[1]])
        if overlapbw > 0:
            elementlist = []
            startelements = []
            #print(DIC[transcript].keys())
            for element in DIC[transcript].keys():
                if element == "start_codon":
                    start_codon = DIC[transcript][element].split(",")
                    start_codon_coor = ";".join([start_codon[0], start_codon[3], start_codon[4]])
                elif element == "stop_codon":
                    stop_codon = DIC[transcript][element].split(",")
                    stop_codon_coor = ";".join([stop_codon[0], stop_codon[3], stop_codon[4]])
                else:
                    elementlist.append(str(element[0]))
                    elementlist.append(str(element[1]))
                    overlapbw2 = getOverlap([start, end], [element[0], element[1]])
                    if overlapbw2 > 0:
                        startelement = DIC[transcript][element]
                        startelements.append(startelement)
            transcriptstart = str(transcript[0])
            transcriptend = str(transcript[1])
            elementlistprint = ",".join(elementlist)
            if startelements:
                for startelement in startelements:
                    coordinateanno = ("^".join([startelement, start_codon_coor, stop_codon_coor, transcriptstart, transcriptend, elementlistprint, overlappercentage]))
                    allannotation.append(coordinateanno)
            else:
                coordinateanno = ("^".join([startelement, start_codon_coor, stop_codon_coor, transcriptstart, transcriptend, elementlistprint, overlappercentage]))
                allannotation.append(coordinateanno)
    
    if elementlist == "None":
        coordinateanno = ("^".join([startelement, start_codon_coor, stop_codon_coor, transcriptstart, transcriptend, elementlist, overlappercentage]))
        allannotation.append(coordinateanno)
        
    # Returns a long string that will be processed with the processannotation function             
        
    
    return [allannotation, start, end]  
          
# Example Annotation:
# chr10,HAVANA,exon,114154676,114154860,.,+,.,gene_id "ENSG00000197142.6"; transcript_id "ENST00000354655.4"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "ACSL5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "ACSL5-002";
# exon_number 2;  exon_id "ENSE00003492611.1";  level 2; tag "alternative_5_UTR"; tag "basic"; tag "appris_principal"; tag "CCDS"; ccdsid "CCDS7573.1"; havana_gene "OTTHUMG00000019060.1"; havana_transcript "OTTHUMT00000050387.1";, transcript_id "ENST00000354655.4",chr10;114154705;114154707,None, .67
# This function takes the string annotation and processes it for the useful information

# It will also decide how to annotate the element based on the following principles:
# 1. If this is an annotation of a downstream exon and an element is on teh same trasncript as the 5' annotation, then those elements will be accepted
# 2. Multiple exons or introns, the most 5' is selected (opposite for minus strand transcripts)
# 3. Exons and Introns 

def processannotation(returnanno, plusminus, previous):
    #print(returnanno)
    astringvec = returnanno[0]
    startexon = returnanno[1]
    endexon = returnanno[2]
    if plusminus == "minus":
        startexon = returnanno[2]
        endexon = returnanno[1]
    
    genetype = []
    genename = []
    chromosome = []
    start = []
    end = []
    exonintron = []
    exonintronnums = []
    startcodon = []
    transcriptstart = []
    transcriptend = []
    elementlist = []
    transcriptid = []
    
    for astring in astringvec:
        splitvec = astring.split('^')
        splitog = splitvec[0]
        if splitog == "None":
            return ["None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None"]
        
        split1 = splitog.split(',')
        chromosome.append(split1[0])
        start.append(int(split1[3]))
        end.append(int(split1[4]))
        exonintron.append(split1[2])
        exonintronstat = split1[2]
        
        splitanno = split1[8].split("\"; ")
        
        tempdesc = split1[8].split("transcript_id \"")[1]
        transcriptid.append(tempdesc.split("\";")[0])
        
        tempdesc = split1[8].split("gene_type \"")[1]
        genetype.append(tempdesc.split("\";")[0])
        
        tempdesc = split1[8].split("gene_name \"")[1]
        genename.append(tempdesc.split("\";")[0])
        
        if exonintronstat == "intron":
            tempdesc = split1[8].split("intron_number ")[1]
            exonintronnum = int(tempdesc.split(";")[0])
        else:
            tempdesc = split1[8].split("exon_number ")[1]
            exonintronnum = int(tempdesc.split(";")[0])
        
        startcodon.append(splitvec[1]) #Check if this is correct
        exonnumber = int(split1[-1])
        
        # I have to change the intron number definitions for the minus strand
        if plusminus == "minus":
            if exonintronstat == "intron":
                 exonintronnum = exonnumber - exonintronnum
        
        exonintronnums.append(str(exonintronnum))
        
        transcriptstart.append(splitvec[3])
        transcriptend.append(splitvec[4])
        elementlist.append(splitvec[5])
    
    #If only one element overlaps, then that will be the annotation
    if len(chromosome) == 1:
        i_r = 0
    else:
        #If the user gives a transcript ID to the function, then only features that are from that transcript ID
        #are then queried to see if an annotation exists. If there are no features from that transcript ID, then
        #it is annotated by features outside of that transcript ID. 
        if previous in transcriptid:
            transindices = [i for i,x in enumerate(transcriptid) if x == previous]
            
            if plusminus == "minus":
                mostupstream = 0
                for transindex in transindices:  
                    if int(start[transindex]) > mostupstream:
                        i_r = transindex
                        mostupstream = int(start[transindex])
            else:
                mostupstream = 9999999999999999
                for transindex in transindices:  
                    if int(start[transindex]) < mostupstream:
                        i_r = transindex
                        mostupstream = int(start[transindex])
        # The most complex condition is when there is no guiding transcript id to decide which annotation to use
        # In this case, it first checks if the annotation is for a single location or for a region.
        #
        # If it is for a single location (endexon="None"), then it will just take all the features and then selects
        # the exon or intron that is the most upstream. (exon given priority to intron)
        #
        # If it is a region, then only elements that include the start of the exon being annotated are included
        # Then the exon or intron that has the lowest number is selected. 
        else:
            if endexon == "None":
                if "exon" in exonintron:
                    exonindices = [i for i,x in enumerate(exonintron) if x == "exon"]
                    if len(exonindices) == 1:
                        i_r = exonindices[0]
                    elif exonindices:
                        minexonnumber = min([ exonintronnums[i] for i in exonindices])
                        for exonindex in exonindices:
                            if exonintronnums[exonindex] == minexonnumber:
                                i_r = exonindex     
                                break
                              
                    
                elif "intron" in exonintron:
                    intronindices = [i for i,x in enumerate(exonintron) if x == "intron"]
                    if len(intronindices) == 1:
                        i_r = intronindices[0]
                    elif intronindices:
                        minintronnumber = min([ exonintronnums[i] for i in intronindices])
                        for intronindex in intronindices:
                            if exonintronnums[intronindex] == minintronnumber:
                                i_r = intronindex     
                                break
            else:
                #Gets the elements that the start index is actually in. Then these elements are 
                startindicesbw = []
                startindicesbw = [i for i,x in enumerate(zip(start,end)) if x[0] <= startexon and x[1] >= startexon]
                #print(startindicesbw)
                   
                if startindicesbw:
                    if "exon" in [exonintron[b] for b in startindicesbw]:
                        exonindices = [i for i,x in enumerate(exonintron) if x == "exon" and i in startindicesbw]
                        if len(exonindices) == 1:
                            i_r = exonindices[0]
                        elif exonindices:
                            minexonnumber = min([ exonintronnums[i] for i in exonindices])
                            for exonindex in exonindices:
                                if exonintronnums[exonindex] == minexonnumber:
                                    i_r = exonindex     
                                    break
                              
                    
                    elif "intron" in [exonintron[b] for b in startindicesbw]:
                        intronindices = [i for i,x in enumerate(exonintron) if x == "intron" and i in startindicesbw]
                        if len(intronindices) == 1:
                            i_r = intronindices[0]
                        elif intronindices:
                            minintronnumber = min([ exonintronnums[i] for i in intronindices])
                            for intronindex in intronindices:
                                if exonintronnums[intronindex] == minintronnumber:
                                    i_r = intronindex     
                                    break
                else:
                    if "exon" in exonintron:
                        exonindices = [i for i,x in enumerate(exonintron) if x == "exon"]
                        if len(exonindices) == 1:
                            i_r = exonindices[0]
                        elif exonindices:
                            minexonnumber = min([ exonintronnums[i] for i in exonindices])
                            for exonindex in exonindices:
                                if exonintronnums[exonindex] == minexonnumber:
                                    i_r = exonindex     
                                    break
                              
                    
                    elif "intron" in exonintron:
                        intronindices = [i for i,x in enumerate(exonintron) if x == "intron"]
                        if len(intronindices) == 1:
                            i_r = intronindices[0]
                        elif intronindices:
                            minintronnumber = min([ exonintronnums[i] for i in intronindices])
                            for intronindex in intronindices:
                                if exonintronnums[intronindex] == minintronnumber:
                                    i_r = intronindex     
                                    break
                                
            
                
    #print(genetype, genename, chromosome, start, end, exonintron, exonintronnums, startcodon, transcriptstart, transcriptend, elementlist, transcriptid)
    return [genetype[i_r], genename[i_r], chromosome[i_r], str(start[i_r]), str(end[i_r]), exonintron[i_r], exonintronnums[i_r], startcodon[i_r], transcriptstart[i_r], transcriptend[i_r], elementlist[i_r], transcriptid[i_r]]

# Array that will store the location of all the lines with a transcript entry. This will let me know the location of the
# transcripts and exons for subsequent analysis
transcriptlines = []

i = 0
with open(sys.argv[1],"r") as gtf:
    next(gtf) # For the two hashtags that are in the file. Should make this a loop where it will parse hashtags till it gets to the actual file itself
    next(gtf) # For the two hashtags that are in the file
    
    # The file name of the peak file has information on the chromosome, start, end, and strand of the peak
    for eachline in gtf:
        #print(eachline)
        
        # GTF File format:
        # chr10     StringTie       transcript      1033329 1034198 1000    .       .       gene_id "STRG.1"; transcript_id "STRG.1.1"; cov "3.585057"; FPKM "1149425.125000"; TPM "500943.937500";
        temp = eachline.strip().split("\t")
        type = temp[2]
        if type == "transcript":
            transcriptlines.append(i)
        i = i + 1
gtf.close()

endoffile = i
basefile = sys.argv[1].split(".")[0]
num_transcripts = len(transcriptlines)

# Load the entire GTF file into an array. This is not the most memory efficient, but
lines = []
with open(sys.argv[1],"r") as gtf:
    lines = gtf.read().splitlines()

num_lines = len(lines)
    
# Load up the genecode dictionaries
plus_dic = {}
with open(argdic['gencodeplusdic'],'r') as DIC:
    plus_dic = pickle.load(DIC)

DIC.close()

minus_dic = {}
with open(argdic['gencodeminusdic'],'r') as DIC:
    minus_dic = pickle.load(DIC)

DIC.close()
    
# Remove the header lines
lines.pop(0)
lines.pop(0)

if filtergenes == True:
    fout=open('{}_annotated_test'.format(sys.argv[1]),'w')
    fout1=open('{}_annotated_filtered_test'.format(sys.argv[1]),'w') 

fout_a=open('{}_annotated_test_all'.format(sys.argv[1]),'w')
fout1_a=open('{}_annotated_filtered_test_all'.format(sys.argv[1]),'w') 

j = 0
# Iterate through the GTF lines and annotate each transcript based on repeatmasker and gencode. Only print those out that begin within a transposable element
for indext in transcriptlines:
    
    temp = lines[indext].strip().split("\t")
    
    # Obtain pertinent transcript information
    chr_trans = temp[0]
    
    
    # Legacy for human. 
    #
    # allowedchr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    
    start_trans = int(temp[3])
    end_trans = int(temp[4])
    
    # Obtain coverage for the transcript
    description = temp[8].split("cov \"")[1]
    coverage = float(description.split("\";")[0])
    
    # Obtain TPM for the transcript
    description = temp[8].split("TPM \"")[1]
    tpmtranscript = float(description.split("\";")[0])
    
    # Obtain the transcriptID
    transcriptinfo = temp[8].split("\"; ")[1]
    transcriptid = transcriptinfo.split(" \"")[1]
    
    # Variable that will let the program know whether to perform the expensive operation of genecode annotation
    skip = False
    
    # If there is no strand, it is not possible to have an accurate annotation
    strand = temp[6]
    
    if strand == ".":
        skip = True
        subfam = []
    
    # If the chromosome is not in the gencode dictionary then an accurate annotation is also not possible
    if chr_trans not in plus_dic.keys():
        skip = True
        subfam = []
    
    if strand == "+":
        res = []
        tmp = ''
        
        try:
                tmp=rmsk.query(chr_trans,start_trans-1,start_trans+1)
        except:
                tmp = ''
        
        for i in tmp:
            res.append(i)
        if res:
            subfam = res
        else:
            subfam = []
            
    elif strand == "-":
        res = []
        tmp = ''
        
        try:
                tmp=rmsk.query(chr_trans,end_trans-1,end_trans+1)
        except:
                tmp = ''
        
        for i in tmp:
            res.append(i)
        if res:
            subfam = res
        else:
            subfam = []
            
    
    # There should only be 1 Subfamily that is returned. If none is returned, then the transcript is not fruther annotate
    if subfam == []:
        transcriptTE = ["None", "None", "None", "None", "None", "None"]
    else:    
        transcriptTE = subfam[0]
        subfamTE = transcriptTE[3]
        
        # I want to specifically look at TE for now, but maybe others later.
        #if subfamTE not in TE.keys():
        #    transcriptTE = ["None", "None", "None", "None", "None", "None"]
             
    # Figure out where the transcript ends    
    if j == num_transcripts - 1:
        indexend = num_lines
    else:
        indexend = transcriptlines[j+1]
    
    if skip != True:
        
        # Extract the the current transcript and annotate it
        cluster = lines[indext:indexend]
        
        if "+" == strand:
            
            # Gencode Dictionary
            plus = plus_dic
            
            # Arrays to store information on the exons AFTER the first one
            exonstarts = []
            exonends = []
            exonannotations = []
            exontypes = []

            i = 0
            for eachline in cluster:
                temp=eachline.strip().split("\t")
                chromosome = temp[0]
                start = int(temp[3])
                end = int(temp[4])
                type = temp[2]
                if type == "transcript":
                    transcriptstart = start
                    transcriptend = end
                else:
                    if i == 0:
                        description = temp[8].split("cov \"")[1]
                        exon1coverage = float(description.split("\";")[0])
                        exon1start = int(temp[3])
                        exon1end = int(temp[4])
                    exonstarts.append(start)
                    exonends.append(end)
                    i = i + 1
            exonstarts.sort()
            exonends.sort()
            tstartannotation=processannotation(annotateregion(transcriptstart, exon1end, plus[chromosome]),"plus", "None")
            splicing = "Yes"
            if i == 1:  # When there is only one item in this array, there is only one exon and no splicing
                splicing = "No"
                tendannotation = processannotation(annotate(transcriptend, plus[chromosome]),"plus", tstartannotation[11])
                firstintronanno = "None"
            elif i == 2: # This means there is only 1 splice site, so that is the one that should be annotated
                tendannotation = processannotation(annotateregion(exonstarts[1], exonends[1], plus[chromosome]),"plus", tstartannotation[11])
                
                if annotateintrons == True:
                    firstintronanno = annotateintron(chromosome, exon1end, exonstarts[1], strand)
                else:
                    firstintronanno = "None"
            else: # This is when there are multiple splice acceptor sites. It will keep iterating through them until it find one that annotated into a transcript
                  # If none of them do, then it will return a list of "None"
                tendannotation = ["None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None"]
                k = 0
                for exon in exonstarts:
                    if k > 0:
                        exonannotation = (processannotation(annotateregion(exon, exonends[k], plus[chromosome]),"plus", tstartannotation[11]))
                        exonannotations.append(exonannotation)
                        exontypes.append(exonannotation[5])                
                    k = k + 1
                
                # I want to find if there are any exons, if not then introns. Finally if its not within any transcript than the none annotation is correct
                if "exon" in exontypes:
                    tendannotation = exonannotations[exontypes.index("exon")]
                elif "intron" in exontypes:
                    tendannotation = exonannotations[exontypes.index("intron")]
                
                if annotateintrons == True:
                    firstintronanno = annotateintron(chromosome, exon1end, exonstarts[1], strand)
                else:
                    firstintronanno = "None"
                
            exonstarts.sort()
            exonends.sort()
            x = 0
            genomiclocations = []
            genomiclocations.append(strand)
            genomiclocations.append(chromosome)
            for exon in exonstarts:
                genomiclocations.append(str(exon))
                genomiclocations.append(str(exonends[x]))
                x = x + 1
            locationreturn = ",".join(genomiclocations)
            annotationreturn = [transcriptid] + tstartannotation + [splicing] + tendannotation + [chromosome] + [str(transcriptstart)] + [str(transcriptend)] + [str(locationreturn)] + [firstintronanno] + transcriptTE + [strand, str(coverage), str(exon1coverage), str(tpmtranscript)]
            stringreturn = "\t".join(annotationreturn)
            
            print >> fout_a, stringreturn
            if splicing == "Yes":
                if transcriptTE[0] != "None":
                    print >> fout1_a, stringreturn
            
            if filtergenes == True:
                if tendannotation[1] in oncogenelist:
                    print >> fout, stringreturn
                    if splicing == "Yes":
                        if transcriptTE[0] != "None":
                            print >> fout1, stringreturn

          
        elif "-" == strand:
            # Gencode Dictionary
            plus = minus_dic
            
            exonstarts = []
            exonends = []
            exonannotations = []
            exontypes = []
            
            i = 0
            for eachline in cluster:
                temp=eachline.strip().split("\t")
                chromosome = temp[0]
                end = int(temp[3])
                start = int(temp[4])
                type = temp[2]
                if type == "transcript":
                    transcriptstart = start
                    transcriptend = end
                else:
                    if i == (len(cluster)-2):
                        description = temp[8].split("cov \"")[1]
                        exon1coverage = float(description.split("\";")[0])
                        exon1start = int(temp[3])
                        exon1end = int(temp[4])
                    exonstarts.append(start)
                    exonends.append(end)
                    i = i + 1
            tstartannotation=processannotation(annotateregion(exon1start, transcriptstart, plus[chromosome]),"minus", "None")
            splicing = "Yes"
            exonstarts.sort(reverse=True)
            exonends.sort(reverse=True)
            if i == 1: # When there is only one item in this array, there is only one exon and no splicing
                splicing = "No"
                tendannotation = processannotation(annotate(transcriptend, plus[chromosome]),"minus", tstartannotation[11])
                firstintronanno = "None"
            elif i == 2: # This means there is only 1 splice site, so that is the one that should be annotated
                tendannotation = processannotation(annotateregion(exonends[1], exonstarts[1], plus[chromosome]),"minus", tstartannotation[11])
                
                if annotateintrons == True:
                    firstintronanno = annotateintron(chromosome, exonstarts[1], exon1start, strand)
                else:
                    firstintronanno = "None"
            else: # This is when there are multiple splice acceptor sites. It will keep iterating through them until it find one that annotated inot a transcript
                  # If none of them do, then it will return a list of "None"
                tendannotation = ["None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None"]
                k = 0
                for exon in exonstarts:
                    if k > 0:  
                        exonannotation = (processannotation(annotateregion(exonends[k], exon, plus[chromosome]),"minus", tstartannotation[11]))
                        exonannotations.append(exonannotation)
                        exontypes.append(exonannotation[5])
                    k = k + 1
                    
                # I want to find if there are any exons, if not then introns. Finally if its not within any transcript than the none annotation is correct
                if "exon" in exontypes:
                    tendannotation = exonannotations[exontypes.index("exon")]
                elif "intron" in exontypes:
                    tendannotation = exonannotations[exontypes.index("intron")]
                
                if annotateintrons == True:
                    firstintronanno = annotateintron(chromosome, exonstarts[1], exon1start, strand)
                else:
                    firstintronanno = "None"    
                
            exonstarts.sort()
            exonends.sort()
            x = 0
            genomiclocations = []
            genomiclocations.append(strand)
            genomiclocations.append(chromosome)
            for exon in exonstarts:
                genomiclocations.append(str(exon))
                genomiclocations.append(str(exonends[x]))
                x = x + 1
            locationreturn = ",".join(genomiclocations)
            annotationreturn = [transcriptid] + tstartannotation + [splicing] + tendannotation + [chromosome] + [str(transcriptstart)] + [str(transcriptend)] + [locationreturn] + [firstintronanno] + transcriptTE + [strand, str(coverage), str(exon1coverage), str(tpmtranscript)]
            stringreturn = "\t".join(annotationreturn)
            
            print >> fout_a, stringreturn
            if splicing == "Yes":
                if transcriptTE[0] != "None":
                    print >> fout1_a, stringreturn
            
            if filtergenes == True:
                if tendannotation[1] in oncogenelist:
                    print >> fout, stringreturn
                    if splicing == "Yes":
                        if transcriptTE[0] != "None":
                            print >> fout1, stringreturn                             
        
    j = j + 1
    #print "{} out of {}".format(str(j), str(num_transcripts))
    
                
                
                
