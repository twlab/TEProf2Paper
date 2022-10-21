#! /usr/bin/env python
# programmer : nshah 
# usage: To be used to get summarized the transcript expression fraction and the TPM calculation into the pipeline

import sys
import os
import glob
import numpy as np

filename = sys.argv[1]
filteredcontent = [i.strip().split('\t') for i in open(filename).readlines()[1:]]
referencecontent = np.loadtxt(filename, dtype=object, delimiter='\t', skiprows=1)
foutfractot=open('{}_frac_tot'.format(sys.argv[1]),'w')
fouttpm=open('{}_tpm'.format(sys.argv[1]),'w')
sumFPKM=sum(referencecontent[:,-1].astype(np.float)) 
for line in filteredcontent:
	geneline = line[9]
	transcriptidline = line[5]
	fpkmline = line[-1]
	subset = referencecontent[np.where(referencecontent[:,9] == geneline)]
	fpkmlist = subset[:,-1].astype(np.float)
	
	#If there is no expression of gene then no division will be performed
	if (max(fpkmlist) == 0.0):
		print >> foutfractot, '\t'.join([transcriptidline] + [str(0.0)])
		print >> fouttpm, '\t'.join([transcriptidline] + [str(0.0)])
	
	#If at least one isoform is expressed then the fraction and tpm will be expressed
	else:
		print >> foutfractot, '\t'.join([transcriptidline] + [str(float(fpkmline)/sum(fpkmlist))])
		print >> fouttpm, '\t'.join([transcriptidline] + [str(float(fpkmline)*1000000.00/sumFPKM)])
	
	