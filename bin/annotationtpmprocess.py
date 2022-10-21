#! /usr/bin/env python
# programmer : nshah
# usage: To be used to process annotation file in order to get statistics

import sys
import os
import glob
import numpy as np

filename = sys.argv[1]
filename_all = filename.replace("filtered_","")
fout = open('{}_c'.format(sys.argv[1]),'w')
filteredcontent = [i.strip().split('\t') for i in open(filename).readlines()]
referencecontent = np.loadtxt(filename_all, dtype=object, delimiter='\t')
for line in filteredcontent:
	geneline = line[15]
	subset = referencecontent[np.where(referencecontent[:,15] == geneline)]
	covlist = subset[:,-3].astype(np.float)
	tpmlist = subset[:,-1].astype(np.float)
	stringout = '\t'.join(line + [str(max(covlist)), str(max(tpmlist)), str(sum(tpmlist))])
	print >> fout, stringout			
