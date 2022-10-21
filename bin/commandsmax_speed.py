#!/usr/bin/python
# programmer : nshah 
# usage: To be used in order to get read information from max files

import sys

distbuffer = 1000

# This is the rest of the command that will create a searchable tabix file for downstream analysis

bamfileloc = sys.argv[2]
fout=open('filterreadcommands.txt','w')

with open(sys.argv[1],"r") as filteredfile:
	for eachline in filteredfile:
		temp = eachline.strip().split("\t")
		number1 = temp[7]
		number2 = temp[20]
		elementsvec = temp[24]
		location = temp[6]
		chromosome = temp[26]
		startlocTE = temp[32]
		endlocTE = temp[33]
		filename = temp[44]
		uniqid = temp[55]
		strand = temp[37]

		if location == 'intron':
			elementvecint = [int(x) for x in elementsvec.split(',')]
			startloc = min(elementvecint) - distbuffer
			endloc = max(elementvecint) + distbuffer
			commandroot = 'samtools view -q 255 -h -f 2 -F 4 ' + bamfileloc + filename + '.bam ' + chromosome + ':' + str(startloc) + '-' + str(endloc) + " > " + uniqid + '--' + filename + ".sam" + ' ; samtools view -q 255 -f 2 -F 4 ' + bamfileloc + filename + '.bam ' + chromosome + ':' + str(startlocTE) + '-' + str(endlocTE) + ' | cut -f1  | sort | uniq > ' + uniqid + '--' + filename + 'IDs.txt ; cat <(samtools view -H ' + uniqid + '--' + filename + '.sam) <(grep -w -F -f ' + uniqid + '--' + filename + 'IDs.txt ' + uniqid + '--' + filename + '.sam)'
		else:
			strand = temp[37]
			elementvecint = [int(x) for x in elementsvec.split(',')]
			if strand == '+':
				startloc = int(temp[32]) - distbuffer
				endloc = max(elementvecint) + distbuffer
			else:
				startloc = min(elementvecint) - distbuffer
				endloc = int(temp[33]) + distbuffer
			commandroot = 'samtools view -q 255 -h -f 2 -F 4 ' + bamfileloc + filename + '.bam ' + chromosome + ':' + str(startloc) + '-' + str(endloc) + " > " + uniqid + '--' + filename + ".sam" + ' ; samtools view -q 255 -f 2 -F 4 ' + bamfileloc + filename + '.bam ' + chromosome + ':' + str(startlocTE) + '-' + str(endlocTE) + ' | cut -f1  | sort | uniq > ' + uniqid + '--' + filename + 'IDs.txt ; cat <(samtools view -H ' + uniqid + '--' + filename + '.sam) <(grep -w -F -f ' + uniqid + '--' + filename + 'IDs.txt ' + uniqid + '--' + filename + '.sam)'


		commandfinish = ' | samtools sort -n -m 1G - | bedtools bamtobed -bedpe -i stdin 2>/dev/null | ' + sys.path[0] + '/rmsk_annotate_bedpe_speed.py ' + chromosome + ',' + startlocTE + ',' + endlocTE + ' ' + elementsvec + ' ' + strand + ' ' + number1 + ' ' + number2 + ' > ./filterreadstats/' +  uniqid + '--' + filename + ".stats ; " + "rm " +  uniqid + '--' + filename + ".sam ; rm " +  uniqid + '--' + filename + 'IDs.txt'
	
		commandall = commandroot + commandfinish
		print >> fout, commandall

			




