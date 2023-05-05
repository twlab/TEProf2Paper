# **TEProF2**  

(**T**ransposable **E**lement **Pro**moter **F**inder **2**) 

TEProf2 is a software pipeline optimized to find TE-gene fusion transcripts from short-read RNA-sequencing data. Upon finding these transcripts, the program will find which ones are enriched in the treatment samples.

In addition, the program can be run in a reference-guided approach based on previous runs or a user's custom transcript set that includes TE-gene fusion genes. 

There are additional modules to facilitate transcript translation.

Please cite the following: https://doi.org/10.1038/s41588-023-01349-3

# Outline

* [<strong>TEProF2</strong>](#teprof2)
* [Outline](#outline)
* [Summary](#summary)
* [Requirements](#requirements)
   * [(1) Software](#1-software)
      * [Conda](#conda)
      * [Add bin folder to PATH](#add-bin-folder-to-path)
   * [(2) Reference Files](#2-reference-files)
      * [(A) Defaults](#a-defaults)
      * [(B) Gencode Dictionary](#b-gencode-dictionary)
      * [(C) Repeatmasker Files](#c-repeatmasker-files)
      * [(D) Repeatmasker Subfamily, Family, Class Mapping](#d-repeatmasker-subfamily-family-class-mapping)
      * [(E) Intron Annotations (Optional)](#e-intron-annotations-optional)
      * [(F) Gene Filter List (Optional)](#f-gene-filter-list-optional)
   * [(3) Input Files](#3-input-files)
* [Usage](#usage)
   * [Quick Run (Not recommended unless using Wang Lab Servers)](#quick-run-not-recommended-unless-using-wang-lab-servers)
   * [(1) Setup arguments.txt](#1-setup-argumentstxt)
      * [Required Arguments](#required-arguments)
      * [Optional Arguments](#optional-arguments)
   * [(2) Run annotation on each GTF file](#2-run-annotation-on-each-gtf-file)
   * [(3) Process annotations files to get rough estimate of relative expression of transcript versus all other transcripts of gene](#3-process-annotations-files-to-get-rough-estimate-of-relative-expression-of-transcript-versus-all-other-transcripts-of-gene)
   * [(4) Aggregate annotation samples across samples](#4-aggregate-annotation-samples-across-samples)
   * [5. Calculate Read Information](#5-calculate-read-information)
      * [(A) Make the directory where the stats will be stored](#a-make-the-directory-where-the-stats-will-be-stored)
      * [(B) Create commands to calculate read information](#b-create-commands-to-calculate-read-information)
      * [(C) Run the commands](#c-run-the-commands)
      * [(D) Combine all the read information files](#d-combine-all-the-read-information-files)
   * [6. Filter Candidates based on read information](#6-filter-candidates-based-on-read-information)
   * [7. Merge with Reference GTF](#7-merge-with-reference-gtf)
   * [8. Annotate Merged GTF](#8-annotate-merged-gtf)
   * [9. Calculate Transcript-Level Expression](#9-calculate-transcript-level-expression)
   * [10. Process and map expression output](#10-process-and-map-expression-output)
      * [(A) Obtain annotations, filter by candidates identified previously, and identify major splicing intron](#a-obtain-annotations-filter-by-candidates-identified-previously-and-identify-major-splicing-intron)
      * [(B) Process stringtie transcript annotation files to get relevant information and aggregate](#b-process-stringtie-transcript-annotation-files-to-get-relevant-information-and-aggregate)
   * [11. Quantification processing, sample identification, and final table creation](#11-quantification-processing-sample-identification-and-final-table-creation)
   * [12. Translating the transcripts with Kozak method and generating FASTA of RNa sequences](#12-translating-the-transcripts-with-kozak-method-and-generating-fasta-of-rna-sequences)
   * [13. CPC2 Translation of the transcripts](#13-cpc2-translation-of-the-transcripts)
   * [14. Finishing translation, classifying antigens, and outputting final list of potential antigens](#14-finishing-translation-classifying-antigens-and-outputting-final-list-of-potential-antigens)
   * [15. Optional further analysis with ballgown](#15-optional-further-analysis-with-ballgown)
* [Usage Reference Guided](#usage-reference-guided)
   * [1. Obtain reference gtf and gff3 file for analysis](#1-obtain-reference-gtf-and-gff3-file-for-analysis)
   * [2. Annotate the gff3 file, or obtain annotations from a previous run](#2-annotate-the-gff3-file-or-obtain-annotations-from-a-previous-run)
   * [3a. Calculate Transcript-Level Expression](#3a-calculate-transcript-level-expression)
   * [3b. Process and map expression output](#3b-process-and-map-expression-output)
      * [(A) Obtain annotations, filter by candidates identified previously, and identify major splicing intron](#a-obtain-annotations-filter-by-candidates-identified-previously-and-identify-major-splicing-intron-1)
      * [(B) Process stringtie transcript annotation files to get relevant information and aggregate](#b-process-stringtie-transcript-annotation-files-to-get-relevant-information-and-aggregate-1)
   * [4. Quantification processing, sample identification, and final table creation](#4-quantification-processing-sample-identification-and-final-table-creation)
   * [5. Translating the transcripts with Kozak method and generating FASTA of RNa sequences](#5-translating-the-transcripts-with-kozak-method-and-generating-fasta-of-rna-sequences)
   * [6. CPC2 Translation of the transcripts](#6-cpc2-translation-of-the-transcripts)
   * [7. Finishing translation, classifying antigens, and outputting final list of potential antigens](#7-finishing-translation-classifying-antigens-and-outputting-final-list-of-potential-antigens)
   * [8. Optional further analysis with ballgown](#8-optional-further-analysis-with-ballgown)
* [Questions?](#questions)

# Summary

TEProF2 v0.1

This pipeline takes assembled RNA-sequencing data (.gtf format) in human or mouse, and then will be able to assemble the data into transcripts, predict the transcripts starting from transposable elements, and calculate expression of the transcripts in comparison to the whole. 

Please cite: x

# Requirements

## (1) Software

stringtie >= 1.3.3

samtools >= 1.3.1

cufflinks >= 2.2.1

python 2.7 (cPickle, pytabix 0.1)

R >= 3.4.1 (ggplot2, bsgenome.hsapiens.ucsc.hg38 (or genome of your choosing), Xmisc, reshape2)

^The above can be installed individually. The pipeline has only been tested for the versions listed and may break with newer versions. 

### Conda

Alternatively, this can be installed with conda. For more information: https://docs.anaconda.com/anaconda/install/

You will create an environment named `teprof2` with the .yaml file provided with the package

```
conda env create -f TEProf2.yml
```

Activate the environment

```
conda activate teprof2
```


Then you will need to install the R package `Xmisc` manually

```
R
```

In the R console install devtools and then install Xmisc from the archive

```
>install.packages('Xmisc')
```

If the above does not work since Xmisc has been recently removed from CRAN, can install it through devtools. Exit R, and run the following. 

```
conda install -c conda-forge r-devtools
```

Start R 

```
R
```

Load devtools and make sure the environment variables are set properly

```
>Sys.setenv(TAR = "/bin/tar")
>install_url("https://cran.r-project.org/src/contrib/Archive/Xmisc/Xmisc_0.2.1.tar.gz")
```

Once this is done, you can exit `R`, and now you have an environment that can run the pipeline

### Add bin folder to PATH

Uncompress and add this folder to your $PATH

## (2) Reference Files

These are the required reference files for the code to work. The `arguments.txt` file that is within the directory of the scripts or that can be specified by the user will define where these files are located. We have default files ready for 3 assemblies (hg19, hg38, and mm10), and we plan on providing other assemblies as well. 

### (A) Defaults

We have created a default set of reference files for hg38, hg19, and mm10. For simple use, download the directory, place it in the rnapipeline folder, and extract. Make sure to update the arguments.txt file with the full path of these files (Explained below in Usage).

Download Link hg38: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefhg38.tar.gz 'Compressed Directory')

Download Link hg19: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefhg19.tar.gz 'Compressed Directory')

Download Link mm10: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefmm10.tar.gz 'Compressed Directory')

Note:
>If you use these default files then you can skip B-F. However, be sure to look at F if you want to use a list of genes different then the default list that is based on the publication. 

### (B) Gencode Dictionary

1. Download Gencode GTF reference file desired here: https://www.gencodegenes.org/

2. Sort the File

```
cat <GENCODE GTF>  | awk '{if($3=="transcript"||$3=="exon"||$3=="start_codon"){print}}' |  awk -F "; " '{print $0"\t"$2}' > <OUTPUT_sorted.gtf>`
```
 
3. Use custom script to create dictionary

```
genecode_to_dic.py <OUTPUT_sorted.gtf>
```
 
This will generate 2 files: (1) `genecode_plus.dic` and (2) `genecode_minus.dic`
 
4. Rename as needed

### (C) Repeatmasker Files

Though this pipeline is optimized for looking at repetitive-element derived transcripts, theoretically any bed file for alternative promoter locations can be used. 

1. Download Repeatmasker Annotation BED. We use the one from the UCSC Table Browser (For this step download it in BED format): [External Link](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=693256623_kh0RR0o6vajdA2WLLTA8OeaAPNB6&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=hg38Patch11&hgta_table=0&hgta_regionType=genome&position=chr12%3A20816734-20825794&hgta_outputType=bed&hgta_outFileName=repeats.bed 'UCSC Table Browser')

Note:
> If you would like to check for overlap with other features such as CpG islands or you need to create a custom file of your features, make sure it has the following format and tab-delimitted:
> 
> chromosome start end label wildcard strand

2. bgzip the file (bgzip comes with samtools)

Note:
> The bed file may need to be sorted
> 
> Command:
> cat <BED> | sort -k1,1 -k2,2n > <SORTED.BED>

2. bgzip the file (bgzip comes with samtools)


```
bgzip rmsk.bed > rmsk.bed.gz
```

3. Create tabix index

```
tabix -p bed rmsk.bed.gz
```

4. Both the bgzipped file and the tabix index must be in the same directory


### (D) Repeatmasker Subfamily, Family, Class Mapping

There needs to be a file that will map the subfamily of the TEs to Class and Family.

This can be done with repeatmasker. Download Repeatmasker Annotation, but this time with All fields NOT the Bed format. This should have 11 columns as well as a header column. ('output format:' field should have 'all fields from selected table')

To extract these annotations from that file:

```
cat repeatmasker.lst | awk 'NR>1{print $11"\t"$12"\t"$13}' | sort | uniq > repeatmasker_hg19_description_uniq.lst
```

### (E) Intron Annotations (Optional)

A useful feature that can be used is that the pipeline will annotate the first intron of each transcript based on a reference file. This can help in deciding whether a candidate has been previously annotated as an alternative transcript. 

Note:
> If you have already created the gencode dictionary by yourself in the 'Gencode Dictionary' step , then you have already done step 1 and 2 and you can skip to step 3

1. Download Gencode GTF reference file desired here: https://www.gencodegenes.org/

2. Sort the File

```
cat <GENCODE GTF>  | awk '{if($3=="transcript"||$3=="exon"||$3=="start_codon"){print}}' |  awk -F "; " '{print $0"\t"$2}' > <OUTPUT_sorted.gtf>`
```
 
3. Use custom script to create intron annotations

```
genecode_introns.py <OUTPUT_sorted.gtf>
```
 
This will generate 2 files: (1) `<OUTPUT_sorted.gtf>_introns_plus` and (2) `<OUTPUT_sorted.gtf>_introns_minus`

4. Sort the intron annotations bed file

```
sort -k1,1 -k2,2n -k3,3n <OUTPUT_sorted.gtf>_introns_plus > <OUTPUT_sorted.gtf>_introns_plus_sorted
sort -k1,1 -k2,2n -k3,3n <OUTPUT_sorted.gtf>_introns_minus > <OUTPUT_sorted.gtf>_introns_minus_sorted
```

5. bgzip the file (bgzip comes with samtools)

```
bgzip <OUTPUT_sorted.gtf>_introns_plus_sorted > <OUTPUT_sorted.gtf>_introns_plus_sorted.gz
bgzip <OUTPUT_sorted.gtf>_introns_minus_sorted > <OUTPUT_sorted.gtf>_introns_minus_sorted.gz
```

6. Create tabix index

```
tabix -p bed <OUTPUT_sorted.gtf>_introns_plus_sorted.gz
tabix -p bed <OUTPUT_sorted.gtf>_introns_minus_sorted.gz
```

7. Both the bgzipped file and the tabix index need to be in the same folder

### (F) Gene Filter List (Optional)

For large sets of analysis, it might be computationally advantageous to limit the analysis to only a small set of genes for subsequent analysis. Create a file with a genesymbol per line. Example:

```
TP63
GAPDH
SYT1
PKIB
...
```

## (3) Input Files

Our pipeline is optimized to work with `stringtie` generated gtf files. We have not tested it with other assembly software such as `cufflinks`. We plan on adding support for other input file styles in future releases. 

With raw FASTQ files, the following steps need to be taken:

1. QC: Using software such as [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ 'fastqc')

2. Adapter Trimming

3. Alignment: We  normally user and prefer [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf 'star').  Make sure that the -XS tag will be outputted since the splice-aware assembly software such as `stringtie` has the information it needs.
[HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml 'hisat2') can also be used, but HISAT uses 60 instead of 255 for uniquelly mapped reads which will interfere with code in future steps.

Note:

> Future steps require sorted and indexed bam files. Thus, outputting sorted BAM files will be beneficial. 

4. Sorting and Indexing: The BAM files need to be sorted by position (samtools sort) and need to be indexed (samtools index)

5. Filtering: Due to problems with mapping to repetitive elements, we usually will filter to only include uniquelly mapped reads. 

6. Assembly: [STRINGTIE](https://ccb.jhu.edu/software/stringtie/ 'stringtie') Default parameters will work if unstranded, but you should look at stringtie documentation to adjust paramaters based on type of sequencing run. For discovering low coverage transcripts we sometimes set the following flags (-c 1 -m 100). Note, future steps will require being able to map the gtf file to the bam file, so the easiest way to do this would be to just change the extension of the files from `*.bam` to `*.gtf` (test1.bam to test1.gtf).

Example of Command for steps 4 & 5:

```
samtools view -q 255 -h EGAR00001163970_636782_1_1.rnaseq.fastqAligned.sortedByCoord.out.bam | stringtie - -o EGAR00001163970_636782_1_1.rnaseq.fastqAligned.sortedByCoord.out.gtf -m 100 -c 1
```
Note:
> Use -q 60 for hisat2 aligned files

# Usage

Once all the gtf files are generated by stringtie for your experiments, the following steps can be used to annotate for alternative promoters from TEs. We have designed this pipeline to be flexible, and thus we provide explanations for each intermediate file that can be used for more custom usage of our code. 

## Quick Run (Not recommended unless using Wang Lab Servers)

If you have set up arguments.txt and have the reference files properly configured, we have a bash script that can run all the steps in the pipeline starting with a directory with the raw FASTQ files.

Run the following command
```
run_pipeline.sh -j <maximum jobs> -r <basename for gtf files> -t <treatment basename> 
```

All the subsequent steps will be run. Note, most should not use this. Also this assumes things about type of sequencing reads. If using this script for analysis for publication, please talk to Nakul first. 

## (1) Setup arguments.txt

The locations of the reference files need to be specified in a tab delimitted file. There should be no extra lines or headers. Here is an example:

```
rmsk /bar/nshah/reference/rmskhg38.bed6.gz
rmskannotationfile /bar/nshah/reference/repeatmasker_description_uniq.lst
gencodeplusdic /bar/nshah/reference/genecode_plus_hg38.dic
gencodeminusdic /bar/nshah/reference/genecode_minus_hg38.dic
focusgenes /bar/nshah/reference/oncogenes_augmented.txt
plusintron /bar/nshah/reference/gencode.v25.annotation.sorted.gtf_introns_plus_sorted.gz
minusintron /bar/nshah/reference/gencode.v25.annotation.sorted.gtf_introns_minus_sorted.gz
```

### Required Arguments

**rmsk:** tabix formatted bed6 files with repeatmasker or other file that user wants to check start location for

**rmskannotationfile:** a tab delimitted with 3 columns: subfamily class family

**gencodeplusdic:** Dictionary of all gencode elements including introns and exons for the plus (+) strand

**gencodeminusdic:** Dictionary of all gencode elements including introns and exons for the minus (-) strand


Note:
> These are all the arguments that are needed. The following arguments.txt file would work:

```
rmsk	/bar/nshah/reference/rmskhg38.bed6.gz
rmskannotationfile /bar/nshah/reference/repeatmasker_description_uniq.lst
gencodeplusdic	/bar/nshah/reference/genecode_plus_hg38.dic
gencodeminusdic	/bar/nshah/reference/genecode_minus_hg38.dic
```

### Optional Arguments

**focusgenes:** The program has two outputs (1) on a focus set of genes (2) with all genes. This file lists the genes that the user wants to filter for originally (Gene Filter List)

**plusintron:** Tabix file of all the plus strand introns (Intron Annotations)

**minusintron:** Tabix file of all the minus strand introns (Intron Annotations)

Note:
> If you do not want to use these options, remove them from the file. There should be no extra lines in the file or it will not work. 

## (2) Run annotation on each GTF file

Run the Program

```
rmskhg38_annotate_gtf_update_test_tpm.py <gtffile> <argumentfile.txt>*
```
>*This is optional. If it is not included then the program will default to the argumentfile.txt within the rnapipeline folder.

**Output File(s):**

(1) \<gtffile\>_annotated_test_all
 
(2) \<gtffile\>_annotated_filtered_test_all
 
(3) \<gtffile\>_annotated_test*
 
(4) \<gtffile\>_annotated_filtered_test*
 
 >*Files (3) and (4) will only be produced if the user would like to focus the analysis on a set of genes such as a list of oncogenes and tumor suppresor genes

Excel Document with Description of Columns: [External Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/Transcript%20Annotation%20Description.xlsx 'GTF Annotation Output')

## (3) Process annotations files to get rough estimate of relative expression of transcript versus all other transcripts of gene

Run the program on the dataset desired. If you would like to stick to the selection of genes then run this on all the \<gtffile\>_annotated_filtered_test files. For all genes, run it on the \<gtffile\>_annotated_filtered_test_all files. 

```
annotationtpmprocess.py \<gtffile\>_annotated_filtered_test_all
```

Output File(s):

(1) \<gtffile\>_annotated_filtered_test_(all)_c

>Three additional columns are made: (1) Maximum coverage of gene of interest (2) Maximum tpm of gene of interest (3) Total of gene of interest.
>Note: For those transcripts that are from a TE but do not splice into a gene, these stats will be compared across all the transcripts like this. Thus the fraction will be very low.  

## (4) Aggregate annotation samples across samples

To aggregate all the data across samples, we have an R script that will be able to aggregate the statistics and summarize an output.

Arguments are positional, and thus making sure they are in correct order is important.

This script can be run standalone with specifying the parameters, or you can do this step in RStudio to add some more custom filters if needed. The columns in the final file must be the same, however for the subsequent steps to work. 

Run the program

```
aggregateProcessedAnnotation.R <options>
```

**Options with Defaults:**

**-e** \<ext_treatment\> (default: ''): The label in the treatment file names that will identify them. If there is not a treatment versus untreated experimental design, the default will call everything as treatment.

**-l** \<exon 1 length max\> (default: 2588): The maximum length of exon 1. We are using the 99th percentile of gencode v25 transcripts

**-s** \<exon skipping max\> (default: 2): Based on genomic contamination, assembly can create transcripts that have intron retention. These could be real, but oftentimes these look like noise. A maximum level of 2 is recommended.

**-n** \<sample total min\> (default: 1): Minimum number of samples that must have a candidate for it to be considered. With low number of samples 1 is recommended. For larger studies this can be increased.

**-t** \<treatment total min\> (default: 0): In case the user would like to only consider candidates that are present in a certain number of treatment samples.

**-x** \<treatment exclusive\> (default: no): In case the user would like a maximum for the untreated samples. This could be if the user wants to only consider treatment-exclusive samples. In that case, the user should put yes instead of no. 

**-k** \<keep none\> (default: no): There can be transcripts from TEs that have splicing, but do not end up splicing into a gene. By default these are removed. This can be changed to yes if the user would like to look at these.

**-f** \<filter for TEs\> (default: yes): Repeatmasker files include many repeats besides TE. This could be an interesting benchmark to compare to, but by default they are removed. The user can specify no if they would like to keep these. NOTE: Currently does not work with downstream steps if filter is set to No.

**-a** \<argument file\> (default: \<directory of pipeline\>/arguments.txt): The arguments.txt file location. By default, the arguments.txt file that is in the rnapipeline directory will be used. If another is desired for use then the fullpath of the file can be specified. 

**Output File(s):**

(1) filter_combined_candidates.tsv: A file with every TE-gene transcript. This file is used for calculating read information in subsequent steps.
 
(2) initial_candidate_list.tsv: A summary of filter_combined_candidates.tsv for each unique candidate. Also lists the treatment and untreated files that the candidate is present in.
 
(3) Step4.RData: Workspace file with data loaded from R session. Subsequent steps load this to save time.

Note:
>This step will give ALL assembled transcripts. Assembly is very noisy, and thus the output of this step will have a high false positive rate.
>Filtering based on tpm or fraction of tpm can help.
>It is reccomended that the user looks at initial_candidate_list.tsv and assures that the candidates are in the format that they desire. 


## 5. Calculate Read Information

It is recommended that candidates are confirmed using read information. We have found this to greatly increase the specificity of the pipeline since assembly can make many errors especially without a reference.

We have developed a fast parser using samtools and bedtools to be able to get read information. There are separate programs for either single and or paired end reads. The only argument needed is the path of the bam files.

Note:
>Bam files should be sorted by position and indexed for this to work.
>Bam file naems should be the same as the gtf files that were the oriignal input for the pipeline except the extention is different.
>It is recommended that paired-end reads are used. Single-end performance is not as robust.

### (A) Make the directory where the stats will be stored
```
mkdir filterreadstats
```

### (B) Create commands to calculate read information

```
commandsmax_speed.py filter_combined_candidates.tsv <bam file location*>
```

If using single-end RNa-seq files, use the following command
```
commandsmax_speed_se.py filter_combined_candidates.tsv <bam file location*>
```

Note:
>*the bam file location should be a full path and should end in a forward slash (e.g. /scratch/nakul/epikate/)



**Output File(s):**

(1) filterreadcommands.txt: A list of all the commands needed to calculate read information. These commands use a combination of samtools, bedtools, and custom scripts. For every combination of candidate and bamfile present in teh dataset, there will be a command.

### (C) Run the commands
It is reccomended that [parallel_GNU](https://www.gnu.org/software/parallel/) or a similar method be used to run the list of commands in parallel. 

```
parallel_GNU -j 4 < filterreadcommands.txt
```
Note:
>(-j) will tell parallel_GNU the number of processes
>The alias for parallel_GNU might just parallel or something else on your server configuration

If using hisat2, do this before running the above commands so that samtools -q 60 is used
```
sed -i 's/ 255 / 60 /g' filterreadcommands.txt
```

### (D) Combine all the read information files

```
find ./filterreadstats -name "*.stats" -type f -maxdepth 1 -print0 | xargs -0 -n128 -P1 grep e > resultgrep_filterreadstatsdone.txt
cat resultgrep_filterreadstatsdone.txt | sed 's/\:/\t/g' > filter_read_stats.txt
```

This will create the file **filter_read_stats.txt** that has all the read information that is used in subsequent analysis.

## 6. Filter Candidates based on read information

Using the read annotation information, only the highest confidence candidates are kept as potential transcripts.

The default is 10 reads starting in the TE in the correct direction and having 2 reads going from the TE to the gene in the file. A limitation is that this can bias the detection against longer transcripts.

Note:
>This requires the Step4.Rdata created in step 4 that has an Rsession saved that will be reloaded for the program. 

Run the program

```
filterReadCandidates.R <options>
```
**Options with Defaults:**

**-r** \<minimum reads in TE\> (default: 10): The number of paired end reads required in a single file that start within the TE.

**-s** \<min start read\> (default: 1): The number of paired-end reads that span distance from TE to gene required across all files.  

**-e** \<exonization read support max percent\> (default: .15): Maximum percentage of files that have reads that start in gene and go to TE of interest, suggesting exonization rather than a promoter event

**-d** \<distance TE\> (default: 2500): Minimum distance between the TE and the start of the trawnscript. Helps remove noise from the promoter. 

**Output File(s):**

(1) read_filtered_candidates.tsv: A file with every TE-gene transcript. This file is used for calculating read information in subsequent steps.
 
(2) candidate_transcripts.gff3: A GFF3 formatted file of all the unique TE-gene transcripts detected. This can be viewed on a genome viewer to see what candidates look like. This is used in subsequent steps. 
 
(3) Step6.RData: Workspace file with data loaded from R session. Subsequent steps load this to save time.


Remove step 4 session data that is no longer needed.

```
rm Step4.RData
```

## 7. Merge with Reference GTF

Using the [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/ 'Cufflinks Github') package, merginging candidate transcripts with the reference. Also, candidates will be merged with each other which can help remove incomplete 5' transcripts that are detected.

Original method (used in paper)
```
gffread -E candidate_transcripts.gff3 -T -o candidate_transcripts.gtf
echo candidate_transcripts.gtf > cuffmergegtf.list
cuffmerge -o ./merged_asm_full -g ~/reference/Gencode/gencode.v25.basic.annotation.gtf cuffmergegtf.list
mv ./merged_asm_full/merged.gtf reference_merged_candidates.gtf
gffread -E reference_merged_candidates.gtf -o- > reference_merged_candidates.gff3
```

## 8. Annotate Merged GTF

Using a lightly modified version of our previous annotation script, the new merged reference can be annotated for a final set of TE-gene candidates. 

```
rmskhg38_annotate_gtf_update_test_tpm_cuff.py reference_merged_candidates.gff3 <argumentfile.txt>*
```
>*This is optional. If it is not included then the program will default to the argumentfile.txt within the rnapipeline folder

**Output File(s):**

(1) reference_merged_candidates.gff3_annotated_test_all
 
(2) reference_merged_candidates.gff3_annotated_filtered_test_all <(Final List of all TE-gene transcripts)
 
(3) reference_merged_candidates.gff3_annotated_test*
 
(4) reference_merged_candidates.gff3_annotated_filtered_test*

Note:
>We recommend users to view the reference_merged_candidates.gff3 on a genome browser to confirm that the merged transcript file has the expected candidates. The user can map transcripts back to their annotation using the transcript ID and confirm proper identification. 

## 9. Calculate Transcript-Level Expression

Transcript-level expression can be calculated to compare the TE-gene transcripts to the canonical gene transcript to find those that are contributing significantly to overall expression of the gene. 

```
samtools view -q 255 -h <bam filename> | stringtie - -o \<bam filename_root\>.gtf -e -b \<bam filename_root\>_stats -p 2 -m 100 -c 1 -G reference_merged_candidates.gtf
```

An example of how to create a file with all the commands across multiple bam files so they can be run in parallel is seen below
```
find <bam_directory> -name "*bam" | while read file ; do xbase=${file##*/}; echo "samtools view -q 255 -h "$file" | stringtie - -o "${xbase%.*}".gtf -e -b "${xbase%.*}"_stats -p 2 -m 100 -c 1 -G reference_merged_candidates.gtf" >> quantificationCommands.txt ; done ;
```
>Subsequently, the commands in quantificationCommands.txt can be run individually or can be run in parallel (e.g. parallel_GNU -j 4 < quantificationCommands.txt)
>
>hisat2 users: 255 needs to be replaced with 60 again

Note:
>The -b flag outputs all the stats to a new folder. This is in a format that is compatible with the [Ballgown](https://github.com/alyssafrazee/ballgown 'Ballgown Github')  downstream pipeline that can be used instead of our own custom methods of transcript identification. 

## 10. Process and map expression output

### (A) Obtain annotations, filter by candidates identified previously, and identify major splicing intron

The following script will aggregate the information from the annotation of the merged transcript and make the final list of candidates.

```
mergeAnnotationProcess.R <options>
```
**Options with Defaults:**

**-f** \<merged gtf annotation file\> (default: reference_merged_candidates.gff3_annotated_filtered_test_all): Reference transcript file to be processed. Default will work if all the previous steps have been done as described with the same names. 

**Output File(s):**

(1) candidate_introns.txt: A large one-line string with all intron locations. This is used in susequent steps to extract intron read information
 
(2) candidate_names.txt: The trnascript names for all candidate transcripts that are left

(3) Step10.RData: Workspace file with data loaded from R session. Subsequent steps load this to save time.

Remove step 6 session data that is no longer needed.

```
rm Step6.RData
```

### (B) Process stringtie transcript annotation files to get relevant information and aggregate

Th stringtie output needs to be formmated into a table, and thus we have a few helper scripts and bash commands to make these tables. 

Obtaining intron coverage information

```
find . -name "*i_data.ctab" > ctab_i.txt

cat ctab_i.txt | while read ID ; do fileid=$(echo "$ID" | awk -F "/" '{print $2}'); cat <(printf 'chr\tstrand\tstart\tend\t'${fileid/_stats/}'\n') <(grep -F -f candidate_introns.txt $ID | awk -F'\t' '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6 }') > ${ID}_cand ; done ;

cat <(find . -name "*i_data.ctab_cand" | head -1 | while read file ; do cat $file | awk '{print $1"\t"$2"\t"$3"\t"$4}' ; done;) > table_i_all

find . -name "*i_data.ctab_cand" | while read file ; do paste -d'\t' <(cat table_i_all) <(cat $file | awk '{print $5}') > table_i_all_temp; mv table_i_all_temp table_i_all; done ;
```


Obtaining the transcript-level expression information for candidates.
```
ls ./*stats/t_data.ctab > ctablist.txt

cat ctablist.txt | while read file ; do echo "stringtieExpressionFrac.py $file" >> stringtieExpressionFracCommands.txt ; done;
```

Run the commands in stringtieExpressionFracCommands.txt (Much faster if done in parallel)
```
parallel_GNU -j <number of jobs> < stringtieExpressionFracCommands.txt
```

Aggregate the stats
```
ls ./*stats/t_data.ctab_frac_tot > ctab_frac_tot_files.txt
ls ./*stats/t_data.ctab_tpm > ctab_tpm_files.txt

cat <(echo "TranscriptID") <(find . -name "*ctab_frac_tot" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_frac_tot
cat ctab_frac_tot_files.txt | while read file ; do fileid=$(echo "$file" | awk -F "/" '{print $2}') ; paste -d'\t' <(cat table_frac_tot) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_frac_tot_temp; mv table_frac_tot_temp table_frac_tot; done ;

cat <(echo "TranscriptID") <(find . -name "*ctab_tpm" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_tpm
cat ctab_tpm_files.txt | while read file ; do fileid=$(echo "$file" | awk -F "/" '{print $2}') ; paste -d'\t' <(cat table_tpm) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_tpm_temp; mv table_tpm_temp table_tpm; done ;

cat <(head -1 table_frac_tot) <(grep -Ff candidate_names.txt table_frac_tot) > table_frac_tot_cand
cat <(head -1 table_tpm) <(grep -Ff candidate_names.txt table_tpm) > table_tpm_cand
```

## 11. Quantification processing, sample identification, and final table creation

The script will aggregate the stringtie annotation and intron coverage information for the candidates. It will predict presence based on parameters set by the user. The tables that are outputted by this step can also be used by the user to do more advanced statistical analysis.

Run the program

```
finalStatisticsOutput.R <options>
```
**Options with Defaults:**

**-e** \<ext_treatment\> (default: ''): The label in the treatment file names that will identify them. If there is not a treatment versus untreated experimental design, the default will call everything as treatment.

**-i** \<minimum reads spanning intron junction\> (default: 1): The minimum number of reads across an intron junction needed to state that a candidate is present. 

**-t** \<minimum gene tpm\> (default: 1): The minimum total TPM of the gene required. (If you plan to do differential gene expression analysis instead, this can be put as 0 to keep everything)

**-a** \<argument file\> (default: \<directory of pipeline\>/arguments.txt): The arguments.txt file location. By default, the arguments.txt file that is in the rnapipeline directory will be used. If another is desired for use then the fullpath of the file can be specified. 


**Output File(s):**

(1) All TE-derived Alternative Isoforms Statistics.xlsx: A file with the final statistics calculated for each candidate. There is also data on the gene expression in both groups (treatment and normal), fraction of gene expression (treatment and normal), the number of reads to main splicing intron (treatment and normal), and treatment enrichment.

Note:
> The **Treatment Count** and **Normal Count** columns are calculated based on the number of files in which the candidate passes the thresholds set on fraction of expression, total gene expression, and number of reads spanning the main splicing intron. The final table has all the
> data used for this, so the user can try using different thresholds to optimize candidate presence based on their data. 
 
(2) allCandidateStatistics.tsv: file with gene expression, fraction expression, transcript expression, and intron junction read information across all the samples for all the candidates. 
 
(3) Step11_FINAL.RData: Workspace file with data loaded from R session. Can be loaded by user to save time to do more advanced analysis. 


Remove old RData, the final one will have all data from previous ones.

```
rm Step10.RData
```

## 12. Translating the transcripts with Kozak method and generating FASTA of RNa sequences

The script will take the data about the coordinates of the transcripts and translate them. Of note, it must know the correct genome to translate based off of. 

Note:
> Should consider filtering transcripts that are in the annotatedcufftranscripts file. There are many ways to filter the data, so one can load the Step11_FINAL.RData into R, filter the annotatedcufftranscripts for the candidates based on filtering, and then save Step11_FINAL.RData and then only the filtered transcripts will be used. In the case of our publication, we filtered based on presence in TCGA normal samples and GTEx normal samples in choosing the candidate transcripts. 

Run the program

```
translatePart1.R <options>
```
**Options with Defaults:**

**-g** \<ext_treatment\> (default: 'BSgenome.Hsapiens.UCSC.hg38'): The genome to use for analysis. Our analysis is based on the hg38 genome. 

**Output File(s):**

(1) candidates.fa: This is the fasta with the RNA sequences of all the candidate transcripts. This can be used in future steps with programs to decide on translational products. 
 
(2) Step12.RData: Workspace file with data loaded from R session. Can be loaded by user to save time to do more advanced analysis. 


Remove old RData, the final one will have all data from previous ones.

```
rm Step11_FINAL.RData
```

## 13. CPC2 Translation of the transcripts 

In addition to the Kozak method that is listed here, we also output the longest ORF that is available. This method also uses the program CPC2 (https://academic.oup.com/nar/article/45/W1/W12/3831091) to assess translation potential. We have modified their script to output the start codon of interest. Install CPC2 from https://github.com/nakul2234/CPC2_Archive/blob/main/bin/CPC2.py and assure 'CPC2.py' is in your $PATH. 

Run the program

```
CPC2.py -i candidates.fa -o candidates_cpcout.fa
```

**Output File(s):**

(1) candidates_cpcout.fa: CPC2 outpout that will parsed in future analysis steps. 

## 14. Finishing translation, classifying antigens, and outputting final list of potential antigens

The script will take the CPC2 output and then incorporate it. It will then output all the potential antigenic peptides from the CPC2 method and the Kozak method that are present. This fasta can be used for searching in mass spectrometry programs. 

Run the program

```
translatePart2.R <options>
```
**Options with Defaults:**

**-g** \<ext_treatment\> (default: 'BSgenome.Hsapiens.UCSC.hg38'): The genome to use for analysis. Our analysis is based on the hg38 genome. 

**Output File(s):**

(1) candidates.fa: This is the fasta with the RNA sequences of all the candidate transcripts. This can be used in future steps with programs to decide on translational products. 
 
(2) Step13.RData: Workspace file with data loaded from R session. Can be loaded by user to save time to do more advanced analysis. 


Remove old RData, the final one will have all data from previous ones.

```
rm Step12.RData
```

## 15. Optional further analysis with ballgown

 [STRINGTIE](https://ccb.jhu.edu/software/stringtie/ 'stringtie') has a downstream analysis pipeline ([Ballgown](https://github.com/alyssafrazee/ballgown 'Ballgown Github')) that allows for transcript level and gene level expression analysis between samples.
 In the folder that has the <*stats> folders run the following commands that will create a folder "ballgown" that will allow for you to do more traditional gene and transcript-level differential expression analysis. 

```
mkdir ballgown
cd ballgown
ls -d ../*_stats | while read file ; do mkdir $(basename $file) ; cd $(basename $file) ; (ls ../${file}/*ctab | while read file2 ; do ln -s $file2 ; done ;) ; cd .. ; done ;
```

There are also options to make this output compatible with more traditional statistical analysis tools. 

# Usage Reference Guided

This is a version of the pipeline where you can skip the custom reference creation. For example, if you are basing this off of our published TCGA analysis, then using the created reference files and annotation from that analysis will suffice and you do not need the other metadata. 

This version assumes that the user has already obtained a reference set of transcripts that they are comfortable with and they want to test across their samples. In this case, they just need to quantify the transcripts.

It is important that you have a valid `arguments.txt` file pointing to the references for the pipeline. Please see step 1 of the previous section to see how to set it up [(1) Setup arguments.txt](#1-setup-argumentstxt). 

## 1. Obtain reference gtf and gff3 file for analysis

The one based on our TCGA 33 tumor analysis can be downloaded here: [External Download Link](https://wangftp.wustl.edu/~nshah/ucsf/TCGA33Download/ 'Reference GTF')

If you are using a custom GTF program from another source, use cufflinks gffread to be able to generate the gff3 file

```
gffread -E <input.gtf> -o- > reference_merged_candidates.gff3
```

## 2. Annotate the gff3 file, or obtain annotations from a previous run

The one based on our TCGA 33 tumor analysis can be downloaded here:

Using a lightly modified version of our previous annotation script, the new merged reference can be annotated for a final set of TE-gene candidates. 

```
rmskhg38_annotate_gtf_update_test_tpm_cuff.py reference_merged_candidates.gff3 <argumentfile.txt>*
```
>*This is optional. If it is not included then the program will default to the argumentfile.txt within the rnapipeline folder

**Output File(s):**

(1) reference_merged_candidates.gff3_annotated_test_all
 
(2) reference_merged_candidates.gff3_annotated_filtered_test_all <(Final List of all TE-gene transcripts)
 
(3) reference_merged_candidates.gff3_annotated_test*
 
(4) reference_merged_candidates.gff3_annotated_filtered_test*

## 3a. Calculate Transcript-Level Expression

Transcript-level expression can be calculated to compare the TE-gene transcripts to the canonical gene transcript to find those that are contributing significantly to overall expression of the gene. 

```
samtools view -q 255 -h <bam filename> | stringtie - -o \<bam filename_root\>.gtf -e -b \<bam filename_root\>_stats -p 2 -m 100 -c 1 -G reference_merged_candidates.gtf
```

An example of how to create a file with all the commands across multiple bam files so they can be run in parallel is seen below
```
find <bam_directory> -name "*bam" | while read file ; do xbase=${file##*/}; echo "samtools view -q 255 -h "$file" | stringtie - -o "${xbase%.*}".gtf -e -b "${xbase%.*}"_stats -p 2 -m 100 -c 1 -G reference_merged_candidates.gtf" >> quantificationCommands.txt ; done ;
```
>Subsequently, the commands in quantificationCommands.txt can be run individually or can be run in parallel (e.g. parallel_GNU -j 4 < quantificationCommands.txt)
>
>hisat2 users: 255 needs to be replaced with 60 again

Note:
>The -b flag outputs all the stats to a new folder. This is in a format that is compatible with the [Ballgown](https://github.com/alyssafrazee/ballgown 'Ballgown Github')  downstream pipeline that can be used instead of our own custom methods of transcript identification. 

## 3b. Process and map expression output

### (A) Obtain annotations, filter by candidates identified previously, and identify major splicing intron

The following script will aggregate the information from the annotation of the merged transcript and make the final list of candidates. It also has additional options for filtering in the reference guided-mode to match the pipeline

```
mergeAnnotationProcess_Ref.R <options>
```
**Options with Defaults:**

**-g** \<merged gtf annotation file\> (default: reference_merged_candidates.gff3_annotated_filtered_test_all): Reference transcript file to be processed. Default will work if all the previous steps have been done as described with the same names. 

**-l** \<exon 1 length max\> (default: 2588): The maximum length of exon 1. We are using the 99th percentile of gencode v25 transcripts

**-s** \<exon skipping max\> (default: 2): Based on genomic contamination, assembly can create transcripts that have intron retention. These could be real, but oftentimes these look like noise. A maximum level of 2 is recommended.

**-k** \<keep none\> (default: no): There can be transcripts from TEs that have splicing, but do not end up splicing into a gene. By default these are removed. This can be changed to yes if the user would like to look at these.

**-f** \<filter for TEs\> (default: yes): Repeatmasker files include many repeats besides TE. This could be an interesting benchmark to compare to, but by default they are removed. The user can specify no if they would like to keep these. NOTE: Currently does not work with downstream steps if filter is set to No.

**-a** \<argument file\> (default: \<directory of pipeline\>/arguments.txt): The arguments.txt file location. By default, the arguments.txt file that is in the rnapipeline directory will be used. If another is desired for use then the fullpath of the file can be specified.  

**Output File(s):**

(1) candidate_introns.txt: A large one-line string with all intron locations. This is used in susequent steps to extract intron read information
 
(2) candidate_names.txt: The trnascript names for all candidate transcripts that are left

(3) Step10.RData: Workspace file with data loaded from R session. Subsequent steps load this to save time.


### (B) Process stringtie transcript annotation files to get relevant information and aggregate

Th stringtie output needs to be formmated into a table, and thus we have a few helper scripts and bash commands to make these tables. 

Obtaining intron coverage information

```
find . -name "*i_data.ctab" > ctab_i.txt

cat ctab_i.txt | while read ID ; do fileid=$(echo "$ID" | awk -F "/" '{print $2}'); cat <(printf 'chr\tstrand\tstart\tend\t'${fileid/_stats/}'\n') <(grep -F -f candidate_introns.txt $ID | awk -F'\t' '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6 }') > ${ID}_cand ; done ;

cat <(find . -name "*i_data.ctab_cand" | head -1 | while read file ; do cat $file | awk '{print $1"\t"$2"\t"$3"\t"$4}' ; done;) > table_i_all

find . -name "*i_data.ctab_cand" | while read file ; do paste -d'\t' <(cat table_i_all) <(cat $file | awk '{print $5}') > table_i_all_temp; mv table_i_all_temp table_i_all; done ;
```


Obtaining the transcript-level expression information for candidates.
```
ls ./*stats/t_data.ctab > ctablist.txt

cat ctablist.txt | while read file ; do echo "stringtieExpressionFrac.py $file" >> stringtieExpressionFracCommands.txt ; done;
```

Run the commands in stringtieExpressionFracCommands.txt (Much faster if done in parallel)
```
parallel_GNU -j <number of jobs> < stringtieExpressionFracCommands.txt
```

Aggregate the stats
```
ls ./*stats/t_data.ctab_frac_tot > ctab_frac_tot_files.txt
ls ./*stats/t_data.ctab_tpm > ctab_tpm_files.txt

cat <(echo "TranscriptID") <(find . -name "*ctab_frac_tot" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_frac_tot
cat ctab_frac_tot_files.txt | while read file ; do fileid=$(echo "$file" | awk -F "/" '{print $2}') ; paste -d'\t' <(cat table_frac_tot) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_frac_tot_temp; mv table_frac_tot_temp table_frac_tot; done ;

cat <(echo "TranscriptID") <(find . -name "*ctab_tpm" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_tpm
cat ctab_tpm_files.txt | while read file ; do fileid=$(echo "$file" | awk -F "/" '{print $2}') ; paste -d'\t' <(cat table_tpm) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_tpm_temp; mv table_tpm_temp table_tpm; done ;

cat <(head -1 table_frac_tot) <(grep -Ff candidate_names.txt table_frac_tot) > table_frac_tot_cand
cat <(head -1 table_tpm) <(grep -Ff candidate_names.txt table_tpm) > table_tpm_cand
```

## 4. Quantification processing, sample identification, and final table creation

The final script will aggregate teh stringtie annotation and intron coverage information for the candidates. It will predict presence based on parameters set by the user. The tables that are outputted by this step can also be used by the user to do more advanced statistical analysis.

Run the program

```
finalStatisticsOutput.R <options>
```
**Options with Defaults:**

**-e** \<ext_treatment\> (default: ''): The label in the treatment file names that will identify them. If there is not a treatment versus untreated experimental design, the default will call everything as treatment.

**-i** \<minimum reads spanning intron junction\> (default: 1): The minimum number of reads across an intron junction needed to state that a candidate is present. 

**-t** \<minimum gene tpm\> (default: 1): The minimum total TPM of the gene required. (If you plan to do differential gene expression analysis instead, this can be put as 0 to keep everything)

**-a** \<argument file\> (default: \<directory of pipeline\>/arguments.txt): The arguments.txt file location. By default, the arguments.txt file that is in the rnapipeline directory will be used. If another is desired for use then the fullpath of the file can be specified. 


**Output File(s):**

(1) All TE-derived Alternative Isoforms Statistics.xlsx: A file with the final statistics calculated for each candidate. There is also data on the gene expression in both groups (treatment and normal), fraction of gene expression (treatment and normal), the number of reads to main splicing intron (treatment and normal), and treatment enrichment.

Note:
> The **Treatment Count** and **Normal Count** columns are calculated based on the number of files in which the candidate passes the thresholds set on fraction of expression, total gene expression, and number of reads spanning the main splicing intron. The final table has all the
> data used for this, so the user can try using different thresholds to optimize candidate presence based on their data. 
 
(2) allCandidateStatistics.tsv: file with gene expression, fraction expression, transcript expression, and intron junction read information across all the samples for all the candidates. 
 
(3) Step11_FINAL.RData: Workspace file with data loaded from R session. Can be loaded by user to save time to do more advanced analysis. 


Remove old RData, the final one will have all data from previous ones.

```
rm Step10.RData
```

## 5. Translating the transcripts with Kozak method and generating FASTA of RNa sequences

The script will take the data about the coordinates of the transcripts and translate them. Of note, it must know the correct genome to translate based off of. 

Note:
> Should consider filtering transcripts that are in the annotatedcufftranscripts file. There are many ways to filter the data, so one can load the Step11_FINAL.RData into R, filter the annotatedcufftranscripts for the candidates based on filtering, and then save Step11_FINAL.RData and then only the filtered transcripts will be used. In the case of our publication, we filtered based on presence in TCGA normal samples and GTEx normal samples in choosing the candidate transcripts. 

Run the program

```
translatePart1.R <options>
```
**Options with Defaults:**

**-g** \<ext_treatment\> (default: 'BSgenome.Hsapiens.UCSC.hg38'): The genome to use for analysis. Our analysis is based on the hg38 genome. 

**Output File(s):**

(1) candidates.fa: This is the fasta with the RNA sequences of all the candidate transcripts. This can be used in future steps with programs to decide on translational products. 
 
(2) Step12.RData: Workspace file with data loaded from R session. Can be loaded by user to save time to do more advanced analysis. 


Remove old RData, the final one will have all data from previous ones.

```
rm Step11_FINAL.RData
```

## 6. CPC2 Translation of the transcripts 

In addition to the Kozak method that is listed here, we also output the longest ORF that is available. This method also uses the program CPC2 (https://academic.oup.com/nar/article/45/W1/W12/3831091) to assess translation potential. We have modified their script to output the start codon of interest. Install CPC2 from https://github.com/nakul2234/CPC2_Archive/blob/main/bin/CPC2.py and assure 'CPC2.py' is in your $PATH. 

Run the program

```
CPC2.py -i candidates.fa -o candidates_cpcout.fa
```

**Output File(s):**

(1) candidates_cpcout.fa: CPC2 outpout that will parsed in future analysis steps. 

## 7. Finishing translation, classifying antigens, and outputting final list of potential antigens

The script will take the CPC2 output and then incorporate it. It will then output all the potential antigenic peptides from the CPC2 method and the Kozak method that are present. This fasta can be used for searching in mass spectrometry programs. 

Run the program

```
translatePart2.R <options>
```
**Options with Defaults:**

**-g** \<ext_treatment\> (default: 'BSgenome.Hsapiens.UCSC.hg38'): The genome to use for analysis. Our analysis is based on the hg38 genome. 

**Output File(s):**

(1) candidates.fa: This is the fasta with the RNA sequences of all the candidate transcripts. This can be used in future steps with programs to decide on translational products. 
 
(2) Step13.RData: Workspace file with data loaded from R session. Can be loaded by user to save time to do more advanced analysis. 


Remove old RData, the final one will have all data from previous ones.

```
rm Step12.RData
```

## 8. Optional further analysis with ballgown

 [STRINGTIE](https://ccb.jhu.edu/software/stringtie/ 'stringtie') has a downstream analysis pipeline ([Ballgown](https://github.com/alyssafrazee/ballgown 'Ballgown Github')) that allows for transcript level and gene level expression analysis between samples.
 In the folder that has the <*stats> folders run the following commands that will create a folder "ballgown" that will allow for you to do more traditional gene and transcript-level differential expression analysis. 

```
mkdir ballgown
cd ballgown
ls -d ../*_stats | while read file ; do mkdir $(basename $file) ; cd $(basename $file) ; (ls ../${file}/*ctab | while read file2 ; do ln -s $file2 ; done ;) ; cd .. ; done ;
```

There are also options to make this output compatible with more traditional statistical analysis tools. 


# Questions?

If you have any questions, please post in the issues of the repository. 


