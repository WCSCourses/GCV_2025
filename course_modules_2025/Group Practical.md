In this group project, you will analyse the genome sequence data of the Hepatitis C Virus (HCV) you generated using Illumina sequencing. Please download the data according to your sample numbers. 


Aims: The aims of the project are to:
1) Check the quality of the data and get the data statistics
2) Clean the data if necessary
3) Find a suitable reference to map the reads
4) Identify the HCV genotype
5) Maps the quality checked reads to the best reference
6) Generate the mapping statistics and the coverage plot
7) Generate the consensus sequence
8) Find the drug resistance mutations, if there are any
9) Any recommended drugs to suggest (drugs to avoid)
10) Check for any other virus(es) in the metagenomic samples (?).

### Download the data

Data is available on 
https://tinyurl.com/GCV2025data

Create a new directory in your home folder and move the downloaded data there. 
```
mkdir ~/Practicals/.
mv ~/Downloads/*.gz ~/Practicals/.
cd ~/Practicals
```


Activate the conda environment prior to your analysis
```
conda activate bioinformatics_env
```

### QC check and filtering
Use `fastqc` program to get the information about "Total Sequences" and Sequence length.

```
fastqc & 
```

Load your sequences and check the "Basic Statistics". Check the "Per base sequence quality" and decide how you want to clean your data.

Use `trim_galore` or `prinseq` or `trimmomatic` program to clean your data. Change the file names as per your sample names. If your data is compressed, you can use it without uncompressing it. Most of these programs can process compressed data. 

Typical usage:
```
trim_galore --length 100 -q 30 --paired file-1.fq file-2.fq
```


### Best reference selection
HCV has eight genotypes and several subtypes. It is crucial to select the right reference to map your reads. Wrong reference selection will lead to an erroneous interpretation of the data. I have compiled and saved all the HCV genotypes in a file. You must map all the reads to these HCV reference genomes and select the best reference. All the HCV references are in the `~/Data-Sreenu/HCV\ genotypes` directory. This file contains HCV reference genomes based on their genotype and GenBank ID.

```
cd ~/Data-Sreenu/HCV\ genotypes/
bwa index hcvGenotypes.fa
```

Use `bwa mem` to map the cleaned reads and save the sorted output in a `bam` file:
```
cd ~/Practicals
bwa mem ~/Data-Sreenu/HCV\ genotypes/hcvGenotypes.fa file-1.fq file-2.fq|samtools sort - -o output.bam
```

Check the mapping of the reads to all the references. Sort the results according to the genome coverage.
```
samtools coverage output.bam |sort -rgk4
```


You have to use the top reference to map the reads. First, extract only the reference genome sequence from the reference genomes file.

```
cd ~/Data-Sreenu/HCV\ genotypes/
samtools faidx hcvGenotypes.fa refID > refID.fa
bwa index refID.fa
```

### Reference mapping

Now follow the reference mapping steps to map all the cleaned reads to this reference.
```
cd ~/Practicals
bwa mem ~/Data-Sreenu/HCV\ genotypes/refID.fa file-1.fq file-2.fq|samtools sort - -o output.bam
samtools index output.bam
```

Use `weeSAM` to get the coverage plot and `samtools coverage` to get all the statistics.

```
samtools coverage output.bam
weeSAM --bam output.bam --html outputName
```


If you prefer, you can generate the coverage plot using other programs by generating the depth file.

```
samtools depth -aa output.bam > depth.txt
```

The `depth.txt` will be a three-column file with reference name, genomic location and read depth.  

### Consensus generation

Generate the consensus sequence using `samtools consensus`

```
samtools consensus --ff UNMAP,SECONDARY,QCFAIL,DUP  --show-del yes --show-ins yes output.bam -o consensusName.fa
```

Please see the `samtools consensus`  help to learn more about the options used here. 

### Finding drug-resistance mutations

Please use the [HCV Glue](https://hcv-glue.cvr.gla.ac.uk/) server to find the drug-resistance mutations in your sequence. 
- Go to the HCV Glue website (https://hcv-glue.cvr.gla.ac.uk)
- Select "Analysis -> Genotyping and Interpretation" from the top menu bar.
- Click "Add files", select your consensus sequence file and Submit.
- Once the analysis is complete, check the "Show response" and view the Full report.


Is the genotype the same as you identified earlier?
What is the best reference suggested?
Are there any drugs to avoid? 



If you have time, check for other viruses in the samples using the steps provided in the metagenomic practicals. 
