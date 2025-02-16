# Command line tools for metagenomics analysis

## Background
One of the main uses of metagenomics in clinical laboratories is in the diagnosis of patients with serious central nervous system infections like encephalitis. Encephalitis can have lots of possible causes, including infection with viruses, bacteria, fungi or parasites, an autoimmune response, or the toxic effects of drugs. We would usually test invasive and precious samples like CSF or brain biopsy, so there isn’t much material to run lots of tests. This makes an untargeted approach like metagenomics useful.

In this tutorial, you will analyse a (not real!) dataset from a CSF sample from a patient with encephalitis. This will help us determine if an infection is causing their disease. This dataset was obtained using metagenomic sequencing on an Illumina NextSeq. You will analyse a small subset of the total dataset to speed up the analysis, although typically millions of reads would be used. A negative control sample has been run alongside the sample, which consists of commercially available human DNA and RNA.

## Before you start

First, activate the conda environment that contains some of the required software. Run the command:
```
conda activate bioinfo_env
```

We’ll run the tutorial in the ~/metagenomics directory (remember the ~ means the home directory). Navigate to this directory and have a look inside. You should find:
- data directory: contains the sample and negative control datasets, as well as the human genome (it’s not the actual human genome due to space constraints, but it will work for this practical)
- kraken_db directory: contains the database that we’ll use for classification later

Throughout this practical, you should process the sample and negative control datasets in the same way. This means you’ll run most commands twice, once for the sample and once for the negative control (don’t worry, there are faster ways to do this if you have lots of samples!). You’ll compare the results at the end.

Try to work out the commands yourself first. If you get stuck, there are some clues within the tutorial, and the answers are in a separate document. Don’t worry if you don’t finish the whole tutorial in the time allowed – you’ll get more out of it if you try to do a few steps yourself and all the commands are available for you to refer to later.

## Quality Control

The first step in most analysis protocols for sequencing data is quality control: removal of sequencing adaptors and any low-quality sequences from the ends of reads. We would also usually check our data with a program like fastqc – we’ll skip this stage today to save time.

The files are available in the folder ~/metagenomics/data

**1.	Write a command to trim adaptors and low quality regions from your data.**

Use your notes from a previous session on this course.

<details>
<summary>Clues</summary>
    
Try using trim_galore from the FileFormats-QC session.    
</details>

## Human removal
We want to find out what microbes are in the sample, so we are not interested in the human reads. We’ll therefore filter them out with alignment before we do anything else.

**2.	Write commands to align the reads to the human genome.**

Again, use your previous notes. Remember to index the genome first.

<details>
<summary>Clues</summary>
    
    Try using bwa mem from the ReferenceAlignment session.    
</details>

**3.	What is the output format of the alignment?**

**4.	How could we use our alignment to get the non-human reads?**

You’ve not yet been given the commands to do this, so they are below. Remember you might need to change the file names if you’ve used a different naming convention. Make sure you understand what they’ve doing before you run them!

```
samtools view -bf 4 -h ~/metagenomics/sample1.sam > ~/metagenomics/sample1_nonhuman.bam

samtools fastq ~/metagenomics/sample1_nonhuman.bam -1 ~/metagenomics/sample1_nonhuman_1.fq -2 ~/metagenomics/sample1_nonhuman_2.fq
```
You can ignore this error message:
```
samtools: /home/manager/miniconda3/envs/bioinfo_env/bin/../lib/libtinfow.so.6: no version information available (required by samtools)
```

## Classification

Now you’ve performed quality control and removed the human reads, we’re ready to run a classifier to determine what species are present. We’re going to use the programs kraken2 and bracken, which are some of the most widely used tools for this purpose.
Kraken2 and Bracken need a database of known reference sequences. This can be found in the ~/metagenomics/kraken_db directory. If you’re trying this yourself on your own computer, you can download prebuilt databases from: https://benlangmead.github.io/aws-indexes/k2
Whenever you run a new bioinformatics tool, it’s a good idea to look at search for the manual online. You can usually find out how to run the tool using the help command, for example:

```
kraken2 --help
kraken2 -h
```

**5.	Write a command to run kraken2 on your data.**
    
Search for the manual online or use the help command. You won't need most of the options - the default parameters are fine for this purpose. Hint: you should use the --report parameter.

<details>
<summary>Clues</summary>
    
   Your command should take the form:

   <pre>
kraken2 --db <i>kraken_database_directory</i> --paired <i>human_filtered_input_file_1 human_filtered_input_file_2</i> --report <i>output_report_filename</i> > <i>read_classifications_filename</i>
    </pre>

    Swap the parts in italics for your file names.
    
</details>

Once you’ve run kraken2, you can perform post-processing with bracken.

**6.	Write a command to run bracken on your kraken2 output.**

Again, use the manual online or the help command to do this. The database you need is the same one as for kraken2. Change the threshold parameter (-t) to 3 for this small test run. You can leave read_len and level as their defaults.

<details>
<summary>Clues</summary>
    
   Your command should take the form:

   <pre>
bracken -d $kraken_db -i <i>kraken_report</i> -o <i>output_file_name</i>
    </pre>

    Swap the parts in italics for your file names.
    
</details>

## Interpretation

**7.	What species are present in your sample?**

<details>
<summary>Clues</summary>
    
    Look at the report produced by Bracken.    
</details>

**8.	Look at the negative control. What can you conclude about what might be causing the disease in the patient?**

<details>
<summary>Clues</summary>
    
    Are there any species that are present in both the sample and negative control?    
</details>

## Extension questions

Save these until the end if you're short on time!

**9.	Using what you've learnt in previous sessions, what further analysis you think might be useful on these datasets?**

<details>
<summary>Clues</summary>
    
    How would you find out where the reads come from in the viral genome? 
</details>


**10.	How could you run the commands in this tutorial for multiple samples at a time?**

<details>
<summary>Clues</summary>
    
    Try writing a for loop in bash.
</details>

**11.	How do Kraken2 and Bracken work?**

<details>
<summary>Clues</summary>
    
    Look at the published articles that describe these protocols.
</details>

#  Online tools for metagenomics analysis
You should now be familiar with running metagenomics analysis on the command line and the steps involved in this process. However, sometimes you might find it easier to use an online tool. Today we're going to test CZID, a freely available online tool for metagenomics analysis.

**If you're using CZID to analyse your own samples, make sure you have permission to upload any human data.**

In future, you can make your own free account, but to save time in the session please use shared teaching account - you'll be provided with login details.

After you've logged in, navigate to the [Medical Detectives dataset](https://czid.org/public?currentDisplay=table&currentTab=samples&mapSidebarTab=summary&projectId=10928&showFilters=true&showStats=true&updatedAt=2024-05-23T14%3A49%3A45.735Z&workflow=short-read-mngs). This is a public dataset that CZID provides for training purposes - it contains simulated samples from a range of locations and tissue types, where metagenomics was performed to identify a pathogen. Have a look at the datasets and familiarise yourself with the interface.

## Quality control

**12.    Are there any samples that look different to the others? What might these be?**

**13.    How many raw reads were there for each sample?**

<details>
<summary>Clues</summary>
    
    Use the + button to the right of the column headings to add additional fields.
</details>

**14.    Roughly what proportion of the reads passed quality control? What about filtering? What do these metrics mean?**

**15. During which part of QC or filtering were most of the reads lost?**
<details>
<summary>Clues</summary>
    
    Use the bar chart button towards the top left of the screen to view some summary plots.
    ERCC reads refer to spike-in controls for quantification - you can ignore these for today!
    
</details>

## Interpretation

Click on the sample for Patient 009 to view the report in more detail.

When using CZID, you should generate a background model using your negative control samples. This takes some time to run so it has been done for you today. Choose the background called Medical_Detectives_v2_WCS at the top of the page.

**16. How can you filter the output to identify species that were present at higher abundances in the sample than in the negative control?**
Hint: look at the descriptions of the various scores [here](https://chanzuckerberg.zendesk.com/hc/en-us/articles/360034790574-Sample-Report-Table).
<details>
<summary>Clues</summary>
    
    Try filtering by Z score.
    
</details>

**17. What is the main viral infection that we might report for this patient? Are there any others?**
<details>
<summary>Clues</summary>
    
    Try using the filters to select just viruses and known pathogens.
    
</details>

**18. Visualise the coverage for the viruses you've found. What do you notice?**
<details>
<summary>Clues</summary>
    
    Hover next to the virus name and click coverage visualisation.
    
</details>

## Extension questions

**19. What other scores are shown on the CZID output? Which ones might be particularly useful?**

**20. Can you generate a consensus sequence for one of the viruses you've found using CZID?**

**21. Can you identify any phylogenetic relationships between viruses found in different samples in the same dataset?**

# Acknowledgements
This tutorial was orginally developed by Sarah Buddle for the Wellcome Connecting Science course in [Genomics for Clinical Virology](https://github.com/WCSCourses/GCV_2025). The second section of the tutorial is adapted from the training materials and documentation produced by CZID at [https://czid.org/](https://czid.org/). 






