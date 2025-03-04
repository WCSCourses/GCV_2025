# Solutions

**1.	Write a command using trim_galore to trim adaptors and low quality regions from your data.**
<details>
<summary><b>Solution</b></summary>
<pre>
 cd ~/metagenomics
 
 trim_galore -q 15 --length 60 --paired ~/metagenomics/data/sample1_1.fq.gz ~/metagenomics/data/sample1_2.fq.gz
 trim_galore -q 15 --length 60 --paired ~/metagenomics/data/neg_control_1.fq.gz ~/metagenomics/data/neg_control_2.fq.gz
</pre>
 
You might have chosen to use different values for -q and --length, but these are the ones used in the previous session and they will work well for our purposes today.
</details>

**2. Write commands to align the reads to the human genome.**
<details>
<summary><b>Solution</b></summary>
<pre>
# Index the human genome
bwa index ~/metagenomics/data/human_genome.fasta
</pre>
<pre>
# Align the reads to the human genome
bwa mem -t 4 ~/metagenomics/data/human_genome.fasta sample1_1_val_1.fq.gz sample1_2_val_2.fq.gz > ~/metagenomics/sample1.sam
bwa mem -t 4 ~/metagenomics/data/human_genome.fasta neg_control_1_val_1.fq.gz neg_control_2_val_2.fq.gz > ~/metagenomics/neg_control.sam
</pre>
</details>

**3. What is the output format of the alignment?**
<details>
<summary><b>Solution</b></summary>
.sam file
</details>

**4. How could we use our alignment to get the non-human reads?**
<details>
<summary><b>Solution</b></summary>
Extract the non aligned reads (in this case to a bam file).

Convert the resulting bam file back to fastq.
</details>

**5. Write a command to run kraken2 on your data.**
<details>
<summary><b>Solution</b></summary>
<pre>
kraken2 --db ~/metagenomics/kraken_db \
--paired ~/metagenomics/sample1_nonhuman_1.fq \
~/metagenomics/sample1_nonhuman_2.fq \
--report ~/metagenomics/sample1_kraken_report.txt \
> ~/metagenomics/sample1_kraken_readclassifications.txt
</pre>
<pre>
kraken2 --db ~/metagenomics/kraken_db \
--paired ~/metagenomics/neg_control_nonhuman_1.fq \
~/metagenomics/neg_control_nonhuman_2.fq \
--report ~/metagenomics/neg_control_kraken_report.txt \
> ~/metagenomics/neg_control_kraken_readclassifications.txt
</pre>
We provide paths to the reference sequence database and the paired-end input fastq files. We output both a summary report and a file that gives the classifications of each read individually.
</details>

**6.	What do the numbers in the kraken2 report mean?**
<details>
<summary><b>Solution</b></summary>
See https://github.com/DerrickWood/kraken2/wiki/Manual#sample-report-output-format.
</details>

**7. Write a command to run bracken on your kraken2 output.**
<details>
<summary><b>Solution</b></summary>
<pre>
bracken -d ~/metagenomics/kraken_db \
-i ~/metagenomics/sample1_kraken_report.txt \
-o ~/metagenomics/sample1_bracken.txt \
-t 3
</pre>
<pre>
bracken -d ~/metagenomics/kraken_db \
-i ~/metagenomics/neg_control_kraken_report.txt \
-o ~/metagenomics/neg_control_bracken.txt \
-t 3
</pre>
We provide the path to the reference database and the kraken2 report as input. We then provide an output file name and set the minimum number of reads required to perform reestimation at 3.
</details>

**8.	What do the numbers in the bracken report mean?**
<details>
<summary><b>Solution</b></summary>
See https://ccb.jhu.edu/software/bracken/index.shtml?t=manual
</details>

**9. What species are present in your sample?**
<details>
<summary><b>Solution</b></summary> 
Sample1 contains human mastadenovirus F and cytomegalovirus.
</details>

**10. Look at the negative control. What can you conclude about what might be causing the disease in the patient?**
<details>
<summary><b>Solution</b></summary>    
The negative control also contains ~5 reads of cytomegalovirus so this is probably a contaminant. Therefore, we would report only the adenovirus clinically.
</details>

**11.    How can you tell which reads were assigned to human mastadenovirus F?**
<details>
<summary><b>Solution</b></summary>
In the ~/metagenomics/sample1_kraken_readclassifications.txt  file, the third column gives the taxon ID of the species that read was assigned to. The second column gives the read ID, which can be found in the read header in the fastq file.
</details>

**12.    Choose a read that was assigned to adenovirus. How could you extract the entry corresponding to this read from the nonhuman fastq file?**
<details>
<summary><b>Solution</b></summary>
<pre>
grep 'A01897:100:HLTLTDRX3:2:2132:8223:9784' ~/metagenomics/sample1_nonhuman_1.fq -A 3 > ~/metagenomics/adenovirus_read.fastq
</pre>
This command searches for a read ID and also extracts the next three lines in the file.
You might have used a different read ID since there are multiple reads classified as adenovirus.
</details>

**13.    Use online blast to analyse this read. What do you notice?**
<details>
<summary><b>Solution</b></summary>
The top BLAST results are all for human adenovirus F (or 40 which is a type of adenovirus F) and the scores such as query cover and percentage identity are good. This gives us more confidence that this read does come from adenovirus and therefore that the virus is in our sample.
</details>

**14.	Using what you've learnt in previous sessions, what further analyses you think might be useful on these datasets?**
<details>
<summary><b>Solution</b></summary>
It could be useful to create genome coverage plots for the viruses we've identified. To do this, you would download a reference genome for the adenovirus and use what you learnt in the alignment session to create a bam file and visualise it.
</details>

**15.	How could you run the commands in this tutorial for multiple samples at a time? Write a suitable script.**
<details>
<summary><b>Solution</b></summary>
For example:
<pre>
for sample1 neg_control; do
    trim_galore -q 15 --length 60 --paired ~/metagenomics/data/${sample}_1.fq.gz ~/metagenomics/data/${sample}_2.fq.gz
done
</pre>
</details>

**16.	How could you adapt your answers to questions 11-12 to extract all the reads that were assigned to human mastadenovirus F?**
<details>
<summary><b>Solution</b></summary>
<pre>
# Select all the reads that were assigned to adenovirus (3rd column is equal to 130309) and print 2nd column (read ID) to a file
awk '$3==130309 {print $2}' ~/metagenomics/sample1_kraken_readclassifications.txt > ~/metagenomics/adenovirus_read_ids.txt
</pre>
<pre>
# Extract the read IDs from the fastq file
grep -F -f ~/metagenomics/adenovirus_read_ids.txt ~/metagenomics/sample1_nonhuman_1.fq -A 3 > ~/metagenomics/adenovirus_reads_1.fq
grep -F -f ~/metagenomics/adenovirus_read_ids.txt ~/metagenomics/sample1_nonhuman_2.fq -A 3 > ~/metagenomics/adenovirus_reads_2.fq
</pre>
</details>

**17.	How do Kraken2 and Bracken work?**
<details>
<summary><b>Solution</b></summary>
See published articles
</details>

**18.    Are there any samples that look different to the others? What might these be?**

<details>
<summary><b>Solution</b></summary>   
H20_1 and H20_2 are negative controls.
</details>

**19.    Roughly how many raw reads were there for each sample?**

<details>
<summary><b>Solution</b></summary>   
Number of raw reads ranges from around 3 million to 150 million (excluding the negative control samples, which have much less).
</details>

**20.    Roughly what proportion of the reads passed quality control? What about filtering? What do these metrics mean?**

<details>
<summary><b>Solution</b></summary>   
During quality control (QC), low quality and complexity and short reads are removed. Filtering happens after QC and is when the human reads are removed. In this dataset, typically 50-90% of reads pass QC and less than 3% pass filtering due to the high human content of the samples.
</details>

**21. During which part of QC or filtering were most of the reads lost?**

<details>
<summary><b>Solution</b></summary>   
Most reads were lost during human filtering, followed by low quality filtering.
</details>

**22. How can you filter the output to identify species that were present at higher abundances in the sample than in the negative control?**

<details>
<summary><b>Solution</b></summary>   
Filter NT Z score >= 0.  (I suggest for the next question you use a Z score filter of 0.1, since filtering with a score of 0 does not always work as expected in this dataset)
</details>

**23. What is the main viral infection that we might report for this patient? Are there any others?**

<details>
<summary><b>Solution</b></summary>   
Chikungunya virus was found at high levels. Human mastadenovirus C, Rotavirus A and Human alphaherpesvirus 2 are found at lower levels.
</details>

**24. Visualise the coverage for the viruses you've found. What do you notice?**

<details>
<summary><b>Solution</b></summary>   
A complete genome with good depth is obtained for Chikungunya virus. For the other viruses, coverage is more patchy.
</details>

**25-26.    Extension questions: ask if you need help!**
