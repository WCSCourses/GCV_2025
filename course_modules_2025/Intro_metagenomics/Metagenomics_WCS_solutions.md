# Solutions

**1.	Write a command to trim adaptors and low quality regions from your data.**
<details>
<summary><b>Solution</b></summary>
<pre>
 trim_galore -q 15 --length 60 --paired ~/metagenomics/data/sample1_1.fq.gz ~/metagenomics/data/sample1_2.fq.gz
 trim_galore -q 15 --length 60 --paired ~/metagenomics/data/neg_control_1.fq.gz ~/metagenomics/data/neg_control_2.fq.gz
</pre>
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
</details>

**6. Write a command to run bracken on your kraken2 output.**
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
</details>

**7. What species are present in your sample?**
<details>
<summary><b>Solution</b></summary> 
Sample1 contains human mastadenovirus F and cytomegalovirus.
</details>

**8. Look at the negative control. What can you conclude about what might be causing the disease in the patient?**
<details>
<summary><b>Solution</b></summary>    
The negative control also contains ~5 reads of cytomegalovirus so this is probably a contaminant. Therefore, we would report only the adenovirus clinically.
</details>

**9-11.    Extension questions: ask if you need help!**

**12.    Are there any samples that look different to the others? What might these be?**

<details>
<summary><b>Solution</b></summary>   
H20_1 and H20_2 are negative controls.
</details>

**13.    Roughly how many raw reads were there for each sample?**

<details>
<summary><b>Solution</b></summary>   
Number of raw reads ranges from around 3 million to 150 million (excluding the negative control samples, which have much less).
</details>

**14.    Roughly what proportion of the reads passed quality control? What about filtering? What do these metrics mean?**

<details>
<summary><b>Solution</b></summary>   
During quality control (QC), low quality and complexity and short reads are removed. Filtering happens after QC and is when the human reads are removed. In this dataset, typically 50-90% of reads pass QC and less than 3% pass filtering due to the high human content of the samples.
</details>

**15. During which part of QC or filtering were most of the reads lost?**

<details>
<summary><b>Solution</b></summary>   
Most reads were lost during human filtering, followed by low quality filtering.
</details>

**16. How can you filter the output to identify species that were present at higher abundances in the sample than in the negative control?**

<details>
<summary><b>Solution</b></summary>   
Filter NT Z score >= 0.  (I suggest for the next question you use a Z score filter of 0.1, since filtering with a score of 0 does not always work as expected in this dataset)
</details>

**17. What is the main viral infection that we might report for this patient? Are there any others?**

<details>
<summary><b>Solution</b></summary>   
Chikungunya virus was found at high levels. Human mastadenovirus C, Rotavirus A and Human alphaherpesvirus 2 are found at lower levels.
</details>

**18. Visualise the coverage for the viruses you've found. What do you notice?**

<details>
<summary><b>Solution</b></summary>   
A complete genome with good depth is obtained for Chikungunya virus. For the other viruses, coverage is more patchy.
</details>

**19-20.    Extension questions: ask if you need help!**
