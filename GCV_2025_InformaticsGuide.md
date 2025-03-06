# <img src="https://coursesandconferences.wellcomeconnectingscience.org/wp-content/themes/wcc_courses_and_conferences/dist/assets/svg/logo.svg" width="300" height="50"> 

# Genomics and Clinical Virology 2025 Informatics Guide

# Mega Fix

```bash
cd ~/Downloads/

wget http://kr.archive.ubuntu.com/ubuntu/pool/universe/g/gconf/libgconf-2-4_3.2.6-6ubuntu1_amd64.deb

wget http://kr.archive.ubuntu.com/ubuntu/pool/universe/g/gconf/gconf2-common_3.2.6-6ubuntu1_all.deb

sudo dpkg -i gconf2-common_3.2.6-6ubuntu1_all.deb

sudo dpkg -i libgconf-2-4_3.2.6-6ubuntu1_amd64.deb
```

ModelTest Fix

```
cd ~/Downloads/

https://github.com/ddarriba/modeltest/files/6192913/modeltest-ng-0.1.7-static.tar.gz

tar -xvzf modeltest-ng-0.1.7-static.tar.gz

./modeltest-ng-static

```

## **Software Used During the Course**

| Software | Version (if not latest) | Module |
|-------------|--------------|----------|
| [bwa](https://github.com/lh3/bwa) | 0.7.17 | Reference alignment, coverage and stats |
| [bowtie2](https://github.com/BenLangmead/bowtie2) | 2.4.4 | Reference alignment, coverage and stats |
| [samtools](https://github.com/samtools/samtools) | 1.16 | Reference alignment, coverage and stats |
| [weeSAM](https://github.com/Centre-for-Infectious-Disease-Genomics-and-One-Health/weeSAM) |  | Reference alignment, coverage and stats |
| [tablet](https://ics.hutton.ac.uk/tablet/) | 1.17.08 | Reference alignment, coverage and stats |
| [ivar](https://github.com/andersen-lab/ivar) | 1.3.1 | Reference alignment, coverage and stats |
| [lofreq](https://csb5.github.io/lofreq/) | 2.1.5 | Reference alignment, coverage and stats |
| [qualimap](https://bitbucket.org/kokonech/qualimap/src/master/) | 2.2.2a-1 | Reference alignment, coverage and stats |
| [prinseq](https://github.com/GenomicsCore/prinseq) | 0.20.4 | Reference alignment, coverage and stats |
| [kraken2](https://github.com/DerrickWood/kraken2) | 2.1.3 | Introduction to metagenomics |
| [bracken](https://github.com/jenniferlu717/Bracken) | 3 | Introduction to metagenomics |
| [mafft](https://mafft.cbrc.jp/alignment/software/) | 7.525 | Phylogenetic Analysis |
| [MEGA-X](https://www.megasoftware.net/) |  | Phylogenetic Analysis |
| [Modeltest-ng](https://github.com/ddarriba/modeltest) | 0.1.7 | Phylogenetic Analysis |
| [IQ-TREE](https://github.com/iqtree/iqtree2) | 2.3.6 | Phylogenetic Analysis |
| [Figtree](https://github.com/rambaut/figtree) | 1.4.4 | Phylogenetic Analysis |
| [Entrez eutilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/) | 2.1.3 | Phylogenetic Analysis |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.12.1 | File Formats & QC |
| [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) | 6.1 | File Formats & QC |
| [unzip](https://linux.die.net/man/1/unzip) | 6 | Variant calling |
| [samtools](https://github.com/samtools/samtools) | 6.00 | Variant calling |
| [iva](https://github.com/sanger-pathogens/iva) | 1.0.11 | De novo |
| [khmer](https://github.com/dib-lab/khmer) | 3.0.0a3 | De novo |
| [numpy](https://github.com/numpy/numpy) | 1.5.1 | De novo |
| [pandas](https://github.com/pandas-dev/pandas) | 1.23.2 | De novo |
| [python](https://www.python.org/) | 3.8.12 | De novo |
| [quast](https://github.com/ablab/quast) | 5.3.0 | De novo |
| [Spades](https://github.com/ablab/spades) | 3.11.1 | De novo |
| [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | 0.39 | De novo |


## **To run the software you will have to activate an environment (bioinformatics_env and samtools_env), What is an Environment?**
An environment is an isolated space where specific software, dependencies, and libraries are installed. It ensures that all required tools run in a controlled and reproducible way, avoiding conflicts with other system applications.

## For this course, we use the **bioinformatics_env** and **samtools_env** environment.

⚠️ **Note:** Software will **not work outside the environment.**

## **How to Activate the bioinformatics_env Environment**
Before using any software, activate the environment with:

```bash
conda activate bioinformatics_env
```

## **How to Activate the samtools_env Environment (has Samtools, BWA, Prinseq, and Tablet installed)**
Before using any software, activate the environment with:

```bash
conda activate samtools_env
```

## **To deactivate the environment when you're done:**

```bash
conda deactivate
```

Link to [bioinformatics_env.yml](https://github.com/WCSCourses/GCV_2025/blob/main/course_data_2025/bioinformatics_env.yml)
Link to [samtools_env.yml](https://github.com/WCSCourses/GCV_2025/blob/main/course_data_2025/samtools_env.yml)

## Informatics Set-Up
We are currently using Oracle VM Virtual Box (https://www.virtualbox.org/) to deliver Informatics, you can find Virtual Box Guides below:
[Virtual Machine SetUp Guide for Intel-Mac and Windows](https://github.com/WCSCourses/index/blob/main/VM%20Guide.pdf). <br />

The Host OS Requirements for Virtual Box <br />
- RAM requirement: 8GB (preferably 12GB) <br />
- Processor requirement: 4 processors (preferably 8) <br />
- Hard disk space: 200GB <br />
- Admin rights to the computer <br />

Note: Please be aware that Virtual Box is currently incompatible with M1/M2/M3 chips on MacBook.
It is exclusively designed for use on Intel-based MacBooks. The current version in use is Virtual Box 7.0.

## Citing and Re-using Course Material

The course data are free to reuse and adapt with appropriate attribution. All course data in these repositories are licensed under the <a rel="license" href="https://creativecommons.org/licenses/by-nc-sa/4.0/">Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)</a>. <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /> 

## Interested in attending a course?

Take a look at what courses are coming up at [Wellcome Connecting Science Courses & Conference Website](https://coursesandconferences.wellcomeconnectingscience.org/our-events/).

---

[Wellcome Connecting Science GitHub Home Page](https://github.com/WCSCourses) 

For more information or queries, feel free to contact us via the [Wellcome Connecting Science website](https://coursesandconferences.wellcomeconnectingscience.org).<br /> 
Find us on socials [Wellcome Connecting Science Linktr](https://linktr.ee/eventswcs)

---
