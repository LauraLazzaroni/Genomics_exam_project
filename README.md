## Genomics_exam_project
This repository contains a README file along with the scripts used to analyse the trios of Lazzaroni Laura and Piovesan Lara. Additionally, all results are provided for each trio, in particular, the cases in which the child is predicted to be affected.


# Trio-Based Exome Analysis for the Identification of Inherited and De Novo Rare Mendelian Disorders

Authors: Lazzaroni Laura and Piovesan Lara

## Abstract
Whole exome sequencing (WES) combined with trio analysis is a powerful clinical approach for identifying causal variants in rare Mendelian disorders. This study evaluated paired-end exome sequencing data from family trios to diagnose suspected rare diseases mapped to chromosome 20. A bioinformatics pipeline was implemented to assess both autosomal dominant and autosomal recessive inheritance patterns.

## Prerequisites
Tools and versions utilized at the time of the analysis:
* FastQC v0.11.9
* Bowtie2 v2.4.4
* Samtools v1.13
* Qualimap v.2.3
* MultiQC v1.34
* Freebayes v1.3.6
* BCFtools v1.13
* Bedtools genomecov v2.30.0

## Input Data
Common files required to perform each trio analysis, that must be stored in the common parent folder `trio_analyses`:
1. **Raw sequencing data (exome):** In FASTQ format (.fq), zipped (.gz) or unzipped.
2. **Reference genome:** Indexed and Bowtie2-indexed (.fa, .fai, .bt2).
3. **Target exonic regions:** Annotations in .bed format.
4. **Samples file:** Plain text file to maintain the order of the BAM files while executing the script.
5. **`analysis_info.txt`:** Fundamental plain text file (provided inside the repository) containing all the information necessary to run the script.

### Format of `analysis_info.txt`
It consists of 7 columns separated by spaces:
* `seq_mode`: `SE` for single-end or `PE` for paired-end.
* `trio_id`: e.g. `trio_1`.
* `child_id`, `father_id`, `mother_id`: e.g. `HG00421`.
* `inheritance`: `AR` for autosomic recessive, `AD_denovo` for autosomic dominant de novo, or `AD_inherited` for autosomic dominant inherited mutation from one of the parents (the affected parent must be specified in the `notes` column).
* `notes`: To be compiled only in case of `AD_inherited` (`father_affected` or `mother_affected`).

*Note: DO NOT add a new empty line at the end of this file, after the column names.*

## Directory Structure
Before running the script, ensure your files are organized as follows:
```text
trio_analyses/
├── analysis_info.txt
└── script.sh
```
After running the script, the files should be organized as follows:
```text
trio_analyses/
├── chr20.fa
├── chr20.fa.fai
├── chr20.1.bt2 (and other bowtie indices)
├── target_regions.bed
├── analysis_info.txt
├── samples.txt
├── script.sh
├── trio_1/
│   ├── HG00421.targets_R1.fq.gz
│   ├── ...
│   └── trio_1_export/
│       └── ...
└── trio_2/ ...
```

## Usage
How to launch the script:
1. Make the script executable: `chmod +x script.sh`
2. Launch the script: `./script.sh`

## Output
At the end of each trio analysis, all the necessary files to visualize the coverage tracks on UCSC Genome Browser are automatically copied inside an export folder (e.g. `trio_1_export`).
This folder will contain:
* `trio_1_multiqc_report.html`
* `trio_1.cand.vcf`
* `trio_1_fatherCov.bg`, `trio_1_motherCov.bg`, `trio_1_childCov.bg`

After exporting and downloading the files locally, remember to delete the export folder.

Enjoy!
