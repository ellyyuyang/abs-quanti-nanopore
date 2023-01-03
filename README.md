# abs-quanti-nanopore
This repo contains the key codes and logic for generating the absolute quantification results in the paper "Rapid Absolute Quantification of Pathogens and ARGs by Nanopore Sequencing" by Yang, Yu et al. 2021. <br>
* This paper is published [here](https://www.sciencedirect.com/science/article/pii/S0048969721072661)
* Raw sequencing data files are available under [BioProject: PRJNA728386](https://dataview.ncbi.nlm.nih.gov/object/PRJNA728386)
* Intended use for research only

***
## Components for reproducible analysis

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Construction of the `Structured Average Genome Size (SAGS)` Database](#1-construction-of-the-structured-average-genome-size-sags-database)
- [End-to-End `Absolute Quantification` workflow](#2-end-to-end-absolute-quantification-workflow)
- [Citations](#if-you-intend-to-use-these-commands-please-cite-these-resources)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## 1. Construction of the `Structured Average Genome Size (SAGS)` Database 
* Tools used: <br>
  * [taxonkit](https://github.com/shenwei356/taxonkit) <br>
  * [MetaPhlAn2: merge_metaphlan_tables.py](https://github.com/biobakery/MetaPhlAn2/tree/master/utils) <br>
* SAGS is built upon the bacterial and archaeal taxonomy and metadata files from [`GTDB_r95`](https://data.gtdb.ecogenomic.org/releases/release95/)

## 2. End-to-End `Absolute Quantification` workflow

#### **A) Tools used:** <br>
  * [seqtk](https://github.com/lh3/seqtk) <br>
  * [seqkit](https://github.com/shenwei356/seqkit) <br>
  * [Kraken2](https://github.com/DerrickWood/kraken2) <br>
  * In case of Illumina metagenomic shotgun reads, [Braken2](https://github.com/jenniferlu717/Bracken) <br>
  * [Minimap2](https://github.com/lh3/minimap2) <br>
  * [filter_fasta_by_list_of_headers.py](https://bioinformatics.stackexchange.com/a/3940) <br>

#### **B) Additional files besides original sequence files required:** (*files bracketed by * should be provided by users*): <br>
  * **Kraken2_gtdb_db**: `*your Kraken2-compatible GTDB index database files*` <br>
  * ***mClover3* fasta file**: `./fasta/mClover3.fa` <br>
  * [**nucleotide ARG database and the structure file**](https://github.com/xiaole99/ARGs-OAP-v2.0-development): `*nucleotide-ARG-DB.fasta*` & `*ARG_structure*` <br>
  * **Structured Avg Genome Size (AGS) database**: `*SAGS*` constructed as [above](#1-construction-of-the-structured-average-genome-size-sags-database) <br>
  * **Nanopore DNA CS fasta file**:  `./fasta/DCS.fasta` <br>
  * **Pathogen list**: `*pathogen.list*` [original list](https://webarchive.nationalarchives.gov.uk/ukgwa/20121212135622/http://www.bis.gov.uk/assets/bispartners/foresight/docs/infectious-diseases/t16.pdf)   <br> *Please refer to our manuscript for details of the conversion to GTDB taxonomy nomenclature* <br>
  * *Original data can be obtained upon request*
#### **C) Logic flow and key codes:**
* Prepare sequencing reads
  * merge reads, convert file types, length filtering by [seqtk](https://github.com/lh3/seqtk) and [seqkit](https://github.com/shenwei356/seqkit); <br>
  * identify ([Minimap2](https://github.com/lh3/minimap2)) and remove DCS reads if DCS is used in ONT library preparation
```
seqtk seq -a input.fq > input.fa
seqkit fx2tab -l input.fa
seqkit seq -m 1000 input.fa > input_1kb.fa
minimap2 -cx map-ont ./fasta/DCS.fasta input.fasta > output_DCS_minimap.paf
```
* Kraken2 for rapid taxonomic classification using GTDB r95 database <br>
  * compile and stratify taxonomic abundance results into different taxonomic resolutions at the number of bases and the number of genome copy levels 
```
kraken2 --db Kraken2_gtdb_db input_1kb.fa  --output kraken2_gtdb_r95 --use-names --report kraken2_report_gtdbr95 --unclassified-out kraken2_gtdb_r95_unclassified --classified-out kraken2_gtdb_r95_classified
```
  
* Spiked marker gene alignment by minimap2 
  * Identify *mClover3* reads by [Minimap2](https://github.com/lh3/minimap2) and filter results with parameters described in our [paper](https://www.sciencedirect.com/science/article/pii/S0048969721072661); <br>
  * calculate *mClover3* gene copy number for a final number of spike cell genome copy number approximation; <br>
```
minimap2 -cx map-ont ./fasta/mClover3.fa input_1kb.fa > minimap_mClover3_algn.paf
```
* ARG identification by Minimap2 against nucleotide ARG database  <br>
  * align reads to nucleotide ARG database by [Minimap2](https://github.com/lh3/minimap2) and filter results with cutoffs from ([here](https://www.nature.com/articles/s41564-019-0626-z))
  * calculate the gene copy number of different ARGs
  * keep those ARG-carrying reads with at least addtional 1kb walkout distance for ARG host tracking <br>
```
minimap2  -cx map-ont nucleotide-ARG-DB.fasta input_1kb.fa > minimap2_ARG_algn.paf
```
* Calculation of the absolute abundance of microbial cells in unit sample volumn <br>
  * refer to our [paper](https://www.sciencedirect.com/science/article/pii/S0048969721072661) for the calculation of scaling factor for converting seqenced genome copy number into cell number per unit sample volumn
  * absolute abundance of pathogens and ARG-carrying hosts can then be extracted
 
#### **D) Running time estimation for major steps:** 
For an input fasta file with size `10 Gb`, an approximated `3-4 hr` data processing time would be expected to generate the final microbial absolute quantification results.
* `Kraken2` for taxonomic classification -- `30 min` with 10 threads and 300 G memory pre-allocated.
* `Minimap2` for *mClover3* (spiked gene) identification -- `2.5 min` with 10 threads and 150 G memory pre-allocated.
* `Processing kraken2 output` to convert the sequenced genome copy numbers to the final absolute cell abundance per unit sample volumn: 
	* Summing bases for all the classified reads to different Kraken2-assigned LCA taxonomic lineages -- `2.5 hr` with 10 threads and 150 G memory pre-allocated.
	* Stratifying the summation results above into different taxonomic levels -- `5 min` with 10 threads and 150 G memory pre-allocated.
	* Convert the sequenced genome copy numbers into the asbolute cell abundance per unit sample volumn -- untimed, but approx. `15 min` with single thread.
		

## If you intend to use these commands, please cite these resources:
[GTDB](https://gtdb.ecogenomic.org/) <br>
[Kraken2](https://github.com/DerrickWood/kraken2) <br>
In case of Illumina metagenomic shotgun reads, [Braken2](https://github.com/jenniferlu717/Bracken) <br>
[Nucleotide ARG database](https://github.com/xiaole99/ARGs-OAP-v2.0-development) <br>
[Minimap2 ](https://github.com/lh3/minimap2) <br>
[taxonkit](https://github.com/shenwei356/taxonkit) <br>
[seqtk](https://github.com/lh3/seqtk) <br>
[seqkit](https://github.com/shenwei356/seqkit) <br>
[MetaPhlAn: merge_metaphlan_tables.py](https://github.com/biobakery/MetaPhlAn/tree/master/metaphlan/utils) <br>
[Pathogen list](https://webarchive.nationalarchives.gov.uk/ukgwa/20121212135622/http://www.bis.gov.uk/assets/bispartners/foresight/docs/infectious-diseases/t16.pdf) <br>

*I try hard to credit all the third-party resources/tools/codes. If any unintentional infringements, please contact elly.yu.yang@gmail.com.*
