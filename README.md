# abs-quanti-nanopore
This repo contains the commands to reproduce the absolute quantification results in the paper "Rapid Absolute Quantification of Pathogens and ARGs by Nanopore Sequencing" by Yang, Yu et al. 2021. <br>
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
***Note:*** <br>
***Different versions of the GTDB database can be used following the same logic below.*** <br>
***Pay attention to the absolute location of the files used under your local device. Modifications to the absolute file paths may apply.***
* Tools used: <br>
  * [taxonkit](https://github.com/shenwei356/taxonkit) <br>
  * [MetaPhlAn: merge_metaphlan_tables.py](https://github.com/biobakery/MetaPhlAn/tree/master/metaphlan/utils) <br>
* SAGS is built upon the bacterial and archaeal taxonomy and metadata files from [`GTDB_r95`](https://data.gtdb.ecogenomic.org/releases/release95/)
* Recommended table heading: **`TaxID|TaxName|Average genome size|Lineage`**

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
  * **Structure Avg Genome Size (AGS) database**: `./files/GTDB_r95_AGS_DB` file constructed above <br>
  * **Nanopore DNA CS fasta file**:  `./fasta/DCS.fasta` <br>
  * **Pathogen list** converted to GTDB taxonomy nomenclature: `./files/pathogen.list`, (for [original list](https://webarchive.nationalarchives.gov.uk/ukgwa/20121212135622/http://www.bis.gov.uk/assets/bispartners/foresight/docs/infectious-diseases/t16.pdf))   <br> *Please refer to the manuscript for details of the conversion to GTDB taxonomy nomenclature* <br>

#### **C) Logic flow and key codes:**
* Prepare sequencing reads
  * merge reads, convert file types, length filtering by [seqtk](https://github.com/lh3/seqtk) and [seqkit](https://github.com/shenwei356/seqkit); <br>
  * identify ([Minimap2](https://github.com/lh3/minimap2)) and remove DCS reads if DCS is used in ONT library preparation
```
minimap2 -cx map-ont ./fasta/DCS.fasta input.fasta > output_DCS_minimap.paf
```
* Kraken2 for rapid taxonomic classification using GTDB r95 database <br>
  * Sum the total number of bases assigned to each taxonomic lineage and calculate the sequenced genome copy number by dividing the sum to the average genome size (prepared [above](#2-construction-of-the-structured-average-genome-size-sags-database)) of the respective taxonomic lineage
  * compile a table with the following information: <br>
  **taxName|taxID|average genome size|sumofbases|cellnumber|lineage of the assigned taxID (`table1`)**
```
kraken2 --db Kraken2_gtdb_db input_1kb.fa  --output kraken2_gtdb_r95 --use-names --report kraken2_report_gtdbr95 --unclassified-out kraken2_gtdb_r95_unclassified --classified-out kraken2_gtdb_r95_classified
```
  
* Spiked marker gene alignment by minimap2 
  * Identify *mClover3* reads by [Minimap2](https://github.com/lh3/minimap2) and filter results with parameters described in our [paper](https://www.sciencedirect.com/science/article/pii/S0048969721072661); <br>
  * calculate *mClover3* gene copy number for a final number of spike cell genome copy number approximation; <br>
  * compile a table with the following info: <br>
  **readID|# of bases aligned to *mClover3*|taxID (`table2`)**
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
***The logic below also applies to the calculation of the absolute abundance of ARGs***
```
With inputs of the following positional parameters:
$1 = int, concentration of spike culture, eg. "3490000" in CFU/uL
$2 = int, volumn of the spike culture dosed, eg. "100" in uL
$3 = int, volumn of the sample to which the spike culture is added, eg. "50" in mL

num_of_mClover3_seqd=$(awk -F"\t" '{sum+=$2} END {print sum/720}' table2);
scaling_factor=$(awk -v mClover3=$num_of_mClover3_seqd -v conc=${1} -v vol=${2}  'BEGIN {print conc*vol/mClover3}' );
awk -F"\t" -v factor=$scaling_factor -v s_vol=${3} '{print (($5*factor)/s_vol)}' table1 |paste table1 - |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($7*1000)"\t"$6}' >table1_actualcellnumberperLsample
```
  ***Note:*** <br>
  ***Absolute abundance of pathogens and ARG-carrying hosts can then be extracted from the **`table1_actualcellnumberperLsample`** generated above***

* **Running time estimation for major steps** <br>
	For an input fastq file with size `10 Gb`, an approximated `3-4 hr` data processing time would be expected to generate the final microbial absolute quantification results.
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
[Metagenomics-Index-Correction: tax_from_gtdb.py](https://github.com/rrwick/Metagenomics-Index-Correction) <br>
[MetaPhlAn: merge_metaphlan_tables.py](https://github.com/biobakery/MetaPhlAn/tree/master/metaphlan/utils) <br>
[filter_fasta_by_list_of_headers.py](https://bioinformatics.stackexchange.com/a/3940) <br>
[Pathogen list](https://webarchive.nationalarchives.gov.uk/ukgwa/20121212135622/http://www.bis.gov.uk/assets/bispartners/foresight/docs/infectious-diseases/t16.pdf) <br>

*I try hard to credit all the third-party resources/tools/codes. If any unintentional infringements, please contact elly.yu.yang@gmail.com.*
