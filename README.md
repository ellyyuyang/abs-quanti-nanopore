# abs-quanti-nanopore
This repo contains the commands to reproduce the absolute quantification results in the paper "Rapid Absolute Quantification of Pathogens and ARGs by Nanopore Sequencing" by Yang, Yu et al. 2021. <br>
* This paper is published [here](https://www.sciencedirect.com/science/article/pii/S0048969721072661)
* Raw sequencing data files are available under [BioProject: PRJNA728386](https://dataview.ncbi.nlm.nih.gov/object/PRJNA728386)
* Intended use for research only
***
## Table of Contents

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Construction of the Kraken2-compatible `GTDB_r95` Database](#1-construction-of-the-kraken2-compatible-gtdb-r95-index-database-files)
- [Construction of the `Structured Average Genome Size (SAGS)` Database](#2-construction-of-the-structured-average-genome-size-sags-database)
- [End-to-End `Absolute Quantification` workflow](#3-end-to-end-absolute-quantification-workflow)
- [Citations](#if-you-intend-to-use-these-commands-please-cite-these-resources)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## Components for reproducible analysis
#### 1. Construction of the `Kraken2-compatible GTDB-r95` index database files
***Note:*** <br>
***Different versions of the GTDB database can be used following the same logic below.*** <br>
***Pay attention to the absolute location of the files used under your local device. Modifications to the absolute file paths may apply.***
* Tools used <br>
  * [Metagenomics-Index-Correction: tax_from_gtdb.py](https://github.com/rrwick/Metagenomics-Index-Correction) <br>
  * [Kraken2](https://github.com/DerrickWood/kraken2) <br>

* Download the [GTDB database](https://data.gtdb.ecogenomic.org/releases/release95/95.0/)
```
wget https://data.gtdb.ecogenomic.org/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz
```
* Download the NCBI Taxonomy
```
kraken2-build  --download-taxonomy --db kraken2_gtdb_r95_20200803_mcloveradded
```
* Convert NCBI taxonomy to GTDB taxonomy
```
tax_from_gtdb.py \
	 --gtdb release95/taxonomy/gtdb_taxonomy.tsv \
	 --assemblies release95/fastani/database/ \
	 --nodes taxonomy/nodes.dmp \
	 --names taxonomy/names.dmp \
 --kraken_dir kraken2_gtdb_r95_20200803_mcloveradded/kraken2_ready_ref
```
* Add GTDB ref to library
```
for file in kraken2_gtdb_r95_20200803_mcloveradded/kraken2_ready_ref/*.fa
do
   kraken2-build --add-to-library $file --db kraken2_gtdb_r95_20200803_mcloveradded/
done
```
* Follow the [official Kraken2 manual](https://github.com/DerrickWood/kraken2) for adding the complete genome of the spike-in into the Kraken2 database. *not neccessary*

#### 2. Construction of the `Structured Average Genome Size (SAGS)` Database 
***Note:*** <br>
***Different versions of the GTDB database can be used following the same logic below.*** <br>
***Pay attention to the absolute location of the files used under your local device. Modifications to the absolute file paths may apply.***
* Tools used: <br>
  * [taxonkit](https://github.com/shenwei356/taxonkit) <br>
  * [MetaPhlAn: merge_metaphlan_tables.py](https://github.com/biobakery/MetaPhlAn/tree/master/metaphlan/utils) <br>
* Download the bacterial and archaeal taxonomy and metadata files from [`GTDB_r95`](https://data.gtdb.ecogenomic.org/releases/release95/):
```
# Taxonomy
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_taxonomy_r95.tsv;
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_taxonomy_r95.tsv;

# Metadata
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz; tar -xf bac120_metadata_r95.tar.gz; rm -rf bac120_metadata_r95.tar.gz;
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/ar122_metadata_r95.tar.gz; tar -xf ar122_metadata_r95.tar.gz;rm -rf ar122_metadata_r95.tar.gz;
```
* Format the downloaded bacterial and archaeal taxonomy and metadata files
```
cat bac120_taxonomy_r95.tsv ar122_taxonomy_r95.tsv|sed 's/\-/_/g' > taxonomy.tsv;
sed '1d' ar122_metadata_r95.tsv >tmp1; sed '1d' bac120_metadata_r95.tsv>tmp2;cat tmp1 tmp2 > metadata.tsv; rm -rf tmp*;

# Select a few useful columns from metadata table and modify the metatable format, obtain NCBI taxid, NCBI species taxid, genome size, GTDB taxonomy
awk -F"\t" '{print $78"\t"$74"\t"$14"\t"$17}' metadata.tsv|sed 's/;/\t/g'|sed '1i ncbi_taxid\tncbi_species_taxid\tgenome_size\tgtdb_taxonomy_k\tgtdb_taxonomy_p\tgtdb_taxonomy_c\tgtdb_taxonomy_o\tgtdb_taxonomy_f\tgtdb_taxonomy_g\tgtdb_taxonomy_s'|sed 's/\-/_/g'> metadata_taxid_speciestaxid_gsize_taxa.tsv;
sed '1d' metadata_taxid_speciestaxid_gsize_taxa.tsv>metadata_taxid_speciestaxid_gsize_taxa_mod.tsv;
```
* Extract the taxonomy names at different taxonomy levels
```
awk -F"\t" '{print$5}' metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|sort|uniq>metadata_taxid_speciestaxid_gsize_taxa_mod_phylum.tsv;
awk -F"\t" '{print$6}' metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|sort|uniq>metadata_taxid_speciestaxid_gsize_taxa_mod_class.tsv;
awk -F"\t" '{print$7}' metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|sort|uniq>metadata_taxid_speciestaxid_gsize_taxa_mod_order.tsv;
awk -F"\t" '{print$8}' metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|sort|uniq>metadata_taxid_speciestaxid_gsize_taxa_mod_family.tsv;
awk -F"\t" '{print$9}' metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|sort|uniq>metadata_taxid_speciestaxid_gsize_taxa_mod_genus.tsv;
awk -F"\t" '{print$10}' metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|sort|uniq>metadata_taxid_speciestaxid_gsize_taxa_mod_species.tsv;
awk -F"\t" '{print$4}' metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|sort|uniq>metadata_taxid_speciestaxid_gsize_taxa_mod_kingdom.tsv;
```
* Calculate the average genome sizes at different taxonomy levels
```
cat metadata_taxid_speciestaxid_gsize_taxa_mod_phylum.tsv|while read line; do grep -w "${line}" metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|awk -F"\t" '{SUM += $3} END {printf("%.1f\n", SUM/NR)}' ; done | paste metadata_taxid_speciestaxid_gsize_taxa_mod_phylum.tsv  - > metadata_phylum_avggsize.txt;
cat metadata_taxid_speciestaxid_gsize_taxa_mod_class.tsv|while read line; do grep -w "${line}" metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|awk -F"\t" '{SUM += $3} END {printf("%.1f\n", SUM/NR)}' ; done | paste metadata_taxid_speciestaxid_gsize_taxa_mod_class.tsv  - > metadata_class_avggsize.txt;
cat metadata_taxid_speciestaxid_gsize_taxa_mod_order.tsv|while read line; do grep -w "${line}" metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|awk -F"\t" '{SUM += $3} END {printf("%.1f\n", SUM/NR)}' ; done | paste metadata_taxid_speciestaxid_gsize_taxa_mod_order.tsv  - > metadata_order_avggsize.txt;
cat metadata_taxid_speciestaxid_gsize_taxa_mod_family.tsv|while read line; do grep -w "${line}" metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|awk -F"\t" '{SUM += $3} END {printf("%.1f\n", SUM/NR)}' ; done | paste metadata_taxid_speciestaxid_gsize_taxa_mod_family.tsv  - > metadata_family_avggsize.txt;
cat metadata_taxid_speciestaxid_gsize_taxa_mod_genus.tsv|while read line; do grep -w "${line}" metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|awk -F"\t" '{SUM += $3} END {printf("%.1f\n", SUM/NR)}' ; done | paste metadata_taxid_speciestaxid_gsize_taxa_mod_genus.tsv  - > metadata_genus_avggsize.txt;
cat metadata_taxid_speciestaxid_gsize_taxa_mod_species.tsv|while read line; do grep -w "${line}" metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|awk -F"\t" '{SUM += $3} END {printf("%.1f\n", SUM/NR)}' ; done | paste metadata_taxid_speciestaxid_gsize_taxa_mod_species.tsv  - > metadata_species_avggsize.txt;
cat metadata_taxid_speciestaxid_gsize_taxa_mod_kingdom.tsv|while read line; do grep -w "${line}" metadata_taxid_speciestaxid_gsize_taxa_mod.tsv|awk -F"\t" '{SUM += $3} END {printf("%.1f\n", SUM/NR)}' ; done | paste metadata_taxid_speciestaxid_gsize_taxa_mod_kingdom.tsv  - > metadata_kingdom_avggsize.txt;
cat *_avggsize.txt>GTDB_r95_structured_genome_size.txt;
```
* Finalizing the database
```
# Combining lineage info from taxid to taxrank file
awk -F"[\t|\t]+" '{print $(NF-1)"\t"$1}' /kraken2_gtdb_r95_20200803_mcloveradded/taxonomy/nodes.dmp |sed '1d' >taxid_taxarank.txt; 
sed 's/\-/_/g' taxid_taxarank.txt >taxid_taxarank_mod.txt; 
awk '{print $NF}' taxid_taxarank.txt >taxid.txt;

# Convert taxid to lineage info and remove the first column  for a file with tab delimiter
taxonkit lineage --data-dir /kraken2_gtdb_r95_20200803_mcloveradded/taxonomy/ --show-status-code taxid.txt | tee lineage.withcode.txt|cut -d$'\t' -f2- >taxid_lineage.txt;

# Combine taxid and taxarank and genome size info
python merge_metaphlan_tables.py taxid_taxarank_mod.txt GTDB_r95_structured_genome_size.txt|sed '1d'|awk -F"\t" '{print $2"\t"$1"\t"$3}' >metadata_taxid_taxarank_avggsize.txt;
awk -F"\t" '{print $2"\t"$1}' taxid_taxarank_mod.txt >temp1; python merge_metaphlan_tables.py temp1 taxid_lineage.txt|sed '1d' |sed 's/\-/_/g'> taxid_taxrank_lineage;rm -rf temp1;

python merge_metaphlan_tables.py metadata_taxid_taxarank_avggsize.txt taxid_taxrank_lineage|awk -F"\t" '$1="tax_"$1 {print $1"\t"$2"\t"$3"\t"$5}'|sed '1d' |sed '1i tax_1\tRoot\t2995590\tRoot' > GTDB_r95_AGS_DB
```

#### 3. End-to-End `Absolute Quantification` workflow

**A) Tools used:** <br>
  * [seqtk](https://github.com/lh3/seqtk) <br>
  * [seqkit](https://github.com/shenwei356/seqkit) <br>
  * [Kraken2](https://github.com/DerrickWood/kraken2) <br>
  * In case of Illumina metagenomic shotgun reads, [Braken2](https://github.com/jenniferlu717/Bracken) <br>
  * [Minimap2](https://github.com/lh3/minimap2) <br>
  * [filter_fasta_by_list_of_headers.py](https://bioinformatics.stackexchange.com/a/3940) <br>

**B) Additional files besides original sequence files required:** (*files bracketed by * should be provided by users*): <br>
  * **Kraken2_gtdb_db**: `*your Kraken2-compatible GTDB index database files*` <br>
  * ***mClover3* fasta file**: `/fasta/mClover3.fa` <br>
  * [**nucleotide ARG database and the structure file**](https://github.com/xiaole99/ARGs-OAP-v2.0-development): `*nucleotide-ARG-DB.fasta*` & `*ARG_structure*` <br>
  * **Structure Avg Genome Size (AGS) database**: `*GTDB_r95_AGS_DB*` file constructed [above](#2-construction-of-the-structured-average-genome-size-sags-database) <br>
  * **Nanopore DNA CS fasta file**:  `/fasta/DCS.fasta` <br>
  * **Pathogen list** converted to GTDB taxonomy nomenclature: `/files/foresight_gtdb_1307`, (for [original list](https://webarchive.nationalarchives.gov.uk/ukgwa/20121212135622/http://www.bis.gov.uk/assets/bispartners/foresight/docs/infectious-diseases/t16.pdf))   *Please refer to the manuscript for details of the conversion to GTDB taxonomy nomenclature* <br>

**C) Logic flow and key codes for reproducibility:**
* Prepare sequencing reads
  * merge, convert file types, filter out reads shorter than 1kb by [seqtk](https://github.com/lh3/seqtk) and [seqkit](https://github.com/shenwei356/seqkit); <br>
  * identify ([Minimap2](https://github.com/lh3/minimap2)) and remove DCS reads (min 85% similarity, min 90% subLcounts, and max 3600bp length) if DCS is used in ONT library preparation
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
  * Identify and filter *mClover3* reads by [Minimap2 with min 75% similarity and min 150 bp aligned bases](https://github.com/lh3/minimap2); <br>
  * calculate *mClover3* gene copy number for a final number of spike cell genome copy number approximation; <br>
  * compile a table with the following info: <br>
  **readID|# of bases aligned to *mClover3*|taxID (`table2`)**
```
minimap2 -cx map-ont ./fasta/mClover3.fa input_1kb.fa > minimap_mClover3_algn.paf
```
* ARG identification by Minimap2 against nucleotide ARG database  <br>
  * align reads to nucleotide ARG database by [Minimap2](https://github.com/lh3/minimap2) and filter results ([min alignment length 200bp, min identity 80%](https://www.nature.com/articles/s41564-019-0626-z))
  * calculate the gene copy number of different ARGs
  * keep those ARG-carrying reads with at least addtional 1kb walkout distance for ARG host tracking <br>
```
minimap2  -cx map-ont nucleotide-ARG-DB.fasta input_1kb.fa > minimap2_ARG_algn.paf
```
* Calculation of the absolute abundance of microbial cells in unit sample volumn <br>
*The logic below also applies to the calculation of the absolute abundance of ARGs*
```
With inputs of the following positional parameters:
$1 = int, concentration of spike culture, eg. "3490000" in CFU/uL
$2 = int, volumn of the spike culture dosed, eg. "100" in uL
$3 = int, volumn of the sample to which the spike culture is added, eg. "50" in mL

for file in *.fastq.fasta
do num_of_mClover3_seqd=$(awk -F"\t" '{sum+=$2} END {print sum/720}' table2);
        scaling_factor=$(awk -v mClover3=$num_of_mClover3_seqd -v conc=${1} -v vol=${2}  'BEGIN {print conc*vol/mClover3}' );
	awk -F"\t" -v factor=$scaling_factor -v s_vol=${3} '{print (($5*factor)/s_vol)}' table1 |paste table1 - |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($7*1000)"\t"$6}' >table1_actualcellnumberperLsample;
        echo "Total number of mClover3 is $num_of_mClover3_seqd";
        echo "Scaling factor is $scaling_factor"
done
```
*Absolute abundance of pathogens and ARG-carrying hosts can then be extracted from the **table1_actualcellnumberperLsample** generated above*

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
