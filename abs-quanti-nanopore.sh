#!/bin/bash

# This script is to process original passed (Q7) nanopore reads from mClover3-spiked samples. Taxonomy classification, ARG annotation, ARG-host classification, as well as potential pathogen identification. This script also subsample raw read file based on the read generation times to compare the community profiling, ARG profile, ARG host profiling, and potential pathogen profileing between different data sizes to the theoretical dataset (all data generated during the entire sequencing run). 

# Written by YANG, Yu, last update on Nov 30 2021.

# Workflow:
# 1. Prepare sequencing reads: 
#	a. convert fastq to fasta;
#	b. calculate read lengths for all reads;
#	c. filter out reads shorter than 1kb;
#	d. convert the read generation time to seconds (a. convert each specific time points to second for each day, and add 24*60*60=86400 secs to the second day and 86400*2 to the third day and so on);
#	e. subsample by generation time for readIDs (3min 5min 10min 15min 30min 1hr 2hr 4hr 8hr 16hr 24hr and All data);
#	f. compile a table with the following information:
#	readID|read length|read generation time (sec) (table1)

# 2. Kraken2 using the mClover-included GTDB r95 database
#	a. compile a table with the following information:
#	taxID|sumofbases|AGS|lineage of the assigned taxID|cellnumber (table2)

# 3. mClover3 alignment by minimap2 
# 	a. align mClover3 reads by minimap2;
#	b. filter alignment results by 75% Identity and min 150bp aligned bases;
#	c. calculate total number of bases aligned to mClover3 and divided by 720bp for a final number of spike cell number approximation;
#	d. readIDs for reads with filtered mClover3 alignment;
#	e. subsample mClover3 minimap2 results according to the generation time of each read algined to mClover3
#	f. compile a table with the following info:
#	readID|# of bases aligned to mClover3|taxID (table3)

# 4. ARG identification by Minimap2 against nucleotide ARG database (min alignment length 200bp, min identity 80% --NM preterm paper)
#	a. align reads to nucleotide ARG database by minimap2;
#	b. filter minimap2 results to keep only a) primary alignments b) min 200bp aligned bases to reference ARG reads c) min identity of 80%;
#	c. calculate the number of ARG bases for each read and keep those with at least addtional 1kb walkout distance for ARG-carrying read taxonomy classification for ARG host tracking
#	d. compile a table with the following information:
#	read ID|num of ARGs|read length|number of bases assigned to ARGs|additional number of bases besides ARGs|taxID (table4)

# 5. Compile a final mothertable containing the following information for each read: 
#	taxID|number of reads|sum of bases|AGS|cell number seq'd|ARG annotation|ARG counts|ARG annotation|potential pathogen by pathogen list|kraken2 LCA classification|lineage of the taxonomy classification

# Other tools used:
# seqtk, seqkit, kraken2, minimap2

# Paths needed:
# Kraken2_gtdb_db: *your Kraken2-compatible GTDB index database files*
# mClover3 fasta file: /fasta/mClover3.fa 
# nucleotide ARG database and the structure file (source: https://github.com/xiaole99/ARGs-OAP-v2.0-development): nucleotide-ARG-DB.fasta, ARG_structure
# structure Avg Genome Size (AGS) database: GTDB_r95_AGS_DB file constructed 
# Nanopore DNA CS fasta file:  /fasta/DCS.fasta
# filter_fasta_by_list_of_headers.py script (source: https://bioinformatics.stackexchange.com/a/3940)


##############################################
# Subsample and kraken2 results taxonomy files
##############################################
#merge nanopore reads from the quality passed fastq folder
cat fastq_pass/*gz > merged.fastq.gz
#convert fq.gz file to fastq file
gunzip merged.fastq.gz > merged.fastq
for file in *.fastq
	#convert fastq to fasta
do seqtk seq -a ${file} > ${file}.fasta;
	#calculate sequence length, compile seq generation time info, subsample seqID according to sequence generation time
	seqkit fx2tab -l ${file}.fasta |awk -F" " '{print $1"\t"$NF"\t"$5}'|sed 's/start_time=//g;s/T/\t/g;s/Z//g' > ${file}.fasta_seqID_length_time.txt;
	cut -d$'\t' -f 3 ${file}.fasta_seqID_length_time.txt|sort|uniq|cat -n|awk -F"\t" '{print $NF"\t"$(NF-1)}'|tr -d " ">${file}.fasta_seqID_length_dates.txt; 
	cut -f 3 ${file}.fasta_seqID_length_time.txt|while read line; do grep -w "${line}" ${file}.fasta_seqID_length_dates.txt ;done|paste ${file}.fasta_seqID_length_time.txt - >${file}.fasta_seqID_length_time.txt_mod;
	cut -f 4 ${file}.fasta_seqID_length_time.txt_mod|awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3}' |paste ${file}.fasta_seqID_length_time.txt_mod - |awk -F"\t" '{print $NF+86400*($(NF-1)-1)}'|paste ${file}.fasta_seqID_length_time.txt_mod - |awk -F"\t" '{print $1"\t"$2"\t"$NF}'|sort -k1 >${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1;
	start_sec=$(cut -f 3 ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort -t$'\t' -n|head -1);
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+3*60)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>1_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_3min_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+5*60)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>2_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_5min_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+10*60)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>3_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_10min_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+15*60)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>4_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_15min_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+30*60)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>5_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_30min_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+60*60)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>6_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_1hr_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+60*60*2)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>7_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_2hr_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+60*60*4)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>8_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_4hr_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+60*60*8)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>9a_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_8hr_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+60*60*16)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>9b_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_16hr_seqID;
	awk -v ini="$start_sec" -F"\t" '$3<=(ini+60*60*24)  {print $1}' ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>9c_${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_24hr_seqID;
	unset start_sec;
	#filter out sequences shorter than 1kb
	seqkit seq -m 1000 ${file}.fasta > ${file}.fasta_1kb_withDCS.fa;
	#filter out DCS reads by blastn with min 85% similarity, min 90% subLcounts, and max 3600bp length (with this cutoff the final % DCS bases is around 1% witch is the dosage percentage, 10ng in 1ug)
	minimap2 -cx map-ont /fasta/DCS.fasta ${file}.fasta_1kb_withDCS.fa > ${file}.fasta_1kb_withDCS_minimapDCS.paf; #path to the DCS fasta file can be modified
	awk '{print 100*($10/$11)}'  ${file}.fasta_1kb_withDCS_minimapDCS.paf|paste ${file}.fasta_1kb_withDCS_minimapDCS.paf - >${file}.fasta_1kb_withDCS_minimapDCS_similarity.paf;
	awk '{print 100*(sqrt(($9-$8)*($9-$8))/$7)}'  ${file}.fasta_1kb_withDCS_minimapDCS.paf|paste ${file}.fasta_1kb_withDCS_minimapDCS_similarity.paf - > ${file}.fasta_1kb_withDCS_minimapDCS_similarity_subLpar.paf;
	grep -w "tp:A:P" ${file}.fasta_1kb_withDCS_minimapDCS_similarity_subLpar.paf|awk -F"\t" '$(NF)>=90 && $(NF-1)>=85 && $2<=3600 {print}' |cut -f 1|sort|uniq > ${file}.fasta_1kb_withDCS_minimapDCS_DCSseqID;
	scripts/filter_fasta_by_list_of_headers.py  ${file}.fasta_1kb_withDCS.fa ${file}.fasta_1kb_withDCS_minimapDCS_DCSseqID > ${file}.fasta_1kb.fa #path to the python script can be modified
	
	#taxonomy classification by kraken2 using GTDB supplemented with mClover genome
	kraken2 --db database/kraken2_gtdb_r95_20200803_mcloveradded/fasta_by_list_of_headers.py ${file}.fasta_1kb.fa  --threads 15 --output ${file}.fasta_1kb_kraken2_gtdb_r95_new --use-names --report ${file}.fasta_1kb_kraken2_report_gtdbr95_new --unclassified-out ${file}.fasta_1kb_kraken2_gtdb_r95_new_unclassified --classified-out ${file}.fasta_1kb_kraken2_gtdb_r95_new_classified; #kraken database path can be replaced here
	seqkit fx2tab -l ${file}.fasta_1kb_kraken2_gtdb_r95_new_classified  |awk -F" " '{print $1"\t"$NF"\t"$(NF-2)}'|sed 's/kraken:taxid|//g'|awk -F"\t" '$3="tax_"$3 {print $1"\t"$2"\t"$3}'|sort -k1 > ${file}.fasta_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid
done

for file in *seqid_length_taxid
	#compile seqID, length, and taxID info
do cut -f 1 ${file}|sort|uniq>${file}_seqid;
	cut -f 3 ${file} |sort|uniq>${file}_uniqtaxids;
	cat ${file}_uniqtaxids |while read line; do grep -w "${line}" ${file} |awk '{ SUM += $2} END{print SUM}'; done | paste ${file}_uniqtaxids - >${file}_uniqtaxids_sumofbases;
	awk -F"\t" '{print $1}'  ${file}_uniqtaxids_sumofbases|while read line; do grep -w "${line}" GTDB_r95_AGS_DB ;done|paste ${file}_uniqtaxids_sumofbases - |awk -F"\t" '{print $1"\t"$4"\t"$2"\t"$5"\t"$6}'>${file}_uniqtaxids_sumofbases_AGS_lineage; #path to the SAGS database can be modified
	awk -F"\t" '{printf("%.15f\n", $3/$4)}' ${file}_uniqtaxids_sumofbases_AGS_lineage |paste ${file}_uniqtaxids_sumofbases_AGS_lineage - |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}'>${file}_uniqtaxids_sumofbases_AGS_lineage_cellnumeber_table2;
	sed 's/;/\t/g' ${file}_uniqtaxids_sumofbases_AGS_lineage_cellnumeber_table2|awk -F"\t" '{$6="d__"$6} {$7="p__"$7} {$8="c__"$8} {$9="o__"$9} {$10="f__"$10} {$11="g__"$11} {$12="s__"$12} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}'|sed 's/[[:alpha:]]__\(\t\|$\)/\t/g'|sed 's/\-/_/g'  > ${file}_uniqtaxids_sumofbases_AGS_lineage_cellnumeber_table2_full
done

# Taxonomy classification results at different taxonomic resolutions
for file in *.fastq.fasta_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid_uniqtaxids_sumofbases_AGS_lineage_cellnumeber_table2_full
do grep "s__" ${file} |sort -k 12n|awk -F"\t" '{print $12}'|sort|uniq > ${file}_uniqspecies;
	grep "g__" ${file} |sort -k 11n|awk -F"\t" '{print $11}'|sort|uniq > ${file}_uniqgenues;
	grep "f__" ${file} |sort -k 10n|awk -F"\t" '{print $10}'|sort|uniq > ${file}_uniqfamilies;
	grep "o__" ${file} |sort -k 9n|awk -F"\t" '{print $9}'|sort|uniq > ${file}_uniqorder;
	grep "c__" ${file} |sort -k 8n|awk -F"\t" '{print $8}'|sort|uniq > ${file}_uniqclass;
	grep "p__" ${file} |sort -k 7n|awk -F"\t" '{print $7}'|sort|uniq > ${file}_uniqphylum;
	grep "d__" ${file} |sort -k 6n|awk -F"\t" '{print $6}'|sort|uniq > ${file}_uniqdomain;
	cat ${file}_uniq* > ${file}_alltaxa;
	cat ${file}_alltaxa|while read line; do grep -w "${line}" ${file} |awk -F"\t" '{SUM+=$3} END {print SUM}'; done | paste ${file}_alltaxa - |sort -k1 -t$'\t' >${file}_alltaxa_sumofbases;
	sed 's/\-/_/g' GTDB_r95_AGS_DB|sort -t$'\t' -k2 >GTDB;#path to the SAGS database can be modified
	awk -F"\t" '{print $1}' ${file}_alltaxa_sumofbases |while read line; do grep -w "${line}" GTDB |head -1;done|paste ${file}_alltaxa_sumofbases - > ${file}_alltaxa_sumofbases_AGS_lineage1;# path to the structured GTDB AGS can be modified
	awk -F"\t" '{print $0, ($1==$4?_:"NOT") "MATCHED"}' ${file}_alltaxa_sumofbases_AGS_lineage1 |awk '{print $NF}'|sort|uniq; # this is to check if all the taxa have been correctly grep (presence of - symbol in grep may cause trouble, better to if all the grep results are correct)
	awk -F"\t" '{print $1"\t"$3"\t"$5"\t"$2"\t"$2/$5"\t"$NF}' ${file}_alltaxa_sumofbases_AGS_lineage1 > ${file}_alltaxa_AGS_sumofbases_cellnumber_lineage;
	rm -rf ${file}_alltaxa_sumofbases_AGS_lineage1
done

###################################################################################
#Minimap2 to find mClover3 genes, identity 75% and min aligned bases to subject 150bps
###################################################################################
for file in *.fasta_1kb.fa
do minimap2 -cx map-ont /fasta/mClover3.fa  ${file} > ${file}_minimap_mClover3_algn.paf; #path to the mClover3 fasta file can be replaced here
        awk '{print 100*($10/$11)}'  ${file}_minimap_mClover3_algn.paf |paste ${file}_minimap_mClover3_algn.paf - >${file}_minimap_mClover3_algn.paf_similarity;
        awk '{print sqrt(($9-$8)*($9-$8))}' ${file}_minimap_mClover3_algn.paf_similarity |paste ${file}_minimap_mClover3_algn.paf_similarity - >${file}_minimap_mClover3_algn.paf_similarity_subLcounts;
        grep -w "tp:A:P" ${file}_minimap_mClover3_algn.paf_similarity_subLcounts|awk -F"\t" '$(NF-1)>=75 && $NF>=150 {print}' >${file}_minimap_mClover3_algn.paf_similarity_subLcounts_primary_S75subL150;
	cut -f 1 ${file}_minimap_mClover3_algn.paf_similarity_subLcounts_primary_S75subL150 |sort|uniq > ${file}_minimap_mClover3_algn.paf_similarity_subLcounts_primary_S75subL150_seqID;
	awk -F"\t" '{print $1"\t"$NF}' ${file}_minimap_mClover3_algn.paf_similarity_subLcounts_primary_S75subL150 > ${file}_minimap_mClover3_algn.paf_similarity_subLcounts_primary_S75subL150_seqID_subLcounts_table3
done

#########################################################################################
#Minimap2 to find ARGs (multidrug removed), min alignment length 200bp, min identity 80%
#########################################################################################
for file in *.fasta_1kb.fa
do minimap2  -cx map-ont nucleotide-ARG-DB.fasta ${file} > ${file}_minimap2_modifiednanoARG_algn.paf; #path to the nucleotide ARG database can be replaced here
	awk '{print 100*($10/$11)}'  ${file}_minimap2_modifiednanoARG_algn.paf |paste ${file}_minimap2_modifiednanoARG_algn.paf - > ${file}_minimap2_modifiednanoARG_algn.paf_similarity;
	awk '{print sqrt(($9-$8)*($9-$8))}' ${file}_minimap2_modifiednanoARG_algn.paf |paste ${file}_minimap2_modifiednanoARG_algn.paf_similarity - > ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts;
	grep -w "tp:A:P" ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts|awk -F"\t" '$(NF-1)>=80 && $NF>=200 {print}' >${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200;
	cut -f 1 ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200|sort|uniq >${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID;
	awk -F"\t" '{print $1"\t"$2"\t"$6"\t"$7"\t"$NF}' ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200|sort -k3|uniq > ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_ARGaccession_subLcounts;
	awk -F"\t" '{print $1"\t"$2}' ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200|sort|uniq > ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length;
	join -t$'\t' -1 3 -2 1 -o 1.1 1.2 1.4 1.5 1.3 2.2 2.3 ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_ARGaccession_subLcounts ARG_structure|sed 's/rifampin/rifamycin/g'|uniq > ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation; #path to the ARG structure file can be modified
	cat ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID|while read line; do grep -w "${line}" ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_ARGaccession_subLcounts|awk -F"\t" '{sum+=$4} END {print sum}' ;done|paste ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID - >${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_ARGbases;
	cat ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID|while read line; do grep -w "${line}" ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_ARGaccession_subLcounts|wc -l ;done|paste ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID - >${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_ARGcounts;
	cat ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID|while read line; do grep -w "${line}" ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200|cut -d$'\t' -f 6|paste -sd ';' ;done|paste ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_ARGcounts - >${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_ARGcounts_ARGaccessions;
	cat ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID|while read line; do grep -w "${line}" ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation|cut -d$'\t' -f 5 |paste -sd ';' ;done|paste ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_ARGcounts_ARGaccessions - >${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_uniqseqID_ARGcounts_ARGaccessions_ARGTypeannotations;
	cat ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID|while read line; do grep -w "${line}" ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation|cut -d$'\t' -f 6 |paste -sd ';' ;done|paste ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_uniqseqID_ARGcounts_ARGaccessions_ARGTypeannotations - >${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_uniqseqID_ARGcounts_ARGaccessions_ARGTypeannotations_ARGsubtypeanno;
	join -t$'\t' -j 1 -o 1.1 1.2 2.2 ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_ARGbases|join -t$'\t' -j 1 -o 1.1 1.2 1.3 2.2 2.3 2.4 2.5 - ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_uniqseqID_ARGcounts_ARGaccessions_ARGTypeannotations_ARGsubtypeanno|awk -F"\t" '{print $1"\t"$4"\t"$2"\t"$3"\t"($2-$3)"\t"$5"\t"$6"\t"$7}' |sort -k1 > ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_uniqseqID_numofARGs_length_ARGbases_extrabases_ARGaccessions_ARGannotation_table4;
	awk -F"\t" '$5>=1000 {print}' ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_uniqseqID_numofARGs_length_ARGbases_extrabases_ARGaccessions_ARGannotation_table4|cut -f 1 |sort|uniq > ${file}_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_uniqseqID_numofARGs_length_ARGbases_extrabases_ARGaccessions_ARGannotation_table4_seqIDwith1kbextra
done

##Count ARG numbers after excluding Escherichia reads###
for file in *.fastq.fasta
do grep -Fwf escherichia_taxID ${file}_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid |cut -f 1 |sort|uniq>escherichia_readID; #path to the escherichia taxID can be modified
	grep -Fvwf escherichia_readID ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation > ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia;
	awk -F"\t" '{print $5"\t"$6"\t"$NF"\t"($4/$3)"\t"$1}' ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia > ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_mod;
	cut -f 5 -d$'\t' ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_mod > ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_mod_seqID;
	cut -f 2 -d$'\t' ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_mod|sort|uniq>${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_uniqARGTypes;
	cut -f 3 -d$'\t' ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_mod|sort|uniq>${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_uniqARGSubtypes;
	cat ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_uniqARGTypes|while read line; do grep -w "${line}" ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_mod |awk -F"\t" '{sum+=$4} END {print sum}';done|paste ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_uniqARGTypes - >${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_uniqARGTypes_counts;
	cat ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_uniqARGSubtypes|while read line; do grep -w "${line}" ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_mod |awk -F"\t" '{sum+=$4} END {print sum}';done|paste ${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_uniqARGSubtypes - >${file}_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_seqID_length_subLcounts_ARGaccession_ARGannotation_noEscherichia_uniqARGSubtypes_counts
done

###########################
# Compile final mothertable
###########################

#for ARG-carrying reads, combine the taxnomy classification and ARG identification results. Count number of ARGs for different taxonomy.
for file in *.fastq
do join -t$'\t' -j 1 -o 1.1 1.2 2.3 1.3 ${file}.fasta_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid ${file}.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1|sort|uniq>tmp1;
	join -t$'\t' -j 1 -o 1.1 1.2 1.3 1.4 2.2 2.4 2.5 2.6 2.7 2.8 tmp1 ${file}.fasta_1kb.fa_minimap2_modifiednanoARG_algn.paf_similarity_subLcounts_primary_S80subL200_uniqseqID_numofARGs_length_ARGbases_extrabases_ARGaccessions_ARGannotation_table4|sort -t$'\t' -k4 -n|uniq>tmp2;
	cut -f 4 -d$'\t' tmp2 |while read line; do grep -w "${line}" GTDB_r95_AGS_DB ;done|paste tmp2 - >ARG-carrying_readID_length_generationsec_taxID_numofARGs_ARGbases_extrabases_ARGaccession_ARGtype_ARGsubtype_taxID_taxon_AGS_lineage;#path to the SAGS database can be modified
	awk -F"\t" '$7>=1000 {print}' tmp2 > tmp2_with1kbextra; #only those reads with at least 2kb walkout distance will be used for ARG host identification
	cut -f 4 -d$'\t' tmp2_with1kbextra |sort|uniq>tmp2_with1kbextra_uniqtaxID;
	cat tmp2_with1kbextra_uniqtaxID |while read line; do grep -w "${line}" tmp2_with1kbextra |awk -F"\t" '{sum+=$5} END {print sum}';done|paste tmp2_with1kbextra_uniqtaxID - >tmp2_with1kbextra_uniqtaxID_ARGnums;
	cat tmp2_with1kbextra_uniqtaxID|while read line; do grep -w "${line}" tmp2_with1kbextra |cut -d$'\t' -f 8|paste -sd ';' ;done|paste tmp2_with1kbextra_uniqtaxID_ARGnums - >tmp2_with1kbextra_uniqtaxID_ARGnums_ARGaccessions;
	cat tmp2_with1kbextra_uniqtaxID|while read line; do grep -w "${line}" tmp2_with1kbextra |cut -d$'\t' -f 9|paste -sd ';' ;done|paste tmp2_with1kbextra_uniqtaxID_ARGnums_ARGaccessions - >tmp2_with1kbextra_uniqtaxID_ARGnums_ARGaccessions_ARGTypes;
	cat tmp2_with1kbextra_uniqtaxID|while read line; do grep -w "${line}" tmp2_with1kbextra |cut -d$'\t' -f 10|paste -sd ';' ;done|paste tmp2_with1kbextra_uniqtaxID_ARGnums_ARGaccessions_ARGTypes - >tmp2_with1kbextra_uniqtaxID_ARGnums_ARGaccessions_ARGTypes_ARGsubtypes;
	cut -f 1 -d$'\t' tmp2_with1kbextra_uniqtaxID_ARGnums_ARGaccessions_ARGTypes_ARGsubtypes|while read line; do grep -w "${line}" GTDB_r95_AGS_DB ;done|paste tmp2_with1kbextra_uniqtaxID_ARGnums_ARGaccessions_ARGTypes_ARGsubtypes - |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$9}' > ARG-carrying_taxID_ARGnums_ARGIDs_ARGTypes_ARGsubtypes_taxon_lineage #path to the structured AGS database can be modified
done

#################
#Results analysis
#################
for file in *.fastq.fasta
do num_of_mClover3_seqd=$(awk -F"\t" '{sum+=$2} END {print sum/720}' ${file}_1kb.fa_minimap_mClover3_algn.paf_similarity_subLcounts_primary_S75subL150_seqID_subLcounts_table3);
        scaling_factor=$(awk -v mClover3=$num_of_mClover3_seqd -v conc=${1} -v vol=${2}  'BEGIN {print conc*vol/mClover3}' );
        all_1kb_bases=$(grep -v ">" ${file}_1kb.fa  | wc | awk '{print $3-$1}');
        all_bases=$(grep -v ">" ${file} | wc | awk '{print $3-$1}');
        all_classified=$(grep -w -v "tax_1" ${file}_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid|awk -F"\t" '{sum+=$2} END {print sum}');
        all_classidied_species=$(grep -w -v "tax_1" ${file}_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid_uniqtaxids_sumofbases_AGS_lineage |grep "s__"|awk -F"\t" '{sum += $3} END {print sum}');
	awk -F"\t" -v factor=$scaling_factor -v s_vol=${3} '{print (($5*factor)/s_vol)}' ${file}_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid_uniqtaxids_sumofbases_AGS_lineage_cellnumeber_table2_full_alltaxa_AGS_sumofbases_cellnumber_lineage |paste ${file}_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid_uniqtaxids_sumofbases_AGS_lineage_cellnumeber_table2_full_alltaxa_AGS_sumofbases_cellnumber_lineage - |awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"($7*1000)"\t"$6}' >${file}_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid_uniqtaxids_sumofbases_AGS_lineage_cellnumeber_table2_full_alltaxa_AGS_sumofbases_cellnumber_lineage_actualcellnumberperLsample;
        echo "Total number of mClover3 is $num_of_mClover3_seqd";
        echo "Scaling factor is $scaling_factor";
        echo "Sum of all bases is $all_bases, sum of all min 1kb bases is $all_1kb_bases, sum of all classified bases (excl. root) is $all_classified, sum of all bases classified to species is $all_classidied_species"
done

##################################################################################################################
#Potential pathogen identification by pathogen list (538 NCBI naming converted to GTDB naming) and pathogenic ARBs
##################################################################################################################
for file in *.fastq.fasta_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid_uniqtaxids_sumofbases_AGS_lineage_cellnumeber_table2_full_alltaxa_AGS_sumofbases_cellnumber_lineage
do sed 's/\-/_/g' ${file}_actualcellnumberperLsample > ${file}_actualcellnumberperLsample_mod;
        sed 's/\-/_/g' /files/foresight_gtdb_1307 > GTDB_foresight_pathogens; #path to the pathogen list can be modified
        grep -Fwf GTDB_foresight_pathogens ${file}_actualcellnumberperLsample_mod |sort -t$'\t' -k2> ${file}_actualcellnumberperLsample_mod_pathogens;
        join -t$'\t' -1 2 -2 1 -o 1.2 1.1 1.3 1.4 1.5 1.6 2.2 2.3 2.4 2.5 2.7 ${file}_actualcellnumberperLsample_mod_pathogens ARG-carrying_taxID_ARGnums_ARGIDs_ARGTypes_ARGsubtypes_taxon_lineage > ARG-carrying_taxID_ARGnums_ARGIDs_ARGTypes_ARGsubtypes_taxon_lineage_pathogens
done

#######################################################################################
#subsample taxonomy classification, and mClover3 alignment results based on subsampled seqID
#######################################################################################
for file in *.fastq.fasta_seqID_length_time.txt_mod_seqID_length_timeinsec_table1_*
do grep -Fwf ${file} *fasta_1kb_kraken2_gtdb_r95_new_classified.seqid_length_taxid_seqid > ${file}_seqid_length_taxid_classified;
        grep -Fwf ${file} *fasta_1kb.fa_minimap_mClover3_algn.paf_similarity_subLcounts_primary_S75subL150 >${file}_minimap2_mClover3_S75subL150
done


