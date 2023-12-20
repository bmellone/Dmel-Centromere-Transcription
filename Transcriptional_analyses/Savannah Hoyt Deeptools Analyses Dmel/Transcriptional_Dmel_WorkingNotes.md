## List of comparisons ##

1.	Compare copies of Jockey-3 in terms of sequence (divergence) and expression (PRO-seq/RNA-seq) with T-tests
    - Full-length vs. truncated (significant difference by length?)
    - Centromeric vs. non-centromeric (significance differences by location?)
    - consider using all three mapping methods (default, k100, k100-filt)

2.	PRO-seq density profiles for Jockey-3 
    - Full-length vs. truncated copies
    - Centromeric vs. non-centromeric (significance differences by location?)
    - consider using all three mapping methods (default, k100, k100-filt)

3. Parallel plots comparing length/divergence, # unique k-mers, and PRO-seq signal
    - with all mapping methods

Notes:
- use all three mapping methods (default, k100, k100-filt)
- look at all developmental stages/tissues

## Bigwig file prep: Ryan D. (Core Lab) way ##

Since the PRO-seq data is paired-end, we need to pull-out just R1 from the sort.bam as this points to the site of engaged RNA polymerase (R2 is really just to help with mapping).
- 0x2: PROPER_PAIR: each segment properly aligned according to the aligner
- $2==99, $2==147: reads come from the 5′->3′ template transcript (positive strand)
- $2==83, $2==163: reads come from the 3′->5′ template transcript (negative strand)

### 1) Start with concordantly aligned minus strand read pairs ###
```
samtools view -hf 0x2 ProSeq_embryo_Iso1_bowtie2_Df_sorted.bam |
  awk '($1 ~/^@/ || $2==99 || $2==147)' |
  samtools sort -n -@ 12 |
  bedtools bamtobed -bedpe -i stdin |
  awk '{OFS="\t"} {print $1,$2,$6,$7,$8,"-"}' > ProSeq_embryo_Iso1_bowtie2_Df_nameSort.bed
```
### 2) Add in concordantly aligned plus strand read pairs ###
```
samtools view -hf 0x2 ProSeq_embryo_Iso1_bowtie2_Df_sorted.bam |
  awk '($1 ~/^@/ || $2==83 || $2==163)' |
  samtools sort -n -@ 12 |
  bedtools bamtobed -bedpe -i stdin |
  awk '{OFS="\t"} {print $1,$2,$6,$7,$8,"+"}' >> ProSeq_embryo_Iso1_bowtie2_Df_nameSort.bed
```
- head of output (62,509,053 lines total): includes start of R1 through end of R2 
```
Contig114	5684	5715	A01204:43:H7CHMDRXY:2:2101:1063:24283	1	-
3L_1	19193141	19193186	A01204:43:H7CHMDRXY:2:2101:1063:24439	42	-
Contig56	13445	13554	A01204:43:H7CHMDRXY:2:2101:1063:25473	1	-
2L_1	15784505	15784634	A01204:43:H7CHMDRXY:2:2101:1063:26443	42	-
3R_28	2244712	2244771	A01204:43:H7CHMDRXY:2:2101:1063:26475	42	-
```
- tail of output 
```
X_1	16260584	16260600	A01204:43:H7CHMDRXY:2:2278:32904:36652	1	+
X_1	18664939	18665129	A01204:43:H7CHMDRXY:2:2278:32913:1438	42	+
X_1	18884762	18884932	A01204:43:H7CHMDRXY:2:2278:32913:2691	42	+
Contig38	153308	153323	A01204:43:H7CHMDRXY:2:2278:32913:3192	1	+
Contig5	438606	438665	A01204:43:H7CHMDRXY:2:2278:32913:3223	1	+
```

### 3) 3' ends of alignments ###
```
cat ProSeq_embryo_Iso1_bowtie2_Df_nameSort.bed |
  awk '{OFS="\t"} $6=="+" {print $1,($3-1),$3,$4,$5"_"($3-$2),$6}; $6=="-" {print $1,$2,($2+1),$4,$5"_"($3-$2),$6}' |
  sort -T /core/labs/Oneill/T2T_working_space/tmp/ -k1,1 -k2,2n > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_sort.bed
```
- head
```
2L_1	820	821	A01204:43:H7CHMDRXY:2:2205:19623:5165	1_115	-
2L_1	1108	1109	A01204:43:H7CHMDRXY:2:2273:24885:12273	1_71	-
2L_1	1133	1134	A01204:43:H7CHMDRXY:2:2136:20338:20682	1_40	+
2L_1	1175	1176	A01204:43:H7CHMDRXY:2:2244:21314:2331	1_86	-
2L_1	1355	1356	A01204:43:H7CHMDRXY:2:2104:16685:30201	1_17	-
```
- tail
```
Y_scaffold7	2069090	2069091	A01204:43:H7CHMDRXY:2:2172:3296:1939	1_46	+
Y_scaffold7	2069230	2069231	A01204:43:H7CHMDRXY:2:2264:26332:8609	1_22	-
Y_scaffold7	2069318	2069319	A01204:43:H7CHMDRXY:2:2173:4173:10285	1_30	-
```
### 4) make bedgraph of 3' ends (position of polymerase) ###
#### sort genome file first to avoid chaos downstream ####
``` 
export LC_COLLATE=C
sort -k1,1 -k2,2n contigs.size > contigs_sort.size
```
#### plus strand ####
```
awk '$6 == "+"' ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_sort.bed | genomeCoverageBed -3 -i /dev/stdin -bg -g contigs_sort.size > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bedgraph
```
- head
```
2L_1	1133	1134	1
2L_1	1674	1675	1
2L_1	3605	3606	1
```
#### minus strand ####
```
awk '$6 == "-"' ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_sort.bed | genomeCoverageBed -3 -i /dev/stdin -bg -g contigs_sort.size > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_m.bedgraph
```
- head
```
2L_1	820	821	1
2L_1	1108	1109	1
2L_1	1175	1176	1
```
#### negate minus strand ####
```
awk '{ $4=$4*-1; print }' ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_m.bedgraph > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bedgraph
```
- head 
```
    2L_1 820 821 -1
    2L_1 1108 1109 -1
    2L_1 1175 1176 -1
```
### 5) sort 3' end bedgraphs ###
```
export LC_COLLATE=C
sort -T /core/labs/Oneill/T2T_working_space/tmp/ -k1,1 -k2,2n ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bedgraph > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus_sort.bedgraph
sort -T /core/labs/Oneill/T2T_working_space/tmp/ -k1,1 -k2,2n ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bedgraph > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus_sort.bedgraph
(non-negated minus strand) sort -T /core/labs/Oneill/T2T_working_space/tmp/ -k1,1 -k2,2n ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_m.bedgraph > ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_m_sort.bedgraph
```
### 6) bedgraph to bigwig ###
```
module load GenomeBrowser/20180626
bedGraphToBigWig ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus_sort.bedgraph contigs_sort.size ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig
bedGraphToBigWig ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus_sort.bedgraph contigs_sort.size ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig
```

## Repeat regions file prep ##

### Jockey-3 only ###

- Jockey-3 in Dmel is labeled as such in the annotation file: Jockey-3_Dmel_08212020
Grab all Jockey-3 in Dmel:
```
grep "Jockey-3_Dmel_08212020" dmel_repeat_annotations.gff > Jockey-3_Dmel_08212020_dmel_repeat_annotations.gff
```
- Rearrange columns to be: chr/contig, start, end, divergence, strand, label
- Manually split the full-length (23) and truncated (306) copies based on the Excel file shared with me: Jockey-3_FL_Copies_Dmel.xlsx
- Edit label to represent whether each copy is full-length or truncated resulting in: FL_Jockey-3_Dmel_08212020 and Trunc_Jockey-3_Dmel_08212020
- Convert to tab-format, for ex:
```
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations.bed > FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab.bed
```
- Sort just to be sure for ex: 
```
export LC_COLLATE=C
sort -T /core/labs/Oneill/T2T_working_space/tmp/ -k1,1 -k2,2n FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab.bed > FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort.bed
```
- Combine and sort: 
```
cat FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab.bed Truncated-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab.bed > FL-Truncated-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab.bed
export LC_COLLATE=C
sort -k1,1 -k2,2n FL-Truncated-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab.bed > FL-Truncated-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort.bed
```

### NOT Jockey-3 ###

- All Dmel repeat annotations: dmel_repeat_annotations.gff (153523 lines) 
- Grab everything except Jockey-3 in Dmel: 
```
grep -v "Jockey-3_Dmel_08212020" dmel_repeat_annotations.gff > NOT-Jockey3_Dmel_08212020_dmel_repeat_annotations.gff
```
- Rearrange columns to be: chr/contig, start, end, divergence, strand, label: 
```
awk 'BEGIN {OFS="\t"}; {print $1,$4,$5,$6,$7,$3}' NOT-Jockey3_Dmel_08212020_dmel_repeat_annotations.gff > NOT-Jockey3_Dmel_08212020_dmel_repeat_annotations_short.gff
```
- Convert to tab-format: 
```
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' NOT-Jockey3_Dmel_08212020_dmel_repeat_annotations_short.gff > NOT-Jockey3_Dmel_08212020_dmel_repeat_annotations_short.bed
```

### Combine Jockey-3 & NOT Jockey-3 = all repeats ###

- combine: 
```
cat FL-Truncated-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort.bed NOT-Jockey3_Dmel_08212020_dmel_repeat_annotations_short.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels.bed
```
- sort:
```
export LC_COLLATE=C
sort -k1,1 -k2,2n ALLrepeats_dmel_repeat_annotations_withJockeyLabels.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_sort.bed
```
### Add CEN labels to repeats file ###

- CEN file: CENcontigs_Dmel.bed
```
tig00057289     1       24561   CEN
3R_5    1       103827  CEN
Contig119       1       93914   CEN
Contig79        1       70181   CEN
Y_Contig26      1       139957  CEN
```
- Print "CEN" to last column if the repeat overlaps the region
```    
bedtools map -c 4 -o collapse -a ALLrepeats_dmel_repeat_annotations_withJockeyLabels_sort.bed -g contigs_sort.size -b CENcontigs_Dmel_sort.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel.bed
```
This resulted in 533 repeats labeled with "CEN" and 152,990 not. 

### Create CEN vs. nonCEN Jockey3 region files ###

- Regardless of FL or Truncated
```
grep "CEN" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly.bed
grep -v "CEN" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN.bed
grep "Jockey-3_Dmel_08212020" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only.bed
grep "Jockey-3_Dmel_08212020" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only.bed
awk '{OFS="\t"} {print $1,$2,$3,$5,$6}' ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short.bed
awk '{OFS="\t"} {print $1,$2,$3,$5,$6}' ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short.bed
grep "+" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_PLUS.bed
grep -v "+" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_MINUS.bed
grep "+" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_PLUS.bed
grep -v "+" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_MINUS.bed
```
- Split into FL and Truncated 
```
grep "Trunc" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_MINUS.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_MINUS_TruncOnly.bed
grep "Trunc" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_PLUS.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_PLUS_TruncOnly.bed
grep -v "Trunc" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_PLUS.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_PLUS_FLonly.bed
grep -v "Trunc" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_MINUS.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_MINUS_FLonly.bed
grep -v "Trunc" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_PLUS.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_PLUS_FLonly.bed
grep -v "Trunc" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_MINUS.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_MINUS_FLonly.bed
grep "Trunc" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_MINUS.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_MINUS_TruncOnly.bed
grep "Trunc" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_PLUS.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_PLUS_TruncOnly.bed 
```

### Add kmer counts to repeats file ###

- add 1 to end of kmer file to be able to count num of unique 21mers: 
```
awk 'NF=NF+1{$NF="1"}1' FS="\t" OFS="\t" dmel_assembly_21mer_SINGLEv1.3.meryl.bed > dmel_assembly_21mer_SINGLEv1.3.meryl_edit.bed
```
- sort and count kmers per repeat
```
export LC_COLLATE=C
sort -T /core/labs/Oneill/T2T_working_space/tmp/ -k1,1 -k2,2n dmel_assembly_21mer_SINGLEv1.3.meryl_edit.bed > dmel_assembly_21mer_SINGLEv1.3.meryl_edit_sort.bed
sort -T /core/labs/Oneill/T2T_working_space/tmp/ -k1,1 -k2,2n dmel_assembly_51mer_SINGLEv1.3.meryl_edit.bed > dmel_assembly_51mer_SINGLEv1.3.meryl_edit_sort.bed
bedtools map -c 4 -o count -a ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel.bed -g contigs_sort.size -b /core/labs/Oneill/T2T_working_space/Dmel_analyses/kmers/dmel_assembly_21mer_SINGLEv1.3.meryl_edit_sort.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel_21merCount.bed
bedtools map -c 4 -o count -a ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel_21merCount.bed -g contigs_sort.size -b /core/labs/Oneill/T2T_working_space/Dmel_analyses/kmers/dmel_assembly_51mer_SINGLEv1.3.meryl_edit_sort.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel_21merCount_51merCount.bed
```
## Deeptools plot preparation ##

### Full-length vs. truncated Jockey-3 ###

- A full-length Jockey-3 is 4,304bp. There are 23 annotated in total. 
- Sort elements by length from top > bottom. 
- Run as 3' anchored since the truncated copies are predominantly 5' truncated.

#### Remove divergence (or other extra columns) from repeat files for deeptools ####
```
awk 'BEGIN {OFS="\t"}; {print $1,$2,$3,$5,$6}' FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort.bed > FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short.bed
```
#### Split FL and Truncated files into elements on the pos/neg strands for deeptools ####
```
grep "+" FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short.bed > FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short_POS.bed
grep -v "+" FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short.bed > FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short_NEG.bed
```
#### Convert unmerged Meryl bed file to bigwig for deeptools ####
```
genome=/core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/contigs_sort.size
#21mer example:
bedtools genomecov -bg -i dmel_assembly_21mer_SINGLEv1.3.meryl.bed -g ${genome} > dmel_assembly_21mer_SINGLEv1.3.meryl.bedgraph
export LC_COLLATE=C
sort -k1,1 -k2,2n dmel_assembly_21mer_SINGLEv1.3.meryl.bedgraph > dmel_assembly_21mer_SINGLEv1.3.meryl_sort.bedgraph
bedGraphToBigWig dmel_assembly_21mer_SINGLEv1.3.meryl_sort.bedgraph ${genome} dmel_assembly_21mer_SINGLEv1.3.meryl.bigwig
```

## Reads (R1) per Repeat ##
```
export LC_COLLATE=C
sort -T /core/labs/Oneill/T2T_working_space/tmp/ -k1,1 -k2,2n /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bam-Bed-Bedgraph/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_sort.bed > /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bam-Bed-Bedgraph/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_sort2.bed

bedtools coverage -counts -sorted -a ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel_21merCount_51merCount.bed -g contigs_sort.size -b /core/labsOneill/T2T_working_space/Dmel_analyses/Bam-Bed-Bedgraph/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_sort2.bed > ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel_21merCount_51merCount_COVcounts_PROseq-embryo-Iso1-bt2-df-R1singlent.bed

grep "Jockey-3_Dmel_08212020" ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENlabel_21merCount_51merCount_COVcounts_PROseq-embryo-Iso1-bt2-df-R1singlent.bed > JOCKEY3_dmel_repeat_annotations_withJockeyLabels_CENlabel_21merCount_51merCount_COVcounts_PROseq-embryo-Iso1-bt2-df-R1singlent.bed
```

## Grab FL fastas of Jockey-3 for alignment (23 total) ##

- Wanted to see if elements were mostly 5' truncated vs. 3' truncated 

```
bedtools getfasta -fi File.S8.Chang_et_al.fasta -bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short.bed > /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/FL-Jockey-3_Dmel_08212020_dmel.fasta
```
