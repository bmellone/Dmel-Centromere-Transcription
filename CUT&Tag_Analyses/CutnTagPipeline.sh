#!/bin/bash
#SBATCH --job-name=cutntag
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=80G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=asna.amjad@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

module load java
module load bwa
module load samtools
module load deeptools
module load picard
module load TrimGalore/0.6.7 cutadapt htseq

##################  REFERENCE #################
REF=/home/FCAM/aamjad/assembly/Het_chromatin_assembly/dmel_scaffold2_plus0310.fasta

######Index #########
bwa index $REF 

################# REMOVE ADAPTERS AND MAP CUT&Tag DATA WITH BOWTIE2 ##############
for i in Iso1_CIDAM_1 Iso1_CIDAM_2 Iso1_CIDAM_3 Iso1_H3K9_1 Iso1_H3K9_3 Iso1_IgG
do

### Remove adapter
mkdir fastq_Trim
trim_galore --paired --nextera --length 35 --phred33 --retain_unpaired -o fastq_Trim "$i"_R1.fastq.gz "$i"_R2.fastq.gz

mkdir Mapping
bowtie2 -p 24 -t -x …/index -1 fastq_Trim/"$i"_R1_val_1.fq.gz -2 fastq_Trim/"$i"_R2_val_2.fq.gz | samtools view -bS -@ 24 -f 0x2 | samtools sort -n -m 8G -O bam -o Mapping/"$i"_Bt2_PairsOnly_NameSort.bam

bedtools bamtobed -bedpe -i Mapping/"$i"_R2_BT2_PairsOnly_NameSort.bam > Mapping/"$i"_R2_BT2_PairsOnly_NameSort.bed

### Remove Duplicate
java -Xmx10g -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar MarkDuplicates \
          REMOVE_DUPLICATES=TRUE \
          I=Mapping/"$i"_Bt2_PairsOnly_NameSort.bam \
          O=Mapping/"$i"_Bt2_PairsOnly_NameSort_noDup.bam \
          M=Mapping/"$i"_Bt2_PairsOnly_NameSort_noDup.txt

### Filtering : keeping uniquely mapping read
samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 2048 -q 30 Mapping/"$i"_TrimGalor_bwa_noDup.bam > Mapping/"$i"_TrimGalor_bwa_noDup_q30.bam
rm Mapping/"$i"_TrimGalor_bwa_noDup.bam
samtools index  Mapping/"$i"_TrimGalor_bwa_noDup_q30.bam


##### Number of read mapping for normalization
COUNT=$(samtools view -c -F 260 "$i"_fil.bam)
NORM=$(echo "scale=5;(1000000/$COUNT)"| bc -l)

mkdir Coverage
##### Coverage
bamCoverage -b Mapping/"$i"_TrimGalor_bwa_fil.bam -o Coverage/"$i"_TrimGalor_bwa_fil_RPM.bw --scaleFactor $NORM -bs 1 --extendReads
bamCoverage -b Mapping/"$i"_TrimGalor_bwa_noDup_q30.bam -o Coverage/"$i"_TrimGalor_bwa_noDup_q30_RPM.bw --scaleFactor $NORM -bs 1 --extendReads


#Need bedgraph for SEACR
bamCoverage -b Mapping/"$i"_TrimGalor_bwa_fil.bam -o Coverage/"$i"_TrimGalor_bwa_fil_RPM.bedgraph --scaleFactor $NORM -bs 1 --extendReads  -of bedgraph
bamCoverage -b Mapping/"$i"_TrimGalor_bwa_noDup_q30.bam -o Coverage/"$i"_TrimGalor_bwa_noDup_q30_RPM.bedgraph --scaleFactor $NORM -bs 1 --extendReads -of bedgraph

### Insert size
java -Xmx10g -jar /isg/shared/apps/picard/picard-tools-2.9.2/picard.jar CollectInsertSizeMetrics \
      I=Mapping/"$i"_TrimGalor_bwa_noDup_q30.bam \
      O=Mapping/"$i"_TrimGalor_bwa_noDup_q30.txt \
      H=Mapping/"$i"_TrimGalor_bwa_noDup_q30_insert_size_histogram.pdf \
      M=0.5
      
done
############################################################
#################### Peak Calling #######################
############################################################

#mkdir SEACR_Peak

#for i in OR_CID_AM_2_S28 OR_CID_AM_1_S27 OR_CID_aff_2_S26 OR_CID_aff_1_S25 OR_H3K27_S23
#do

#bash /scratch/alarracu_lab/Cecile/SEACR-master/SEACR_1.3.sh Coverage/"$i"_TrimGalor_bwa_fil_RPM.bedgraph Coverage/OR_IgG_S24_TrimGalor_bwa_fil_RPM.bedgraph norm stringent SEACR_Peak/"$i"_TrimGalor_bwa_fil_RPM
#bash /scratch/alarracu_lab/Cecile/SEACR-master/SEACR_1.3.sh Coverage/"$i"_TrimGalor_bwa_noDup_q30_RPM.bedgraph Coverage/OR_IgG_S24_TrimGalor_bwa_noDup_q30_RPM.bedgraph norm stringent SEACR_Peak/"$i"_TrimGalor_bwa_noDup_q30_RPM

#done 
