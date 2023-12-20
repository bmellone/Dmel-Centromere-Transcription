## Savannah Hoyt Code

#!/bin/bash
#SBATCH --job-name=Summarize_Repeats
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=20G
#SBATCH --mail-user=usernamen@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

#Define input file
INfile=SAMPLE1_PROseq_COV_repeats.bed

#count repeats without read overlap
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' LINE_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' SINE_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' LTR_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' Retroposon_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' DNA_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' Satellite_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' SR_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' LC_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' RC_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' scRNA_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' srpRNA_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' snRNA_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' tRNA_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' rRNA_${INfile}
awk '{OFS="\t"}{count[$7]++}END{print count["0"]}' Unknown_${INfile}
