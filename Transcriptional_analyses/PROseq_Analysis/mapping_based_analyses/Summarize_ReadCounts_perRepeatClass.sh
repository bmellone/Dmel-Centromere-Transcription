## Savannah Hoyt Code

#!/bin/bash
#SBATCH --job-name=summarizeReads
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=20G
#SBATCH --mail-user=username@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#Define input file
INfile=SAMPLE1_PROseq_COV_repeats.bed
#Pull out each repeat class into new file:
grep "LINE" ${INfile} > LINE_${INfile}
grep "SINE" ${INfile} > SINE_${INfile}
grep "LTR" ${INfile} > LTR_${INfile}
grep "Retroposon" ${INfile} > Retroposon_${INfile}
grep "DNA" ${INfile} > DNA_${INfile}
grep "Satellite" ${INfile} > Satellite_${INfile}
grep "Simple_repeat" ${INfile} > SR_${INfile}
grep "Low_complexity" ${INfile} > LC_${INfile}
grep "RC" ${INfile} > RC_${INfile}
grep "scRNA" ${INfile} > scRNA_${INfile}
grep "srpRNA" ${INfile} > srpRNA_${INfile}
grep "snRNA" ${INfile} > snRNA_${INfile}
grep "tRNA" ${INfile} > tRNA2_${INfile}
#Because some SINE families have "tRNA" in the name, need to remove those
grep -v "SINE" tRNA2_${INfile} > tRNA_${INfile}
#And then delete the original tRNA file
rm tRNA2_${INfile}
grep "rRNA" ${INfile} > rRNA_${INfile}
grep "Unknown" ${INfile} > Unknown_${INfile}
#Print number of lines per repeat class file to get total repeats per class:
wc -l LINE_${INfile}
wc -l SINE_${INfile}
wc -l LTR_${INfile}
wc -l Retroposon_${INfile}
wc -l DNA_${INfile}
wc -l Satellite_${INfile}
wc -l SR_${INfile}
wc -l LC_${INfile}
wc -l RC_${INfile}
wc -l scRNA_${INfile}
wc -l srpRNA_${INfile}
wc -l snRNA_${INfile}
wc -l tRNA_${INfile}
wc -l rRNA_${INfile}
wc -l Unknown_${INfile}
#Sum number of reads overlapping repeat class, and write to out
awk '{SUM+=$7}END{print SUM}' LINE_${INfile}
awk '{SUM+=$7}END{print SUM}' SINE_${INfile}
awk '{SUM+=$7}END{print SUM}' LTR_${INfile}
awk '{SUM+=$7}END{print SUM}' Retroposon_${INfile}
awk '{SUM+=$7}END{print SUM}' DNA_${INfile}
awk '{SUM+=$7}END{print SUM}' Satellite_${INfile}
awk '{SUM+=$7}END{print SUM}' SR_${INfile}
awk '{SUM+=$7}END{print SUM}' LC_${INfile}
awk '{SUM+=$7}END{print SUM}' RC_${INfile}
awk '{SUM+=$7}END{print SUM}' scRNA_${INfile}
awk '{SUM+=$7}END{print SUM}' srpRNA_${INfile}
awk '{SUM+=$7}END{print SUM}' snRNA_${INfile}
awk '{SUM+=$7}END{print SUM}' tRNA_${INfile}
awk '{SUM+=$7}END{print SUM}' rRNA_${INfile}
awk '{SUM+=$7}END{print SUM}' Unknown_${INfile}
