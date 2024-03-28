# Dmel-Centromere-Transcription

This repository is for all analyses related to Drosophila melanogaster centromere transcription datasets (PROseq, RNAseq, and CUT&Tag).
All analyses are mainly done on the University of Connecticut HPC cluster under SLURM software.
Other tools and codes used are listed here:

Each code block is uploaded as a basic set of commands, and was run on the University of Connecticut HPC cluster under the SLURM management software. This code block is included here:
#!/bin/bash 
#SBATCH --job-name=X
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 15
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=100g
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
This block was modified for memory and performance requirements accordingly. 

Meryl-associated code and unique kmer:
https://github.com/bmellone/Dmel-Centromere-Transcription/tree/main/Unique_Kmers

Drosophila heterochromatin assembly (PLoS Bio paper 2019)
Drosophila heterochromatin repeat annotations (PLoS Bio paper 2019)
