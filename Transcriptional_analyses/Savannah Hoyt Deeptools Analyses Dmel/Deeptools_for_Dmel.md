# Deeptools for D. melanogaster #

- Load these modules: deeptools/3.5.0, python/3.8.1
- For file prep see: Transcriptional_Dmel_WorkingNotes.md

## For regions files ##
- Should be 5 columns (tab-delimited): 
```
chr start  end strand   label
```

## Deeptools code ##

### Plot PROseq signal across genes algined at TSS ###
- Also ran with: [--averageTypeSummaryPlot sum], [--averageTypeBins sum]
```
# plus strand regions
computeMatrix reference-point --referencePoint TSS -p 12 \
-R /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/File.S10.Chang_et_al_genes_short_plus.bed \
-S /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig \
-o Dmel_PROseq_PLUS_Genes_mycode.mat.gz \
-b 5000 -a 50000 -bs 100 \
--averageTypeBins max \
--missingDataAsZero
#relabel:
computeMatrixOperations relabel -m Dmel_PROseq_PLUS_Genes_mycode.mat.gz -o Dmel_PROseq_PLUS_Genes_mycode_rename.mat.gz --groupLabels "Genes"
# minus strand regions
computeMatrix reference-point --referencePoint TSS -p 12 \
-R /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/File.S10.Chang_et_al_genes_short_neg.bed \
-S /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig \
 -o Dmel_PROseq_MINUS_Genes_mycode.mat.gz \
 -b 5000 -a 50000 -bs 100 \
 --averageTypeBins max \
 --missingDataAsZero
 # relabel:
 computeMatrixOperations relabel -m Dmel_PROseq_MINUS_Genes_mycode.mat.gz -o Dmel_PROseq_MINUS_Genes_mycode_rename.mat.gz --groupLabels "Genes"
# combine plus and minus strand region matrices
computeMatrixOperations rbind -m Dmel_PROseq_PLUS_Genes_mycode_rename.mat.gz Dmel_PROseq_MINUS_Genes_mycode_rename.mat.gz -o Dmel_PROseq_OUT_Genes_mycode_rename.mat.gz
# create heatmap
plotHeatmap -m Dmel_PROseq_OUT_Genes_mycode_rename.mat.gz --sortUsing region_length \
-o Dmel_PROseq_OUT_Genes_mycode_SE.htmp.pdf \
--dpi 800 --plotType se \
--colorMap Blues Reds \
--averageTypeSummaryPlot mean \
--samplesLabel "Antisense" "Sense" \
--heatmapHeight 15 --heatmapWidth 4 \
--zMin 0 --yMin 0
```
### Plot unique k-mer density across Jockey3 aligned at TES ###
- Also ran with: unique 51mers
```
computeMatrix reference-point --referencePoint TES -p 12 \
-R /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/Truncated-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short.bed \
-S /core/labs/Oneill/T2T_working_space/Dmel_analyses/kmers/dmel_assembly_21mer_SINGLEv1.3.meryl.bigwig \
-o /core/labs/Oneill/T2T_working_space/Dmel_analyses/deeptools/FLvsTruncated-Jockey3_Dmel_08212020_Meryl.single.k21_OUT.mat.gz \
-b 5000 -a 200 -bs 20 \
--averageTypeBins max \
--missingDataAsZero
#relabel:
computeMatrixOperations relabel -m /core/labs/Oneill/T2T_working_space/Dmel_analyses/deeptools/FLvsTruncated-Jockey3_Dmel_08212020_Meryl.single.k21_OUT.mat.gz -o /core/labs/Oneill/T2T_working_space/Dmel_analyses/deeptools/FLvsTruncated-Jockey3_Dmel_08212020_Meryl.single.k21_OUT_rename.mat.gz --groupLabels "FL_Jockey-3_Dmel_08212020" "Trunc_Jockey-3_Dmel_08212020"
# create heatmap
plotHeatmap -m /core/labs/Oneill/T2T_working_space/Dmel_analyses/deeptools/FLvsTruncated-Jockey3_Dmel_08212020_Meryl.single.k21_OUT_rename.mat.gz --sortUsing region_length \
-o /core/labs/Oneill/T2T_working_space/Dmel_analyses/deeptools/FLvsTruncated-Jockey3_Dmel_08212020_Meryl.single.k21_Avg-3anchor-SE-Greys.htmp.pdf \
--dpi 800 --plotType se \
--colorMap Greys \
--colorList white,lightgrey \
--averageTypeSummaryPlot mean \
--samplesLabel "unique_21mer" \
--heatmapHeight 15 --heatmapWidth 4 \
--zMin 0 --yMin 0
```

### Plot PROseq signal across Jockey3 (FL vs. Trunc) aligned at TES ###
- Also ran with: [-b 5000 -a 100 -bs 10], [-b 4500 -a 200 -bs 20], [-b 4300 -a 200 -bs 20]
```
#plus strand regions
computeMatrix reference-point --referencePoint TES -p 12 \
-R /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short_POS.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/Truncated-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short_POS.bed \
-S /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig \
-o Dmel_PROseq-k100_PLUS_Jockey3-FLvsTrunc.mat.gz \
-b 5000 -a 200 -bs 20 \
--averageTypeBins max \
--missingDataAsZero
#relabel:
computeMatrixOperations relabel -m Dmel_PROseq-k100_PLUS_Jockey3-FLvsTrunc.mat.gz -o Dmel_PROseq-k100_PLUS_Jockey3-FLvsTrunc_rename.mat.gz --groupLabels "FL_Jockey-3_Dmel_08212020" "Trunc_Jockey-3_Dmel_08212020"
# minus strand regions
computeMatrix reference-point --referencePoint TES -p 12 \
-R /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/FL-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short_NEG.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/Truncated-Jockey-3_Dmel_08212020_dmel_repeat_annotations_tab_sort_short_NEG.bed \
-S /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig \
-o Dmel_PROseq-k100_MINUS_Jockey3-FLvsTrunc.mat.gz \
-b 5000 -a 200 -bs 20 \
--averageTypeBins max \
--missingDataAsZero
# relabel:
computeMatrixOperations relabel -m Dmel_PROseq-k100_MINUS_Jockey3-FLvsTrunc.mat.gz -o Dmel_PROseq-k100_MINUS_Jockey3-FLvsTrunc_rename.mat.gz --groupLabels "FL_Jockey-3_Dmel_08212020" "Trunc_Jockey-3_Dmel_08212020"
# combine plus and minus strand region matrices
computeMatrixOperations rbind -m Dmel_PROseq-k100_PLUS_Jockey3-FLvsTrunc_rename.mat.gz Dmel_PROseq-k100_MINUS_Jockey3-FLvsTrunc_rename.mat.gz \
-o Dmel_PROseq-k100_Jockey3-FLvsTrunc_rename_OUT.mat.gz
# create heatmap
plotHeatmap -m Dmel_PROseq-k100_Jockey3-FLvsTrunc_rename_OUT.mat.gz --sortUsing region_length \
-o Dmel_PROseq-k100_Jockey3-FLvsTrunc_rename_OUT_SE_b5000-a200-bs20.htmp.pdf \
--dpi 800 --plotType se \
--colorMap Blues Reds \
--averageTypeSummaryPlot mean \
--samplesLabel "Antisense" "Sense" \
--heatmapHeight 15 --heatmapWidth 4 \
--zMin 0 --yMin 0
```

### Plot PROseq signal across Jockey3 (CEN vs. NonCEN) aligned at TES ###
```
#plus strand regions
computeMatrix reference-point --referencePoint TES -p 12 \
-R /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_PLUS.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_PLUS.bed \
-S /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig \
-o Dmel_PROseq_PLUS_Jockey3-CENvsNot.mat.gz \
-b 5000 -a 200 -bs 20 \
--averageTypeBins max \
--missingDataAsZero
#relabel:
computeMatrixOperations relabel -m Dmel_PROseq_PLUS_Jockey3-CENvsNot.mat.gz -o Dmel_PROseq_PLUS_Jockey3-CENvsNot_rename.mat.gz --groupLabels "CEN" "nonCEN"
# minus strand regions
computeMatrix reference-point --referencePoint TES -p 12 \
-R /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_MINUS.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_MINUS.bed \
-S /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig \
-o Dmel_PROseq_MINUS_Jockey3-CENvsNot.mat.gz \
-b 5000 -a 200 -bs 20 \
--averageTypeBins max \
--missingDataAsZero
# relabel:
computeMatrixOperations relabel -m Dmel_PROseq_MINUS_Jockey3-CENvsNot.mat.gz -o Dmel_PROseq_MINUS_Jockey3-CENvsNot_rename.mat.gz --groupLabels "CEN" "nonCEN"
# combine plus and minus strand region matrices
computeMatrixOperations rbind -m Dmel_PROseq_PLUS_Jockey3-CENvsNot_rename.mat.gz Dmel_PROseq_MINUS_Jockey3-CENvsNot_rename.mat.gz \
-o Dmel_PROseq_Jockey3-CENvsNot_rename_OUT.mat.gz
# create heatmap
plotHeatmap -m Dmel_PROseq_Jockey3-CENvsNot_rename_OUT.mat.gz --sortUsing region_length \
-o Dmel_PROseq_Jockey3-CENvsNot_rename_OUT_SE_b5000-a200-bs20.htmp.pdf \
--dpi 800 --plotType se \
--colorMap Blues Reds \
--averageTypeSummaryPlot mean \
--samplesLabel "Antisense" "Sense" \
--heatmapHeight 15 --heatmapWidth 4 \
--zMin 0 --yMin 0
```
### Plot PROseq signal across Jockey3 (CEN FL vs. NonCEN FL vs. CEN Trunc vs. NonCEN Trunc) aligned at TES ###
```
#plus strand regions
computeMatrix reference-point --referencePoint TES -p 12 \
-R /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_PLUS_FLonly.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_PLUS_TruncOnly.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_PLUS_FLonly.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_PLUS_TruncOnly.bed \
-S /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig \
-o Dmel_PROseq_PLUS_Jockey3-CENvsNot-FLvsTrunc.mat.gz \
-b 5000 -a 200 -bs 20 \
--averageTypeBins max \
--missingDataAsZero
#relabel:
computeMatrixOperations relabel -m Dmel_PROseq_PLUS_Jockey3-CENvsNot-FLvsTrunc.mat.gz -o Dmel_PROseq_PLUS_Jockey3-CENvsNot-FLvsTrunc_rename.mat.gz --groupLabels "CEN_FL" "CEN_Trunc" "nonCEN_FL" "nonCEN_Trunc"
# minus strand regions
computeMatrix reference-point --referencePoint TES -p 12 \
-R /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_MINUS_FLonly.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_CENonly_Jockey3only_short_MINUS_TruncOnly.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_MINUS_FLonly.bed /core/labs/Oneill/T2T_working_space/Dmel_analyses/Annotations/Jockey3_only/ALLrepeats_dmel_repeat_annotations_withJockeyLabels_notCEN_Jockey3only_short_MINUS_TruncOnly.bed \
-S /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_plus.bigwig /core/labs/Oneill/T2T_working_space/Dmel_analyses/Bigwigs/ProSeq_embryo_Iso1_bowtie2_Df_nameSort_3ends_minus.bigwig \
-o Dmel_PROseq_MINUS_Jockey3-CENvsNot-FLvsTrunc.mat.gz \
-b 5000 -a 200 -bs 20 \
--averageTypeBins max \
--missingDataAsZero
# relabel:
computeMatrixOperations relabel -m Dmel_PROseq_MINUS_Jockey3-CENvsNot-FLvsTrunc.mat.gz -o Dmel_PROseq_MINUS_Jockey3-CENvsNot-FLvsTrunc_rename.mat.gz --groupLabels "CEN_FL" "CEN_Trunc" "nonCEN_FL" "nonCEN_Trunc"
# combine plus and minus strand region matrices
computeMatrixOperations rbind -m Dmel_PROseq_PLUS_Jockey3-CENvsNot-FLvsTrunc_rename.mat.gz Dmel_PROseq_MINUS_Jockey3-CENvsNot-FLvsTrunc_rename.mat.gz \
-o Dmel_PROseq_Jockey3-CENvsNot-FLvsTrunc_rename_OUT.mat.gz
# create heatmap
plotHeatmap -m Dmel_PROseq_Jockey3-CENvsNot-FLvsTrunc_rename_OUT.mat.gz --sortUsing region_length \
-o Dmel_PROseq_Jockey3-CENvsNot-FLvsTrunc_rename_OUT_SE_b5000-a200-bs20.htmp.pdf \
--dpi 800 --plotType se \
--colorMap Blues Reds \
--averageTypeSummaryPlot mean \
--samplesLabel "Antisense" "Sense" \
--heatmapHeight 15 --heatmapWidth 4 \
--zMin 0 --yMin 0
```
