MACS2=MACS2
OUTPUT=ABC

######Computing QC
mkdir -p $OUTPUT/SKMEL-5/qc
python /gpfs/data/HernandoLab/home/bericp01/ABC-Enhancer-Gene-Prediction/workflow/scripts/grabMetrics.py \
--macs_peaks $MACS2/SKMEL-5-ATAC_narrow_Peaks.bed \
--preds_file $OUTPUT/SKMEL-5/predictions_filtered/output.tsv \
--neighborhood_outdir $OUTPUT/SKMEL-5 \
--chrom_sizes Hg38_chrom_size.tsv \
--outdir $OUTPUT/SKMEL-5/qc \
--output_qc_summary $OUTPUT/SKMEL-5/qc/qc.tsv \
--output_qc_plots $OUTPUT/SKMEL-5/qc/qc.pdf \
--hic_gamma 1.0528504866340533 \
--hic_scale 6.455244375638351 

mkdir -p $OUTPUT/SKMEL-147/qc
python /gpfs/data/HernandoLab/home/bericp01/ABC-Enhancer-Gene-Prediction/workflow/scripts/grabMetrics.py \
--macs_peaks $MACS2/SKMEL-147-ATAC_narrow_Peaks.bed \
--preds_file $OUTPUT/SKMEL-147/predictions_filtered/output.tsv \
--neighborhood_outdir $OUTPUT/SKMEL-147 \
--chrom_sizes Hg38_chrom_size.tsv \
--outdir $OUTPUT/SKMEL-147/qc \
--output_qc_summary $OUTPUT/SKMEL-147/qc/qc.tsv \
--output_qc_plots $OUTPUT/SKMEL-147/qc/qc.pdf \
--hic_gamma 0.9239335729836557 \
--hic_scale 4.846981848320076

mkdir -p $OUTPUT/501mel/qc
python ABC-Enhancer-Gene-Prediction/workflow/scripts/grabMetrics.py \
--macs_peaks $MACS2/501mel-ATAC_narrow_Peaks.bed \
--preds_file $OUTPUT/501mel/predictions_filtered/output.tsv \
--neighborhood_outdir $OUTPUT/501mel \
--chrom_sizes Hg38_chrom_size.tsv \
--outdir $OUTPUT/501mel/qc \
--output_qc_summary $OUTPUT/501mel/qc/qc.tsv \
--output_qc_plots $OUTPUT/501mel/qc/qc.pdf \
--hic_gamma 0.9456060921860431 \
--hic_scale 5.081208553261949
