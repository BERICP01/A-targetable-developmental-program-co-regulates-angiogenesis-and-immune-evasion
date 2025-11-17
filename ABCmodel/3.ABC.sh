OUTPUT=ABC
SAMPLE=("SKMEL-5" "SKMEL-147" "501mel")
CURRENT_SAMPLE=${SAMPLE[$SLURM_ARRAY_TASK_ID-1]}

####Filter prediction
mkdir -p $OUTPUT/${CURRENT_SAMPLE}/predictions_filtered
python ABC-Enhancer-Gene-Prediction/workflow/scripts/filter_predictions.py \
--pred_file $OUTPUT/${CURRENT_SAMPLE}/predictions/EnhancerPredictionsAllPutative.tsv.gz \
--pred_nonexpressed_file $OUTPUT/${CURRENT_SAMPLE}/predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz \
--only_expressed_genes true \
--include_self_promoter false \
--score_column ABC.Score \
--threshold 0.02 \
--output_tsv_file $OUTPUT/${CURRENT_SAMPLE}/predictions_filtered/output.tsv \
--output_slim_tsv_file $OUTPUT/${CURRENT_SAMPLE}/predictions_filtered/output_slim.tsv \
--output_bed_file $OUTPUT/${CURRENT_SAMPLE}/predictions_filtered/output.bed \
--output_gene_stats_file $OUTPUT/${CURRENT_SAMPLE}/predictions_filtered/output_gene_stats.tsv
