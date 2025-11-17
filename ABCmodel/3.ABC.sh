OUTPUT=ABC
SAMPLE=("SKMEL-5" "SKMEL-147" "501mel")

####Filter prediction
mkdir -p $OUTPUT/${SAMPLE}/predictions_filtered
python ABC-Enhancer-Gene-Prediction/workflow/scripts/filter_predictions.py \
--pred_file $OUTPUT/${SAMPLE}/predictions/EnhancerPredictionsAllPutative.tsv.gz \
--pred_nonexpressed_file $OUTPUT/${SAMPLE}/predictions/EnhancerPredictionsAllPutativeNonExpressedGenes.tsv.gz \
--only_expressed_genes true \
--include_self_promoter false \
--score_column ABC.Score \
--threshold 0.02 \
--output_tsv_file $OUTPUT/${SAMPLE}/predictions_filtered/output.tsv \
--output_slim_tsv_file $OUTPUT/${SAMPLE}/predictions_filtered/output_slim.tsv \
--output_bed_file $OUTPUT/${SAMPLE}/predictions_filtered/output.bed \
--output_gene_stats_file $OUTPUT/${SAMPLE}/predictions_filtered/output_gene_stats.tsv
