INPUT=HiC_matrix_folder
OUTPUT=ABC

######Computing the ABC Score
python ABC-Enhancer-Gene-Prediction/workflow/scripts/predict.py \
--enhancers $OUTPUT/SKMEL-5/EnhancerList.txt \
--accessibility_feature 'ATAC' \
--genes $OUTPUT/SKMEL-5/GeneList.txt \
--hic_file $INPUT/SKMEL-5.allValidPairs.hic \
--hic_pseudocount_distance 5000 \
--hic_resolution 5000 \
--chrom_sizes Hg38_chrom_size.tsv \
--cellType SKMEL-5  \
--outdir $OUTPUT/SKMEL-5/predictions \
--make_all_putative \
--scale_hic_using_powerlaw \
--hic_gamma 1.0528504866340533 \
--hic_scale 6.455244375638351 \
--hic_gamma_reference 0.87 \
--hic_pseudocount_distance 5000

python ABC-Enhancer-Gene-Prediction/workflow/scripts/predict.py \
--enhancers $OUTPUT/SKMEL-147/EnhancerList.txt \
--accessibility_feature 'ATAC' \
--genes $OUTPUT/SKMEL-147/GeneList.txt \
--hic_file $INPUT/SKMEL-147.allValidPairs.hic \
--hic_pseudocount_distance 5000 \
--hic_resolution 5000 \
--chrom_sizes Hg38_chrom_size.tsv \
--cellType SKMEL-147  \
--outdir $OUTPUT/SKMEL-147/predictions \
--make_all_putative \
--scale_hic_using_powerlaw \
--hic_gamma 0.9239335729836557 \
--hic_scale 4.846981848320076 \
--hic_gamma_reference 0.87 \
--hic_pseudocount_distance 5000


python ABC-Enhancer-Gene-Prediction/workflow/scripts/predict.py \
--enhancers $OUTPUT/501mel/EnhancerList.txt \
--accessibility_feature 'ATAC' \
--genes $OUTPUT/501mel/GeneList.txt \
--hic_file $INPUT/501mel.allValidPairs.hic \
--hic_pseudocount_distance 5000 \
--hic_resolution 5000 \
--chrom_sizes Hg38_chrom_size.tsv \
--cellType 501mel \
--outdir $OUTPUT/501mel/predictions \
--make_all_putative \
--scale_hic_using_powerlaw \
--hic_gamma 0.9456060921860431 \
--hic_scale 5.081208553261949 \
--hic_gamma_reference 0.87 \
--hic_pseudocount_distance 5000





