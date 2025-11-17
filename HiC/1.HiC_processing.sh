DATA_DIR=fastq_directory
OUTPUT_DIR=HiC-Pro
CONFIG=config-hicpro.txt
CHR_SIZE=Hg38_chrom_size.bed #hg38 bed file containing two column chromosome name and bp size
TMP=tmp_folder
OUTPUT_DIR=output_matrix
RESFRAG=hg38_arima.bed #hg38 genome digested with arima enzymes
sample_name=("501mel" "SKMEL-147" "SKMEL-5" "A375-IgG-HiChIP" "A375-shH2-K27ac-HiChIP" "A375-shH3-K27ac-HiChIP" "A375-shNTC-K27ac-HiChIP-rep1" "A375-shNTC-K27ac-HiChIP-rep2")
# Run HiC-Pro

HiC-Pro_3.1.0/bin/HiC-Pro -i $DATA_DIR -o $OUTPUT_DIR -c $CONFIG -p

# Generate HiC matrix

HiC-Pro_3.1.0/bin/utils/hicpro2higlass.sh -i $OUTPUT_DIR/hic_results/data/${sample_name}/${sample_name}.allValidPairs \
-r 5000 -c $CHR_SIZE -n -t $TMP -o $OUTPUT_DIR2
cooler zoomify \
    -r 5000,10000,25000,40000,50000,100000,250000,500000,1000000,2500000,5000000 \
    -o $OUTPUT_DIR2/${sample_name}.mcool \
    $OUTPUT_DIR2/${sample_name}.cool

for res in 5000 10000 25000 40000 50000 100000 250000 500000 1000000 2500000 5000000
do
    echo "Balancing resolution $res ..."
    cooler balance $OUTPUT_DIR2/${sample_name}.mcool::resolutions/$res
done

HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh -i $OUTPUT_DIR/hic_results/data/${sample_name}/${sample_name}.allValidPairs \
-g $CHR_SIZE -j $JUICER -r $RESFRAG -t $TMP -o $OUTPUT
