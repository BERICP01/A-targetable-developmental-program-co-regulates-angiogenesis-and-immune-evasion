FASTQ='fastq'
REFERENCE='PIPseeker_mm10'
SAMPLE=("NTC-CD45-NEG" "NTC-CD45-POS" "SH-HOX-CD45-NEG" "SH-HOX-CD45-POS")
CURRENT_SAMPLE=${SAMPLE[$SLURM_ARRAY_TASK_ID-1]}

# pipseeker
mkdir -p pipseeker
mkdir -p pipseeker/B16F10-${CURRENT_SAMPLE}
pipseeker full --fastq $FASTQ/B16F10-${CURRENT_SAMPLE} \
--star-index-path $REFERENCE --output-path pipseeker/B16F10-${CURRENT_SAMPLE} \
--chemistry V
