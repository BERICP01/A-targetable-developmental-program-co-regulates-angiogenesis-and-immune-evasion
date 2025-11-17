FASTQ='fastq'
REFERENCE='PIPseeker_mm10'
SAMPLE=("NTC-CD45-NEG" "NTC-CD45-POS" "SH-HOX-CD45-NEG" "SH-HOX-CD45-POS")

# pipseeker
mkdir -p pipseeker
mkdir -p pipseeker/B16F10-${SAMPLE}
pipseeker full --fastq $FASTQ/B16F10-${SAMPLE} \
--star-index-path $REFERENCE --output-path pipseeker/B16F10-${SAMPLE} \
--chemistry V
