FASTQ=fastq
PIP=pipseeker
REFERENCE=genome_references/TRUST4/mouse
OUTPUT=B16_scRNAseq/trust4
SAMPLE=("NTC-CD45-POS" "SH-HOX-CD45-POS")
CURRENT_SAMPLE=${SAMPLE[$SLURM_ARRAY_TASK_ID-1]}
mkdir -p $OUTPUT/${CURRENT_SAMPLE}

cat $FASTQ/${CURRENT_SAMPLE}_S*_L*_R1_001.fastq.gz > $FASTQ/${CURRENT_SAMPLE}_R1_merged.fastq.gz
cat $FASTQ/${CURRENT_SAMPLE}_S*_L*_R2_001.fastq.gz > $FASTQ/${CURRENT_SAMPLE}_R2_merged.fastq.gz

pipspeak \
-c pipspeak/data/config_v3.yaml \
-i $FASTQ/${CURRENT_SAMPLE}_R1_merged.fastq.gz \
-I $FASTQ/${CURRENT_SAMPLE}_R2_merged.fastq.gz

run-trust4 -f $REFERENCE/bcrtcr.fa --ref $REFERENCE/IMGT+C.fa \
-u $OUTPUT/${CURRENT_SAMPLE}/pipspeak_R2.fq.gz \
--barcode $OUTPUT/${CURRENT_SAMPLE}/pipspeak_R1.fq.gz --readFormat bc:0:15 \
--barcodeWhitelist $OUTPUT/${CURRENT_SAMPLE}/pipspeak_whitelist.txt \
-o ${CURRENT_SAMPLE} --od $OUTPUT/${CURRENT_SAMPLE} -t 16 
