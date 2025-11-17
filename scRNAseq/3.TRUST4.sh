FASTQ=fastq
PIP=pipseeker
REFERENCE=genome_references/TRUST4/mouse
OUTPUT=B16_scRNAseq/trust4
SAMPLE=("NTC-CD45-POS" "SH-HOX-CD45-POS")
mkdir -p $OUTPUT/${SAMPLE}

cat $FASTQ/${SAMPLE}_S*_L*_R1_001.fastq.gz > $FASTQ/${SAMPLE}_R1_merged.fastq.gz
cat $FASTQ/${SAMPLE}_S*_L*_R2_001.fastq.gz > $FASTQ/${SAMPLE}_R2_merged.fastq.gz

pipspeak \
-c pipspeak/data/config_v3.yaml \
-i $FASTQ/${SAMPLE}_R1_merged.fastq.gz \
-I $FASTQ/${SAMPLE}_R2_merged.fastq.gz

run-trust4 -f $REFERENCE/bcrtcr.fa --ref $REFERENCE/IMGT+C.fa \
-u $OUTPUT/${SAMPLE}/pipspeak_R2.fq.gz \
--barcode $OUTPUT/${SAMPLE}/pipspeak_R1.fq.gz --readFormat bc:0:15 \
--barcodeWhitelist $OUTPUT/${SAMPLE}/pipspeak_whitelist.txt \
-o ${SAMPLE} --od $OUTPUT/${SAMPLE} -t 16 
