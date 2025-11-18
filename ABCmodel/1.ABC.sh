SAMPLE=("501mel" "SKMEL-5" "SKMEL-147")
MACS2=MACS2_folder
BAM=bowtie2_folder
RNA=RNAseq/cpm_folder
INPUT=HiC_matrix_folder
OUTPUT=ABC_folder
JUICER=juicer_tools_1.19.02.jar

# Call candidate regions
python ABC-Enhancer-Gene-Prediction/workflow/scripts/makeCandidateRegions.py \
--narrowPeak $MACS2/${SAMPLE}-ATAC_narrow_Peaks.bed \
--accessibility $BAMATAC/${SAMPLE}-ATAC.bam \
--outDir $OUTPUT/${SAMPLE} \
--chrom_sizes Hg38_chrom_size.tsv \
--chrom_sizes_bed Hg38_chrom_size.bed \
--regions_blocklist hg38-blacklist.v2.bed \
--regions_includelist GENCODE_v44_Hg38_wholegene.sorted.bed \
--peakExtendFromSummit 250 \
--nStrongestPeaks 150000

# Quantifying Enhancer Activity
python ABC-Enhancer-Gene-Prediction/workflow/scripts/run.neighborhoods.py \
--candidate_enhancer_regions $OUTPUT/${SAMPLE}/${SAMPLE}-ATAC_narrow_Peaks.candidateRegions.bed \
--primary_gene_identifier symbol \
--genes GENCODE_v44_Hg38_wholegene.sorted.bed \
--H3K27ac $BAMK27/${SAMPLE}-K27ac.bam \
--ATAC $BAMATAC/${SAMPLE}-ATAC.bam \
--expression_table $RNA/${SAMPLE}-cpm.csv \
--chrom_sizes Hg38_chrom_size.tsv \
--chrom_sizes_bed Hg38_chrom_size.bed \
--ubiquitously_expressed_genes UbiquitouslyExpressedGenes.txt \
--cellType ${SAMPLE} \
--outdir $OUTPUT/${SAMPLE}

###IMPORTANT!!! KR normalization from juicer is not ideal and some chromosomes may not have BP 5000, modify "juicebox_dump.py" by replacing all KR with VC. 
#I have included a modified version of the juicebox_dump.py script you can overwrite with the original one.
python ABC-Enhancer-Gene-Prediction/workflow/scripts/juicebox_dump.py \
--hic_file $INPUT/${SAMPLE}.allValidPairs.hic \
--juicebox "java -jar $JUICER" \
--outdir $OUTPUT/${SAMPLE}/HiC \
--resolution 5000 \
--include_raw \

###IMPORTANT!!! Before running, remove "interpolate_nan=False" and turn "allow_vc=True" from "compute_powerlaw_fit_from_hic.py".
#I have included a modified version of the compute_powerlaw_fit_from_hic.py script you can overwrite with the original one.
#Fit HiC data to powerlaw model and extract parameters
python ABC-Enhancer-Gene-Prediction/workflow/scripts/compute_powerlaw_fit_from_hic.py \
--hic_dir $OUTPUT/${SAMPLE}/HiC \
--outDir $OUTPUT/${SAMPLE} \
--maxWindow 1000000 \
--minWindow 5000 \
--hic_resolution 5000 \
