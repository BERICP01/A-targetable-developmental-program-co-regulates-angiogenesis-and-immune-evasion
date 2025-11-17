SAMPLE=("A375-shH2-K27ac" "A375-shH3-K27ac" "A375-shNTC-K27ac-rep1" "A375-shNTC-K27ac-rep2")

hichipper --out hichipper_${SAMPLE}  ${SAMPLE}_hichipper.yaml --keep-samples $CURRENT_SAMPLE --keep-temp-files --min-dist 1000 --read-length 150
