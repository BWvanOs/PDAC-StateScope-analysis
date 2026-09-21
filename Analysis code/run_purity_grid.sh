INPUT_DIR="/INPUT/DIR/"
OUTPUT_DIR="${INPUT_DIR}/purity_grid"

mkdir -p "${OUTPUT_DIR}"

# 0.30, 0.35, ... 1.00
PURITIES=$(seq 0.30 0.05 1.00)

for CNS in "${INPUT_DIR}"/*.final.cns; do

    SAMPLE=$(basename "${CNS}" .final.cns)
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"

    mkdir -p "${SAMPLE_OUT}"

    for PURITY in ${PURITIES}; do

        OUT="${SAMPLE_OUT}/purity_${PURITY}.call.cns"

        echo "  purity = ${PURITY}"

        cnvkit.py call \
            "${CNS}" \
            -m clonal \
            --purity "${PURITY}" \
            --ploidy 2 \
            -o "${OUT}"

    done

done

echo
echo "Finished purity grid."