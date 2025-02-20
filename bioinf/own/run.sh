nextflow run nf-core/sarek \
    -profile docker \
    --input samplesheet.csv \
    --outdir ./sarek_results \
    --genome GATK.GRCh38
