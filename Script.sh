
#!/bin/bash


# Trio-Based Exome Analysis for the Identification of Inherited and De Novo Rare Mendelian Disorders
# Authors: Lazzaroni Laura and Piovesan Lara


### START ###

# Start inside the 'trio_analyses' parent folder
# Ensure to have the analysis_info.txt file and the script.sh file inside this folder
# After running the script, this folder will contain all the files required for each trio analysis (aka 'common files')
# Inside this folder, a sub-folder will be automatically created for each trio
# The sub-folder will contain the files required for each specific trio analysis (aka 'trio-specific files')

# Remember to change the path to the directories at the beginning!


### IMPORTING COMMON FILES ###

echo "Creating links to the common files"

# Reference genome (indexed and Bowtie2-indexed)
ln -s /home/BCG2026_exam/chr20.* .

# Annotations
ln -s /home/BCG2026_exam/chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed target_regions.bed

# Samples file
echo $'child\nfather\nmother' > samples.txt


### START THE ITERATION ###

# Now, based on the information provided in the analysis_info.txt file, the script will run accordingly

while read SEQ_MODE TRIO_NUM CHILD_ID FATHER_ID MOTHER_ID INHERITANCE NOTES; do


    ### NEW DIRECTORY AND IMPORTING TRIO-SPECIFIC FILES ###

    echo "Creating a new directory and importing the trio-specific files"
    mkdir -p "${TRIO_NUM}"

    echo "Creating links to the trio-specific files"
    cd ${TRIO_NUM}

    # Raw sequencing data (exome)
    ln -s /home/BCG2026_exam/BCG2026_Lazzaroni_L/${TRIO_NUM}/* .


    ### PRE-ALIGNMENT QUALITY CONTROL ###

    echo "Pre-alignment quality control"

    # FastQC
    fastqc *.fq.gz


    ### ALIGNMENT TO THE REFERENCE GENOME ###

    echo "Alignment"

    # If the sequencing mode is single-end
    if [ "$SEQ_MODE" == 'SE' ]; then
        bowtie2 -x ../chr20 -U ${CHILD_ID}.fq.gz --rg-id "child" --rg "SM:child" | samtools view -Sb | samtools sort -o child.bam &
        bowtie2 -x ../chr20 -U ${FATHER_ID}.fq.gz --rg-id "father" --rg "SM:father" | samtools view -Sb | samtools sort -o father.bam &
        bowtie2 -x ../chr20 -U ${MOTHER_ID}.fq.gz --rg-id "mother" --rg "SM:mother" | samtools view -Sb | samtools sort -o mother.bam &
        wait
    # If the sequencing mode is paired-end
    elif [ "$SEQ_MODE" == 'PE' ]; then
        bowtie2 -x ../chr20 -1 ${CHILD_ID}.targets_R1.fq.gz -2 ${CHILD_ID}.targets_R2.fq.gz --rg-id "child" --rg "SM:child" | samtools view -Sb | samtools sort -o child.bam &
        bowtie2 -x ../chr20 -1 ${FATHER_ID}.targets_R1.fq.gz -2 ${FATHER_ID}.targets_R2.fq.gz --rg-id "father" --rg "SM:father" | samtools view -Sb | samtools sort -o father.bam &
        bowtie2 -x ../chr20 -1 ${MOTHER_ID}.targets_R1.fq.gz -2 ${MOTHER_ID}.targets_R2.fq.gz --rg-id "mother" --rg "SM:mother" | samtools view -Sb | samtools sort -o mother.bam &
        wait
    fi

    # Indexing
    samtools index child.bam
    samtools index father.bam
    samtools index mother.bam


    ### POST-ALIGNMENT QUALITY CONTROL ###

    echo "Post-alignment quality control"

    # Qualimap
    qualimap bamqc -bam child.bam --feature-file ../target_regions.bed --outdir child &
    qualimap bamqc -bam father.bam --feature-file ../target_regions.bed --outdir father &
    qualimap bamqc -bam mother.bam --feature-file ../target_regions.bed --outdir mother &
    wait

    # MultiQC
    multiqc .
    mv multiqc_report.html ${TRIO_NUM}_multiqc_report.html


    ### VARIANT CALLING AND FILTERING ###

    echo "Variant calling"

    # Variant calling
    freebayes -f ../chr20.fa -m 20 -C 5 -Q 10 -q 10 --min-coverage 10 child.bam father.bam mother.bam > ${TRIO_NUM}.vcf

    # Compress VCF file and allow random access
    bgzip ${TRIO_NUM}.vcf

    # Indexing
    bcftools index ${TRIO_NUM}.vcf.gz

    echo "Filtering"

    # Filtering
    if [ "$INHERITANCE" == 'AR' ]; then
        echo "Recessive"
        bcftools view -R ../target_regions.bed ${TRIO_NUM}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0]="AA" && GT[1]="RA" && GT[2]="RA"' | bcftools filter -i 'QUAL>20' -Ov -o ${TRIO_NUM}.cand.vcf
    elif [ "$INHERITANCE" == 'AD_inherited' ]; then
        if [ "$NOTES" == 'father_affected' ]; then
            echo "Dominant, father affected"
            bcftools view -R ../target_regions.bed ${TRIO_NUM}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0]="RA" && GT[1]="RA" && GT[2]="RR"' | bcftools filter -i 'QUAL>20' -Ov -o ${TRIO_NUM}.cand.vcf
        elif [ "$NOTES" == 'mother_affected' ]; then
            echo "Dominant, mother affected"
            bcftools view -R ../target_regions.bed ${TRIO_NUM}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0]="RA" && GT[1]="RR" && GT[2]="RA"' | bcftools filter -i 'QUAL>20' -Ov -o ${TRIO_NUM}.cand.vcf
        fi
    elif [ "$INHERITANCE" == 'AD_denovo' ]; then
        echo "Dominant, de novo"
        bcftools view -R ../target_regions.bed ${TRIO_NUM}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0]="RA" && GT[1]="RR" && GT[2]="RR"' | bcftools filter -i 'QUAL>20' -Ov -o ${TRIO_NUM}.cand.vcf
    fi


    ### COVERAGE TRACKS ###

    echo "Coverage tracks"

    bedtools genomecov -ibam father.bam -bg -trackline -trackopts 'name="father"' -max 100 > ${TRIO_NUM}_fatherCov.bg &
    bedtools genomecov -ibam mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > ${TRIO_NUM}_motherCov.bg &
    bedtools genomecov -ibam child.bam -bg -trackline -trackopts 'name="child"' -max 100 > ${TRIO_NUM}_childCov.bg &
    wait


    ### FILE EXPORT ###

    echo "File export"

    EXPORT_DIR="${TRIO_NUM}_export"
    mkdir -p "$EXPORT_DIR"

    cp ${TRIO_NUM}_multiqc_report.html "$EXPORT_DIR/"
    cp ${TRIO_NUM}.cand.vcf "$EXPORT_DIR/"
    cp *.bg "$EXPORT_DIR/"

    cd ..


done < analysis_info.txt

