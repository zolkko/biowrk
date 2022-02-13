#!/bin/bash


set -uex


# Directory for the genome assembly and related files
GENDIR=hg38

# Directories in which reads will be stored
REFSDIR=refs

# Directories related to quality control
QCDIR=qc
MQCDIR=mqc

# Output directory for the generate sorted BAM
BAMDIR=bams

# Output file name (GSE id) for the BAM and the output-log file.
OUTF=GSE63577

# Directory to store counts dataset
CNTDIR=counts


# ===============================================
# Create hg38 genome directory
# ===============================================

echo "=========================================="
echo "Create ${GENDIR} for hg38 genome assembly"
echo "=========================================="
mkdir -p ${GENDIR}


# ===============================================
# Download hg38 genome
# ===============================================

echo "=========================================="
echo "Download hg38 genome assembly to ${GENDIR}"
echo "=========================================="
wget -P ${GENDIR} http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip ${GENDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz


# ===============================================
# Fix chromosome names inside the genome
# ===============================================

echo "=========================================="
echo "Fix chromosome names in the hg38 genome"
echo "=========================================="
mv ${GENDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${GENDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.origin
cut ${GENDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.origin |\
    sed 's/^>\([[:digit:]]\+\)/>chr\1/g' |\
    sed 's/^>X/>chrX/g' |\
    sed 's/^>Y/>chrY/g' > ${GENDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa


# ===============================================
# Create the index for the genome
# ===============================================

echo "=========================================="
echo "Create the index for the genome"
echo "=========================================="
samtools faidx ${GENDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa


# ===============================================
# Download annotations
# ===============================================

echo "=========================================="
echo "Download the annotation for the hg38 gnome assemly"
echo "=========================================="
wget -P ${GENDIR} https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
gunzip ${GENDIR}/gencode.v39.annotation.gtf.gz


# ===============================================
# Create supplementary files for indexing the genome assembly
#
# For this assey the exon regions are not nessessary.
# ===============================================

echo "=========================================="
echo "Create splice-site file"
echo "=========================================="
extract_splice_sites.py ${GENDIR}/gencode.v39.annotation.gtf > ${GENDIR}/gencode.v39.annotation.ss


echo "=========================================="
echo "Create exom file"
echo "=========================================="
hisat2_extract_exons.py ${GENDIR}/gencode.v39.annotation.gtf > ${GENDIR}/gencode.v39.annotation.exon


# ===============================================
# Create HiSat2 index
#
# If your system unable to generate the index
# with SS and EXON provided, then just comment out
# corresponding parameters.
# ===============================================

echo "=========================================="
echo "Create HiSat2 index"
echo "=========================================="

hisat2-build \
	-p 10 \
	--seed 100\
	--ss ${GENDIR}/gencode.v39.annotation.ss \
	--exon ${GENDIR}/gencode.v39.annotation.exon \
	${GENDIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	${GENDIR}/GRCh38


# ===============================================
# Download reads.
#
# Fist three reads are Young then go Old reads.
# ===============================================

echo "=========================================="
echo "Download reads"
echo "=========================================="
RUNS=("SRR1660555" "SRR1660556" "SRR1660557" "SRR1660558" "SRR1660559" "SRR1660560")

mkdir -p ${REFSDIR}

for RUN in ${RUNS[@]}
do
	echo "=============================="
	echo "Downloading $RUN..."
	echo "=============================="
	fastq-dump --gzip --outdir ${REFSDIR} $RUN &
done
wait


# ===============================================
# Quality control
# ===============================================

echo "=========================================="
echo "Run quality control"
echo "=========================================="
mkdir -p ${QCDIR}
mkdir -p ${MQCDIR}

fastqc -o ${QCDIR} ${REFSDIR}/*.fastq.gz
multiqc -o ${MQCDIR} ${QCDIR}

read -n 1 -s -r -p "Review the QC reports and Press any key to continue"


# ===============================================
# Align reads
#
# Here I am aligning all the reads at once
# ===============================================

READS=$(find ${REFSDIR} -iname "*.fastq.gz" -type f | tr '\n' ',')

echo "=========================================="
echo "Map read: ${READS}"
echo "Save BAM files to sorted: ${BAMDIR}/${OUTF}.bam"
echo "=========================================="

mkdir -p ${BAMDIR}

hisat2 \
	--no-softclip \
	--no-unal \
	--new-summary \
	--time \
	-x ${GENDIR}/GRCh38 \
	--known-splicesite-infile ${GENDIR}/gencode.v39.annotation.ss \
	-U ${READS} \
	--threads 10 \
	--seed 100 \
	--summary-file ${BAMDIR}/${OUTF}.log | samtools view -Sb > ${BAMDIR}/${OUTF}.bam


# ===============================================
# Now I must split one big BAM file into multiple
# BAM files and sort them
# ===============================================

echo "=========================================="
echo "Split ${OUTF}.bam by qname (RUN id)"
echo "=========================================="

for RUN in ${RUNS[@]}
do
    (samtools view -H ${BAMDIR}/${OUTF}.bam &&\
        (samtools view ${BAMDIR}/${OUTF}.bam | grep "^${RUN}")) |\
        samtools sort > ${BAMDIR}/${RUN}.bam &&\
        samtools index ${BAMDIR}/${RUN}.bam &
done
wait


# ===============================================
# Once BAM files are split and their
# indexes are ready, I can count gene hits.
# ===============================================

echo "=========================================="
echo "Count gene Hits"
echo "=========================================="

mkdir -p ${CNTDIR}

for RUN in ${RUNS[@]}
do
    htseq-count -f bam \
        -t gene \
        -s no \
        -m union \
        ${BAMDIR}/${RUN}.bam \
        ${GENDIR}/gencode.v39.annotation.gtf | tee ${CNTDIR}/${RUN}.txt &
done
wait


# ===============================================
# Combine multiple counts-files into one dataset.
# ===============================================

echo "=========================================="
echo "Combine the counts dataset"
echo "=========================================="

# head_runs=$(printf "\t%s" "${RUNS[@]}")
# head_runs=$(printf "ENSEMBLID\t%s" "${head_runs:1}")

{ echo -e "ENSEMBLID\tSRR1660555\tSRR1660556\tSRR1660557\tSRR1660558\tSRR1660559\tSRR1660560" & join -t $'\t' ${CNTDIR}/SRR1660555.txt $CNTDIR/SRR1660556.txt |\
	join -t $'\t' $CNTDIR/SRR1660557.txt - |\
	join -t $'\t' $CNTDIR/SRR1660558.txt - |\
	join -t $'\t' $CNTDIR/SRR1660559.txt - |\
	join -t $'\t' $CNTDIR/SRR1660560.txt - | head -n -5; } > $CNTDIR/counts_data.tsv

