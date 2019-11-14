#!/bin/bash

SCRIPT=`basename $0`
# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 6 ]; then
cat << EOF 
Usage: qsub Platypus_callVariants_fromVCF_parallel.sh -F \"bamDir refPath bamSuffix outDir numCores VCF\"
This script implements Platypus' callVariants function. All SGE run parameters 
should be specified at submission.  

bamDir: Directory containing bam files

refPath: Full path to reference sequence fasta

bamSuffix: bam suffix that links bam to a reference (e.g. B73v4 will identify
bams ending in B73v4.bam)

outDir: Directory to write output

numCores: number of threads to parallelize across.  Ensure that the requested number
is appropriate for the requested queue.

VCF: vcf file that you want the new samples to be genotyped at  

EOF
  exit 1
fi

# check that input directory exists and is a directory
if ! [ -d "$1" ]; then
  echo "$1 is not a directory" >&2
  exit 1
fi
if ! [ -e "$1" ]; then
  echo "$1 not found" >&2
  exit 1
fi

# Check that reference file exists
if ! [ -e "$2" ]; then
  echo "$2 not found" >&2
  exit 1
fi


source activate py27
module load picard-tools/2.18.16
module load liblzma
module load samtools/1.7

export C_INCLUDE_PATH=/panfs/roc/msisoft/htslib/1.6/include
export LD_LIBRARY_PATH=//panfs/roc/msisoft/htslib/1.6/lib

#	Path to directory containing BAMs
BAM_DIR="$1"

#	Path to the reference
REF="$2"
REF_DICT="${REF%.*}.dict"

#	 Build the sample list
SUF="${3}"
SAMPLE_LIST=($(find ${BAM_DIR} -name "*${SUF}.bam"))

#Set output directory
OUTPUT_DIR="${4}"
mkdir -p ${OUTPUT_DIR}

#VCF containing variants to be genotyped
VCF="${6}"

FILE="${SAMPLE_LIST[${PBS_ARRAYID}]}"
SAMPLE="$(basename -- $FILE)"



#Make temporary copy of template VCF

cp ${VCF} ${VCF%.vcf.gz}.${SAMPLE}.vcf.gz
cp ${VCF}.tbi ${VCF%.vcf.gz}.${SAMPLE}.vcf.gz.tbi

#Platypus will use the sample name in the read group to label samples in the VCF, which was merging duplicate fastqs.  We want to keep fastqs separate to make sure these dupliciates are grouping as sister on the tree.  This command replaces the sample name read group with the fastq name.
S=${FILE##*/}
NEWFILE="${FILE%.bam}.SM.bam"
java -jar /panfs/roc/msisoft/picard/2.18.16/picard.jar AddOrReplaceReadGroups I="$FILE" O="$NEWFILE" RGSM="${S%.bam}" RGLB="${S%.bam}" RGID="${S%.bam}" RGPL=illumina RGPU="${S%.bam}"
samtools index $NEWFILE

#Run platypus
python /home/hirschc1/pmonnaha/software/Platypus_0.8.1/Platypus.py callVariants \
                --nCPU=${5} \
                --refFile=${REF} \
                --minMapQual=30 \
                --bamFiles=${NEWFILE} \
                --source=${VCF%.vcf.gz}.${SAMPLE}.vcf.gz \
                --getVariantsFromBAMs=0 \
                --minPosterior=0 \
                --genIndels=0 \
                --skipDifficultWindows=1 \
                -o ${OUTPUT_DIR}/${SAMPLE}_${SUF}_Platy.knwn.vcf

rm ${VCF%.vcf.gz}.${SAMPLE}.vcf.gz
rm ${VCF%.vcf.gz}.${SAMPLE}.vcf.gz.tbi
rm $NEWFILE
bgzip ${OUTPUT_DIR}/${SAMPLE}_${SUF}_Platy.knwn.vcf
tabix -p vcf ${OUTPUT_DIR}/${SAMPLE}_${SUF}_Platy.knwn.vcf.gz

exit
