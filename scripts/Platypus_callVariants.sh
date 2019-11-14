#!/bin/bash

SCRIPT=`basename $0`
# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 5 ]; then
cat << EOF 
Usage: qsub -t 1-10 Platypus_callVariants.sh -F \"bamDir refPath bamSuffix outDir numCores\"
This script implements Platypus' callVariants function. All SGE run parameters 
should be specified at submission.  

bamDir: Directory containing bam files

refPath: Full path to reference sequence fasta

bamSuffix: bam suffix that links bam to a reference (e.g. B73v4 will identify
bams ending in B73v4.bam)

outDir: Directory to write output

numCores: number of threads to parallelize across.  Ensure that the requested number
is appropriate for the requested queue.  

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

#	Put them into a format that will be accepted by the Platypus command line
BAM_STR=()
for s in "${SAMPLE_LIST[@]}"
do
	BAM_STR+=("${s},")
done
BAM_STR=($(echo ${BAM_STR[@]} | sed 's/.$//' | sed 's/ //g'))

#   Get which sequence we are analyzing with the task ID
SEQNAMES=($(cut -f 2 ${REF_DICT} | grep -E '^SN' | cut -f 2 -d ':'))
CURRENT=${SEQNAMES[${PBS_ARRAYID}]}

#Set output directory
OUTPUT_DIR="${4}"
mkdir -p ${OUTPUT_DIR}

# Make sure output file does not already exist
if [ "$(ls -A ${OUTPUT_DIR}/${SUF}_${CURRENT}_Platy.vcf )" ]; then
     echo "Output VCF already exists in output directory" >&2
     exit 1
fi

python /home/hirschc1/pmonnaha/software/Platypus_0.8.1/Platypus.py callVariants \
	--nCPU=${5} \
	--refFile=${REF} \
  --regions=${CURRENT} \
	--bamFiles=${BAM_STR[@]}\
	-o ${OUTPUT_DIR}/${SUF}_${CURRENT}_Platy.vcf
