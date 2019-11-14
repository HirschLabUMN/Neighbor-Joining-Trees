# NJTrees
Make neighbor-joining trees.

The original purpose of this repository was to make simple neighbor-joining trees to confirm that samples with multiple fastq files were correctly labelled.  That is, these fastq files should group as sister in the NJ tree.  Variants are called using Platypus and trees are built using a combination of R packages (vcfR, adegenet, and StAMPP).

Although Platypus requires much less memory than other variant callers, it is still largely unable to handle cohorts > ~300.  For this reason, Platypus was run initially on the first 100 WiDiv samples, and then, newer samples were genotyped individually (again, using Platypus) at the SNPs identified from the initial run.  Following the individual-genotyping, VCF files are merged, filtered, and downsampled to ~50k sites.  This number of sites is typically sufficient for building an NJ tree, yet small enough not to crash R.

# Requirements
Platypus (v 0.8.1) \
bcftools (v1.9) \
cowplot (v1.0.0) \
lattice (v0.20-38) \
StAMPP  (v1.5.1) \
pegas (v0.12) \
ape (v5.3) \
adegraphics (v1.0-15) \
adegenet (v2.1.1) \
ade4 (v1.7-13) \
vcfR (v1.8.0)

# Initial variant calling
Variant calling for the initial 100 samples was performed via *scripts/Platypus_callVariants.sh* 

    Usage: qsub -t 1-10 Platypus_callVariants.sh -F "bamDir refPath bamSuffix outDir numCores"
    This script implements Platypus' callVariants function. All SGE run parameters
    should be specified at submission.

    bamDir: Directory containing bam files

    refPath: Full path to reference sequence fasta

    bamSuffix: bam suffix that links bam to a reference (e.g. B73v4 will identify
    bams ending in B73v4.bam)

    outDir: Directory to write output

    numCores: number of threads to parallelize across.  Ensure that the requested number
    is appropriate for the requested queue.

This shell script is meant to be submitted as a task array (i.e. -t 1-10), which will do variant calling separately for the 10 maize chromosomes.  However, be sure to specify all SGE run parameters at submission.  For example:

    qsub -t 1-10 Platypus_callVariants.sh -F "/panfs/roc/scratch/pmonnaha/Maize/widiv_bams/merged/ /home/hirschc1/pmonnaha/references/W22_chr1-10.fasta W22v12 /panfs/roc/scratch/pmonnaha/Maize/platy/w22 16" -e ~/OandE/Platypus_callVariants_W22.e -o ~/OandE/Platypus_callVariants_W22.o -A mcgaughs -m abe -M pmonnaha@umn.edu -l mem=62gb,nodes=1:ppn=16,walltime=48:00:00 -q sb
    
I have found that this amount of resources was sufficient for all jobs to complete.

# Individual genotyping from pre-existing VCF
Genotyping of individual samples based on a pre-existing VCF is accomplished via *scripts/Platypus_callVariants_fromVCF_parallel.sh* . Given that ~50k SNPs is sufficient to build a tree, you can probably get away with just doing the genotyping on a VCF from a single chromosome.

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
    
This script is also meant to be run as a task array, although the user will need to determine the total number of BAM files in the input directory.  For example, if there are 320 BAM files, then do *qsub -t 1-320*.  SGE parameters, again, must be specified at the command line upon submission.  I have found the following resources to be sufficient for moderately sized BAM files (20 - 50 Gb). 

    -l mem=16gb,nodes=1:ppn=1,walltime=06:00:00
    
# Merging, filtering, and downsampling
Once all samples have been genotyped, we can proceed to producing the VCF that will be used for building the NJ tree.  Much of this is done via *bcftools*

## Merge samples
Navigate to the output folder containing all of the individually-genotyped VCFs (make sure additional VCFs are not present in the directory).

    bcftools merge *vcf.gz -o <output_file> -Oz
    
## Filter samples
I leave this as a separate step, so that the user can adjust the filtration criteria as they see fit.

    bcftools filter -i 'F_PASS(GQ>10 & GT!="mis") > 0.5 && && (N_PASS(GT="het") / N_SAMPLES) < 0.1' <input_file> | bgzip > <output_file>
    
This call will keep only sites in which at least 50% of the samples have genotype quality greater than 10 and a non-missing genotype field. It will also filter sites with more than 20 heterozygous genotypes.  See http://samtools.github.io/bcftools/bcftools.html for additional filtration options.

## Downsampling
Using the full VCF file will likely cause R to crash, so we downsample to a more manageable amount of data.  Downsampling is accomplished via the *shuf* command.  However, we must strip off the VCF header prior to running this command.

    #Strip header
    zgrep "#" <filtered_vcf> | bgzip > head.vcf.gz
    #Randomly downsample non-header lines
    zgrep -v "#" <filtered_vcf> | shuf -n 50000 | bgzip > tmp.vcf.gz
    #Reassemble and sort VCF
    zcat head.vcf.gz tmp.vcf.gz | bcftools sort | bgzip > <final_filtered_downsampled_vcf>
    tabix -p vcf <final_filtered_downsampled_vcf>

Change the -n flag in the *shuf* command to sample more or less sites.

# Make NJ tree
Making the NJ tree is accomplished within R.  Instead of creating a standard script with input/output, I have created an R Notebook (/Markdown) file with the necessary functions to prepare and plot the data (*analysis/NJ_Tree.Rmd*).  This interactive setting situates the user to more easily implement small changes to the process without having to rerun an entire script.  

There are a number of accessory files, which facilitate labelling and coloring of samples (e.g. according to heterotic group).  I have included these files in *accessories*.

# Specific Notes

