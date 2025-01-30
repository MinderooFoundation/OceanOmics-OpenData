# OceanOmics Public Datasets on Amazon Web Services

This document describes the format and gives simple examples for getting started with the OceanOmics data stored on AWS.

# Ocean Genomes

The OceanOmics Ocean Genomes data aims to sequence all marine vertebrates. Its main product is single individual-genome sequencing data and genome assemblies.

Data is stored in FASTA, BAM, and FASTQ format in the following buckets:

(to fill once we have bucket names)


The structure of the data follows the [GenomeArk structure](https://genomeark.s3.amazonaws.com/index.html). 

One folder per species with the accepted Latin name,  
this folder contains samples in the form of Tree of Life IDs (ToLID), 
then several folders containing genome assembly and sequencing data:

| Folder name | Folder contents | Subfolder | Subfolder contents |
| ------------- | ------------- | ------------- | ------------- |
| genomic_data/ | raw read data split by data type. | NA | NA |
| " | " | illumina/ | raw Illumin reads in fastq.gz format |
| " | " | pacbio/ | raw PacBio reads in bam format |
| " | " | 10x/ | raw HiC data in fastq.gz format |
| assembly_hifiasm/ | First assembly step, PacBio data assembled using hifiasm. Contains haplotype 1 (hap1), haplotype 2 (hap2), and alternative contigs. | NA | NA |
| assembly_yahs/ | Second assembly step, the hifiasm assembly scaffolded using yahs and HiC data. Contains hap1 and hap2. | NA | NA |
| assembly_tiara/ | Third assembly step, the scaffolded assembly decontaminated using tiara. Contains hap1 and hap2, but also the combined scaffolds of hap1 and hap2. | NA | NA |
| assembly_curated/ | 'Final' genome assembly after manual curation using HiC contact matrices. Contains fasta.gz files of pseudo-chromosomes or scaffolds. | NA | NA |

Example: Assembly OG88, *Enoplosus armatus*

    Enoplosus_armatus/fEnoArm2/genomic_data/illumina/ - contains the raw Illumina reads OG88.ilmn.230814.R1.fastq.gz and OG88.ilmn.230814.R2.fastq.gz
    Enoplosus_armatus/fEnoArm2/assembly_curated/ - contains the final curated assembly. OG88G_v230728.hic1.3.curated.hap1.chr_level.fa and OG88G_v230728.hic1.3.curated.hap2.chr_level.fa

# OceanOmics eDNA

The OceanOmics eDNA data aims to learn how to use environmental DNA (eDNA) to assess ecosystem health. Its main product is metabarcoding raw reads and metagenomics raw reads.

The data is structured by OceanOmics voyage, several FASTQ files per OceanOmics expedition (voyage). The S3 bucket contains one folder per voyage, with each voyage folder containing paired end reads of every sample of the voyage.
