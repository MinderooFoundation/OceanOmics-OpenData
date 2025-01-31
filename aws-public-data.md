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
| " | " | hic/ | raw HiC data in fastq.gz format |
| assembly_hifiasm/ | First assembly step, PacBio data assembled using hifiasm. Contains haplotype 1 (hap1), haplotype 2 (hap2), and alternative contigs. | NA | NA |
| assembly_yahs/ | Second assembly step, the hifiasm assembly scaffolded using yahs and HiC data. Contains hap1 and hap2. | NA | NA |
| assembly_tiara/ | Third assembly step, the scaffolded assembly decontaminated using tiara. Contains hap1 and hap2, but also the combined scaffolds of hap1 and hap2. | NA | NA |
| assembly_curated/ | 'Final' genome assembly after manual curation using HiC contact matrices. Contains fasta.gz files of pseudo-chromosomes or scaffolds. | NA | NA |

Example: Assembly OG88, *Enoplosus armatus*

```
.
└── Enoplosus_armatus
    └── fEnoArm2
        ├── assembly_curated
        │   ├── OG88G_v230728.hic1.3.curated.hap1.chr_level.fa
        │   └── OG88G_v230728.hic1.3.curated.hap2.chr_level.fa
        ├── assembly_hifiasm
        │   ├── OG88G_v230728.hic1.0.hifiasm.alt_ctg.fasta
        │   ├── OG88G_v230728.hic1.0.hifiasm.hap1.p_ctg.fasta
        │   ├── OG88G_v230728.hic1.0.hifiasm.hap1.p_ctg.gfa
        │   ├── OG88G_v230728.hic1.0.hifiasm.hap2.p_ctg.fasta
        │   ├── OG88G_v230728.hic1.0.hifiasm.hap2.p_ctg.gfa
        │   ├── OG88G_v230728.hic1.0.hifiasm.p_ctg.gfa
        │   └── OG88G_v230728.hic1.0.hifiasm.pri_ctg.fasta
        ├── assembly_tiara
        │   ├── OG88G_v230728.hic1.2.tiara.hap1.hap2_combined_scaffolds.fa
        │   ├── OG88G_v230728.hic1.2.tiara.hap1_scaffolds.fa
        │   └── OG88G_v230728.hic1.2.tiara.hap2_scaffolds.fa
        ├── assembly_yahs
        │   ├── OG88G_v230728.hic1.1.yahs.hap1_scaffolds.fa
        │   ├── OG88G_v230728.hic1.1.yahs.hap2_scaffolds.fa
        │   └── OG88G_v230728.hic1.1.yahs.pri_scaffolds.fa
        └── genomic_data
            ├── hic
            │   ├── OG88L.230801.hic.R1.fastq.gz
            │   └── OG88L.230801.hic.R2.fastq.gz
            ├── illumina
            │   ├── OG88.ilmn.230814.R1.fastq.gz
            │   └── OG88.ilmn.230814.R2.fastq.gz
            └── pacbio
                ├── OG88G_m64497e_230726_221540.hifi_reads.bc2017--bc2017.bam
                ├── OG88G_m64497e_230726_221540.unbarcoded.hifi_reads.bam
                ├── OG88G_m64497e_230728_091325.hifi_reads.bc2017--bc2017.bam
                └── OG88G_m64497e_230728_091325.unbarcoded.hifi_reads.bam
```

# OceanOmics eDNA

The OceanOmics eDNA data aims to learn how to use environmental DNA (eDNA) to assess ecosystem health. Its main product is metabarcoding raw reads and metagenomics raw reads.

The data is structured by OceanOmics voyage, several FASTQ files per OceanOmics expedition (voyage). The S3 bucket contains one folder per voyage, with each voyage folder containing paired end reads of every sample of the voyage. Each voyage is split into several folders, one folder per chosen assay (12S, 16S, etc. pp.). Each assay contains one folder named `Unknown` which contains demultiplexed reads that could not be assigned to a known barcode-pair.

For example,

```
.
└── V10_CKI_P1
    ├── 12S
    │   ├── V10_CKI_V_6_1.R1.fq.gz
    │   ├── V10_CKI_V_6_1.R2.fq.gz
    │   ├── V10_CKI_V_6_2.R1.fq.gz
    │   ├── V10_CKI_V_6_2.R2.fq.gz
    │   ├── Unknown
    │   │   ├── F10-F10.R1.fq.gz
    │   │   ├── F10-F10.R2.fq.gz
    │   │   ├── F10-F11.R1.fq.gz
    │   │   ├── F10-F11.R2.fq.gz
    |   |   ├── ..... (more libraries)
    |   ..... (more libraries)
    ├── 16S
    │   ├── V10_CKI_V_6_1.R1.fq.gz
    │   ├── V10_CKI_V_6_1.R2.fq.gz
    │   ├── V10_CKI_V_6_2.R1.fq.gz
    │   ├── V10_CKI_V_6_2.R2.fq.gz
    │   ├── Unknown
    │   │   ├── F10-F10.R1.fq.gz
    │   │   ├── F10-F10.R2.fq.gz
    │   │   ├── F10-F11.R1.fq.gz
    │   │   ├── F10-F11.R2.fq.gz
    |   |   ├── ..... (more libraries)
    |   ..... (more libraries)
```


