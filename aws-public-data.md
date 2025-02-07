# OceanOmics Public Datasets on Amazon Web Services

This document describes the format and gives simple examples for getting started with the OceanOmics data stored on AWS.

# Ocean Genomes

The OceanOmics Ocean Genomes data aims to sequence all marine vertebrates. Its main product is single individual-genome sequencing data and genome assemblies.

Data is stored in FASTA, BAM, and FASTQ format in the following buckets:

(to fill once we have bucket names)

The structure of the data follows the [GenomeArk structure](https://genomeark.s3.amazonaws.com/index.html). 

One folder per species with the accepted scientific name,  
this folder contains samples in the form of Tree of Life IDs (ToLID), 
then several folders containing genome assembly and sequencing data:

| Folder name | Folder contents |
| ------------- | ------------- | 
| assembly_curated/ | 'Final' genome assembly after manual curation using HiC contact matrices. Contains fasta.gz files of pseudo-chromosomes or scaffolds. | 
| assembly_oceangenomes_v1 | contains evaluation data and pre-curation assembly steps | 
| assembly_oceangenomes_v1/evaluation/busco | Contains BUSCO genome completeness results |
| assembly_oceangenomes_v1/evaluation/genomescope | Contains GenomeScope genome size estimate results | 
| assembly_oceangenomes_v1/evaluation/{hap1,hap2} | Contains BUSCO results per haplotype | 
| assembly_oceangenomes_v1/evaluation/merqury | Contains Merqury genome quality evaluation results |
| assembly_oceangenomes_v1/intermediates/ | Contains hifiasm assembled contigs and contig graphs for hap1 and hap2 | 
| genomic_data/ | contains raw reads of all data types |
| genomic_data/illumina | raw Illumina paired ends reads, fastq.gz | 
| gneomic_data/hic | raw HiC data in fastq.gz format |
| genomic_data/pacbio_hifi | raw PacBio HiFi reads in bam and fastq.ggz format |

Example: Assembly fArrgeo1, *Arripis georgianus*

```
Arripis_georgianus/
└── fArrGeo1
    ├── assembly_MT
    │   ├── fArrGeo1.20230529.mitohifi.fasta
    |   └── fArrGeo1.20230529.mitohifi.gb
    ├── assembly_curated
    │   ├── fArrGeo1.hap1.cur.20230529.agp
    │   ├── fArrGeo1.hap1.cur.20230529.fasta.gz
    │   ├── fArrGeo1.hap1.cur.20230529.gfastats.txt
    │   ├── fArrGeo1.hap1.cur.map.pretext
    │   ├── fArrGeo1.hap1.cur.pretext_snapshot.png
    │   ├── fArrGeo1.hap2.cur.20230529.agp
    │   ├── fArrGeo1.hap2.cur.20230529.fasta.gz
    │   ├── fArrGeo1.hap2.cur.20230529.gfastats.txt
    │   ├── fArrGeo1.hap2.cur.map.pretext
    │   └── fArrGeo1.hap2.cur.pretext_snapshot.png
    ├── assembly_oceangenomes_v1
    │   ├── evaluation
    │   │   ├── busco
    │   │   │   ├── fArrGeo1.hap1.cur.busco_full_table.tsv
    │   │   │   ├── fArrGeo1.hap1.cur.busco_short_summary.txt
    │   │   │   ├── fArrGeo1.hap1.cur.missing_busco_list.tsv
    │   │   │   ├── fArrGeo1.hap1.hap2.cur.busco_figure.png
    │   │   │   ├── fArrGeo1.hap2.cur.busco_full_table.tsv
    │   │   │   ├── fArrGeo1.hap2.cur.busco_short_summary.txt
    │   │   │   └── fArrGeo1.hap2.cur.missing_busco_list.tsv
    │   │   ├── genomescope
    │   │   │   ├── fArrGeo1_genomescope_linear_plot.png
    │   │   │   ├── fArrGeo1_genomescope_log_plot.png
    │   │   │   ├── fArrGeo1_genomescope_model.txt
    │   │   │   ├── fArrGeo1_genomescope_progress.txt
    │   │   │   ├── fArrGeo1_genomescope_summary.txt
    │   │   │   ├── fArrGeo1_genomescope_transformed_linear_plot.png
    │   │   │   └── fArrGeo1_genomescope_transformed_log_plot.png
    │   │   ├── hap1
    │   │   │   ├── busco
    │   │   │   │   ├── fArrGeo1.hap1.hap2.s2.busco_figure.png
    │   │   │   │   ├── fArrGeo1.hap1.s2.busco_full_table.tsv
    │   │   │   │   ├── fArrGeo1.hap1.s2.busco_short_summary.txt
    │   │   │   │   └── fArrGeo1.hap1.s2.missing_busco_list.tsv
    │   │   │   ├── gfastats
    │   │   │   │   └── fArrGeo1.hap1.s2.gfastats.txt
    │   │   │   └── pretext
    │   │   │       ├── fArrGeo1.hap1.s2.map.pretext
    │   │   │       └── fArrGeo1.hap1.s2.pretext_snapshot.png
    │   │   ├── hap2
    │   │   │   ├── busco
    │   │   │   │   ├── fArrGeo1.hap2.s2.busco_full_table.tsv
    │   │   │   │   ├── fArrGeo1.hap2.s2.busco_short_summary.txt
    │   │   │   │   └── fArrGeo1.hap2.s2.missing_busco_list.tsv
    │   │   │   ├── gfastats
    │   │   │   │   └── fArrGeo1.hap2.s2.gfastats.txt
    │   │   │   └── pretext
    │   │   │       ├── fArrGeo1.hap2.s2.map.pretext
    │   │   │       └── fArrGeo1.hap2.s2.pretext_snapshot.png
    │   │   └── merqury
    │   │       ├── fArrGeo1_png
    │   │       │   ├── fArrGeo1.hap1.chr_level.spectra-cn.fl.png
    │   │       │   ├── fArrGeo1.hap1.chr_level.spectra-cn.ln.png
    │   │       │   ├── fArrGeo1.hap1.chr_level.spectra-cn.st.png
    │   │       │   ├── fArrGeo1.hap2.chr_level.spectra-cn.fl.png
    │   │       │   ├── fArrGeo1.hap2.chr_level.spectra-cn.ln.png
    │   │       │   ├── fArrGeo1.hap2.chr_level.spectra-cn.st.png
    │   │       │   ├── fArrGeo1.spectra-asm.fl.png
    │   │       │   ├── fArrGeo1.spectra-asm.ln.png
    │   │       │   ├── fArrGeo1.spectra-asm.st.png
    │   │       │   ├── fArrGeo1.spectra-cn.fl.png
    │   │       │   ├── fArrGeo1.spectra-cn.ln.png
    │   │       │   └── fArrGeo1.spectra-cn.st.png
    │   │       ├── fArrGeo1_qv
    │   │       │   ├── fArrGeo1.hap1.chr_level.qv
    │   │       │   ├── fArrGeo1.hap2.chr_level.qv
    │   │       │   └── fArrGeo1.summary.qv
    │   │       └── fArrGeo1_stats
    │   │           └── fArrGeo1.completeness.stats
    │   ├── fArrGeo1.hap1.hap2.s2.fasta.gz
    │   └── intermediates
    │       ├── fArrGeo1_hap1_contigs.fasta.gz
    │       ├── fArrGeo1_hap2_contigs.fasta.gz
    │       ├── hifiasm
    │       │   ├── fArrGeo1_hap1_contig_graph.gfa.gz
    │       │   └── fArrGeo1_hap2_contig_graph.gfa.gz
    │       └── meryl
    │           └── fArrGeo1_v230525_meryldb.tar.gz
    └── genomic_data
        ├── hic
        │   ├── OG95H_D_HICL_S4_L001_R1_001.fastq.gz
        │   ├── OG95H_D_HICL_S4_L001_R2_001.fastq.gz
        │   └── fArrGeo1_hic_files.md5
        └── pacbio_hifi
            ├── OG95H_m64497e_230528_010643.hifi_reads.bc2010--bc2010.bam
            ├── OG95H_m64497e_230528_010643.hifi_reads.bc2010--bc2010.filt.fastq.gz
            ├── OG95H_m64497e_230528_010643.unbarcoded.hifi_reads.bam
            ├── OG95H_m64497e_230528_010643.unbarcoded.hifi_reads.filt.fastq.gz
            ├── OG95H_m64497e_230529_102918.hifi_reads.bc2010--bc2010.bam
            ├── OG95H_m64497e_230529_102918.hifi_reads.bc2010--bc2010.filt.fastq.gz
            ├── OG95H_m64497e_230529_102918.unbarcoded.hifi_reads.bam
            ├── OG95H_m64497e_230529_102918.unbarcoded.hifi_reads.filt.fastq.gz
            └── fArrGeo1_hifi_files.md5

```

# OceanOmics eDNA

The OceanOmics eDNA data aims to learn how to use environmental DNA (eDNA) to assess ecosystem health. Its main product is metabarcoding raw reads and metagenomics raw reads.

The data is structured by OceanOmics project (often but not always a voyage), several FASTQ files per OceanOmics expedition (voyage). The S3 bucket contains one folder per project code, with each project folder containing paired end reads of every sample of the project. Each project is split into several folders, one folder per chosen assay (12S, 16S, etc. pp.). Each assay contains one folder named `Unknown` which contains demultiplexed reads that could not be assigned to a known barcode-pair.

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


