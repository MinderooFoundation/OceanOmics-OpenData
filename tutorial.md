This file describes several potential use cases for the data in the OceanOmics AWS Open Data account, ranging from assemblies to eDNA.

- [Ocean Genomes - working with assemblies](#ocean-genomes---working-with-assemblies)
- [eDNA - working with environmental DNA data](#eDNA---working-with-environmental-DNA-data)

# Ocean Genomes - working with assemblies

Here are several examples on what one can do with the genome assembly data.

## Raw reads

The most common use case for raw reads is to assemble them into a genome assembly.

### Genome assembly

Here's a simple example to assemble the raw reads into a genome assembly using [MEGAHIT](https://github.com/voutcn/megahit):

    megahit -1 OG88.ilmn.230814.R1.fastq.gz -2 OG88.ilmn.230814.R2.fastq.gz -o OG88_assembly
    
This will generate a folder named `OG88_assembly/` and within that folder, `final_contigs.fa` will contain the genome assembly.

### Sourmash-based diversity studies

You can add these reads as a reference to a sourmash database - that will let you search your metagenomics datasets while including the genome as a potential reference. You can answer questions like, 'does my environmental sample contain DNA from *Enoplosus armatus*?' See the [sourmash manual](https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html) for more information.

    sourmash sketch dna -p k=31 OG88.ilmn.230814.R1.fastq.gz OG88.ilmn.230814.R2.fastq.gz \
        --name "sample" -o fEnoArm2.sig

That will result in a 'sketch' of the DNA we have sequenced for OG88, *E. armatus*.  You can query this sketch using your metagenomics library:

First, we sketch your metagenome sequencing reads:

    sourmash sketch dna -p scaled=10000,k=31 your_metagenome.fastq.gz -o your_metagenome.sig
    
    sourmash search your_metagenome.sig fEnoArm2.sig --containment

The --containment flag calculates what percentage of your metagenome is contained in the *E. armatus* genome assembly. Chances are it's tiny! Fish DNA always makes up a tiny fraction of a marine water sample (<<1%).

## Final curated assembly

The most common use case for a final curated genome assembly is to use it as a reference for downstream tasks such as taxonomic assignment.

### Use as BLAST reference database

You can either add the assembly to an existing fasta file, or you can make a BLAST database directly:

    makeblastdb -dbtype nucl -in ./fEnoArm2/assembly_curated/OG88G_v230728.hic1.3.curated.hap1.chr_level.fa

This will generate a BLAST reference database. Here's a way to query that:

    blastn -query your_query.fasta -d OG88G_v230728.hic1.3.curated.hap1.chr_level.fa

You should see some hits appear.

Other taxonomic classifiers or sequence comparison tools can be used with this assembly - MMSeqs2 and Kraken2 are commonly used for this task.

### Search for genes

We can also annotate the genome assembly to find genes in the genome. Let's use [Tiberius](https://github.com/Gaius-Augustus/Tiberius), a Deep Learning-based gene predicter.

     python bin/tiberius.py --genome OG88G_v230728.hic1.3.curated.hap1.chr_level.fa --out OG88G_v230728.hic1.3.curated.hap1.chr_level.annotation.gtf --no_softmasking

After a few hours this will result in a gtf file detailing where gene models sit in the genome assembly. These gene models can then be extracted for further functional annotation, for example, by assigning Gene Ontology terms via [PANNZER2](http://ekhidna2.biocenter.helsinki.fi/sanspanz/) or similar tools. You can then look for your gene families of interest - for example, by searching for the GO term GO:0042288 'MHC class I protein binding' you will find genes implicated in disease resistance.

### Search for relatives

We can grab a random piece of DNA and search for it in public databases.

For example, here's a small region from the second pseudomolecule of *Enoplosus armatus*:

<pre>
CACACTGTTCATTATTTGAGTCATTCATATTTTAGCTATATTCTTTCTGACCGAGTCTGA
CTGAGCTCCAAACGTACAGTGAAAGCTGAAGTTGGAATTTACCAGCTGGAGACTAGCATC
TTTATCAGTGCCCCTCCGAGAGCTGCATGTGTCGTCACAATAAGATGTGACAGGCGACTG
GAAAGAGATTTACACCTGGTTAATAAGCTTTCCTACTAACTAATCATGAAACTTTATCGC
TCATCTCTGCTCCTGTGTTCAACTATTACCACGTTAACATCTTTCCAATCAGTGCTTTCA
CTTCTCCTCCTCTCAACACCCCCTTCACTTGCTTCAGTTGTCTCTGTCCCTTCACCTGTA
ACAGTTTTCATAAATAAATAAAAAAAGCATTGCACAGTCACAGGTCAGTGTATCTACCCA
CACAGAATATTCAACCACTGAAAAGAAATGACAAAAATACAGAGACATGCCAAACACCTT
</pre>

You can put that into online [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to see what it is:
![BLAST results](https://github.com/user-attachments/assets/11ebc286-cc07-4fa2-9088-5a6efe9ff5e1)

That's all tumor suppressing genes in various fish, almost all distant relatives of *Enoplosus* showing how little diversity of fish is present in public genomic databases. The best hit, *Siniperca*, is a Chinese freshwater fish from a different family, but the same order.

You can put this DNA into [Logan](https://logan-search.org/), which searches all 'raw' sequencing data and plots where the samples are from. This particular sequence hits nothing! No metagenomics project has sampled a similar fish (fish DNA makes up very little of metagenomics data, the chance to also hit this particular piece of DNA is very, very low).

# eDNA - working with environmental DNA data

We can run similar analyses with the raw reads in the eDNA data. There are three types of data in the eDNA data: 12S metabarcoding data, which is a specific region of the 12S gene in the mitogenome, 16S metabarcoding, the same but for 16S, and shotgun metabarcoding, which is 'raw' sequencing data from everything found in a water sample.

Let' pull out one of the 12S reads, the first read from V10_CKI_N_1_5_R1:

<pre>
@VH00640:1:AACCK23M5:1:1101:55591:1057 1:N:0:1
GTCGGTAAAACTCGTGCCAGCCACCGCGGTTATACGAGAGGCCCAAGTCGATAGGCATACGGCGTAAAGAGTGGTTAAGGAACATTTACACTAGGGTCGAATATTTTCAAAGCCGTCATACGCCTATGAGAATGAGAAACTCA
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CC;CCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCC;CCCC-CC;C-CCCCCC-CCCCCCCC-C-CCC-C-CCC;C;
</pre>

The second row is the DNA sequence. Let's BLAST that:

![image](https://github.com/user-attachments/assets/c393cd2e-6ef7-4a68-833d-2c223a252724)

The first hit is for *Plectrogenium kamoharai*, but only with an identity of 89%, too low for most researchers in this space. We prefer more than 97%. So this does look like some kind of fish, but we don't know which one. 

Let's put this read into Logan as well:

![image](https://github.com/user-attachments/assets/83438ce9-399d-41b6-90a3-e227b7552008)

Interesting! Many hits, most of them decent. You can see all of those hits next to Australia - all of these hits are OceanOmics samples we have previously uploaded to SRA. The read we searched is within this dataset, too. Nobody else has detected this fish. You can see how little marine metagenomics data there is!

There's an interesting hit next to the US, but that seems to be an error somewhere as the metadata itself indicates this to be a OceanOmics sample from Australia.

At this point you would take these entire datasets and run them through an amplicon/metabarcoding analysis pipeline. OceanOmics has a pipeline, too: [OceanOmics-amplicon-nf](https://github.com/MinderooFoundation/OceanOmics-amplicon-nf/). These pipelines merge duplicated reads into amplicon sequence variants (ASVs), assign taxonomic labels, and tell you which species are present in entire sample.
