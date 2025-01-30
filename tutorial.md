# Ocean Genomes

Here are several examples on what one can do with the genome assembly data.

## Raw reads

The most common use case for raw reads is to assembly them into a genome assembly.

### Genome assembly

Here's a simple example to assemble the raw reads into a genome assembly using [MEGAHIT](https://github.com/voutcn/megahit):

    megahit -1 OG88.ilmn.230814.R1.fastq.gz -2 OG88.ilmn.230814.R2.fastq.gz -o OG88_assembly
    
This will generate a folder named `OG88_assembly/` and within that folder, `final_contigs.fa` will contain the genome assembly.

### Sourmash-based diversity studies



## Final curated assembly

The most common use case for a final curated genome assembly is to use it as a reference for downstream tasks such as taxonomic assignment.

### Use as BLAST reference database

### Add to an existing Kraken2 database

### Search for genes

We can also annotate the genome assembly to find genes in the genome. Let's use [Tiberius](https://github.com/Gaius-Augustus/Tiberius), a Deep Learning-based gene predicter.

     python bin/tiberius.py --genome OG88G_v230728.hic1.3.curated.hap1.chr_level.fa --out OG88G_v230728.hic1.3.curated.hap1.chr_level.annotation.gtf --no_softmasking

This will result in a gtf file detailing where gene models sit in the genome assembly. These gene models can then be extracted for further funcational annotation, for example, by assigning Gene Ontology terms via [PANNZER2](http://ekhidna2.biocenter.helsinki.fi/sanspanz/) or similar tools.
