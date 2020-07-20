# Process FIMO output

This pipelines filter motifs identified with the tool FIMO to keep those located at specific positions in the genome. In this case, we will select motifs identified in the genome of Botrytis cinerea to select those located within the intergenic-1000bp region upstream of Botrytis genes.

## Create conda environment to run the programs needed

```
conda create -n bed
conda activate bed
conda install -c bioconda bedops bedtools samtools=1.9
```  
## Botrytis GFF file is transformed to BED file

```
convert2bed -i gff < Botrytis_cinerea.ASM83294v1.47.gff3  > Botrytis_cinerea.ASM83294v1.47.bed
```

## Generate forward and reverse Botrytis genes files

```
awk '$6=="-" && $8 == "gene"' Botrytis_cinerea.ASM83294v1.47.bed  > Botrytis_cinerea.ASM83294v1.47.genes.rev.bed
awk '$6=="+" && $8 == "gene"' Botrytis_cinerea.ASM83294v1.47.bed  > Botrytis_cinerea.ASM83294v1.47.genes.fwd.bed
```

## Create file with length of each chromosome of Botrytis

```
samtools faidx genome_fasta
cut -f1,2 genome_fasta.fai > genome_fasta.contig.size
```

## Sort Botrytis genes BED file according to genomes name

```
sortBed -g  genome_fasta.contig.size -i Botrytis_cinerea.ASM83294v1.47.genes.rev.bed > Botrytis_cinerea.ASM83294v1.47.genes.rev.sorted.bed
sortBed -g  genome_fasta.contig.size -i Botrytis_cinerea.ASM83294v1.47.genes.fwd.bed > Botrytis_cinerea.ASM83294v1.47.genes.fwd.sorted.bed
```

## Create file with intergenic regions in Botrytis

```
complementBed -i Botrytis_cinerea.ASM83294v1.47.genes.rev.sorted.bed -g genome_fasta.contig.size > Botrytis_cinerea.ASM83294v1.47.rev.intergenic.bed
complementBed -i Botrytis_cinerea.ASM83294v1.47.genes.fwd.sorted.bed -g genome_fasta.contig.size > Botrytis_cinerea.ASM83294v1.47.fwd.intergenic.bed
```

## FIMO output GFF files are transformed to BED file

```
mkdir Default_bed/
for folder in Default/finalmemeuni*.txt/; do base=${folder##*/finalmemeuni}; convert2bed -i gff < $folder/fimo.gff >  Default_bed/${base%.*}.bed; done
```

## Generate forward and reverse motif files

```
mkdir StrandedMotifs/
for file in Default_bed/*.bed; do i=${file##*/}; awk '$6=="-"' $file  > StrandedMotifs/${i%.*}.rev.bed; awk '$6=="+"' $file  > StrandedMotifs/${i%.*}.fwd.be$
```

## Filter motifs to the ones located in intergenic regions

```
mkdir MotifsIntergenic/
for file in StrandedMotifs/*rev.bed; do i=${file##*/}; bedtools intersect -wa -a $file -b Botrytis_cinerea.ASM83294v1.47.rev.intergenic.bed > MotifsIntergen$
for file in StrandedMotifs/*fwd.bed; do i=${file##*/}; bedtools intersect -wa -a $file -b Botrytis_cinerea.ASM83294v1.47.fwd.intergenic.bed > MotifsIntergen$
```

## Find genes with motifs located up to 1000bp of each gene

```
mkdir MotifsUpstream1000bp
for file in MotifsIntergenic/*fwd.intergenic.bed; do i=${file##*/}; bedops --range -1000:0 --everything Botrytis_cinerea.ASM83294v1.47.genes.fwd.bed | bedma$
for file in MotifsIntergenic/*rev.intergenic.bed; do i=${file##*/}; bedops --range -1000:0 --everything Botrytis_cinerea.ASM83294v1.47.genes.rev.bed | bedma$
```

## Concatenate results

```
mkdir Results_Default/
for file in MotifsUpstream1000bp/*fwd*; do i=${file##*/}; cat $file ${file%fwd*}rev.intergenic.1000bp.bed > Results_Default/${i%.*.*.*.*}.intergenic.1000bp.b$
```
