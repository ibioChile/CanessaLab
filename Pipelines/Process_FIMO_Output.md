# Process FIMO output

This pipelines filter motifs identified with the tool FIMO to keep those located at specific positions in the genome. In this case, we will select motifs identified in the genome of Botrytis cinerea to select those located within the intergenic-1000bp region upstream of Botrytis genes.

## 1. Process Botrytis GFF file

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

## #Select 'genes' from Botrytis bed file:
```
awk '$8 == "gene"' Botrytis_cinerea.ASM83294v1.47.bed  > Botrytis_cinerea.ASM83294v1.47.genes.bed
```
## Generate forward and reverse Botrytis genes files

    awk '$6=="-"' Botrytis_cinerea.ASM83294v1.47.genes.bed > Botrytis_cinerea.ASM83294v1.47.genes.rev.bed
    awk '$6=="+"' Botrytis_cinerea.ASM83294v1.47.genes.bed	> Botrytis_cinerea.ASM83294v1.47.genes.fwd.bed

## Create file with length of each chromosome of Botrytis

```
samtools faidx genome_fasta
cut -f1,2 genome_fasta.fai > genome_fasta.contig.size
```

## Sort Botrytis genes BED file according to genomes name

```
sortBed -g  genome_fasta.contig.size -i Botrytis_cinerea.ASM83294v1.47.genes.bed > Botrytis_cinerea.ASM83294v1.47.genes.sorted.bed```
```

## Create file with intergenic regions in Botrytis

```
complementBed -i Botrytis_cinerea.ASM83294v1.47.genes.sorted.bed -g genome_fasta.contig.size > Botrytis_cinerea.ASM83294v1.47.intergenic.bed
```

## 2. Process FIMO output

## FIMO output GFF files are transformed to BED file

```
mkdir Default_bed/
for folder in Default/finalmemeuni*.txt/; do base=${folder##*/finalmemeuni}; convert2bed -i gff < $folder/fimo.gff >  Default_bed/${base%.*}.bed; done
```

## 3. Option1: Search motifs located on the same strand than gene

## Generate forward and reverse motif files

```
mkdir StrandedMotifs/
for file in Default_bed/*.bed; do i=${file##*/}; awk '$6=="-"' $file  > StrandedMotifs/${i%.*}.rev.bed; awk '$6=="+"' $file  > StrandedMotifs/${i%.*}.fwd.bed; done
```

## Filter motifs to the ones located in intergenic regions

```
mkdir MotifsIntergenic_ss/
for file in StrandedMotifs/*rev.bed; do i=${file##*/}; bedtools intersect -wa -a $file -b Botrytis_cinerea.ASM83294v1.47.intergenic.bed > MotifsIntergenic_ss/${i%.*}.intergenic.bed; done
for file in StrandedMotifs/*fwd.bed; do i=${file##*/}; bedtools intersect -wa -a $file -b Botrytis_cinerea.ASM83294v1.47.intergenic.bed > MotifsIntergenic_ss/${i%.*}.intergenic.bed; done
```

## Find genes with motifs located up to 1000bp of each gene

```
mkdir MotifsUpstream1000bp_ss
for file in MotifsIntergenic_ss/*fwd.intergenic.bed; do i=${file##*/}; bedops --range -1000:0 --everything Botrytis_cinerea.ASM83294v1.47.genes.fwd.bed | bedmap --echo --echo-map-id --echo-map-range - $file | bedops --range 1000:0 --everything - > MotifsUpstream1000bp_ss/${i%.*}.1000bp.bed; done
for file in MotifsIntergenic_ss/*rev.intergenic.bed; do i=${file##*/}; bedops --range -0:1000 --everything Botrytis_cinerea.ASM83294v1.47.genes.rev.bed | bedmap --echo --echo-map-id --echo-map-range - $file | bedops --range 0:-1000 --everything - > MotifsUpstream1000bp_ss/${i%.*}.1000bp.bed; done
```

## Concatenate results

```
mkdir Results_Default_SameStrand/
for file in MotifsUpstream1000bp_ss/*fwd*; do i=${file##*/}; cat $file ${file%fwd*}rev.intergenic.1000bp.bed > Results_Default_SameStrand/${i%.*.*.*.*}.intergenic.1000bp.ss.bed; done

```

## 4. Option2: Search motifs located on both strands, upstream of genes

## Filter motifs to the ones located in intergenic regions

    mkdir MotifsIntergenic_bs/
    for file in Default_bed/*.bed; do i=${file##*/}; bedtools intersect -wa -a $file -b Botrytis_cinerea.ASM83294v1.47.intergenic.bed > MotifsIntergenic_bs/${i%.*}.intergenic.bed; done

## Find genes with motifs located up to 1000bp of each gene

    mkdir MotifsUpstream1000bp_bs
    for file in MotifsIntergenic_bs/*.intergenic.bed; do i=${file##*/}; bedops --range -1000:0 --everything Botrytis_cinerea.ASM83294v1.47.genes.fwd.bed | bedmap --echo --echo-map-id --echo-map-range - $file | bedops --range 1000:0 --everything - > MotifsUpstream1000bp_bs/${i%.*}.fwd.1000bp.bed; done
    for file in MotifsIntergenic_bs/*.intergenic.bed; do i=${file##*/}; bedops --range -0:1000 --everything Botrytis_cinerea.ASM83294v1.47.genes.rev.bed | bedmap --echo --echo-map-id --echo-map-range - $file | bedops --range 0:-1000 --everything - > MotifsUpstream1000bp_bs/${i%.*}.rev.1000bp.bed; done

## Concatenate results

    mkdir Results_Default_BothStrands/
    for file in MotifsUpstream1000bp_bs/*fwd*; do i=${file##*/}; cat $file ${file%fwd*}rev.1000bp.bed > Results_Default_BothStrands/${i%.*.*.*.*}.intergenic.1000bp.bs.bed; done


## Create a summary of results
Use [this python script](https://github.com/ibioChile/CanessaLab/blob/master/Scripts/Summary_Filtered_FIMO.py) to create a summary of resulting bed files. Before running, modify input folder (contining bed files) and output file in script.
