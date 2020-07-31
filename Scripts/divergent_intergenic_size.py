from statistics import mean, stdev
import pandas as pd

prev_chrom = ''
prev_strand = ''
prev_gene_id = ''
inter_size_list = []
inter_size_df = []

gff_bot = open("Botrytis_cinerea.ASM83294v1.47.genes.bed",'r')

for current_line in gff_bot:
    L = current_line.strip().split()
    gene_id = str(L[3]).split(":")[1]
    chromosome = L[0]
    strand = L[5]
    start = L[1]
    end = L[2]
    if chromosome == prev_chrom:
        # If genes are divergent:
        if prev_strand == "-" and strand == "+":
            start_inter = prev_end
            stop_inter = start
            inter_size = float(start) - float(prev_end)
            inter_size_list.append(inter_size)
            inter_size_df.append([prev_gene_id, gene_id, inter_size])
    prev_chrom = chromosome
    prev_strand = strand
    prev_gene_id = gene_id
    prev_start = start
    prev_end = end

print("Average of Divergent intergenic size: " + str(round(mean(inter_size_list))) + "+-" + str(round(stdev(inter_size_list))) + "bp" )

df = pd.DataFrame(inter_size_df, columns=['Gene1', 'Gene2', 'Intergenic size (bp)'])
df.to_csv(r'/Users/pamelacamejo/Documents/IBIO/Paulo_Canessa/FIMO_Botrytis/Results/Intergenic_Divergent_Genes_Bot.csv', index=False)
