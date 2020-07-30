from statistics import mean, stdev

prev_chrom = ''
prev_strand = ''
inter_size_list = []

gff_bot = open("Botrytis_cinerea.ASM83294v1.47.genes.bed",'r')

for current_line in gff_bot:
    L = current_line.strip().split()
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

    prev_chrom = chromosome
    prev_strand = strand
    prev_start = start
    prev_end = end

print("Average of Divergent intergenic size: " + str(round(mean(inter_size_list))) + "+-" + str(round(stdev(inter_size_list))) + "bp" )
