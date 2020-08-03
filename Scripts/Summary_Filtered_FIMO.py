import glob
import os
import pandas as pd

arr = glob.glob("/FIMO_filtered_Results/*.bed")

motifs_info = []
for bed_file in arr:
    motif_name = os.path.split(bed_file)[1].split(".")[0]
    with open(bed_file)as f:
        for line in f:
            L = line.strip().split()
            Genome = L[0]
            Gene_ID = str(L[3]).split(":")[1]
            Gene_start = L[1]
            Gene_end = L[2]
            Strand = L[5]
            motif_ids = str(L[9]).split("|")[1]
            if motif_ids:
                motif_id = str(motif_ids).split(";")
                for mid in motif_id:
                    motifs_info.append([motif_name, Genome, Gene_ID, Gene_start, Gene_end, Strand])
            else:
                motifs_info.append(["Not Motif", Genome, Gene_ID, Gene_start, Gene_end, Strand])

df = pd.DataFrame(motifs_info, columns=['motif_name', 'Genome', 'Gene_ID', 'Gene_start', 'Gene_end', 'Strand'])
df_pivot = df.pivot_table(index=['Genome', 'Gene_ID', 'Gene_start', 'Gene_end', 'Strand'], columns=['motif_name'],
                                  aggfunc=len).fillna(0)
df_pivot_filt = df_pivot.drop(labels='Not Motif', axis=1)
df_pivot_filt.to_csv(r'FIMO_Bot_intergenic_1000bp_summary.csv', index=True)
