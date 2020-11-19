import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO

TFs = "/Users/pamelacamejo/Documents/IBIO/Paulo_Canessa/Projects/Botrytis_PWM/To_PWMs/IDs_Tatroviride_TFs.txt"
InterProScan = "/Users/pamelacamejo/Documents/IBIO/Paulo_Canessa/Projects/Botrytis_PWM/To_PWMs/Tatroviridev2_FrozenGeneCatalog_20100319.proteins.fasta.tsv"
IDs_Interpro_PFAM = "/Users/pamelacamejo/Documents/IBIO/Paulo_Canessa/Projects/Botrytis_PWM/To_PWMs/IDs_Interpro_PFAM.txt"
Fasta_TF = "/Users/pamelacamejo/Documents/IBIO/Paulo_Canessa/Projects/Botrytis_PWM/To_PWMs/Total_Triat2_TFs.fasta"

# Read TF list
with open(TFs) as TF:
    TF_list = TF.read().splitlines()

# Read INTERPRO and PFAM IDs and convert to lists
df_IPRO_PFAM = pd.read_csv(IDs_Interpro_PFAM, sep="\t")
rows = df_IPRO_PFAM.apply(lambda x: x.tolist(), axis=0)
PFAM_ID = rows['Pfam'].tolist()
PFAM_ID = [x for x in PFAM_ID if x != 'nan']
InterPRO_ID = rows['InterPRO'].tolist()
InterPRO_ID.extend(rows['Interpro alternativo'].tolist())
InterPRO_ID = [x for x in InterPRO_ID if x != 'nan']

# Read InterPro Scan results and filter by TF, INTERPRO and PFAM Ids
DBD_list = {k: [] for k in ['GeneID','DB','DB_ID','start','end']}
DBD_info = []
with open(InterProScan)as f:
    for line in f:
        L = line.strip().split("\t")
        if L[0] in TF_list:
            if L[4] in PFAM_ID:
                info = [L[0],L[6],L[7]]
                if info not in DBD_info:
                    DBD_info.append(info)
                    DBD_list["GeneID"].append(L[0])
                    DBD_list["DB"].append(L[3])
                    DBD_list["DB_ID"].append(L[4])
                    DBD_list["start"].append(L[6])
                    DBD_list["end"].append(L[7])
            elif (len(L) > 11 and L[11] in InterPRO_ID):
                info = [L[0], L[6], L[7]]
                if info not in DBD_info:
                    DBD_info.append(info)
                    DBD_list["GeneID"].append(L[0])
                    DBD_list["DB"].append(L[3])
                    DBD_list["DB_ID"].append(L[11])
                    DBD_list["start"].append(L[6])
                    DBD_list["end"].append(L[7])

# Create table with TFs information of DBD
DBD_df = pd.DataFrame.from_dict(DBD_list)
DBD_df.to_csv(r'/Users/pamelacamejo/Documents/IBIO/Paulo_Canessa/Projects/Botrytis_PWM/To_PWMs/DBD_Tricho_df.csv', index=False)

# Extract DBD from aminoacids fasta file
DBD_frag_list = []

for i in range(int(len(DBD_list["GeneID"]))):
    Gene = DBD_list["GeneID"][i]
    start = DBD_list["start"][i]
    end = DBD_list["end"][i]
    ID = DBD_list["DB_ID"][i]
    TF_seqs = SeqIO.parse(Fasta_TF, "fasta")
    for record in TF_seqs:
        if Gene == record.id:
            DBD = record.seq[int(start) - 1:int(end)]
            header = record.id + "_" + ID + "_" + start + "_" + end
            fragment = SeqRecord(DBD, header, "", "")
            DBD_frag_list.append(fragment)

# Create fasta with DBDs
fasta_out = FastaIO.FastaWriter("/Users/pamelacamejo/Documents/IBIO/Paulo_Canessa/Projects/Botrytis_PWM/To_PWMs/DBD_Tricho.fasta", wrap=None)
fasta_out.write_file(DBD_frag_list)