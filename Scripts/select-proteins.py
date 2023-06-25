import re
import os
import sys
from Bio import SeqIO
import pandas as pd

threshold = int(snakemake.params[0])
mmseqs_res_path = snakemake.input[0]

tsv_read = pd.read_csv(mmseqs_res_path, header = None, delim_whitespace=True)
representative_proteins = tsv_read[0]
all_proteins = tsv_read[1]

selected_protein_list = [] 
counter = 1

for i in range(0, all_proteins.size-1):
    if (representative_proteins[i] != "-"):
        each_cluster = []
        each_cluster.append(representative_proteins[i])
        for j in range(i+1,representative_proteins.size-1):
            if (representative_proteins[j] == representative_proteins[i]):
                each_cluster.append(all_proteins[j])
                representative_proteins[j] = "-"
            else:
                break

        clades_list = []
        for k in range(0, len(each_cluster)-1):
            each_protein = each_cluster[k]
            
            # Parse protein name to get clade name                
            if (each_protein[0:4] == "HIPS"):
                each_clade = re.sub("_BOGUS_[0-9]{2,}", "", each_protein)
                each_clade = re.sub("_ctg.{2,}$", "", each_clade)
            elif (each_protein[0:2] == "SB"):
                each_clade = re.sub("_[0-9]{2,}", "", each_protein)
            elif (each_protein[0:4] == "NOSO"):
                each_clade = re.sub("_[0-9]{2,}", "", each_protein)
            else:
                each_clade = re.sub(".gb.{2,}", "", each_protein)
                each_clade = re.sub("_[a-z | A-Z].*.", "", each_clade) 
            
            clades_list.append(each_clade)
            #print(each_clade)
            print("{0:.0%}".format(counter/len(all_proteins)), end = "\r")
            counter += 1

        clades_set = set(clades_list) 
        if (len(clades_set) >= int(threshold)):
            for each_protein in each_cluster:
                selected_protein_list.append(each_protein)

# Output Selected Proteins as Fasta file 
input_protein_folder = snakemake.input[1]

counter = 1
selected_prot = ""

for each_protein in selected_protein_list:
    # Parse protein name to get clade name                
    if (each_protein[0:4] == "HIPS"):
        each_clade = re.sub("_BOGUS_[0-9]{2,}", "", each_protein)
        each_clade = re.sub("_ctg.{2,}$", "", each_clade)
        each_clade = each_clade + ".fasta"
    elif (each_protein[0:2] == "SB"):
        each_clade = re.sub("_[0-9]{2,}", "", each_protein)
        each_clade = each_clade + ".fasta"
    elif (each_protein[0:4] == "NOSO"):
        each_clade = re.sub("_[0-9]{2,}", "", each_protein)
        each_clade = each_clade + ".fasta"
    else:
        each_clade = re.sub(".gb.{2,}", "", each_protein)
        each_clade = re.sub("_[a-z | A-Z].*.", "", each_clade) 
        each_clade = each_clade + ".fasta"
    
    #print(each_clade)
    input_protein_fasta = os.path.join(input_protein_folder, each_clade) 
    for seq_record in SeqIO.parse(input_protein_fasta, "fasta"):
        if (seq_record.id == each_protein):
            record_protein = ">" + seq_record.id + "\n" + seq_record.seq + "\n"
            selected_prot += record_protein
            print("{0:.0%}".format(counter/len(selected_protein_list)), end = "\r")
            counter += 1
            break 
        else:
            pass


    

#out_f_name = "Selected_Prot.fasta"
out_f_name = snakemake.output[0]
print("Outputting the selected proteins ...")
with open(out_f_name, "w") as file:
    file.write(str(selected_prot))
print("Finished.")

