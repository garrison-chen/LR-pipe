import re
import os
from Bio import SeqIO
import subprocess

#cmd_str = "mkdir toSeparate"
#subprocess.run(cmd_str, shell=True)

subprocess.run(["mkdir", snakemake.output[0]]) 

#clade_name = ""
#clade_records = ""

all_clade_names = []
input_protein_fasta = snakemake.input[0]
#input_protein_fasta = "Selected_Prot.fasta"
output_dir = snakemake.output[0]

print("--------------------")

for seq_record in SeqIO.parse(input_protein_fasta, "fasta"):
    
    # Parse protein name to get clade name                
    if (seq_record.id[0:4] == "HIPS"):
        each_clade = re.sub("_BOGUS_[0-9]{2,}", "", seq_record.id)
        each_clade = re.sub("_ctg.{2,}$", "", each_clade)
    elif (seq_record.id[0:2] == "SB"):
        each_clade = re.sub("_[0-9]{2,}", "", seq_record.id)
    elif (seq_record.id[0:4] == "NOSO"):
        each_clade = re.sub("_[0-9]{2,}", "", seq_record.id)
    else:
        each_clade = re.sub(".gb.{2,}", "", seq_record.id)
        each_clade = re.sub("_[a-z | A-Z].*.", "", each_clade) 

    
    all_clade_names.append(each_clade)

unique_clade_name = set(all_clade_names)
print(len(unique_clade_name))
#print(unique_clade_name)

counter = 1
for clade_name in unique_clade_name:
    
    output_proteins = ""
    
    for seq_record in SeqIO.parse(input_protein_fasta, "fasta"):
        
        # Parse protein name to get clade name                
        if (seq_record.id[0:4] == "HIPS"):
            each_clade = re.sub("_BOGUS_[0-9]{2,}", "", seq_record.id)
            each_clade = re.sub("_ctg.{2,}$", "", each_clade)
        elif (seq_record.id[0:2] == "SB"):
            each_clade = re.sub("_[0-9]{2,}", "", seq_record.id)
        elif (seq_record.id[0:4] == "NOSO"):
            each_clade = re.sub("_[0-9]{2,}", "", seq_record.id)
        else:
            each_clade = re.sub(".gb.{2,}", "", seq_record.id)
            each_clade = re.sub("_[a-z | A-Z].*.", "", each_clade) 

        record_protein = ">" + seq_record.id + "\n" + seq_record.seq + "\n"
        
        if each_clade  == clade_name:
            output_proteins += record_protein
            
    out_f_name = clade_name + ".fasta"
    
    with open(os.path.join(output_dir, out_f_name), "w") as file:
        file.write(str(output_proteins))
    print("{0:.0%}".format(counter/len(unique_clade_name)), end = "\r")
    counter += 1

