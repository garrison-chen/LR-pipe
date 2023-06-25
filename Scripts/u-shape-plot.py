import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

mmseqs_res_path = snakemake.input[0]
tsv_read = pd.read_csv(mmseqs_res_path, header = None, delim_whitespace=True)
representative_proteins = tsv_read[0]
all_proteins = tsv_read[1]


# Loop through file to store each cluster as list of list
clades_counter_list = []
all_clades_list = []
counter = 1
for i in range(0,all_proteins.size-1):
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
        for k in range(0, len(each_cluster)):
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
            all_clades_list.append(each_clade)
            print("{0:.0%}".format(counter/len(all_proteins)), end = "\r")
            counter += 1
            
        clades_set = set(clades_list)
        clades_counter_list.append(len(clades_set))

all_clades_set = set(all_clades_list)


# plot the U-shaped Figure
figure(figsize=(10, 8), dpi=300)
max_cluster_size = max(clades_counter_list)

u_shape_bar_list = []
for i in range(1,max_cluster_size):
    counter = 1;
    for j in clades_counter_list:
        if (i == j):
            counter += 1
    u_shape_bar_list.append(counter)
    
br = np.arange(len(u_shape_bar_list))
plt.bar(br, u_shape_bar_list, width = 0.9, color = (0.3,0.1,0.4,0.6)) 
plt.title("U-shaped plot (SI 0.7)", fontweight ='bold', fontsize = 13) 
plt.xlabel("unique clade number per cluster", fontweight ='bold', fontsize = 13) 
plt.ylabel("Count", fontweight ='bold', fontsize = 13)
plt.savefig(snakemake.output[0])
