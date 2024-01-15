import pandas as pd

import os
from pathlib import Path
import shutil

# path to the directory that contains "combined_count_table.txt", one of the outputs that GraftM generates for marker recovery, for all samples
# need sample name as each file name 
txt_dir = "EDITHERE"

# path to output txt file 
output_path_filename = "EDITHERE"

os.chdir(txt_dir)

ITS2_hit_num = []
sample_ITS2_dict = {}
its_hit_count = 0
for file in os.listdir():
    with open(file) as f:
        ITS2_dict = {}
        sample_name = file.split(".")[0]
        for line in f.readlines():
            items = line.strip()
            if items.startswith("#"):
                header = items.split("\t")
            else:
                hit_ITS2 = items.split("\t")
                tax_ITS2 = hit_ITS2[-1].split("; ")
                if len(tax_ITS2) > 5:
                    ITS2_name = tax_ITS2[-1].split("__")[-1]
                    ITS2_reads = hit_ITS2[1]
                    ITS2_dict[ITS2_name]= [ITS2_reads]
                else:
                    ITS2_name = tax_ITS2[-1].split("__")[-1] +"_sp"
                    ITS2_reads = hit_ITS2[1]
                    ITS2_dict[ITS2_name]= [ITS2_reads]
        sample_ITS2_dict[sample_name]=ITS2_dict

        # num of its2 found in each sample 
        ITS2_hit_num.append(len(ITS2_dict.keys()))
        
        # num of sample with no its hits 
        if len(ITS2_dict.keys()) ==0:
            its_hit_count +=1   
          
print("mean num of its2 in each sample ",sum(ITS2_hit_num)/len(ITS2_hit_num))
print("min num of its2 in each sample", min(ITS2_hit_num))
print("max num of its2 in each sample", max(ITS2_hit_num))
print("num of sample where no its2 found", its_hit_count, "out of 965")

its_list = []
for sample, its2_info in sample_ITS2_dict.items():
    for its2, info in its2_info.items():
        if its2 not in its_list:
            its_list.append(its2)
print(its_list)

print(len(its_list))
its2_read_list = []
hit_count = 0
miss_count = 0
for sample, its2_info in sample_ITS2_dict.items():
    sample_its2_read_dict = {}
    sample_its2_read_dict["0_sample"] = sample
    for its2 in its_list:
        if its2 in its2_info.keys():
            sample_its2_read_dict[its2]=its2_info[its2][0]
        else:
            sample_its2_read_dict[its2]=0
    sorted_dict = dict(sorted(sample_its2_read_dict.items()))
    its2_read_list.append(sorted_dict)

# convert into dataframe
df = pd.DataFrame(data=its2_read_list)

#convert into excel
df.to_csv("output_path_filename", index=False)




