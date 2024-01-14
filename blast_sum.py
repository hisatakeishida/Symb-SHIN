# python script to summarise blast hits 

import pandas as pd
import os
from pathlib import Path
import shutil

# directory that contains blast hits across samples
txt_dir = "EDITHERE"

# output text path
out_txt = "EDITHERE"

os.chdir(txt_dir)
marker_hit_num = []
seq_hit_list = []
for file in os.listdir():
    with open(file) as f:
        sample_list = []
        seq_hit = []
        if str(file).endswith("_contigs"):
            pass
        else: #looking for contigs.sorted
            print(str(file))
            sample_name = file.split("_contigs")[0]
            print(sample_name)
            line = f.readline().strip('\n')
            items = line.split('	')
            sample_list.append(sample_name)
            if len(items) > 1:
                seq_hit.append(sample_name)
                seq_hit.append(items[1])
                seq_hit_list.append(seq_hit)
            sample_list = sample_list + items
            marker_hit_num.append(sample_list)

#convert into dataframe
df = pd.DataFrame(data=marker_hit_num)
print(df)

#convert into excel
df.to_csv(out_txt, index=False)




