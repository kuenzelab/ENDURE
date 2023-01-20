import pandas as pd
import shutil

for i in range(1, 3):
    data = pd.read_fwf(f'data/scores/t6d{i}.txt', skiprows=1)
    lowest = data.at[data['total_score'].idxmin(), 'description']
    shutil.copy(f'data/results/{lowest}.pdb', f'../PROSS/data/pdb/adjusted/{lowest[0:5]}.pdb')
