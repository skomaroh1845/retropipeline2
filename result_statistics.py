import parameters as p
import pandas as pd
import os
from pathlib import Path
import gzip

# 0 step

directory = os.path.abspath(p.INPUTDIR) + '/'
files = Path(directory).glob('*.gz')
cols = []
i = 0
for file in files:
    if i % 2 == 0:
        f_name = file.name
        cols.append(file.name.split('_')[0] + '_' + file.name.split('_')[1])
    i += 1

df = pd.DataFrame(columns=cols)

df.loc['before (num reads)'] = pd.Series(index=df.columns).fillna(0)

i = 0
files = Path(directory).glob('*.gz')
for file in files:
    if i % 2 == 0:
        with gzip.open(file, 'r') as f:
            df.loc['before (num reads)'][i] = sum(1 for _ in f)
    i += 1


# 1 step
directory = os.path.abspath(p.OUTPUTDIR) + '/preprocessing/'
files = Path(directory).glob('*.txt')
i = 0
df.loc['trimm'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        df.loc['trimm'][i] = sum(1 for _ in f)
        i += 1

# 2 step
directory = os.path.abspath(p.OUTPUTDIR) + '/mapping/'
files = Path(directory).glob('*.sam')
i = 0
df.loc['mapping'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        df.loc['mapping'][i] = int(sum(1 if _[0] != '@' else 0 for _ in f) / 2)
        i += 1

# 3 step
directory = os.path.abspath(p.OUTPUTDIR) + '/table/'
files = Path(directory).glob('*.txt')
i = 0
df.loc['filter_sam_files'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        df.loc['filter_sam_files'][i] = sum(1 for _ in f) - 1
        i += 1

# 4 step
directory = os.path.abspath(p.OUTPUTDIR) + '/collapse_table/'
files = Path(directory).glob('*humanread.txt')
i = 0
df.loc['grouping'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        df.loc['grouping'][i] = sum(1 for _ in f) - 1
        i += 1

# 6 step
directory = os.path.abspath(p.OUTPUTDIR) + '/notfpALU_table/'
files = Path(directory).glob('*.txt')
i = 0
df.loc['no_fpALU'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        df.loc['no_fpALU'][i] = sum(1 for _ in f) - 1
        i += 1

# 7 step
directory = os.path.abspath(p.OUTPUTDIR) + '/filter_rsite/'
files = Path(directory).glob('*.txt')
i = 0
df.loc['no_rsite'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        df.loc['no_rsite'][i] = sum(1 for _ in f) - 1
        #print(file, df.loc['step7'][i])
        i += 1

# 8 step
directory = os.path.abspath(p.OUTPUTDIR) + '/filter_primer/'
files = Path(directory).glob('*.txt')
i = 0
df.loc['no_flank_primer'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        df.loc['no_flank_primer'][i] = sum(1 for _ in f) - 1
        #print(file, df.loc['step8'][i])
        i += 1

# 9 step
directory = os.path.abspath(p.OUTPUTDIR) + '/filter_re/'
files = Path(directory).glob('*.txt')
i = 0
df.loc['no_flank_re'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        df.loc['no_flank_re'][i] = sum(1 for _ in f) - 1
        i += 1

# 10 step
directory = os.path.abspath(p.OUTPUTDIR) + '/pre_metatable/'
files = Path(directory).glob('*.txt')
i = 0
df.loc['no_f_repeats'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        df.loc['no_f_repeats'][i] = sum(1 for _ in f) - 1
        i += 1


# 11 step
directory = os.path.abspath(p.OUTPUTDIR) + '/metatable_somatic/'
files = Path(directory).glob('*humanread.txt')

df.loc['metaclustering'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        line = f.readline()
        while line != '':
            line = f.readline()
            if line == '':
                break
            line = line.split('\t')[1]
            line = line.split('_')[0] + '_' + line.split('_')[1]
            df.loc['metaclustering'][line] += 1


# 12 step
directory = os.path.abspath(p.OUTPUTDIR) + '/result_somatic/'
files = Path(directory).glob('*_ts.txt')

df.loc['template_switch'] = pd.Series(index=df.columns).fillna(0)
for file in files:
    with open(file, 'r') as f:
        line = f.readline()
        while line != '':
            line = f.readline()
            if line == '':
                break
            line = line.split('\t')[1]
            line = line.split('_')[0] + '_' + line.split('_')[1]
            df.loc['template_switch'][line] += 1

df.to_csv('general_stats.csv')
