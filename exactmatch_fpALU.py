import pandas as pd
import intervaltree as it
import sys, os, re
from os import listdir
from os.path import isfile, join
from tqdm import tqdm_notebook
import numpy as np
from datetime import datetime
from joblib import Parallel, delayed
import pysam


def create_db(x, restrict_site, max_dist, min_dist, ref):
    x.columns = ['CHR', 'START', 'END', 'STRAND', 'NAME']
    x['IDX'] = range(x.shape[0])
    chromosome_name = ['chr{}'.format(i) for i in range(1, 22+1)] + ['chrX', 'chrY']
    x = x[x['CHR'].isin(chromosome_name)]
    r_site_control = []
    r_site_recompile = list(map(re.compile, restrict_site.keys()))
    for i, row in tqdm_notebook(x.iterrows(), total=x.shape[0], desc='create database'):
        row_info = [row['CHR'], row['START'], row['END'], row['STRAND'], row['NAME'], row['IDX']]
        pos_left = int(row['START']) - max_dist
        try:
            flank_left = ref.fetch(row['CHR'], pos_left-1, int(row['START'])-min_dist-1).upper()
        except ValueError:
            flank_left = -1
        if flank_left != -1:
            fiter_left_dict = {}
            for r in r_site_recompile:
                fiter = list(r.finditer(flank_left))
                if fiter:
                    fiter_left_dict[fiter[0].group(0)] = fiter
            if fiter_left_dict:
                fiter_pos = {}
                for r_site, pos in fiter_left_dict.items(): 
                    fiter_pos[r_site] = pos[-1].start()
                clos_pos = max(fiter_pos, key=fiter_pos.get)
                r_site_left_pos = int(row['START']) - min_dist - (len(flank_left) - fiter_pos[clos_pos])
                row_info.append(clos_pos)
                row_info.append(r_site_left_pos)
            else:
                row_info.extend([np.nan, np.nan])
        pos_right = int(row['END']) + max_dist
        try:
            flank_right = ref.fetch(row['CHR'], int(row['END'])+min_dist, pos_right).upper()
        except ValueError:
            flank_right = -1
        if flank_right != -1:
            fiter_right_dict = {}
            for r in r_site_recompile:
                fiter = list(r.finditer(flank_right))
                if fiter:
                    fiter_right_dict[fiter[0].group(0)] = fiter
            if fiter_right_dict:
                fiter_pos = {}
                for r_site, pos in fiter_right_dict.items(): 
                    fiter_pos[r_site] = pos[0].start()
                clos_pos = min(fiter_pos, key=fiter_pos.get)
                r_site_right_pos = int(row['END']) + min_dist + fiter_pos[clos_pos] + 1
                row_info.append(clos_pos)
                row_info.append(r_site_right_pos)
            else:
                row_info.extend([np.nan, np.nan])
        r_site_control.append(row_info)
    colnames = ['CHR', 'START', 'END', 'STRAND', 'NAME', 'IDX', 'R_SITE_LEFT', 'R_SITE_POS_LEFT', 'R_SITE_RIGHT', 'R_SITE_POS_RIGHT']
    r_site_control = pd.DataFrame(r_site_control, columns=colnames)
    return r_site_control

def correct_prime(x, prime):
    r_site_pos = 'R_SITE_POS'
    r_site = 'R_SITE'
    if prime == 5:
        correct_r_site = np.where(x['STRAND'] == '+', x['{}_LEFT'.format(r_site)], x['{}_RIGHT'.format(r_site)])
        correct_pos = np.where(x['STRAND'] == '+', x['{}_LEFT'.format(r_site_pos)], x['{}_RIGHT'.format(r_site_pos)])
    elif prime == 3:
        correct_r_site = np.where(x['STRAND'] == '+', x['{}_RIGHT'.format(r_site)], x['{}_LEFT'.format(r_site)])
        correct_pos = np.where(x['STRAND'] == '+', x['{}_RIGHT'.format(r_site_pos)], x['{}_LEFT'.format(r_site_pos)])
    else:
        sys.exit("prime must by 5 or 3")
    x[r_site] = correct_r_site
    x[r_site_pos] = correct_pos
    x = x.drop(['{}_LEFT'.format(r_site), '{}_LEFT'.format(r_site_pos), '{}_RIGHT'.format(r_site), '{}_RIGHT'.format(r_site_pos)], axis=1)
    x = x.dropna()
    x[r_site_pos] = x[r_site_pos].astype(int)
    x = x[['CHR', 'START', 'END', 'STRAND', 'NAME', r_site, r_site_pos, 'IDX']]
    return x

def intersection(filename, inputdir, outputdir, replibrary, min_read, window):
    rep = replibrary.copy()
    replib_group = rep.groupby(['CHR', 'STRAND'])
    replib_tree = {}
    for name, group in tqdm_notebook(replib_group, desc='replib'):
        if name[1] == '+':
            start_group = [pos-window for pos in list(group['START'])]
            end_group = [pos+1+window for pos in list(group['START'])]
        else:
            end_group = [pos+2+window for pos in list(group['END'])]
            start_group = [pos+1-window for pos in list(group['END'])]
        group_list = map(list, group.values)
        group_list = ['\t'.join([x[0], str(x[1]), str(x[2]), x[3], x[4], x[5], str(x[6])]) for x in group_list]
        group_list = [(x,y) for x,y in zip(group_list, list(group.index))]
        replib_tree[name[0] + name[1]] = it.IntervalTree(it.Interval(start, end, data)
         for start, end, data in zip(start_group, end_group, group_list))

    readsname = os.path.splitext(filename)[0].split('_humanread')[0] 

    replib_out = open(outputdir + readsname + '_ematch.txt', 'w')
    replib_out.write('\t'.join(['CHR',
                                'START',
                                'END',
                                'STRAND',
                                'NAME',
                                'R_SITE',
                                'R_SITE_POS',
                                'READNAME',
                                'IDX']) + '\n')
    df = pd.read_table(inputdir + filename)
    for i in tqdm_notebook(range(np.shape(df)[0]), desc=readsname):
        row = df.iloc[i, ]
        if row['NUM_READS'] >= min_read:
            if row['CHR'] + row['INS_STRAND'] in replib_tree:
                finter = replib_tree[row['CHR'] + row['INS_STRAND']][row['POS']]
                if len(finter) > 0:
                    for x in finter:
                        replib_out.write(x.data[0] + '\t' + row['READNAME'] + '\t' + str(x.data[1]) + '\n')
                        #replib_tree[row['CHR'] + row['INS_STRAND']].removei(x.begin, x.end, x.data)
    replib_out.close()


def main(inputdir, outputdir, replibrary, refway, restrict_site, max_dist, min_dist, min_read, inswindow, direction, n_core):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    ref = pysam.Fastafile(refway)
    replibrary_data = correct_prime(create_db(pd.read_table(replibrary), restrict_site, max_dist, min_dist, ref=ref), direction)
    print(replibrary_data.shape)
    print(replibrary_data.head())

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]
    onlyfiles = [f for f in onlyfiles if re.search('humanread', f)]

    if len(onlyfiles) == 1:
        filename = onlyfiles[0]
        stat_series = intersection(filename,
                              inputdir, outputdir, replibrary_data, min_read, inswindow)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core, prefer="threads")(delayed(intersection)(filename,
                                            inputdir, outputdir, replibrary_data, min_read, inswindow)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()
