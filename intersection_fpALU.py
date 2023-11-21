import pandas as pd
import intervaltree as it
import os, re
from os import listdir
from os.path import isfile, join
# from tqdm import tqdm_notebook
import numpy as np
from datetime import datetime
from joblib import Parallel, delayed




def intersection(filename, inputdir, outputdir, outputdir_fix, replib_inputdir, inswindow, fix_ins, prime):

    readsname = os.path.splitext(filename)[0].split('_humanread')[0]

    rep = pd.read_table(replib_inputdir + readsname + '_ematch.txt')
    replib_group = rep.groupby(['CHR', 'STRAND'])
    replib_tree = {}
    # for name, group in tqdm_notebook(replib_group, desc='rep: '+readsname):
    for name, group in replib_group:
        if prime == 5:
            if name[1] == '+':
                start_group = [pos-int((pos-r_site)*0.9)
                                 for r_site, pos in zip(list(group['R_SITE_POS']), list(group['START']))]
                end_group = [pos+inswindow+1 for pos in list(group['START'])]
            else:
                end_group = [pos+int((r_site-pos)*0.9)+1
                                 for r_site, pos in zip(list(group['R_SITE_POS']), list(group['END']))]
                start_group = [pos-inswindow for pos in list(group['END'])]
        elif prime == 3:
            if name[1] == '+':
                end_group = [pos+int((r_site-pos)*0.9)+1
                                for r_site, pos in zip(list(group['R_SITE_POS']), list(group['END']))]
                start_group = [pos-inswindow for pos in list(group['END'])]
            else:
                start_group = [pos-int((pos-r_site)*0.9)
                                for r_site, pos in zip(list(group['R_SITE_POS']), list(group['START']))]
                end_group = [pos+inswindow+1 for pos in list(group['START'])]
        replib_tree[name[0] + name[1]] = it.IntervalTree(it.Interval(start, end, ins_name)
         for start, end, ins_name in zip(start_group, end_group, list(group['NAME'])))

    df = pd.read_table(inputdir + filename)
    if df.shape[0] == 0:
        pass
    else:
        df_out = open(outputdir + readsname + '.txt', 'w')
        df_fix = open(outputdir_fix + readsname + '.txt', 'w')
        df_out.write('\t'.join(['CLUSTER_ID',
                                    'READNAME',
                                    'CHR',
                                    'POS',
                                    'INS_STRAND',
                                    'PRIMER',
                                    'RE',
                                    'RE_AMOUNT',
                                    'RE_HAMMING',
                                    'R1',
                                    'R2',
                                    'TLEN',
                                    'CIGAR_R1',
                                    'MDFLAG_R1',
                                    'MD_SUM',
                                    'MISMATCH',
                                    'INSERTION',
                                    'DELETION',
                                    'NUM_READS',
                                    'NUM_BC',
                                    'NUM_TLEN']) + '\n')
        df_fix.write('\t'.join(['CLUSTER_ID',
                                    'READNAME',
                                    'CHR',
                                    'POS',
                                    'INS_STRAND',
                                    'PRIMER',
                                    'RE',
                                    'RE_AMOUNT',
                                    'RE_HAMMING',
                                    'R1',
                                    'R2',
                                    'TLEN',
                                    'CIGAR_R1',
                                    'MDFLAG_R1',
                                    'MD_SUM',
                                    'MISMATCH',
                                    'INSERTION',
                                    'DELETION',
                                    'NUM_READS',
                                    'NUM_BC',
                                    'NUM_TLEN',
                                    'NAME']) + '\n')
        
        # for i in tqdm_notebook(range(np.shape(df)[0]), desc=readsname):
        for i in range(np.shape(df)[0]):
            row = df.iloc[i, ]
            if row['CHR'] + row['INS_STRAND'] in replib_tree:
                find_iter = replib_tree[row['CHR'] + row['INS_STRAND']][row['POS']]
                if len(find_iter) == 0:
                    pd.DataFrame(row).T.to_csv(df_out, sep='\t', mode='a', header=False, index=False)
                else:
                    ins_names = [ins.data for ins in find_iter]
                    if fix_ins is None:
                        fix_row = pd.concat([row, pd.Series([', '.join(ins_names)], index=['NAME'])])
                        pd.DataFrame(fix_row).T.to_csv(df_fix, sep='\t', mode='a', header=False, index=False)
                    else:
                        if set(fix_ins).intersection(set(ins_names)):
                            fix_row = pd.concat([row, pd.Series([', '.join(ins_names)], index=['NAME'])])
                            pd.DataFrame(fix_row).T.to_csv(df_fix, sep='\t', mode='a', header=False, index=False)
            else:
                pd.DataFrame(row).T.to_csv(df_out, sep='\t', mode='a', header=False, index=False)
        df_out.close()
        df_fix.close()


def main(inputdir, outputdir, outputdir_fix, replib_inputdir, inswindow, fix_ins, n_core, prime):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'
    outputdir_fix = os.path.abspath(outputdir_fix) + '/'
    replib_inputdir = os.path.abspath(replib_inputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if not os.path.exists(outputdir_fix):
        os.makedirs(outputdir_fix)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]
    onlyfiles = [f for f in onlyfiles if re.search('humanread', f)]

    if len(onlyfiles) == 1:
        filename = onlyfiles[0]
        stat_series = intersection(filename,
                              inputdir, outputdir, outputdir_fix, replib_inputdir, inswindow, fix_ins, prime)
        #stat_df = stat_series.to_frame().transpose()
    else:
        stat_series = Parallel(n_jobs=n_core, prefer="threads")(delayed(intersection)(filename,
                                            inputdir, outputdir, outputdir_fix, replib_inputdir, inswindow, fix_ins, prime)
                                                for filename in onlyfiles)
        #stat_df = pd.concat(stat_series, axis=1).transpose()