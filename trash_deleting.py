import pandas as pd
import os
from os import listdir
from os.path import isfile, join
from datetime import datetime


def main(inputdir, outputdir,
         re_hamming,
         flank_errors,
         rsite,
         repeat,
         m_primer, primer_name,
         m_re, re_name):

    before = datetime.now()
    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    onlyfiles = [f for f in listdir(inputdir) if (isfile(join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]

    for filename in onlyfiles:
        df = pd.read_table(inputdir + filename)
        #print(filename, 'before_filters', df.shape)
        if re_hamming:
            #print(filename, 're_hamming', df[df['RE_HAMMING'] < re_hamming].shape)
            df = df[df['RE_HAMMING'] < re_hamming]
        if flank_errors:
            df['TMP'] = df['MISMATCH'].values + df['INSERTION'].values + df['DELETION'].values
            #print(filename, 'flank_errors', df[df['TMP'] < flank_errors].shape)
            df = df[df['TMP'] < flank_errors]
            df = df.drop('TMP', axis=1)
        if rsite:
            #print(filename, 'rsite', df[df[rsite].isnull()].shape)
            df = df[df[rsite].isnull()]
        if repeat:
            #print(filename, 'repeat', df[df['REPEAT_FUNC'] < repeat].shape)
            df = df[df['REPEAT_FUNC'] < repeat]
        if m_primer:
            #print(filename, 'm_primer', df[df['MISS_'+primer_name+'_HAMMING'] > m_primer].shape)
            df = df[df['MISS_'+primer_name+'_HAMMING'] > m_primer]
        if m_re:
            #print(filename, 'm_re', df[df['MISS_'+re_name+'_HAMMING'] > m_re].shape)
            df = df[df['MISS_'+re_name+'_HAMMING'] > m_re]
        if df.shape[0] == 0:
            continue
        df.to_csv(outputdir + filename, index=False, sep='\t')
