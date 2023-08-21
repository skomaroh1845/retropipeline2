import parameters as p

import os

import misseq

misseq.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/filter_rsite/',
            outputdir = os.path.abspath(p.OUTPUTDIR) + '/filter_primer/',
            refway = p.REF_GENOME,
            mseq = p.PRIMER,
            mname = 'MISS_P_HAMMING',
            shift = p.SHIFT_MISS_PRIMER,
            n_core = p.NCORE)