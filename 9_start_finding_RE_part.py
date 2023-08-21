import parameters as p

import os

import misseq

misseq.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/filter_primer/',
            outputdir = os.path.abspath(p.OUTPUTDIR) + '/filter_re/',
            refway = p.REF_GENOME,
            mseq = None,
            mname = 'MISS_RE_HAMMING',
            shift = p.SHIFT_MISS_RE,
            n_core = p.NCORE)