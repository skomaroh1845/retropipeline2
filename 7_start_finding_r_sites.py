import parameters as p

import os

import misseq

misseq.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/notfpALU_table/',
            outputdir = os.path.abspath(p.OUTPUTDIR) + '/filter_rsite/',
            refway = p.REF_GENOME,
            mseq = p.RESTRICT_SITE,
            mname = "R_SITE",
            shift = p.SHIFT_RESTRICT_SITE,
            n_core = p.NCORE)