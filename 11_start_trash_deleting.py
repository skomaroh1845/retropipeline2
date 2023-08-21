import parameters as p

import os

import trash_deleting

trash_deleting.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/filter_replib/',
                    outputdir = os.path.abspath(p.OUTPUTDIR) + '/pre_metatable/',
                    re_hamming = p.RE_HAMMING,
                    flank_errors = p.FLANK_ERRORS,
                    rsite = 'R_SITE',
                    repeat = p.REPEAT,
                    m_primer = p.MISS_PRIMER,
                    primer_name = 'P',
                    m_re = p.MISS_RE,
                    re_name = 'RE')