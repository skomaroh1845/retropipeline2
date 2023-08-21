import parameters as p

import os

import metacluster

metacluster.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/pre_metatable/',
                 outputdir = os.path.abspath(p.OUTPUTDIR) + '/metatable_somatic/',
                 pcdir = os.path.abspath(p.OUTPUTDIR) + '/collapse_table/',
                 target_re = p.RE,
                 window = p.WINDOW,
                 blen = p.BLEN)