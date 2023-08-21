import parameters as p

import os

import collapser

collapser.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/table/',
               outputdir = os.path.abspath(p.OUTPUTDIR) + '/collapse_table/',
               target_re = p.RE,
               n_core = p.NCORE)