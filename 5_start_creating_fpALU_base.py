import parameters as p

import os

import exactmatch_fpALU

exactmatch_fpALU.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/collapse_table/',
                      outputdir = os.path.abspath(p.OUTPUTDIR) + '/ematch_table/',
                      replibrary = p.ALU_REF_LIB,
                      refway = p.REF_GENOME,
                      restrict_site = p.RESTRICT_SITE,
                      max_dist = p.MAX_DIST,
                      min_dist = 100,
                      min_read = p.MIN_READ,
                      inswindow = p.FIX_WINDOW,
                      direction = 5,
                      n_core = p.NCORE)