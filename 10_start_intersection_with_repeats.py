import parameters as p

import os

import intersection_replib

intersection_replib.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/filter_re/',
                         outputdir = os.path.abspath(p.OUTPUTDIR) + '/filter_replib/',
                         repeatway = p.REPEATS_REF_LIB,
                         n_core = p.NCORE)