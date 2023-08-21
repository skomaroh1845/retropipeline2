import parameters as p

import os

import samfilter_pysam

samfilter_pysam.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/mapping/',
                     outputdir = os.path.abspath(p.OUTPUTDIR) + '/table/',
                     flags = p.FLAGS,
                     n_core = p.NCORE)