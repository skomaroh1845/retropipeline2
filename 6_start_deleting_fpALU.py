import parameters as p

import os

import intersection_fpALU

intersection_fpALU.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/collapse_table/',
                        outputdir = os.path.abspath(p.OUTPUTDIR) + '/notfpALU_table/',
                        outputdir_fix = os.path.abspath(p.OUTPUTDIR) + '/fix_ins_table/',
                        replib_inputdir = os.path.abspath(p.OUTPUTDIR) + '/ematch_table/',
                        inswindow = p.INSWINDOW_FIX,
                        fix_ins = None,
                        n_core = p.NCORE,
                        prime = p.WHICH_TAIL)