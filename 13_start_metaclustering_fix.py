import parameters as p

import os

import metacluster

metacluster.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/fix_ins_table/',
                 outputdir = os.path.abspath(p.OUTPUTDIR) + '/metatable_fix/',
                 pcdir = os.path.abspath(p.OUTPUTDIR) + '/collapse_table/',
                 target_re = p.RE,
                 window = p.WINDOW,
                 blen = p.BLEN)