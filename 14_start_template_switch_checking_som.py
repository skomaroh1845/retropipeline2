import parameters as p

import os

import template_switch

template_switch.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/metatable_somatic/',
                     outputdir = os.path.abspath(p.OUTPUTDIR) + '/result_somatic/',
                     primer = p.PRIMER,
                     main_flank_len = p.MAIN_FLANK_LEN,
                     template_switch_md = p.TEMPLATE_SWITCH_MD,
                     refway = p.BWA_INDEX,
                     bwaway = 'bwa')