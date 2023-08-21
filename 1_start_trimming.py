import parameters as p

import os
import trimmRE

trimmRE.main(inputdir = p.INPUTDIR,
             outputdir = os.path.abspath(p.OUTPUTDIR) + '/preprocessing/',
             primer = p.PRIMER,
             ad1 = p.AD1,
             ad2 = p.AD2,
             blen = p.BLEN,
             shift = p.SHIFT,
             mist = p.MIST,
             restrict_site = p.RESTRICT_SITE,
             re_part = p.RE,
             is_short_flank = p.IS_SHORT_FLANK,
             chaos = p.CHAOS,
             n_core = p.NCORE,
             trimm_n = p.TRIMM_N,
             trimm_poly_N = p.TRIMM_POLY_N,
             poly_n_r1 = p.POLY_N_R1,
             poly_n_r2 = p.POLY_N_R2,
             skip_short_reads = p.SKIP_SHORT_READS,
             mid_mist_short_reads = p.MID_MIST_SHORT_READS,
             end_mist_short_reads = p.END_MIST_SHORT_READS,
             place_of_search_tail = p.PLACE_OF_SEARCH_TAIL,
             min_seq_len_after_trimm = p.MIN_SEQ_LEN_AFTER_TRIMM)