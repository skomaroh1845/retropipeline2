import sys

#PATHWAYS
INPUTDIR = sys.argv[1]  # '/home/common/komarov.na/data/artificial_multi'
OUTPUTDIR = sys.argv[2]  #  '/home/common/komarov.na/results/artificial_multi'

MAPPER = 'bowtie2' #bwa or bowtie2
BOWTIE_INDEX = '/home/common/komarov.na/reference/UCSC_hg38/ref_bowtie2_index'
BWA_INDEX = '/home/common/komarov.na/reference/UCSC_hg38/ref_bwa_index'

ALU_REF_LIB = '/home/common/komarov.na/reference/Alu_replibrary_hg38.txt'  # fixed and polymorphic Alu insertions
REF_GENOME = '/home/common/komarov.na/reference/UCSC_hg38/genome.fa'
REPEATS_REF_LIB = '/home/common/komarov.na/reference/repeats_hg38.tabular'  # all known fixed repeats

#SYS PARAMETRS
NCORE = 16  # threads

#PREPROCESSING
PRIMER = 'GAGCCACCGCGC'
SHIFT = 5
MIST = 1
AD1 = 'GCGTGCTGCGG'
AD2 = 'AGGGCGGT'
BLEN = 9
RE = 'CCGGCC'
RESTRICT_SITE = {'AGCT': 'CT'} #, 'TCGA': 'GA'} # {'GTAC': 'TAC', 'CTAG': 'TAG'} # {'restrict site': 'r2_start'}
IS_SHORT_FLANK = False
CHAOS = True
TRIMM_N = 0
TRIMM_POLY_N = False
POLY_N_R1 = (15, 0.8, 7) # poly_n_win_size_r1, poly_n_th_r1, poly_n_shift_r1
POLY_N_R2 = (15, 0.8, 7) # poly_n_win_size_r2, poly_n_th_r2, poly_n_shift_r2
SKIP_SHORT_READS = 50
MID_MIST_SHORT_READS = (1, 2) # (r1 ,r2) max mismatches in AD2 or PRIMER which will remove from reads
END_MIST_SHORT_READS = (5, 5) # (r1 ,r2) min length of AD2 or PRIMER at the reads ends
PLACE_OF_SEARCH_TAIL = (None, None) # (r1 ,r2)
MIN_SEQ_LEN_AFTER_TRIMM = (0, 0) # (r1 ,r2)

#SAMFILTER
FLAGS = [99, 83, 147, 163] # + [355, 339, 403, 419]

R_SITE_DISTANCE = 1000

#FILTERS
FIX_WINDOW = 1
MIN_READ = 1
MAX_DIST = 1000
MIN_DIST = 0
INSWINDOW_FIX = 15 # 40
SHIFT_RESTRICT_SITE = 3
SHIFT_MISS_PRIMER = 20 # 12
SHIFT_MISS_RE = 8

#THRESHOLDS
RE_HAMMING = 2
FLANK_ERRORS = 5
REPEAT = 0.8
MISS_PRIMER = 4
MISS_RE = 2
TEMPLATE_SWITCH_MD = 30

#METACLUSTERING
WINDOW = 100
MAIN_FLANK_LEN = 11 # template switch?
