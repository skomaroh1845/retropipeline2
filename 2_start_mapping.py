import parameters as p

import os

import bwamem

if p.MAPPER == 'bowtie2':
    mapper_execline = f'bowtie2 -p {p.NCORE} -I 25 -X 1000 --dovetail'
    refway = p.BOWTIE_INDEX
else:
    if p.MAPPER != 'bwa':
        print('BWA')
    mapper_execline = f'bwa mem -t {p.NCORE}'
    refway = p.BWA_INDEX

bwamem.main(inputdir = os.path.abspath(p.OUTPUTDIR) + '/preprocessing/',
            outputdir = os.path.abspath(p.OUTPUTDIR) + '/mapping/',
            refway = refway,
            bwaline = mapper_execline)