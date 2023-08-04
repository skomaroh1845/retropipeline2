import os
import re
import pandas as pd
import numpy as np
import intervaltree as it
import tqdm

from intersection_re import RsiteRetroDB, AbstractIntersectionClass


class IntersectionLINE(AbstractIntersectionClass):
    
    def _create_tree(self, strand, group):
        
        if strand == '+':
            start = [x - 10 for x in group['RSITE_POS']]
            end = [x + 10 + 1 for x in group['RSITE_POS']]
        else:
            start = [x - 1s0 for x in group['RSITE_POS']]
            end = [x + 10 + 1 for x in group['RSITE_POS']]
        
        data = [{'name': x} for x in group['NAME']]
        
        return it.IntervalTree(it.Interval(s, e, d) for s, e, d in zip(start, end, data))



def intersection(filename, inputdir, somatic_dir, fix_dir, intersector, min_read):
    
    df = pd.read_csv(os.path.join(inputdir, filename), '\t')
    readsname = os.path.splitext(filename)[0].split('_humanread')[0]
    
    somatic_file = open(os.path.join(somatic_dir, readsname) + '.txt', 'w')
    fix_file = open(os.path.join(fix_dir, readsname) + '.txt', 'w')
    
    fix_columns = list(df.columns) + ['NAME']
    somatic_file.write('\t'.join(list(df.columns)) + '\n')
    fix_file.write('\t'.join(fix_columns) + '\n')
    
    
    for i in tqdm.tqdm_notebook(range(df.shape[0]), desc=readsname):
        row = df.iloc[i, ]
        if row['NUM_READS'] >= min_read:
            try:
                matches = intersector.intersect(row['CHR'], row['INS_STRAND'], row['POS'])
                if len(matches) == 0:
                    pd.DataFrame(row).T.to_csv(somatic_file, sep='\t', mode='a', header=None, index=None)
                else:
                    new_row = row.copy()
                    new_row.set_value('NAME', matches[0].data['name'])
                    pd.DataFrame(new_row).T.to_csv(fix_file, sep='\t', mode='a', header=None, index=None)
            except:
                pd.DataFrame(row).T.to_csv(somatic_file, sep='\t', mode='a', header=None, index=None)
        else:
            pd.DataFrame(row).T.to_csv(somatic_file, sep='\t', mode='a', header=None, index=None)
    somatic_file.close()
    fix_file.close()



def main(inputdir, outputdir, outputdir_fix, replibrary, refway, restrict_site, max_dist, min_dist, min_read, inswindow):

    inputdir = os.path.abspath(inputdir) + '/'
    outputdir = os.path.abspath(outputdir) + '/'
    outputdir_fix = os.path.abspath(outputdir_fix) + '/'

    db_creator = RsiteRetroDB(refway, restrict_site=restrict_site)
    line_db = db_creator.create_db(pd.read_csv(replibrary, sep='\t'), prime=3)
    dist = np.minimum(np.abs(line_db['RSITE_POS'].values - line_db['START'].values),
                      np.abs(line_db['RSITE_POS'].values - line_db['END'].values))
    line_db = line_db[(dist >= min_dist) & (dist <= max_dist)]

    intersector = IntersectionLINE(line_db)
    intersector.create_tree_dict(line_db)

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    if not os.path.exists(outputdir_fix):
        os.makedirs(outputdir_fix)

    onlyfiles = [f for f in os.listdir(inputdir) if (os.path.isfile(os.path.join(inputdir, f))
                                                     and os.path.splitext(f)[1] == '.txt')]
    onlyfiles = [f for f in onlyfiles if re.search('humanread', f)]


    for filename in onlyfiles:
        intersection(filename, inputdir, outputdir, outputdir_fix, intersector, min_read)
