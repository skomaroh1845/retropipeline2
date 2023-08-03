import os
import re
import sys
import pandas as pd
import numpy as np
import intervaltree as it
import tqdm
import pysam


def reverse_complement(seq):
    bases = str.maketrans('ACGTRYKMBDHVacgtrykmbdhv', 'TGCAYRMKVHDBtgcayrmkvhdb')
    return seq[::-1].translate(bases)


class RsiteHunter:
    
    
    def __init__(self, refway, restrict_site):
        self.ref = pysam.FastaFile(refway)
        self.rsite = restrict_site
        self.rsite_recompile = re.compile('|'.join(restrict_site.keys()))
    
    
    def find_rsite_in_seq(self, seq, keep='first'):
        if keep == 'first':
            idx = 0
        elif keep == 'last':
            idx = -1
        else:
            raise ValueError("Invalid value for keep. Expecting 'first' or 'last'")
        
        if seq is not None:
            rsite_found = [x for x in self.rsite_recompile.finditer(seq)]
            if rsite_found:
                s = rsite_found[0].start()
                e = rsite_found[0].end()
                rsite_name, closest_pos = seq[s:e], s
                return rsite_name, closest_pos

        return None, None
    
        
    def find_rsite(self, chrom, pos, side):
        if side not in ['left', 'right']:
            raise ValueError("Invalid value for keep. Expecting 'left' or 'right'")

        is_found_rsite = False
        batch_100_nt = 0
        max_len_rsite = max([len(x) for x in self.rsite.keys()])
        while not is_found_rsite:
            batch_100_nt += 1
            try:
                if side == 'left': 
                    start = pos - 1 - batch_100_nt * 100
                    end = pos - 1 - (batch_100_nt - 1) * 100 + max_len_rsite
                    seq = reverse_complement(self.ref.fetch(chrom, start, end).upper())
                elif side == 'right':
                    start = pos + (batch_100_nt - 1) * 100
                    end = pos + batch_100_nt * 100 + max_len_rsite
                    seq = self.ref.fetch(chrom, start, end).upper()
            except ValueError:
                seq = None
                break
            rsite_name, rsite_pos = self.find_rsite_in_seq(seq, keep='first')
            if rsite_name is not None:
                is_found_rsite = True
            
        if is_found_rsite:
            rsite_fragment = len(self.rsite[rsite_name]) - 1
            if side == 'left':
                rsite_pos = start + (len(seq) - rsite_pos) - rsite_fragment
            elif side == 'right':
                rsite_pos = start + rsite_pos + 1 + rsite_fragment

        return rsite_name, rsite_pos


class RsiteRetroDB(RsiteHunter):

    
    def create_db(self, db, prime):
        if prime not in [3, 5]:
            raise ValueError("Invalid value for prime. Expecting 3 or 5")
        
        db_corrected = db.copy()
        # db_corrected['IDX'] = range(db_corrected.shape[0])
        db_corrected = db_corrected[db_corrected['CHR'].isin([f'chr{i}' for i in range(1, 22+1)] + ['chrX', 'chrY'])]

        new_row_holder = []
        for i, row in tqdm.tqdm_notebook(db_corrected.iterrows(), total=db_corrected.shape[0], desc='create database'):
            new_row = [row['CHR'], row['START'], row['END'], row['STRAND'], row['NAME']]

            name_left, pos_left = self.find_rsite(row['CHR'], row['START'], side='left')
            name_right, pos_right = self.find_rsite(row['CHR'], row['END'], side='right')

            new_row.extend([name_left, pos_left, name_right, pos_right])
            new_row_holder.append(new_row)
        
        colnames = list(db_corrected.columns)
        colnames_rsite_temp = ['RSITE_NAME_LEFT', 'RSITE_POS_LEFT', 'RSITE_NAME_RIGHT', 'RSITE_POS_RIGHT']
        colnames.extend(colnames_rsite_temp)
                    
        
        db_corrected = pd.DataFrame(new_row_holder, columns=colnames)
        prime_rule = {
            5: {'+': 'LEFT', '-': 'RIGHT'},
            3: {'+': 'RIGHT', '-': 'LEFT'}
        }
        rsite_name = np.where(db_corrected['STRAND'] == '+',
                              db_corrected['RSITE_NAME_' + prime_rule[prime]['+']],
                              db_corrected['RSITE_NAME_' + prime_rule[prime]['-']])
        rsite_pos = np.where(db_corrected['STRAND'] == '+',
                             db_corrected['RSITE_POS_' + prime_rule[prime]['+']],
                             db_corrected['RSITE_POS_' + prime_rule[prime]['-']])
        
        db_corrected['RSITE_NAME'] = rsite_name
        db_corrected['RSITE_POS'] = rsite_pos
        db_corrected = db_corrected.drop(colnames_rsite_temp, axis=1)
        db_corrected = db_corrected.dropna()

        return db_corrected


class AbstractIntersectionClass:
    
    
    def __init__(self, data):
        self.data = data
        self.tree_dict = {}
    
    
    def create_tree_dict(self, data=None, progress=True):
        if data is None:
            data = self.data
        
        data_groupby = data.groupby(['CHR', 'STRAND'])
        if progress:
            pbar = tqdm.tqdm_notebook(total=len(data_groupby))

        for (chrom, strand), group in data_groupby:
            self.tree_dict[(chrom, strand)] = self._create_tree(strand, group)
            if progress:
                pbar.update(1)
        
        if progress:
            pbar.close()
        
        return self.tree_dict
    
    
    def _create_tree(self, strand, group):
        raise NotImplementedError("Method '_create_tree' not implemented")
        
        
    def _modify_interval(self, interval):
        if isinstance(interval, int):
            x1, x2 = interval, interval + 1
        else:
            x1, x2 = sorted(interval)
        return x1, x2
    
    
    def intersect(self, chrom, strand, start, end=None):
        """Intersect point or interval with interval tree
        
        Parameters
        ---------
        chrom : str
            chromosome name
        
        strand : str
            insertion strand ('+' or '-')
        
        start : int
            position on chromosome (inclusive)
        
        end : int, optional
            end position on chromosome (exclusive) for interval
        
        Returns
        -------
        matches : list
            list of intersections (it.Interval)
        """
        if end is None:
            end = start + 1
        matches = list(self.tree_dict[(chrom, strand)][start : end])
        
        return matches
        