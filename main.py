'''
The entry point for the phenotype simulator
'''
import sys
import numpy as np
import random
import string
import pandas 
import json

from itertools import starmap
from pathlib import Path
'''
int, float -> list of factors
factors{}
'''
def createFactors(n_factors, p_regulatory):

    def genID():
        letters = string.ascii_lowercase
        length = 5
        return ''.join([random.choice(letters) for i in range(length)])

    def createBaseFactor(n):
        return {'id': genID(), 
                'type': 'base'}
    
    def createRegulatoryFactor(target):
        return {'id': genID(),
                'type': 'reg',
                'target': target['id']}

    n_regulatory = int(np.ceil(n_factors * p_regulatory))
    n_base = int(n_factors - n_regulatory)

    base_factors = list(map(createBaseFactor, np.arange(n_base)))
    reg_factors = list(map(createRegulatoryFactor, [random.choice(base_factors) for x in range(n_regulatory)]))

    return base_factors, reg_factors

'''
the input file is a linear chromosome painting file. Should be really well formatted, right?
'''
def read_input(in_file):
    df = pandas.read_csv(in_file, sep='\t')
    #we're dropping 3D7 as well.
    df = df.drop(df.columns[:3], axis=1)

    df = df.transpose()    #right now we drop the chr names and positions
    return df

'''
paintings is a pandas df. cols are positions. 
'''
def annotateFactors(base_factors, regulatory_factors, paintings):

    pos_list = random.sample(list(paintings.columns), len(base_factors) + len(regulatory_factors))
    itr = iter(pos_list)
    #do base factors first
    base_eff_min = -5
    base_eff_max = 5
    #value assignment scheme can be changed
    #current scheme is for the second most common allele to be active and the rest to be inactive.
    for factor in base_factors:
        pos = next(itr)
        data = list(paintings[pos])
        bases = sorted(list(set(data)), key=lambda x: data.count(x), reverse=True)
        eff = random.randint(base_eff_min, base_eff_max)
        vals = [0]*len(bases)
        vals[1] = eff

        factor['pos'] = pos
        factor['alleles'] = list(bases)
        factor['eff_sizes'] = vals

    #then regulatory factors
    reg_eff_min = -1
    reg_eff_max = 1
    for factor in regulatory_factors:
        pos = next(itr)
        data = list(paintings[pos])
        bases = sorted(list(set(data)), key=lambda x: data.count(x), reverse=True)
        eff = random.uniform(reg_eff_min, reg_eff_max)
        vals = [0]*len(bases)
        vals[1] = eff

        factor['pos'] = pos
        factor['alleles'] = list(bases)
        factor['eff_sizes'] = vals
    
    return base_factors, regulatory_factors

def calculatePhenotype(base_factors, regulatory_factors, painting):

    #the base days to clear is 10.
    base_val = 10
    reg_stack = {}
    for reg_factor in regulatory_factors:
        allele = painting[reg_factor['pos']]
        eff = reg_factor['eff_sizes'][reg_factor['alleles'].index(allele)]
        
        #if the target isn't on the stack, init it to 1 and then apply eff
        try:
            reg_stack[reg_factor['target']] += eff 
        except KeyError:
            reg_stack[reg_factor['target']] = 1 + eff

    result = base_val
    for factor in base_factors:
        allele = painting[factor['pos']]
        eff = factor['eff_sizes'][factor['alleles'].index(allele)]

        try:
            result += eff * reg_stack[factor['id']]
        except KeyError:
            result += eff
        
    return result 

def write_factors(base_factors, regulatory_factors, outfile):
    
    with open(outfile, 'w') as output:
        output.write(json.dumps(base_factors + regulatory_factors))

def write_phenotypes(phenotypes, sample_list, outfile):
    #so we're gonna format according to the plasmo metadata
    with open(outfile, 'w') as output:
        #header
        cols = ['Sample_Name', 'clearance half-life']
        output.write('{0},{1}\n'.format(cols[0], cols[1]))


        for sample, val in zip(sample_list, phenotypes):
            output.write('{0},{1}\n'.format(sample, val))

#convert the names to the ones present in meta.csv
def convert_to_plasmo(sample_list):
       
    key_path = Path('/d/data/plasmo/filtered_runfile.txt')    
    with open(key_path) as input:
        key_df = pandas.read_csv(input, sep='\t', header=0)
    
    key = {}
    for row in key_df.iterrows():
        data = row[1]
        run = data['Run']
        name = data['Sample_Name']
        key[run] = name

    return [key[x] for x in sample_list]


if __name__ == '__main__':
    args = sys.argv
    print('start')
    #arguments
    n_factors = int(args[1])
    p_regulatory = float(args[2])
    in_file = args[3]

    in_file = Path(in_file)
    parent = in_file.parent
    factor_outfile = parent / 'factors.json'
    phenotype_outfile = parent / 'sim_meta.csv'

    #create factors
    base_factors, regulatory_factors = createFactors(n_factors, p_regulatory)

    #annotate factors
    chr_paintings = read_input(in_file)
    base_factors, regulatory_factors = annotateFactors(base_factors, regulatory_factors, chr_paintings)
    #calculate phenotypes
    phenotypes = starmap(calculatePhenotype, [(base_factors, regulatory_factors, chr_paintings.loc[x]) for x in chr_paintings.index])
    
    #we need to convert the run IDs to names
    names = convert_to_plasmo(list(chr_paintings.index))
    #some sort of output
    write_factors(base_factors, regulatory_factors, factor_outfile)
    write_phenotypes(phenotypes, names, phenotype_outfile)
    print('finish')