#! /usr/bin/env python3

desc = '''Identification of compositional tree-heterogeneous sequences using a Chi square test (Foster, 2004) combined with different p-value correction methods.
Optionally, it removes the composition-heterogeneous sequences from the alignment file.

Usage:
$ python chiTest.py alignment.fasta -t tree.phy -m model.dat -c bonferroni -r

João Brazão version 1
'''

import os
import sys
import subprocess
import argparse
import textwrap
from p4 import *
import statsmodels.stats.multitest as smt

var.doCheckForDuplicateSequences = False

#################################################

def matrixLowerTriangleToUpperTriangle(file, dim = 20):
    # Convert a PAML substitution rate matrix to P4 format
    with open (file,'r') as f:
        matrix_values = [float(n) for n in f.read().split()]
    if len(matrix_values) == 210:
        compositions = (matrix_values[191:])
        rates = matrix_values[0:190]
    elif len(matrix_values) == 190:
        rates = matrix_values
        print("Compositions will be optimize,")
    else:
        sys.exit(f'{len(matrix_values)}: Exit. Your matrix has not the correct number of parameters')

    assert isinstance(dim, int)
    expectedLen = int(((dim * dim) - dim)/2)
    assert len(rates) == expectedLen, f"Lower triangle len is {len(rates)}, expected {expectedLen}"

    bigM = []
    for rNum in range(dim):
        bigM.append(['0.0'] * dim)

    i = 0
    for rNum in range(1,dim):
        for cNum in range(rNum):
            bigM[cNum][rNum] = rates[i]
            i += 1
    upper = []
    for rNum in range(dim-1):
        for cNum in range(rNum+1, dim):
            upper.append(bigM[rNum][cNum])
    return upper

def do_test(alignment, tree, rates):

    # Read in the tree and give it a name
    t = func.readAndPop(tree)
    # Read in the data and give it a name
    read(alignment)
    d = Data()

    # Attach the data and a model to the tree.
    t.data = d
    t.newComp(free=1, spec='empirical')
    t.newRMatrix(free=0, spec='specified', val=rates)
    t.setNGammaCat(nGammaCat=4)
    t.newGdasrv(free=1, val=0.5)
    t.setPInvar(free=0, val=0.0)
    t.optLogLike()

    sys.stdout = open(f'{alignment}_compoTestP4.txt','wt')
    t.compoTestUsingSimulations(nSims=100, doIndividualSequences=1, doChiSquare=0, verbose=1)

def correct_pvalue(correction_pvalue,alignment):

    with open(f'{(alignment.fName).split()[0]}_compoTestP4.txt','r') as f:
        taxa_and_pvalues = [n.split() for n in f.readlines()[8:] if n != '\n' if not n.startswith('\tP-value')]

    if ['=================================================='] in taxa_and_pvalues: #why I put this here?
        taxa_and_pvalues = taxa_and_pvalues[0:taxa_and_pvalues.index(['=================================================='])]

    assert len(alignment.taxNames) == len(taxa_and_pvalues), "Error the number of lines in the *compoTestP4 do not match the taxa number. Check the file."

    if correction_pvalue:
        if correction_pvalue in ['bonferroni','sidak','holm','holm-sidak','fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']:
            pvalues = [float(n[1]) for n in taxa_and_pvalues]
            rejected, pvalues_corrected, alphacSidak, alphacBonf = smt.multipletests(pvalues, method=correction_pvalue)
        else:
            try:
                my_pvalue = float(correction_pvalue)
                rejected = [(n < my_pvalue) for n in [float(n[1]) for n in taxa_and_pvalues]]
            except ValueError:
                sys.exit('Method do not exit and is neither a float. Choose a number or from bonferroni, sidak, holm,holm-sidak,fdr_bh, fdr_by, fdr_tsbh, or fdr_tsbky.')

    else:
        ## P-values not corrected. It will be used a alfa of 0.05
        rejected = [(n < 0.05) for n in [float(n[1]) for n in taxa_and_pvalues]]

    print(f'Method to correct pvalue: {correction_pvalue}\nOffensive taxa: {len([n for n in rejected if n == True])}')

    offensive_taxa = []
    for assumption_rejected,taxon in zip(rejected, [n[0] for n in taxa_and_pvalues]):
        if assumption_rejected == True:
            offensive_taxa.append(taxon)

    return (offensive_taxa)

#######################################################################################
def main(datafile, tree, matrix, correction_pvalue, remove_taxa, overwrite):

    if not datafile:
        print("Error: No datafile specified.")
        sys.exit()
    else:
        if not os.path.exists(datafile):
            print("Error: cannot find datafile %s" % datafile)
            sys.exit()

    align = func.readAndPop(datafile)
    if 1:
        rates_matrix = matrixLowerTriangleToUpperTriangle(matrix, dim = 20)
        do_test(datafile, tree, rates_matrix)

    if correction_pvalue:
        offensive_taxa = correct_pvalue(correction_pvalue,align)

    #Write file report
    if len(offensive_taxa) != 0:
        with open(f"{datafile}_{correction_pvalue}_offensive_taxa.txt",'w') as f:
            for each in offensive_taxa:
                f.write(each+'\n')

    with open('resume.txt','a') as f:
        f.write(f"{datafile};{correction_pvalue};{len(offensive_taxa)}\n")

    if remove_taxa:

        for taxon in offensive_taxa:
            align.sequences.pop(align.taxNames.index(taxon))

        if overwrite:
            print('Alignment file will be overwrite in phylip format')
            align.writePhylip(datafile)
        else:
            align.writeFasta(f'{datafile.split(".")[0]}_TH_homogeneous.fasta')

######################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(desc),
            )
    parser.add_argument("datafile",
                        help="Data matrix filename.",
                        )
    parser.add_argument("-t", "--tree",
                        help="Tree file for compositional test.",
                        )
    parser.add_argument("-m", "--matrix",
                        help="Rate matrix for the compositional test.",
                        )
    parser.add_argument("-c", "--correction_pvalue",
                        help="Choose a pvalue threshold or a method to correct the pvalue. Choose from bonferroni, sidak, holm, holm-sidak, fdr_bh, fdr_by, fdr_tsbh, or fdr_tsbky. Default = 0.05"
                        )
    parser.add_argument("-r", "--remove_taxa",
                        help="Activate removal of TH taxon sequences.",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-i", "--overwrite",
                        help="Overwrite the current alignment file. Default=False",
                        default=False,
                        action='store_true'
                        )

    args = parser.parse_args()
    main(args.datafile, args.tree, args.matrix, args.correction_pvalue, args.remove_taxa, args.overwrite)
