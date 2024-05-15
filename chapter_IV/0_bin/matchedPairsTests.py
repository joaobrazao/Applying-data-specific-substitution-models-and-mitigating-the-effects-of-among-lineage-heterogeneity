#! /usr/bin/env python3

desc = """
Identification of tree-heterogeneous sequences using the matched-pairs tests combined with different p-value correction methods.
Optionally, it removes identified heterogeneous sequences from the alignment file.

Usage:
$ python matchedPairsTests.py alignment.fasta -c bonferroni -r -o results.csv

João Brazão version 1
"""

import os
import sys
import subprocess
import argparse
import textwrap
from p4 import *
import numpy as np
from Bio import SeqIO
import statsmodels.stats.multitest as smt

var.doCheckForDuplicateSequences = False
var.doCheckForAllGapColumns = False

def upper_triangle_of_distance_matrix(matrix):
    array = numpy.array(matrix, dtype=float)
    return array[np.triu_indices(array.shape[0], k = 1)]

def do_tests(align, protein_alig_name, MPTH_test, alpha, pvalue_correction_method, output):

    #get taxon pairs tested
    taxon_pairs = []
    for taxon in align.taxNames:
        taxon_pairs.append([(taxon,taxon2) for taxon2 in align.taxNames])
    taxa_array = np.array(taxon_pairs)
    taxon_pairs = taxa_array[np.triu_indices(taxa_array.shape[0], k = 1)]

    #MPTH - Calculate p-values
    QB, QS, QR, PB, PS, PR = align.matchedPairsTests()
    bow = upper_triangle_of_distance_matrix(PB.matrix)
    stu = upper_triangle_of_distance_matrix(PS.matrix)
    abab = upper_triangle_of_distance_matrix(PR.matrix)

    #Do the three tests - default: 'all'
    test = ['bowker','stuart','ababneh']
    matrices = [bow, stu, abab]
    # or just one
    if MPTH_test != 'all':
        matrices = [a for a,b in zip(matrices,test) if b == MPTH_test]
        test = MPTH_test.split()

    offensive_taxa = []
    for MPTH,pvalues in zip(test,matrices):
        if pvalue_correction_method:
            try:
                print("Converting nan's: ",np.count_nonzero(np.isnan(pvalues)))
                np.nan_to_num(pvalues, copy=False, nan=0.5)
                assumption_status, pvalues_corrected, alphacSidak, alphacBonf = smt.multipletests(pvalues, method=pvalue_correction_method)
            except ValueError:
                sys.exit('Method does not exit. Choose from bonferroni, sidak, holm, holm-sidak, fdr_bh, fdr_by, fdr_tsbh, or fdr_tsbky.')

        else:
            my_pvalue_threshold = float(alpha)
            assumption_status = [(n < my_pvalue_threshold) for n in pvalues]
            pvalues_corrected = pvalues
            alpha = my_pvalue_threshold

        print(f"\n{30*'#'} {MPTH}'s Test {30*'#'}")
        print(f"Test of {protein_alig_name} -> Count of Failed Tests under a alpha ({alpha}) {pvalue_correction_method} = {list(assumption_status).count(True)} Percentage = {list(assumption_status).count(True)/len(pvalues)}")

        pvalues_taxa = list()
        TH_pvalues_taxa = list()
        for assumption_rejected,pvalue,taxon_pair in zip(assumption_status, pvalues_corrected, taxon_pairs):
            pvalues_taxa.append((pvalue,list(taxon_pair)))
            if assumption_rejected == True:
                TH_pvalues_taxa.append((pvalue,list(taxon_pair)))
                ##################################################################
                with open(f'pvalues_resume.txt','a') as f:
                    #MPTH_test   pvalue_adjust_method    status  pvalue
                    f.write(f'{protein_alig_name}\t{taxon_pair}\t{MPTH_test}\t{pvalue_correction_method}\trejected\t{str(pvalue)}')
                    f.write('\n')
            else:
                with open(f'pvalues_resume.txt','a') as f:
                    f.write(f'{protein_alig_name}\t{taxon_pair}\t{MPTH_test}\t{pvalue_correction_method}\taccepted\t{str(pvalue)}')
                    f.write('\n')

                ##################################################################

        if 1:
            if len(TH_pvalues_taxa) == 0:
                with open(output,'a') as f:
                    f.write(f"{protein_alig_name};{MPTH};0;{len(assumption_status)}\n")
                print(f"No time-heterogeneous taxon sequences detected using the {MPTH} test.")

            else:
                def get_taxa(TH_list):
                    taxon1=sorted_list[0][1][0]
                    taxon2=sorted_list[0][1][1]
                    taxon_pairs = [n[1] for n in sorted_list]
                    count1 = sum(1 for each in taxon_pairs if taxon1 in each)
                    count2 = sum(1 for each in taxon_pairs if taxon2 in each)
                    print(f"The lowest p-value is {sorted_list[0][0]} from the pairwise comparison between {taxon1} and {taxon2}.")
                    if count1 > count2:
                        return taxon1
                    elif count2 > count1:
                        return taxon2
                    else:
                        # Occasionally two taxon sequences has the same number of pairwise rejections
                        # The taxon with the lowest p-values summed (overall) is discarded.
                        pvalue_sum_1 = sum([n[0] for n in pvalues_taxa if taxon1 in n[1]])
                        pvalue_sum_2 = sum([n[0] for n in pvalues_taxa if taxon2 in n[1]])
                        if pvalue_sum_1 < pvalue_sum_2:
                            print(f'{taxon1} and {taxon2} had the same number of occurences ({count1}). But {taxon1} p-values (all summed) were lower.')
                            return taxon1
                        elif pvalue_sum_2 < pvalue_sum_1:
                            print(f'{taxon2} and {taxon1} had the same number of occurences ({count1}). But {taxon2} p-values (all summed) were lower.')
                            return taxon2
                        else:
                            print("The two taxa failed the same number of pairwise comparasions and have the same p-value summed.")
                            return (taxon1,taxon2)

                TH_taxon_offenders = []
                while len(TH_pvalues_taxa) !=0:
                    sorted_list = sorted(TH_pvalues_taxa, key=lambda x: x[0])
                    TH_taxon = get_taxa(sorted_list)
                    if type(TH_taxon) == tuple:
                        print(f"{TH_taxon[0]} and {TH_taxon[1]} had the same number of occurences.\n")
                        TH_taxon_offenders.append(TH_taxon[0])
                        TH_taxon_offenders.append(TH_taxon[1])
                        with open(f"{protein_alig_name}_TH_taxa_{MPTH_test}_{pvalue_correction_method}.txt",'a') as f:
                            f.write(f"{MPTH}'s {pvalue_correction_method} of {protein_alig_name}: p-value = {sorted_list[0][0]} - offender {TH_taxon[0]}\n")
                            f.write(f"{MPTH}'s {pvalue_correction_method} of {protein_alig_name}: p-value = {sorted_list[0][0]} - offender {TH_taxon[1]}\n")
                        TH_pvalues_taxa =[n for n in TH_pvalues_taxa if TH_taxon[0] not in n[1]]
                        TH_pvalues_taxa =[n for n in TH_pvalues_taxa if TH_taxon[1] not in n[1]]

                    else:
                        print(f"{TH_taxon} had more occurences.\n")
                        TH_taxon_offenders.append(TH_taxon)
                        with open(f"{protein_alig_name}_TH_taxa_{MPTH_test}_{pvalue_correction_method}.txt",'a') as f:
                            f.write(f"{MPTH}'s {pvalue_correction_method} of {protein_alig_name}: p-value = {sorted_list[0][0]} - offender {TH_taxon}\n")
                        TH_pvalues_taxa =[n for n in TH_pvalues_taxa if TH_taxon not in n[1]]

                with open(output,'a') as f:
                    f.write(f"{protein_alig_name};{MPTH};{len([n for n in assumption_status if n == True])};{len(pvalues_corrected)}\n")

                for each in TH_taxon_offenders:
                    if each not in offensive_taxa:
                        offensive_taxa.append(each)

    return(offensive_taxa)

#######################################################################################
def main(datafile, remove_taxa, test, pvalue_correction_method, alpha, overwrite, output):

    if not datafile:
        print("Error: No datafile specified.")
        sys.exit()
    else:
        if not os.path.exists(datafile):
            print("Error: cannot find datafile %s" % datafile)
            sys.exit()

    align = func.readAndPop(datafile)
    print(f"\n{40*'#'} {align.fName}'s Test {40*'#'}")

    #get alignment file name
    filename = os.path.basename(align.fName)
    protein_alig_name = os.path.splitext(filename)[0]
    #do test
    offenders = do_tests(align,protein_alig_name, test, alpha, pvalue_correction_method, output)

    if remove_taxa:
        if offenders != None:
            print("Removing tree-heterogeneous sequences...")

            itaxa = align.taxNames
            total_offensive_taxa = []

            for taxon in offenders:
                align.sequences.pop(align.taxNames.index(taxon))
                total_offensive_taxa.append(taxon)

            if overwrite:
                print('Alignment file will be overwrite in nexus format')
                align.writeNexus(f'{protein_alig_name}.nex')
            else:
                align.writeFasta(f'{protein_alig_name}_{test}_{pvalue_correction_method}_THfree.fasta')
        else:
            align.writeFasta(f'{protein_alig_name}_{test}_{pvalue_correction_method}_THfree.fasta')

######################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(desc),
            )
    parser.add_argument("datafile",
                        help="Data matrix filename.",
                        )
    parser.add_argument("-r", "--remove_taxa",
                        help="Remove tree-heterogeneous sequences.",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-t", "--test",
                        help="Choose one of the MPTH (bowker, stuart or ababneh). Default = Do the three.",
                        default = 'all',
                        type = str
                        )
    parser.add_argument("-c", "--pvalue_correction_method",
                        help="Method to correct p-values. Options:'bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky' ",
                        )
    parser.add_argument("-a", "--alpha",
                        help="P-value threshold. Default = 0.05",
                        default = 0.05,
                        type = float
                        )
    parser.add_argument("-i", "--overwrite",
                        help="Overwrite the current alignment file. Only works for nexus format now. Default=False",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-o", "--output",
                        help="A output line with file, test, failed tests, and total tests. It is written in an append-way. Default=results.csv",
                        default="results.csv"
                        )

    args = parser.parse_args()
    main(args.datafile, args.remove_taxa, args.test, args.pvalue_correction_method, args.alpha, args.overwrite, args.output)
