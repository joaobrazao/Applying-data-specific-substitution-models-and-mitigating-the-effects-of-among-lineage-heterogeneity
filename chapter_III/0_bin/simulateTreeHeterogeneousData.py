#! /usr/bin/env python3

desc="""Alignment data simulator - Compositional and/or rate heterogeneous among lineages.

João Brazão version 1
"""
import os
import re
from p4 import *
import argparse
import textwrap
import random

def main(length, number, model, rates, TH_parameters, TH_compositions, output):

    if output:
        if number != None:
            sys.exit("ERROR. If you set output name, it does just one replicate. Turn off replicate number option.")
    else:
        if not number:
            number = 1

    # Lineage-specific (model) settings
    ## TH_models
    TH_models = TH_parameters.split()

    # 30-taxa tree with equal tip-branch lengths (0.3) and equal internal branches (0.2)
    tree = '(((A1:0.4, (A2:0.4, A3:0.4):0.3):0.3, (A4:0.4, (A5:0.4, A6:0.4):0.3):0.1):0.3, ((((B1:0.4, (B2:0.4, B3:0.4):0.3):0.3, (B4:0.4, (B5:0.4, B6:0.4):0.3):0.1):0.3, ((C1:0.4, (C2:0.4, C3:0.4):0.3):0.3, (C4:0.4, (C5:0.4, C6:0.4):0.3):0.3):0.3):0.1, (((D1:0.4, (D2:0.4, D3:0.4):0.3):0.3, (D4:0.4, (D5:0.4, D6:0.4):0.3):0.1):0.3, ((E1:0.4, (E2:0.4, E3:0.4):0.3):0.3, (E4:0.4, (E5:0.4, E6:0.4):0.3):0.3):0.3):0.3):0.3);'
    tree0 = func.readAndPop(tree)

    # Preparing tree
    node_tips = [1, 14, 25, 37, 48]
    random_tips = random.sample(node_tips,len(TH_models))
    combined_list = list(zip(random_tips, TH_models))
    #Re-do tree with tips name changed (using TH_models model name): important for downstrean analyses.
    for each in combined_list:
        for mynode in tree0.getAllLeafNames(tree0.node(each[0])):
            tree = tree.replace(tree0.node(mynode).name,tree0.node(mynode).name + '_' + each[1])
    mytree = func.readAndPop(tree)
    taxNames = mytree.getAllLeafNames(0)

    # Make data
    a = func.newEmptyAlignment(dataType='protein', taxNames=taxNames, length=length)
    d = Data([a])
    mytree.data = d

    # General composition and rates
    ## Rates
    if rates:
        r0 = mytree.newRMatrix(free=0, spec = model)
        mytree.setModelComponentOnNode(r0, node=mytree.root, clade=1)
    else:
        mytree.newRMatrix(free=0, spec = model)
        print('Non-stationary compositions.')

    if TH_compositions:
        c0 = mytree.newComp(partNum=0, free=0, spec = model)
        mytree.setModelComponentOnNode(c0, node=mytree.root, clade=1)
    else:
        mytree.newComp(partNum=0, free=0, spec = model)

    # Lineage specific sequences/branches
    for each in combined_list:

        if TH_compositions:
            c1 = mytree.newComp(partNum=0, free=0, spec=each[1])
            mytree.setModelComponentOnNode(c1, node=each[0], clade=1) #compositions
        if rates == 'free':
            print('Lineage-specific rates are free.')
            r1 = mytree.newRMatrix(free=1)
            mytree.setModelComponentOnNode(r1, node=each[0], clade=1) #rates
        elif rates == 'fixed':
            print('Lineage-specific rates are fixed.')
            r1 = mytree.newRMatrix(free=0, spec = each[1])
            mytree.setModelComponentOnNode(r1, node=each[0], clade=1) #rates
        else:
            print('Rates-homogeneous simulation.')

    # No among-site rate variation, to keep it simple
    if 1:
        mytree.setNGammaCat(nGammaCat=1)
        mytree.setPInvar(free=0, val=0.0)
    else:
        mytree.setNGammaCat(nGammaCat=4)
        mytree.newGdasrv(free=0, val=0.75)
        mytree.setPInvar(free=0, val=0.0)

    if output:
        if not var.gsl_rng:
            var.gsl_rng = pf.gsl_rng_get()
        mySeed = int(output)    # your chosen int seed
        pf.gsl_rng_set(var.gsl_rng, mySeed)

        mytree.simulate()
        #d.writeNexus(f"alignment_{length}_{output}.nex", writeDataBlock=True)
        d.alignments[0].writePhylip(f"alignment_{length}_{output}.phy")
    else:
        func.reseedCRandomizer(os.getpid())
        for n in range(1,number+1):
            mytree.simulate()
            d.writeNexus(f"alignment_{length}_{n}.nex", writeDataBlock=True)
            a.writePhylip(f"alignment_{length}_{n}.phy")
        print ("Simulated alignments created.")

#################################################################

if __name__ == "__main__":

    parser=argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(desc),
            )

    parser.add_argument("length",type=int,
            help="Simulated alignment length."
            )
    parser.add_argument("-n","--number",type=int,
            help="Number of replicates."
            )
    parser.add_argument("-m","--model",
            help="Rates and compositions for the default model.\nP4 has 'equal', 'empirical', 'specified','cpREV', 'd78', 'jtt', 'mtREV24', 'mtmam', 'wag', 'blosum62', 'rtRev', 'tmjtt94', 'tmlg99', 'lg', 'hivb', 'mtart', 'mtzoa', 'gcpREV', 'stmtREV', 'vt', 'pmb'.\nDefault = wag",
            default='wag'
            )
    parser.add_argument("-r","--rates",
            help="'free' or 'fixed' to the model used for instantaneous rates.",
            )
    parser.add_argument("-t","--TH_parameters",
            help="Substitution models used to feed the TH rate parameters. Set the rates to fixed. Default: 'jtt cpREV'",
            default='jtt cpREV'
            )
    parser.add_argument("-c", "--TH_compositions",
                        help="Tree-heterogeneous sequences. These sequences have the same compositions of the model indicated in 'TH_parameters'. Default = False.",
                        default=False,
                        action='store_true'
                        )
    parser.add_argument("-o","--output",
            help="Output file name.",
            )

    args=parser.parse_args()
    main(args.length, args.number, args.model, args.rates, args.TH_parameters, args.TH_compositions, args.output)
