#!/usr/bin/env python
# coding: utf-8
#-------------------------------------------------------------------------------
import pandas as pd
import argparse
#-------------------------------------------------------------------------------
# Author: Aidan Shands
# USAGE
# python Deliminate_Genome.py -i Firs.csv -l l-value
# Outputs:
#
# References:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3091767/
# https://www.frontiersin.org/articles/10.3389/fgene.2020.00579/full
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input file')
parser.add_argument('-l', '--value', help='L-value')
args = parser.parse_args()
input_FIR= args.input
L = int(args.value)
#-------------------------------------------------------------------------------
FIRs = pd.read_csv(input_FIR)
GSR = FIRs.loc[(FIRs['fiveprime'] >= L) & (FIRs['threeprime'] >= L)].copy()
GDR = FIRs.loc[(FIRs['fiveprime'] < L) & (FIRs['threeprime'] < L)].copy()
BTW = FIRs.loc[((FIRs['fiveprime'] >= L)) & ((FIRs['threeprime'] < L)).copy() | ((FIRs['fiveprime'] < L)) & ((FIRs['threeprime'] >= L))].copy()
ND = FIRs[FIRs.isna().any(axis=1)].copy()

GSR["Type"] = "GSR"
GDR["Type"] = "GDR"
BTW["Type"] = "Inbetween"
ND["Type"] = "ND"

Final_FIRs = pd.concat([GSR, GDR, BTW, ND])

Final_FIRs.to_csv("Delimited_FIRs.csv", index = None)
