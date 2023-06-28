#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import argparse, re
# Aidan Shands
# USAGE:
# python Simulate_L.py -i FIRS.csv -sco Single_Copy_Orthologs.csv -s 100 -e 5100 -b 100
# Single_Copy_Orthologs.csv is a csv with 1 column (ID), and a gene id one per line
# FIRS.csv is the output file from Calculate_FIR_length.pl script. 
# Citation: Saunders, D. G. et al. (2014). (https://figshare.com/authors/Sylvain_Raffaele/415674)
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='FIRs.csv')
parser.add_argument('-sco', '--sco', help='Single-copy orthologs.csv')
parser.add_argument('-s', '--start', help='range start')
parser.add_argument('-e', '--end', help='range end')
parser.add_argument('-b', '--by', help='range iterations')
args = parser.parse_args()
FIR_input = args.input
SCO_input = args.sco
Range_start = int(args.start)
Range_end = int(args.end)
Range_by = int(args.by)
#-------------------------------------------------------------------------------
FIRs = pd.read_csv(FIR_input)
SCOs = pd.read_csv(SCO_input)
# If the IDs have a "." then replace the "." and values after with empty string
SCOs['ID'] = SCOs['ID'].map(lambda x: re.sub(r'\.\d+$', '', x) if '.' in x else x)
SCOs_L = SCOs["ID"].to_list() 
L_list = [*range(Range_start, Range_end, Range_by)]# Define the L-values to simulate
# Function that simulates L
def Lsim(PcFIRs, Ortho_List, Type):
    Summary = pd.DataFrame(columns = ['L-value', 'Length_of_Dataset', 'Ortholog_Count'])
    Gene_dict = {}
    for name in L_list:
        Gene_dict[name] = pd.DataFrame()
    keys_values = Gene_dict.items()
    Gene_dict = {str(key): value for key, value in keys_values}
    for (name, df), L in zip(Gene_dict.items(), L_list):
        # Gene-sparse
        if Type == "GSR":
            df = FIRs.loc[(FIRs['fiveprime'] >= L) & (FIRs['threeprime'] >= L)].copy()
        # In-between
        elif Type == "BTWN":
            df = FIRs.loc[((FIRs['fiveprime'] >= L)) & ((FIRs['threeprime'] < L)).copy() | ((FIRs['fiveprime'] < L)) & ((FIRs['threeprime'] >= L))].copy()
        # Gene-dense
        elif Type == "GDR":
            df = FIRs.loc[(FIRs['fiveprime'] < L) & (FIRs['threeprime'] < L)].copy()
        df.loc[df.geneid.isin(Ortho_List), "Ortholog"] = "TRUE" 
        # Writing results to summary
        Summary.loc[name, 'L-value'] = str(name) 
        Summary.loc[name, 'Length_of_Dataset'] = len(df.index) 
        Summary.loc[name, 'Ortholog_Count'] = df.loc[df.Ortholog == 'TRUE', 'Ortholog'].count()
        Summary.to_csv(Type +"_Simulated_L.csv", index = None)
#-------------------------------------------------------------------------------
# Run function
Lsim(FIRs, SCOs_L, "GDR")
Lsim(FIRs, SCOs_L, "GSR")
Lsim(FIRs, SCOs_L, "BTWN")
print('Success!')
