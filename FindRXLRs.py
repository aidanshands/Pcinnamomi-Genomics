#!/usr/bin/env python
#-------------------------------------------------------------------------------
# Aidan Shands
# This script will search protein fasta files and will export a .fasta file with
# proteins that match the regular expression and a .txt file that contains the
# protein ID and the motif position.
# Usage: python FindRXLRs.py -i input.fasta
#-------------------------------------------------------------------------------
import sys, re, argparse
from Bio import SeqIO
#-------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input file')
args = parser.parse_args()
Fasta= args.input
#-------------------------------------------------------------------------------
def find_RXLRs(fasta, motif):
    if motif == "without -EER":
        motifs = {"RXLR": r"^[a-yA-Y]{10,110}?R[A-Y]LR",
                  "QXLR": r"^[a-yA-Y]{10,110}?Q[A-Y]LR",
                  "RXL": r"^[a-yA-Y]{10,110}?R[A-Y]L[^R]",
                  "XLR": r"^[a-yA-Y]{10,110}?[^QR][A-Y]LR"}
    if motif == "with -EER":
        motifs = {'QXLR-EER': r"^[a-yA-Y]{10,110}Q[A-Y]LR[A-Y]*EER",
                  'RXL-EER': r"^[a-yA-Y]{10,110}R[A-Y]L[A-Y]*EER",
                  'RXLR-EER': r"^[a-yA-Y]{10,110}R[A-Y]LR[A-Y]*EER",
                  'XLR-EER': r"^[a-yA-Y]{10,110}[A-Y][A-Y]LR[A-Y]*EER",
                  "EER": r"^[a-yA-Y]{34,150}?EER"}

    for motif, regex in motifs.items():
        fasta_sequences = SeqIO.parse(open(fasta),'fasta')
        ID_File = motif + "_ID_and_Motifs.txt"
        result_file = motif + '_Positives.fasta'
        with open(result_file, 'w') as f, open(ID_File, 'w') as ID:
            for sequence in fasta_sequences:
                matches = re.search(regex, str(sequence.seq))
                if matches:
                    SeqIO.write([sequence], f, "fasta")
                    ID_Loc = matches.end()
                    if motif == "without -EER":
                        ID_Pos = "%s %s %d" % (sequence.id + '\t', motif + 'motif begins at:''\t', (ID_Loc - 3))
                    else:
                        ID_Pos = "%s %s %d" % (sequence.id + '\t', 'EER motif begins at:''\t', (ID_Loc - 2))
                    ID.write(str(ID_Pos) + "\n")
#-------------------------------------------------------------------------------
find_RXLRs(Fasta, "without -EER")
find_RXLRs(Fasta, "with -EER")
