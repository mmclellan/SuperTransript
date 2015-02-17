__author__ = 'michaelmclellan'

from Bio import *
from Bio import SeqIO

fasta = "ESR1_reconstructed_transcripts.fasta"

for seq_record in SeqIO.parse(fasta, "fasta"):

    print(seq_record.id)
    print(len(seq_record.seq))