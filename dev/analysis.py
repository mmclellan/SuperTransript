__author__ = 'michaelmclellan'

from Bio import *
from Bio import SeqIO

fasta = "sorted.fasta"

for seq_record in SeqIO.parse(fasta, "fasta"):
    seq = seq_record.seq

print(seq[0:292])
print(seq[1571:1571+292])