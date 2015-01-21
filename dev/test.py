import re
import SeqBlock
from igraph import *
import sys
from Bio import SeqIO




# Takes blat .psl file and returns the genomic co-ordinates for the blocks
def psl_parse(filename):
    f = open(filename, 'r')
    blat_blocks = []

    for i, line in enumerate(f):
        if i > 4 and line and line != '\n':
            split = re.split('\t', line)
            for b, bl in enumerate(split):
                if not bl:
                    split.pop(b)

            blocks = [split[9], split[13], split[18], split[19]]

            j = 2
            while j < 4:
                blocks[j] = [x.strip() for x in blocks[j].split(',')]
                blocks[j] = [x for x in blocks[j] if x != '']
                j += 1

            blat_blocks.append(blocks)

    return blat_blocks

output = psl_parse('output_2.psl')
# for x in output:
#     print x


# Retrieves the sequence data for the blocks determined by BLAT
# Need to determine what to do with sequence data and how to store it and then
# Translate it to a graph.
def get_sequence(blat_blocks, fasta_file):
    output_records = []
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        for x in blat_blocks:
            if x[0] == seq_record.id:
                for i, z in enumerate(x[3]):
                    block_start = int(z)
                    block_size = int(x[2][i])
                    rec = seq_record[block_start:(block_start+block_size)]
                    output_records.append(rec)

    # print output_records
    return output_records


# Writes list of block sequences to fasta file.
# To Do: need to rework how to manage naming of blocks as they currently share names from query sequence
def write_fasta(list_records, filename):
    SeqIO.write(list_records, filename, 'fasta')

fasta = "ESR1.fasta"
blocks = get_sequence(output, fasta)
write_fasta(get_sequence(output, fasta), "final_success")


# for i, block in enumerate(blocks):
#     print i
#     print block
#     adjacency_list[SeqBlock.Block.num_blocks]


# Creates hash table using number of blocks as keys and uses block objects as values
def adjacency_list():

    adjacency_list = {}

    for b in blocks:

        print "Current block", b.id
        # check if substring is found in current values in the hash table
        for x in adjacency_list.keys():

            print adjacency_list[x].query_id, " -> ", adjacency_list[x].target_id

            if adjacency_list[x].sequence.find(b.seq) == -1:
                print "Not found"
            else:
                print "String found at: ", adjacency_list[x].sequence.find(b.seq)
            print "\n"

        adjacency_list[SeqBlock.Block.num_blocks] = SeqBlock.Block(b.seq, b.id, 'id2')

    # print adjacency_list.keys()
    # print adjacency_list.values()

seq = []

for block in blocks:
    seq.append(str(block.seq))
# print seq

k = 11 # k-mer value may need revising


kmer_list = {}

for j in seq:
    x = 9
    for c in j:
        kmer = (j[x:k+x])

        if len(kmer) == k:
            if not (kmer in kmer_list):
                kmer_list[kmer] = 0
            kmer_list[kmer] += 1
        x += 1

#print kmer_list

# for x in kmer_list:
#     print x


dbg_edge_elements = set()

edges = []

for kmer in kmer_list:
    for x in kmer_list:
        if kmer[0:k-1] == x[1:k]:
            e = [x, kmer]
            edges.append(e)


# for kmer in kmer_list:
#     dbg_edge_elements.add(kmer)
#
# edge = lambda element: '('+element[0:k-1]+ ', '+element[1:k]+')'
#
# dbg_edge_elements = [edge(element) for element in dbg_edge_elements]
# print dbg_edge_elements

g = Graph()

for x in kmer_list:
    g.add_vertices(x)
for x in edges:
    g.add_edge(x[0], x[1])

# print(g)

# Create a PNG plot of the graph - currently doesn't work
# layout = g.layout("kamada_kawai")
# plot(g, layout=layout)


print(g.get_edgelist())



# cut up blocks into k-mers
# use dictionary to construct the de bruijn graph
# keys are k-mers and values are the count for k-mers
# see what happens, it might work fine and produce
# using k-mers create de bruijn graph
# sort de bruijn graph to produce a topological order



