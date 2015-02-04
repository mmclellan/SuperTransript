from igraph import *
from Bio import SeqIO
import re

# Main function for working
def main():

    # fasta file for corresponding genes
    fasta_filename = "ESR1.fasta"

    print('\nSuper transcript magic algorith!... \n')

    print("\n>> Retrieving transcript IDs...")
    transcripts = get_transcript_ids(fasta_filename)

    print("\n>> Constructing empty graph structure... ")
    global graph
    graph = Graph(directed=True)

    # Sets graph attributes, so each node has an empty hash table with transcript as key
    for x in transcripts:
        graph.vs[x] = ""
    graph.vs["Base"] = ""

    output = psl_parse('output_3.psl') # Parse the blat file to retrieve genomic information about transcripts
    blat_blocks = get_sequence(output, "ESR1.fasta") # Retrieves sequences for genomic blocks

    print(graph.summary())

    f = open('graph_output.txt', 'w')
    for vertex in graph.vs():
        f.write(str(vertex))
        f.write("\n")
    f.close()

    # Constructs a graph using hard coded edges and nodes for ESR1 gene and sorts it
    # graph = construct_graph()
    # draw_svg_graph(graph)
    # print(khan_sort(graph))


# Constructs graph object with nodes being unique blocks
# Uses igraph to handle objects
def construct_graph_hardCode():

    print(">> Constructing graph....")

    # Hard coded graph for sorting algorithm
    g = Graph(directed=True)

    n = 12
    i = 0
    while i != 12:
        g.add_vertex(name=i, sequence="abcd")
        i += 1

    g.add_edge(0, 1)
    g.add_edge(1, 3)
    g.add_edge(2, 3)
    g.add_edge(3, 4)
    g.add_edge(3, 8)
    g.add_edge(4, 5)
    g.add_edge(5, 6)
    g.add_edge(5, 7)
    g.add_edge(6, 7)
    g.add_edge(7, 8)
    g.add_edge(8, 9)
    g.add_edge(9, 10)
    g.add_edge(10, 11)

    # g.write_svg("ESR1.svg", layout = "drl")

    return g


# Method to retrieve list of transcript IDs from a fasta
def get_transcript_ids(fasta):

    transcript_ids = []
    filename = fasta

    for seq_record in SeqIO.parse(filename, "fasta"):
        transcript_ids.append(seq_record.id)

    return transcript_ids

# TODO Need to re-write khan algorithm to work with base focused graph instead of large blocks of sequence
# - Depends on how I simplify and compact the graph
# Topologically sorts a graph using the Khan algorithm
# seems to work without worrying about deleting edges, simply takes care of itself? Sweet
def khan_sort(graph):

    print(">> Sorting graph using khan algorithm....")
    # creates list of nodes in the graph and number of nodes
    node_list = graph.vs()
    num_nodes = len(graph.vs())

    S = [] # List for storing final sorted elements
    L = [] # list for storing nodes with no incoming edges # seems to work only using L
    sorted_sequence = ''

    while num_nodes > 0:
        # Find nodes with no incoming edges - add them to L
        for x in node_list:
            if x.degree(type="in") == 0:
                # print("Node with no incoming edges name:  ", x.__getitem__("name"))
                L.append(x.__getitem__("name"))
                sorted_sequence += x.__getitem__("sequence")
                graph.delete_vertices(x)
                num_nodes -= 1
    print(">> Graph sorted, returning ordered nodes....")
    return L


# Takes blat .psl file and returns the genomic co-ordinates for the blocks determined by BLAT
def psl_parse(filename):

    print(">> Opening psl file....")
    f = open(filename, 'r')
    blat_blocks = []

    for i, line in enumerate(f):
        if i > 4 and line and line != '\n':
            split = re.split('\t', line)
            for b, bl in enumerate(split):
                if not bl:
                    split.pop(b)

            blocks = [split[9], split[10], split[13], split[14], split[18], split[19], split[20]]

            j = 4
            while j < 7:
                blocks[j] = [x.strip() for x in blocks[j].split(',')]
                blocks[j] = [x for x in blocks[j] if x != '']
                j += 1

            blat_blocks.append(blocks)

    print(">> Blat blocks found.... \n")
    return blat_blocks


# Retrieves the sequence data for the blocks determined by BLAT and adds found nodes to graph structure
# May need to add +1 to block starts and ends to account for 0-indexing
def get_sequence(blat_blocks, fasta_file):

    global graph

    print(">> Finding genomic sequence from from fasta file....")
    print(">> Opening fasta file...")

    output_records = []

    for block in blat_blocks:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if blat_blocks[0][0] == seq_record.id:  # Finds the correlating sequence from the fast file

                first_base = True  # True if this is the first base of the transcript

                for i, z in enumerate(block[4]):

                    print("\nPrinting i: ", i)

                    num_blocks = len(block[4])

                    block_size = int(block[4][i])
                    query_id = block[0]
                    query_block_start = int(block[5][i])
                    query_block_end = int(block[5][i])+block_size
                    query_length = int(block[1])

                    target_id = str(block[2])
                    target_block_start = int(block[6][i])
                    target_block_end = int(block[6][i])+block_size
                    target_length = int(block[3])
                    print(query_id, " ", target_id)

                    # while query ID and target ID do not match

                    if query_id != target_id:

                        # adds unaligned section before aligned blocks
                        if i == 0 and query_block_start != 0 and first_base is True:
                            for b, base in enumerate(seq_record[0:query_block_start]):
                                add_node_graph(base, [query_id], [b])
                            print("Added block: ", query_id, " ", 0, "-", query_block_start)
                            first_base = False

                        if i == 0 and target_block_start != 0 and first_base is True:
                            for b, base in enumerate(seq_record[0:target_block_start]):
                                add_node_graph(base, [target_id], [b])
                            print("Added block: ", target_id, " ", 0, "-", target_block_start)
                            first_base = False

                         # add basses from qstart and tstart length of block size
                        for b, base in enumerate(seq_record[query_block_start:query_block_end]):
                            add_node_graph(base, [query_id, target_id], [query_block_start+b, target_block_start+b])
                        print("Added block: ", query_id, " ", query_block_start, "-", query_block_end)
                        print("Added block: ", target_id, " ", target_block_start, "-", target_block_end)

                        # Adds bases after blocks to end of query or target sequence
                        if i == num_blocks-1 and query_block_end != query_length:
                            for b, base in enumerate(seq_record[query_block_end:query_length]):
                                add_node_graph(base, [query_id], [(query_block_end+b)])
                            print("Added block: ", query_id, " ", query_block_end, "-", query_length)

                        if i == num_blocks-1 and target_block_start != target_length:
                            for b, base in enumerate(seq_record[target_block_end:target_length]):
                                add_node_graph(base, [target_id], [(target_block_end+b)])
                            print("Added block: ", target_id, " ", target_block_end, "-", target_length)

    print(">> Returning block sequence....\n")
    return output_records


# Writes list of block sequences to fasta file.
# To Do: need to rework how to manage naming of blocks as they currently share names from query sequence
def write_fasta(list_records, filename):
    SeqIO.write(list_records, filename, 'fasta')


# Method for adding nodes to the graph; updates graph nodes
def add_node_graph(base, list_transcripts, coordinate):

    global graph

    # If this is the first base added to the graph add the node
    if len(graph.vs()) == 0:
        graph.add_vertex()
        graph.vs[len(graph.vs())-1]["Base"] = base

    # TODO need to check if nodes already exist
    # - search if base in transcript[0] with co-ordinate exists:
    # - if it does - do nothing
    # - if it doesn't - update node by adding co-ordinate for that transcript

    base_exist = []
    [base_exist.append(False) for t in list_transcripts]

    # TODO find a more efficient way to do this, this is going to get out of hand when working on larger data sets
    # Sub Graph of main graph containing all the bases of one currently needing to be added
    base_sub_graph = graph.vs.select(Base=base)

    for v in graph.vs.select(Base=base):
        for t, trans in enumerate(list_transcripts):
            if v.__getitem__(list_transcripts[t]) == coordinate[t]:
                base_exist[t] = True

    # No transcripts and co-ordinates in there, so add new node
    # May need to re-arrange the logic in adding new nodes - but think its good here
    if True not in base_exist:

        # Adding new nodes to graph
        graph.add_vertex()
        graph.vs[len(graph.vs())-1]["Base"] = base
        graph.add_edge(len(graph.vs())-2, len(graph.vs())-1)

        for z, ba in enumerate(list_transcripts):
            graph.vs[len(graph.vs())-1][list_transcripts[z]] = coordinate[z]

    else:

        # TODO Fix this method of checking if nodes exist
        # - If one of the transcripts returns true then do not add new nodes
        # - Update the existing node with the new transcripts
        # - Fix, so that if base_exist is completely try do not add at all.

        if True and False in base_exist:
            # get index that gave true
            true_indices = [i for i, x in enumerate(base_exist) if x is True]

            # Get index that are false
            false_indices = [i for i, x in enumerate(base_exist) if x is False]

            for x in true_indices:
                # gets vertex from sub graph and checks that the node exists
                for vertex in base_sub_graph:
                    if vertex.__getitem__(list_transcripts[x]) == coordinate[x]:
                        # Updates the truth node with the nodes that are not in the graph
                        for j in false_indices:
                            vertex[list_transcripts[j]] = coordinate[j]

# TODO Need to write how to add edges
# - Adding edges between two most recent nodes
# - Look at if working with a different transcript? Not sure


def check_base_in_graph(base, transcripts, coords):

    global graph

    for x in graph.vs.select(Base=base):
        if x.__getitem__(transcripts[0]) == coords[0]:
            print("Item is in the graph")

# TODO Implement graph simplification
# - Compress nodes

if __name__ == "__main__":
    main()



