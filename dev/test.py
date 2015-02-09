from igraph import *
from Bio import SeqIO
import re
import time


# Main function for working
def main():

    # fasta file for corresponding genes
    fasta_filename = "ESR1.fasta"

    print('\nSuper transcript magic algorith!... \n')
    print("\n>> Retrieving transcript IDs...")
    transcripts = get_transcript_ids(fasta_filename)

    print("\n>> Constructing empty graph structure and hash table... ")
    global graph
    graph = Graph(directed=True)

    global transcript_hash_table
    transcript_hash_table = {}

    # Sets graph attributes, so each node has an empty hash table with transcript as key
    for x in transcripts:
        graph.vs[x] = ""
        transcript_hash_table[x] = []
    graph.vs["Base"] = ""

    # Get the block co-ordinates from BLAT output
    output = psl_parse('output_2.psl')

    # Gets the sequence from fasta file. Also adds the nodes to the graph as sequences are retrieved
    blat_blocks = get_sequence(output, "ESR1.fasta")

    # Simplify graph conactenating nodes that do not fork
    print("Simplifying graph")
    simplify_graph(graph)

    # Graph summary
    print(graph.summary())

    # Sorts graph using khan algorithm
    sorted_sequence = khan_sort(graph)

    # Write graph information to text file for checking
    f = open('graph_output.txt', 'w')
    for vertex in graph.vs():
        f.write(str(vertex))
        f.write("\n")
    f.write("\n")
    for key in transcript_hash_table.keys():
        f.write(key)
        f.write(" : ")
        for i in transcript_hash_table[key]:
            f.write(str(i))
            f.write(" ")
        f.write("\n")

    f.write("sorted sequecne \n")
    f.write(sorted_sequence)
    f.close()

    # Write graph structure to SVG file
    graph.write_svg("small_graph.svg", layout="kk")

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

                sorted_sequence += x.__getitem__("Base")
                graph.delete_vertices(x)
                num_nodes -= 1

    print(">> Graph sorted, returning ordered nodes....")

    return sorted_sequence


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
    global transcript_hash_table

    # print(list_transcripts)

    # Checks if this is the first base being added to the graph
    # If this is the first base added to the graph add the node
    if len(graph.vs()) == 0:
        graph.add_vertex()
        graph.vs[len(graph.vs())-1]["Base"] = base

        for z, ba in enumerate(list_transcripts):
            graph.vs[len(graph.vs())-1][list_transcripts[z]] = coordinate[z]
            # Added transcript coordinates to hash table
            transcript_hash_table[list_transcripts[z]].append(coordinate[z])

    # Creates list of booleans for evaluating if transcript co-ordinates exist in the graph
    base_exist = []
    [base_exist.append(False) for t in list_transcripts]

    '''
    # TODO Sub graph method - keeping just in case decide but hash table is much faster
    # Sub Graph of main graph containing all the bases of one currently needing to be added
    base_sub_graph = graph.vs.select(Base=base)
    # Updates base_exist list with True at same index of transcripts if they exist in graph
    for v in graph.vs.select(Base=base):
        for t, trans in enumerate(list_transcripts):
            if v.__getitem__(list_transcripts[t]) == coordinate[t]:
                base_exist[t] = True
    '''

    for t, trans in enumerate(list_transcripts):
        list_coords = transcript_hash_table[trans]
        if coordinate[t] in list_coords:
            base_exist[t] = True

    # No transcripts and co-ordinates in there, so add new node - base_exist is still all false
    # May need to re-arrange the logic in adding new nodes - but think its good here
    if True not in base_exist:
        # Adding new nodes to graph
        graph.add_vertex()
        graph.vs[len(graph.vs())-1]["Base"] = base
        # graph.add_edge(len(graph.vs())-2, len(graph.vs())-1)

        for z, ba in enumerate(list_transcripts):
            graph.vs[len(graph.vs())-1][list_transcripts[z]] = coordinate[z]

            # Added transcript coordinates to hash table
            transcript_hash_table[list_transcripts[z]].append(coordinate[z])

    # - If one of the transcripts returns true then do not add new nodes
    # - Update the existing node with the new transcripts
    elif True and False in base_exist:
        true_indices = [i for i, x in enumerate(base_exist) if x is True]
        false_indices = [i for i, x in enumerate(base_exist) if x is False]

        # Updates true index transcripts at nodes with false index coordinates
        for vertex in graph.vs.select(Base=base):
            for t in true_indices:
                if vertex[list_transcripts[t]] == coordinate[t]:
                    for j in false_indices:
                        vertex[list_transcripts[j]] = coordinate[j]
                        transcript_hash_table[list_transcripts[j]].append(coordinate[j])

    # # Add edges for each transcript currently not working,
    # for g, transc in enumerate(list_transcripts):
    #     add_edge(base, transc, coordinate[g])


# TODO look at re-write so that it only runs through the graph once,
# - So look through list of transcripts - add edges the currenty way but do only need to iterate through the graph once
# Add an edge between most recently added edges in the one transcript - current method implemented
# OR
# After nodes are added - iterate through graph along each transcript sequence and co-ordinates added edges,
# if edge exist between nodes already then no need to add them and can simplify that node?
def add_edge(base, transcript, coordinate):

    global graph

    # Don't work very fast - very slow
    # TODO re write add_edge method very slow - need to implement a method to make this faster

    if coordinate > 0:

        for vertex in graph.vs():

            if vertex[transcript] == (coordinate-1):
                node_a = vertex
            if vertex[transcript] == coordinate:
                node_b = vertex

        try:
            if len(graph.es.select(_within=[node_a.index, node_b.index])) == 0:
                graph.add_edge(node_a, node_b)

        except UnboundLocalError:
            print("Didn't find nodes or something")

        simplify_graph(graph)


# TODO Implement graph simplification
# - Travel along each edge, and check the number of in and out degrees
# - If both are 1, then compress the nodes to 1, and change the base value to A+B
# - Coordinates will be the smaller value, can find the co-ordinates by the length of
#   the bases and the starting point of the node
def simplify_graph(graph):

    i = 0

    while i != len(graph.vs()):
        print("Lenght of graph : ", len(graph.vs()))
        print("i : ", i)

    # for each node in the graph
    # for vertex in graph.vs():

        vertex = graph.vs()[i]

        # vertex out degree is int but vertex indegree is list - because neighbour is of vertexseq object
        if vertex.outdegree() == 1:
            neighbor_id = graph.neighbors(vertex, mode="out")
            neighbor = graph.vs()[neighbor_id[0]]

            # there should on be one value in the list
            if neighbor.indegree() == 1:

                v_base = vertex["Base"]
                n_base = neighbor["Base"]
                new_base = v_base + n_base
                vertex["Base"] = new_base

                # Edges are added from vertex to its new neighbours
                for index in graph.neighbors(neighbor_id[0], mode="out"):
                    graph.add_edge(vertex, graph.vs()[index])

                # Neighbour vertex is deleted from the graph
                graph.delete_vertices(neighbor_id)
                print("deleted a node")
            else:
                i += 1
        else:
            i += 1


if __name__ == "__main__":
    start_time = time.clock()
    main()
    print("Running time: ", time.clock()-start_time, "s")


