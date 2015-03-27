__author__ = 'michaelmclellan'

import networkx as nx
import re
import time
import bisect
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

## TO-DO LIST ##
# TODO - Look at validation method, once graph is sorted,
# - check each node that the adjacent nodes transcript co-ordinate is not lower then its parent
# TODO - Retrieve transcript sequences from graph, make sure the graph is correctly holding data


# Main function for working
def main():

    parameters = sys.argv

    input_filename = parameters[1]

    # fasta file for corresponding genes
    fasta_filename = "ESR1_transcripts_for_Michael.fasta"

    print('\nSuper transcript magic algorithm!... \n')
    print("\n>> Retrieving transcript IDs...")

    transcripts = get_transcript_ids(fasta_filename)

    print("\n>> Constructing empty graph structure and hash table... ")

    # Global variables #
    global graph
    global transcript_hash_table

    global adjacency_list

    adjacency_list = {}

    graph = nx.DiGraph(name=input_filename)
    transcript_hash_table = {}

    # Sets graph attributes, so each node has an empty hash table with transcript as key
    for tran in transcripts:
        graph.graph[tran] = ''
        transcript_hash_table[tran] = []
    graph.graph['Base'] = ''

    # Get the block co-ordinates from BLAT output
    output = psl_parse(input_filename)

    # Gets the sequence from fasta file. Also adds the nodes to the graph as sequences are retrieved
    get_sequence(output, "ESR1_transcripts_for_Michael.fasta")

    # Hash table containing the length of transcripts
    transcript_lengths = {}
    for out in output:
        transcript_lengths[out[0]] = out[1]
    for out in output:
        transcript_lengths[out[2]] = out[3]

    add_edge_entire_graph(transcript_lengths)

    graph.remove_node(None)

    # Print info on graph
    print(nx.info(graph))

    # fh = open("test.edgelist", 'wb')
    # print(nx.write_edgelist(graph, fh, data=True))

    # print(graph.summary())

    # print(">>> Combinging edges...")
    # graph.simplify()

    # print(graph.summary())

    # Simplify graph conactenating nodes that do not fork
    # simplify_graph(graph)

    # Writing out each node sequence for blatting
    # write_vertex_to_fasta(graph)

    # Reconstruct transcripts
    # reconstruct_transcript_sequence(graph, transcripts)

    # Sorts graph using khan algorithm
    # sorted = khan_sort(graph)

    # print(">>> Writing supertranscript to fasta file...")
    # sorted_sequence = Seq(sorted)
    # sorted_seq_rec = SeqRecord(sorted_sequence)
    # sorted_seq_rec.id = "SORT00000000000.1"
    # SeqIO.write(sorted_seq_rec, "ESR1_sorted_singlefunction.fasta", 'fasta')

    # Write graph information to text file for checking
    f = open('graph_output.txt', 'w')
    for node in graph.nodes(data=True):
        f.write(str(node))
        f.write("\n")
    f.write("\n")
    f.close()

    # print(">>> Drawing graph... ")
    # graph.write_svg("graph.svg", layout="kk")


def get_transcript_ids(fasta):

    # Method to retrieve list of transcript IDs from a fasta #

    transcript_ids = []
    filename = fasta

    for seq_record in SeqIO.parse(filename, "fasta"):
        transcript_ids.append(seq_record.id)

    return transcript_ids


# TODO - Method, so that when choosing a new no incoming node, checks that has a larger coordinate,
# - Use the genomic coordinates to help the sorting progress, to get the correct order


# Topologically sorts a graph using the Khan algorithm
def khan_sort(graph, transcripts):

    print(">> Sorting graph using khan algorithm....")

    sorted_sequence = ''

    L = []
    S = []

    print("Initial length of graph: ", len(graph.vs()))

    # Get initial list of nodes with no incoming edges
    print("get inital no incoming nodes")
    for vertex in graph.vs():
        if vertex.degree(mode="IN") == 0:
            L.append(vertex)

    print(len(L))

    # Sorting using the lists
    while L:

        for tran in transcripts:

            # For each transcript, store the gene co-ordinates if any to see how the transcripts will line up

            list_of_coordinates = []

            for l in L:
                list_of_coordinates.append(l.__getitem__(tran))

        node = L.pop(0)

        S.append(node)

        print(node.index)

        neighbours = graph.neighbors(node)

        for m in neighbours:
            print(node)
            print(graph.vs[m])

            graph.delete_edges(graph.get_eid(node.index, m))

            if graph.vs[m].degree(mode="IN") == 0:
                L.append(graph.vs()[m])

    for s in S:
        sorted_sequence += s.__getitem__("Base")

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

    print(">> Finding genomic sequence from from fasta file....")
    print(">> Opening fasta file...")

    # Hash table storing all the sequence record objects for the transcripts
    fasta_sequence = {}
    for block in blat_blocks:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            fasta_sequence[seq_record.id] = seq_record

    global graph

    print(">> Building nodes... ")

    output_records = []

    # TODO - Rewrite function to use the hash table to get sequence record objects

    for block in blat_blocks:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if blat_blocks[0][0] == seq_record.id:  # Gets query sequence from the fasta file

                first_base = True  # True if this is the first base of the transcript

                for i, z in enumerate(block[4]):

                    # print("\nPrinting i: ", i)

                    num_blocks = len(block[4])

                    block_size = int(block[4][i])
                    query_id = block[0]
                    query_block_start = int(block[5][i])
                    query_block_end = int(block[5][i])+block_size
                    query_length = int(block[1])
                    query = fasta_sequence[query_id]

                    target_id = str(block[2])
                    target_block_start = int(block[6][i])
                    target_block_end = int(block[6][i])+block_size
                    target_length = int(block[3])
                    target = fasta_sequence[target_id]

                    # TODO Add exception handling here if there is missing transcript from the fasta file

                    # print("Number blocks: ", num_blocks)
                    # print(query_id, query_block_start, query_block_end, query_length)
                    # print(target_id, target_block_start, target_block_end, target_length)

                    # grab the next block if there is one to find gap if there is one
                    query_gap_present = False
                    query_gap_size = 0

                    target_gap_present = False
                    target_gap_size = 0

                    # Determines if there are gaps or not in the transcripts
                    if num_blocks > 1:
                        try:
                            query_next_block_start = int(block[5][i+1])
                            if query_block_end != query_next_block_start:
                                query_gap_present = True
                                query_gap_size = query_next_block_start - query_block_end

                        except IndexError: # exception if there is no block
                            pass

                        try:
                            target_next_block_start = int(block[6][i+1])
                            if target_block_end != target_next_block_start:
                                target_gap_present = True
                                target_gap_size = target_next_block_start - target_block_end

                        except IndexError: # exception if there is no block
                            pass

                    # Ignores the BLATs between the same sequences
                    if query_id != target_id:

                        # ADDS UNALIGNED SECTION FOR TRANSCRIPTS 
                        if i == 0 and query_block_start != 0 and first_base is True:
                            for b1, base1 in enumerate(seq_record[0:query_block_start]):
                                add_node_graph(base1, [query_id], [b1])
                            # print("Added block: ", query_id, " ", 0, "-", query_block_start)
                            first_base = False

                        if i == 0 and target_block_start != 0 and first_base is True:
                            for b2, base2 in enumerate(seq_record[0:target_block_start]):
                                add_node_graph(base2, [target_id], [b2])
                            # print("Added block: ", target_id, " ", 0, "-", target_block_start)
                            first_base = False

                        # ADD BASSES FROM ALIGNED SECTION
                        for b3, base3 in enumerate(seq_record[query_block_start:query_block_end]):
                            add_node_graph(base3, [query_id, target_id], [query_block_start+b3, target_block_start+b3])
                            # print("Added block: ", query_id, " ", query_block_start, "-", query_block_end)
                            # print("Added block: ", target_id, " ", target_block_start, "-", target_block_end)

                        #TODO These don't seem to work, not sure why
                        # ADDS BLOCKS AFTER ALIGNMENT SECTIONG
                        if i == num_blocks-1:
                            if query_block_end != query_length:
                                for b, base in enumerate(seq_record[query_block_end:query_length]):
                                    add_node_graph(base, [query_id], [query_block_end+b])
                                # print("Added block: ", query_id, " ", query_block_end, "-", query_length)

                        if i == num_blocks - 1:
                            if target_block_end != target_length:
                                for b, base in enumerate(target[target_block_end:target_length]):
                                    add_node_graph(base, [target_id], [target_block_end+b])
                                # print("Added block: ", target_id, " ", target_block_end, "-", target_length)

                        # ADDS GAPS IF PRESENT
                        if query_gap_present is True:
                            query_gap_end = query_block_end+query_gap_size
                            for b4, base4 in enumerate(seq_record[query_block_end:query_gap_end]):
                                add_node_graph(base4, [query_id], [query_block_end+b4])
                            # print("Added block: ", query_id, " ", query_block_end, "-", query_gap_end)

                        if target_gap_present is True:
                            target_gap_end = target_block_end+target_gap_size
                            for b5, base5 in enumerate(seq_record[target_block_end:target_gap_end]):
                                add_node_graph(base5, [target_id], [target_block_end+b5])
                            # print("Added block: ", target_id, " ", target_block_start, "-", target_gap_end)

    print(">> Returning block sequence....\n")
    return output_records


# Writes list of block sequences to fasta file.
# To Do: need to rework how to manage naming of blocks as they currently share names from query sequence
def write_fasta(list_records, filename):
    SeqIO.write(list_records, filename, 'fasta')


# TODO VALIDATE
# Method for adding nodes to the graph; updates graph nodes
def add_node_graph(given_base, list_transcripts, coordinate):

    global graph
    global transcript_hash_table

    base = given_base.upper()
    graph_size = graph.number_of_nodes()

    # Checks if this is the first base being added to the graph
    # If this is the first base added to the graph add the node
    if graph_size == 0:

        # Adds empty node to end of graph - grabs teh recently added node
        graph.add_node(graph.number_of_nodes()+1, Base=base)

        newest_node = graph.node[graph.number_of_nodes()]

        for trans_hash in transcript_hash_table.keys():
            newest_node[trans_hash] = None

        for z, ba in enumerate(list_transcripts):
            newest_node[list_transcripts[z]] = coordinate[z]
            transcript_hash_table[list_transcripts[z]].append(coordinate[z])  # Adds coordinates to hash table

    # Creates list of booleans for evaluating if transcript co-ordinates exist in the graph
    base_exist = [False] * len(list_transcripts)

    for t, trans in enumerate(list_transcripts):
        list_coords = transcript_hash_table[trans]
        if coordinate[t] in list_coords:
            base_exist[t] = True

    # No transcripts and co-ordinates in there, so add new node - base_exist is still all false
    if True not in base_exist:

        # Adding new nodes to graph - grabs newest node
        graph.add_node(graph.number_of_nodes()+1, Base=base)
        newest_node = graph.node[graph.number_of_nodes()]

        for trans_hash2 in transcript_hash_table.keys():
            newest_node[trans_hash2] = None

        for z, ba in enumerate(list_transcripts):
            newest_node[list_transcripts[z]] = coordinate[z]
            transcript_hash_table[list_transcripts[z]].append(coordinate[z])  # Added coordinates to hash table


    # If one transcript returns true - do not add new nodes - Update the existing node with the new transcripts
    elif True and False in base_exist:
        true_indices = [i for i, x in enumerate(base_exist) if x is True]
        false_indices = [i for i, x in enumerate(base_exist) if x is False]

        # This seems to work ok, will need to test on larger data set

        length_graph = graph.__len__()

        # TODO - Double check where graph is starting, if needs to start at 0 or 1
        x = 1

        while x != length_graph:
            # print(graph.node[x])

            for tr in true_indices:
                if graph.node[x][list_transcripts[tr]] == coordinate[tr]:
                    for j in false_indices:
                        graph[x][list_transcripts[t]] = coordinate[t]
                        transcript_hash_table[list_transcripts[j]].append(coordinate[j])
            break
            x += 1


# TODO - VALIDATE
# TODO - NOT WORKING CORRECTLY
def add_edge_entire_graph(dict_lengths):

    print(">> Adding edges...")

    global graph

    # Creates a dictionary to contain the index information for each transcript
    # Each value for the dictionary will be a list of vertex indexes from the graph, the index in the list will be the
    # coordinate of the base and the value stored will be the index of the graph it is stored in
    edge_index = {}

    for tran in dict_lengths:
        length = int(dict_lengths[tran])
        edge_index[tran] = [None] * length

    # Printing insanity checking here - Creating correct size lists
    # print("Insanity checking")
    # for edge in edge_index:
    #     print(edge, len(edge_index[edge]))

    # Visits each node
    for nod in graph.nodes():

        # Get transcript base locations for each transcript
        for tran in edge_index:

            # Check if the base is in transcripts
            if graph.node[nod][tran] is not None:

                # Add node ID to dictionary at index of genomic coordinate
                # print("Coordinate : ", graph.node[nod][tran], " in ", tran, " located in node : ", nod)
                edge_index[tran][graph.node[nod][tran]] = nod

    for edge in edge_index:
        index_list = edge_index[edge]
        i = len(index_list)
        # i = 0

        # while i < len(index_list)-1:
        #     graph.add_edge(index_list[i], index_list[i+1])
        #     i += 1

        while i > 1:
            graph.add_edge(index_list[i-2], index_list[i-1])
            i -= 1



# TODO CONVERT TO NETWORKX
# TODO - Need to optimize this method, it is very slow - not sure how
# Simplifies the graph by compacting nodes that have only one in degree and one out degree
def simplify_graph(graph):

    print(">>> Simplifying graph...")
    i = 0
    while i < len(graph.vs()):

        vertex = graph.vs()[i]

        nod = graph.node[i]

        # Check if current node has only a single outgoing edge
        if vertex.outdegree() == 1:
            neighbor_id = graph.neighbors(vertex, mode="out")
            neighbor = graph.vs()[neighbor_id[0]][0]

        if graph.out_degree(nod) == 1:

            # get neighbour id
            # get neighbour node

            neighbor = graph.neighbors(nod)

            # Checks if the adjacent outgoing edge of current node has only one incoming edge
            if neighbor.indegree() == 1:  # there should on be one value in the list

                # Gets base values for the current and neighbour nodes and combines them
                # v_base = vertex["Base"]
                # n_base = neighbor["Base"]
                # new_base = v_base + n_base
                # vertex["Base"] = new_base

                vertex["Base"] = vertex["Base"] + neighbor["Base"]

                # Edges are added from vertex to its new neighbours
                for index in graph.neighbors(neighbor_id[0], mode="out"):
                    graph.add_edge(vertex, graph.vs()[index])

                # TODO - Double check what deleting nodes does to its edges,
                # edges may still exist, should delete them too
                # Neighbour vertex is deleted from the graph
                graph.delete_vertices(neighbor_id)

            else:
                i += 1  # increases iteration if there is more then one incoming edge on neighbour node
        else:
            i += 1  # increases iteration if there is more then one outgoing edge on node


# Writes out each node to a fasta file format, with each node being a sequence in the file
# Used for checking sequences of individual nodes pre-sort or when ever
def write_vertex_to_fasta(graph):

    print(">>> Printing out vertex sequences into fasta file... ")

    records = []

    for vert in graph.vs():
        print(vert.__getitem__("Base"))

    i = 1
    for vertex in graph.vs():
        vert_seq = Seq(vertex.__getitem__("Base"))
        vert_seq_rec = SeqRecord(vert_seq)
        vert_seq_rec.id = ("ENST00000000000." + str(i))
        records.append(vert_seq_rec)
        i+=1

    SeqIO.write(records, "vertex_sequences.fasta", 'fasta')


# Method retrieves the sequence for each transcript
def reconstruct_transcript_sequence(graph, transcript_list):

    transcripts = ["ENST00000456483.2", "ENST00000338799.5"]

    print(">>> Reconstructing transcripts from graph... ")

    records = []

    for x in transcript_list:
        # print("X : ", x)

        seq_list = []
        list_coords = []

        for vertex in graph.vs():

            if vertex.__getitem__(x) is not None:

                ind = bisect.bisect(list_coords, vertex.__getitem__(x))
                seq_list.insert(ind, vertex.__getitem__("Base"))

                bisect.insort(list_coords, vertex.__getitem__(x))


        sequence = ""
        for s in seq_list:
            sequence += str(s)

        seq = Seq(sequence)
        seq_rec = SeqRecord(seq)
        seq_rec.id = (x)

        records.append((seq_rec))

    SeqIO.write(records, "ESR1_reconstructed_transcripts.fasta", "fasta")


if __name__ == "__main__":
    start_time = time.clock()
    main()
    print("Running time: ", time.clock()-start_time, "s")








