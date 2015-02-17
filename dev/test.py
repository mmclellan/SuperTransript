from igraph import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import  SeqRecord
import re
import time

from collections import deque


# Main function for working
def main():

    # fasta file for corresponding genes
    fasta_filename = "ESR1_transcripts_for_Michael.fasta"

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
    get_sequence(output, "ESR1_transcripts_for_Michael.fasta")

    # Hash table containing the length of transcripts
    transcript_lengths = {}
    for out in output:
        transcript_lengths[out[0]] = out[1]
    for out in output:
        transcript_lengths[out[2]] = out[3]

    add_edge_entire_graph(transcript_lengths)

    print(graph.summary())

    print(">>> Combinging edges...")
    graph.simplify()

    print(graph.summary())

    # Simplify graph conactenating nodes that do not fork
    simplify_graph(graph)

    # print(">>> Drawing graph... ")
    # graph.write_svg("graph.svg", layout="kk")

    # Graph summary
    print(graph.summary())

    # Writing out each node sequence for blatting
    # write_vertex_to_fasta(graph)

    # Sorts graph using khan algorithm
    sorted = khan_sort(graph)

    print(">>> Writing supertranscript to fasta file...")
    sorted_sequence = Seq(sorted)
    sorted_seq_rec = SeqRecord(sorted_sequence)
    sorted_seq_rec.id = "SORT00000000000.1"
    SeqIO.write(sorted_seq_rec, "ESR1_sorted.fasta", 'fasta')


    '''
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

    f.write("sorted sequence \n")
    f.write(sorted_sequence)
    f.close()
    '''

# Method to retrieve list of transcript IDs from a fasta
def get_transcript_ids(fasta):

    transcript_ids = []
    filename = fasta

    for seq_record in SeqIO.parse(filename, "fasta"):
        transcript_ids.append(seq_record.id)

    return transcript_ids


#TODO BUG - Think there is a bug here adding multiple nodes to the output
# Topologically sorts a graph using the Khan algorithm
# seems to work without worrying about deleting edges, looks to take care of its self
def khan_sort(graph):

    print(">> Sorting graph using khan algorithm....")

    # creates list of nodes in the graph and number of nodes
    node_list = graph.vs()

    sorted_sequence = ''

    L = []
    S = []

    i = 0
    graph_length = len(graph.vs())
    num_nodes = graph_length

    # Add initial nodes with no incoming edges to s for sorting

    print("Initial length of graph: ", len(graph.vs()))

    while len(graph.vs()) > 0:

        if len(graph.vs()) == 1:
            L.append(graph.vs()[0].__getitem__("Base"))
            break

        print("The length of graph: ", len(graph.vs))

        for x in graph.vs():
                if x.degree(mode="IN") == 0:
                    print("Vertex \n", x, " has no incoming edges")
                    S.append(x)

                    for s in S:

                        # need to get out going neighbour from s
                        neighbours = graph.neighbors(s, mode="OUT")

                        for j in neighbours:
                            graph.delete_edges(graph.get_eid(s.index, j))
                            print("Delete edge: ", s.index, j)

                        L.append(s.__getitem__("Base"))
                        S.remove(s)

                        print("Delete vertex: ", s.index)
                        graph.delete_vertices(s.index)

                    break  # break out of for loop and start through the graph from scratch

                else:
                    print("Cycle must be present")

    print(len(L))
    print(L)

    for l in L:
        sorted_sequence += l

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

    # Hash table storing all the sequence record objects for the transcripts
    fasta_sequence = {}
    for block in blat_blocks:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            fasta_sequence[seq_record.id] = seq_record

    global graph

    print(">> Finding genomic sequence from from fasta file....")
    print(">> Opening fasta file...")

    output_records = []

    # TODO - Rewrite function to use the hash table to get sequence record objects

    for block in blat_blocks:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if blat_blocks[0][0] == seq_record.id:  # Gets query sequence from the fasta file

                first_base = True  # True if this is the first base of the transcript

                for i, z in enumerate(block[4]):

                    print("\nPrinting i: ", i)

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

                    # TODO Add exception handling here if there is missing transcript from the fasta file
                    target = fasta_sequence[target_id]

                    # print("Number blocks: ", num_blocks)
                    # print(query_id, query_block_start, query_block_end, query_length)
                    # print(target_id, target_block_start, target_block_end, target_length)

                    # grab the next block if there is one to find gap if there is one
                    query_gap_present = False
                    query_gap_size = 0
                    target_gap_present = False
                    target_gap_size = 0

                    if num_blocks > 1:
                        try:
                            query_next_block_start = int(block[5][i+1])
                            if query_block_end != query_next_block_start:
                                query_gap_present = True
                                query_gap_size = query_next_block_start - query_block_end

                        except IndexError:
                            print("No next query block")

                        try:
                            target_next_block_start = int(block[6][i+1])
                            if target_block_end != target_next_block_start:
                                target_gap_present = True
                                target_gap_size = target_next_block_start - target_block_end

                        except IndexError:
                            print("No next target block")

                    # while query ID and target ID do not match
                    if query_id != target_id:

                        # adds unaligned section before aligned blocks
                        if i == 0 and query_block_start != 0 and first_base is True:

                            for b1, base1 in enumerate(seq_record[0:query_block_start]):
                                add_node_graph(base1, [query_id], [b1])
                            print("Added block: ", query_id, " ", 0, "-", query_block_start)
                            first_base = False

                        if i == 0 and target_block_start != 0 and first_base is True:

                            for b2, base2 in enumerate(seq_record[0:target_block_start]):
                                add_node_graph(base2, [target_id], [b2])
                            print("Added block: ", target_id, " ", 0, "-", target_block_start)
                            first_base = False

                        # Add gaps found in target or query sequence
                        if query_gap_present is True:
                            query_gap_end = query_block_end+query_gap_size
                            for b4, base4 in enumerate(seq_record[query_block_end:query_gap_end]):
                                add_node_graph(base4, [query_id], [query_block_end+b4])
                            print("Added block: ", query_id, " ", query_block_end, "-", query_gap_end)

                        if target_gap_present is True:
                            target_gap_end = target_block_end+target_gap_size
                            for b5, base5 in enumerate(seq_record[target_block_end:target_gap_end]):
                                add_node_graph(base5, [target_id], [target_block_end+b5])
                            print("Added block: ", target_id, " ", target_block_start, "-", target_gap_end)


                        #TODO These don't seem to work, not sure why
                        # Adds bases after blocks to end of query or target sequence
                        if i == num_blocks-1:
                            if query_block_end != query_length:
                                for b, base in enumerate(seq_record[query_block_end:query_length]):
                                    add_node_graph(base, [query_id], [query_block_end+b])
                                print("Added block: ", query_id, " ", query_block_end, "-", query_length)

                        if i == num_blocks - 1:
                            if target_block_end != target_length:
                                for b, base in enumerate(target[target_block_end:target_length]):
                                    add_node_graph(base, [target_id], [target_block_end+b])
                                print("Added block: ", target_id, " ", target_block_end, "-", target_length)



                            # add basses from qstart and tstart to the length of block size
                        for b3, base3 in enumerate(seq_record[query_block_start:query_block_end]):
                            add_node_graph(base3, [query_id, target_id], [query_block_start+b3, target_block_start+b3])
                        print("Added block: ", query_id, " ", query_block_start, "-", query_block_end)
                        print("Added block: ", target_id, " ", target_block_start, "-", target_block_end)

    print(">> Returning block sequence....\n")
    return output_records


# Writes list of block sequences to fasta file.
# To Do: need to rework how to manage naming of blocks as they currently share names from query sequence
def write_fasta(list_records, filename):
    SeqIO.write(list_records, filename, 'fasta')


# Method for adding nodes to the graph; updates graph nodes
def add_node_graph(given_base, list_transcripts, coordinate):

    base = given_base.upper()

    global graph
    global transcript_hash_table

    # print("adding node.. ")


    # Checks if this is the first base being added to the graph
    # If this is the first base added to the graph add the node
    if len(graph.vs()) == 0:
        graph.add_vertex()
        graph.vs[len(graph.vs())-1]["Base"] = base

        # print("Creating first node for the graph... ")

        for z, ba in enumerate(list_transcripts):
            graph.vs[len(graph.vs())-1][list_transcripts[z]] = coordinate[z]
            # Added transcript coordinates to hash table
            transcript_hash_table[list_transcripts[z]].append(coordinate[z])

    # Creates list of booleans for evaluating if transcript co-ordinates exist in the graph
    base_exist = []
    [base_exist.append(False) for t in list_transcripts]

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

        for z, ba in enumerate(list_transcripts):
            graph.vs[len(graph.vs())-1][list_transcripts[z]] = coordinate[z]

            # print("Node does not exist -> Adding ", list_transcripts[z], "at coord", coordinate[z])

            # Added transcript coordinates to hash table
            transcript_hash_table[list_transcripts[z]].append(coordinate[z])

    # - If one of the transcripts returns true then do not add new nodes
    # - Update the existing node with the new transcripts
    elif True and False in base_exist:
        true_indices = [i for i, x in enumerate(base_exist) if x is True]
        false_indices = [i for i, x in enumerate(base_exist) if x is False]

        # TODO this is quite slow this method here
        # Updates true index transcripts at nodes with false index coordinates
        for vertex in graph.vs.select(Base=base):
            for t in true_indices:
                if vertex[list_transcripts[t]] == coordinate[t]:
                    for j in false_indices:
                        vertex[list_transcripts[j]] = coordinate[j]
                        transcript_hash_table[list_transcripts[j]].append(coordinate[j])
                        # print("Update node: ", list_transcripts[t], "-", coordinate[t], " with ", list_transcripts[j], "-", coordinate[j])

    # elif False not in base_exist:
    #     print("Node already exists: not adding node")

# TODO - Delete this method as it is no longer in use
def add_edge(transcript, coordinate):

    global graph

    # pass through transcript length and transcript stopping coordinate
    # Use counter i or something and set it to the length of
    # While i is greater then 0
        # count back from i to 0
        # add edge from coordinate i to i-1, then i -= 1


    # search graph and add ids to list
        # iterate through the list of vertex indexes and add nodes between them?

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
            print("Did not find node")


# TODO - Rework method, it is adding redundant nodes that are already present in graph
def add_edge_entire_graph(dict_lengths):

    print(">>> Adding edges...")

    global graph

    # Creates a dictionary to contain the index information for each transcript
    # Each value for the dictionary will be a list of vertex indexes from the graph, the index in the list will be the
    # coordinate of the base and the value stored will be the index of the graph it is stored in
    edge_index = {}

    for tran in dict_lengths:
        length = int(dict_lengths[tran])
        edge_index[tran] = [None] * length

    # Printing insanity checking here - looks like it is making the right size empty lists
    # print("Insanity checking")
    # for edge in edge_index:
    #     print(edge, len(edge_index[edge]))

    for vertex in graph.vs():
        for tran in edge_index:
            if vertex.__getitem__(tran) is not None:
                edge_index[tran][vertex.__getitem__(tran)] = vertex.index

    # for e in edge_index:
    #     print(e, "\n")
    #     print(edge_index[e])

    for edge in edge_index:
        index_list = edge_index[edge]
        i = len(index_list)

        while i > 1:
            graph.add_edge(index_list[i-2], index_list[i-1])
            # print("edge from: ", index_list[i-2], " -> ", index_list[i-1], " in ", edge)
            i -= 1



#TODO - Need to optimize this method, it is very slow
# Simplifies the graph by compacting nodes that have only one in degree and one out degree
def simplify_graph(graph):

    print(">>> Simplifying graph...")

    counter = 0

    i = 0

    while i != len(graph.vs()):

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

            else:
                i += 1
        else:
            i += 1


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


if __name__ == "__main__":
    start_time = time.clock()
    main()
    print("Running time: ", time.clock()-start_time, "s")







