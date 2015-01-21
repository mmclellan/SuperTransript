import re
import sys
import getopt
from Bio import SeqIO

def main(argv):

    inputfile = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile", "ofile="])
    except getopt.GetoptError:
        print('test. py -i <butt> <butt>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print 'test.py -i <inputfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
    print 'Input file is "', inputfile
    print psl_parse(inputfile)


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

            print blocks

            blat_blocks.append(blocks)

    return blat_blocks


# Returns the sequence data for the blocks from the blat data
def get_sequence(blat_blocks, fasta_file):

    # for each list of co-ordinates from the blat_blocks list
    # create a

    # Reads in the fasta file
    record = SeqIO.read('ESR1.fasta', 'fasta')
    print record


if __name__ == "__main__":
    main(sys.argv[1:])