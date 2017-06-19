import numpy as np
# from numpy import mean, median, std, nansum
import sys
# from string import maketrans
# import re

from cpmodule import fickett
from cpmodule import orf_extraction
# from cpmodule import fasta
from cpmodule import FrameKmer
from cpmodule import kozak
from cpmodule import gtf_exons
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-g", "--gene", action="store", dest="gene_file",
                  help="Transcripts either in BED format or mRNA sequences in FASTA format: If this is BED format file, '-r' must be specified; if this is mRNA sequence file in FASTA format, ignore the '-r' option. The input BED or FASTA file could be regular text file or compressed file (*.gz, *.bz2) or accessible url.")
parser.add_option("-o", "--outfile", action="store", dest="out_file",
                  help="output file. Tab separated text file: geneID <tab> mRNA size <tab> ORF size <tab> Fickett Score <tab> Hexamer Score<tab>Coding Probability.")
parser.add_option("-x", "--hex", action="store", dest="hexamer_dat",
                  help="Prebuilt hexamer frequency table (Human, Mouse, Fly, Zebrafish). Run 'make_hexamer_tab.py' to make this table out of your own training dataset.")
parser.add_option("-s", "--start", action="store", dest="start_codons", default='ATG',
                  help="Start codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). default=%default")
parser.add_option("-t", "--stop", action="store", dest="stop_codons", default='TAG,TAA,TGA',
                  help="Stop codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). Multiple stop codons should be separated by ','. default=%default")
parser.add_option("--gtf", action="store", dest="gtf",
                  help="Please, provide gtf/gff file for extraction exon-related features.")
parser.add_option("--fformat", action="store", dest="fformat", default='gtf',
                  help="Please, provide file format: -ff gtf or -ff gff. default gtf.")
parser.add_option("-l", "--lines-drop", action="store", dest="lines_drop", default=1,
                  help="Provide how many lines drop from the start of the gtf/gff file. Default = 1")

(options, args) = parser.parse_args()


# build hexamer table from hexamer frequency file

def count_gc(seq):
    seq = seq.upper()
    dict_nucl = {"A": 0, "T": 0, "G": 1, "C": 1}
    GC = np.sum([dict_nucl.get(c, 0) for c in seq]) / (len(seq) * 1.0)
    return GC


def extract_feature_from_seq(seq, c_tab, g_tab):
    '''extract features of sequence from fasta entry'''

    mRNA_seq = seq.upper()
    mRNA_size = len(seq)

    orf_finder = orf_extraction.ORFFinder(mRNA_seq)
    tmp = orf_finder.find_longest()

    ''' in the case if start codon have not been found '''
    if tmp==-1:
        return [0] * 9

    starts, orf_seq, orf_size, mean_orf_length, orf_coverage = tmp

    fickett_score = fickett.fickett_value(orf_seq)

    k34, k21, k6 = kozak.find_kozak_feat(mRNA_seq, starts)
    hexamer = FrameKmer.kmer_ratio(orf_seq, 6, 3, c_tab, g_tab)

    return (mRNA_size, orf_size, mean_orf_length, orf_coverage, fickett_score, hexamer, k34, k21, k6)


coding = {}
noncoding = {}
for line in open(options.hexamer_dat):
    line = line.strip()
    fields = line.split()
    if fields[0] == 'hexamer': continue
    coding[fields[0]] = float(fields[1])
    noncoding[fields[0]] = float(fields[2])

exon_max, exon_mean, exon_num = gtf_exons.gtf_parser(options.gtf, int(options.lines_drop), options.fformat)

TMP = open(options.out_file + '.txt', 'w')
TMP.write('\t'.join(("sname", "mRNA_size", "ORF_size", "mean_orf_length", "orf_coverage", "fickett_score", "hexamer", "gc_content", "kozak34", "kozak21",
                     "kozak6", "exon_max", "exon_mean", "exon_num")) + '\n')

count = 0

for sname, seq in FrameKmer.seq_generator(options.gene_file):
    
    gc_content = count_gc(seq)
    count+=1

    mRNA_size, orf_size, mean_orf_length, orf_coverage, \
    fickett_score, hexamer, k34, k21, k6  = extract_feature_from_seq(seq=seq,
                                                                   c_tab=coding,
                                                                   g_tab=noncoding)
    TMP.write('\t'.join(str(i) for i in (
        sname, mRNA_size, orf_size, mean_orf_length, orf_coverage, fickett_score, hexamer, gc_content, k34, k21, k6, exon_max.get(sname, 0),
        exon_mean.get(sname, 0), exon_num.get(sname, 0))) + '\n')

    if count % 100 == 0:
        print(sys.stderr, "%d genes finished\r" % count,)

TMP.close()
