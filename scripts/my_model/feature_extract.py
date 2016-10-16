import numpy as np
from numpy import mean, median, std, nansum
import os, sys
from string import maketrans
import re

from cpmodule import fickett
from cpmodule  import orf
from cpmodule  import fasta
# from cpmodule  import annoGene
from cpmodule  import FrameKmer
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-g","--gene",action="store",dest="gene_file",help="Transcripts either in BED format or mRNA sequences in FASTA format: If this is BED format file, '-r' must be specified; if this is mRNA sequence file in FASTA format, ignore the '-r' option. The input BED or FASTA file could be regular text file or compressed file (*.gz, *.bz2) or accessible url.")
parser.add_option("-o","--outfile",action="store",dest="out_file",help="output file. Tab separated text file: geneID <tab> mRNA size <tab> ORF size <tab> Fickett Score <tab> Hexamer Score<tab>Coding Probability.")
parser.add_option("-x","--hex",action="store",dest="hexamer_dat",help="Prebuilt hexamer frequency table (Human, Mouse, Fly, Zebrafish). Run 'make_hexamer_tab.py' to make this table out of your own training dataset.")
parser.add_option("-s","--start",action="store",dest="start_codons",default='ATG',help="Start codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). default=%default")
parser.add_option("-t","--stop",action="store",dest="stop_codons",default='TAG,TAA,TGA',help="Stop codon (DNA sequence, so use 'T' instead of 'U') used to define open reading frame (ORF). Multiple stop codons should be separated by ','. default=%default")

(options,args)=parser.parse_args()

#build hexamer table from hexamer frequency file

def extract_feature_from_seq(seq,stt,stp,c_tab,g_tab):
	'''extract features of sequence from fasta entry'''
	
	stt_coden = stt.strip().split(',')
	stp_coden = stp.strip().split(',')
	transtab = maketrans("ACGTNX","TGCANX")
	mRNA_seq = seq.upper()
	mRNA_size = len(seq)
	tmp = orf.ORFFinder(mRNA_seq)
	(CDS_size1, CDS_frame1, CDS_seq1) = tmp.longest_orf(direction="+",start_coden=stt_coden, stop_coden=stp_coden)
	fickett_score1 = fickett.fickett_value(CDS_seq1)
	hexamer = FrameKmer.kmer_ratio(CDS_seq1,6,3,c_tab,g_tab)
	return (mRNA_size, CDS_size1, fickett_score1,hexamer)



def count_gc(seq):
        
        seq = seq.upper()
        dict_nucl = {"A":0, "T":0, "G":1, "C":1}
        GC = np.sum([dict_nucl.get(c,0) for c in seq])/(len(seq)*1.0)
        return GC


def kozak_feat(seq, find_num):
    
    # first_feature {-3,+4}
    # second -- {-2, -1}
    # third -- {-6}
    # method returns 3 boolean values, 
    # that mean find or not find for every feature
    
    
    if find_num < 6:
        return -1
    
    else:
        feat34 = 0
        feat21 = 0
        feat6 = 0

        num6 = seq[find_num-6]
        if len(seq) >= find_num+4:
            num4 = seq[find_num+3]
        else:
            num4 = -1

        num3 = seq[find_num-3]
        num2 = seq[find_num-2]
        num1 = seq[find_num-1]
        
    if num6 == "G":
        feat6 = 1
    if num4 == "G" and (num3 == "A" or num3 == "G"):
        feat34 = 1
    if num2 == "C" and num1 == "C":
        feat21 = 1
    
    return feat34, feat21, feat6

def find_kozak_feat(seq):
    
    atg_addr = [m.start() for m in re.finditer("ATG", seq)]
    atg_feat = []
    feat_all = False
    
    if len(atg_addr) <1:
        return (0,0,0)
    
    for atg in atg_addr:
        feat_all = kozak_feat(seq, atg)
        atg_feat.append(feat_all)
        if feat_all == -1:
            return (0,0,0)
        if feat_all == (1, 1, 1):
            feat_all = 1
            break
   
    if feat_all==1:
        return (1, 1, 1)
    else:
        return atg_feat[np.argmax([sum(i) for i in atg_feat])]


coding={}
noncoding={}	
for line in open(options.hexamer_dat):
	line = line.strip()
	fields = line.split()
	if fields[0] == 'hexamer':continue
	coding[fields[0]] = float(fields[1])
	noncoding[fields[0]] =  float(fields[2])

count = 0
TMP = open(options.out_file + '.txt', 'w')
TMP.write('\t'.join(("sname", "mRNA_size", "ORF_size", "fickett_score", "hexamer", "gc_content", "kozak34", "kozak21", "kozak6"))+'\n')
for sname,seq in FrameKmer.seq_generator(options.gene_file):	
	count +=1
	print(sname)
	gc_content = count_gc(seq)
	k34, k21, k6 = find_kozak_feat(seq)

	(mRNA_size, CDS_size, fickett_score,hexamer) = extract_feature_from_seq(seq = seq, stt = options.start_codons,stp = options.stop_codons,c_tab=coding,g_tab=noncoding)
	TMP.write('\t'.join(str(i) for i in (sname, mRNA_size, CDS_size, fickett_score, hexamer, gc_content, k34, k21, k6))+'\n')
	if count % 30 == 0:
		print sys.stderr, "%d genes finished\r" % count,
TMP.close()


