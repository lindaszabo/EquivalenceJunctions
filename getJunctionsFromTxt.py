# Determines equivalent junction sequences from txt file containing chr, donor, and acceptor coordinates.
# Handling coordinates and strand is based on how data were presented in circBase txt downloads, so
# may need to be modified for other data sources.

import sys
import argparse
from Bio import SeqIO
from BCBio import GFF

# usage: python getJunctionsFromTxt.py -f /scratch/PI/horence/linda/EquivJuncPaper/hg19_genome.fa
#                                           -t /scratch/PI/horence/linda/EquivJuncPaper/hsa_hg19_circRNA.txt
#                                           -o /scratch/PI/horence/linda/EquivJuncPaper/hsa_hg19_circRNA_EquivSeqs.txt
#                                           -v

# output: tab delimited file with equivJuncSeq, chr, donor, acceptor 

# TO GENERATE REPORT FROM THIS OUTPUT: cat hsa_hg19_circRNA_EquivSeqs.txt Sim1deletionEquivSeqs.txt | sort -k1 | uniq | cut -f1 | sort | uniq -c | sort -k1nr
# this removes double-counting cases if the same junction is listed multiple times.

CHECK_EQUIV_LEN = 30

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def getRevComp(inStr):
    revStr = inStr[::-1]
    revComp = []
    for a in revStr:
        revComp.append(complement[a])
    return "".join(revComp)


def getEquivSeqs(txtFile, chrCol, donorCol, acceptorCol, strandCol):
    t_handle = open(txtFile, "rU")
    for line in t_handle:
        try:
            splits = line.strip().split()
            juncChr = splits[chrCol]
            juncStrand = splits[strandCol]
            if juncStrand == "+":
                donorLoc = int(splits[acceptorCol])
                acceptorLoc = int(splits[donorCol])
            else:
                donorLoc = int(splits[donorCol])
                acceptorLoc = int(splits[acceptorCol])
            if args.verbose:
                print juncChr, juncStrand, donorLoc, acceptorLoc

            if juncChr in f_dict:
                if juncStrand == "+":
                    dseq_prev = str(f_dict[juncChr][donorLoc-CHECK_EQUIV_LEN:donorLoc].seq).upper()
                    dseq_next = str(f_dict[juncChr][donorLoc:donorLoc+CHECK_EQUIV_LEN].seq).upper()
                    
                    aseq_prev = str(f_dict[juncChr][acceptorLoc-CHECK_EQUIV_LEN:acceptorLoc].seq).upper()
                    aseq_next = str(f_dict[juncChr][acceptorLoc:acceptorLoc+CHECK_EQUIV_LEN].seq).upper()
                else:
                    dseq = str(f_dict[juncChr][donorLoc-CHECK_EQUIV_LEN:donorLoc+CHECK_EQUIV_LEN].seq).upper()
                    dseq = getRevComp(dseq)
                    dseq_prev = dseq[:CHECK_EQUIV_LEN]
                    dseq_next = dseq[CHECK_EQUIV_LEN:]
                    
                    aseq = str(f_dict[juncChr][acceptorLoc-CHECK_EQUIV_LEN:acceptorLoc+CHECK_EQUIV_LEN].seq).upper()
                    aseq = getRevComp(aseq)
                    aseq_prev = aseq[:CHECK_EQUIV_LEN]
                    aseq_next = aseq[CHECK_EQUIV_LEN:]
                
                if args.verbose:    
                    print "dseq_prev:", dseq_prev
                    print "dseq_next:", dseq_next
                    print "aseq_prev:", aseq_prev
                    print "aseq_next:", aseq_next
                    
                    
                # generate equiv junc seq
                equivJuncSeq = ""
                
                # check moving breakpoint more upstream in gene
                e_ind = CHECK_EQUIV_LEN - 1
                while e_ind != -1:
                    if dseq_prev[e_ind] != aseq_prev[e_ind]:
                        break
                    equivJuncSeq = equivJuncSeq + dseq_prev[e_ind]
                    e_ind = e_ind - 1
                
                equivJuncSeq = equivJuncSeq[::-1]  # since it was added from end 1 by 1, need to reverse to get coorrect order
                
                if args.verbose:    
                    print "equiv from moving upstream:", equivJuncSeq
                    
                # check moving breakpoint more downstream in gene
                e_ind = 0
                while e_ind < CHECK_EQUIV_LEN:
                    if dseq_next[e_ind] != aseq_next[e_ind]:
                        break
                    equivJuncSeq = equivJuncSeq + dseq_next[e_ind]
                    e_ind = e_ind + 1
                    
                if args.verbose:    
                    print "equiv from checking both directions:", equivJuncSeq
                
                o_handle.write("\t".join([equivJuncSeq, juncChr, str(donorLoc), str(acceptorLoc)])) 
                o_handle.write("\n")
        except Exception as e:
            print "Exception"
            print e
            print "error:", sys.exc_info()[0]
            print "parsing line", line
    t_handle.close()
            
if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastaFile', required=True, help='path to fasta file with chromosome sequences')
    parser.add_argument('-t', '--txtFile', required=True, help='txt file containing splice junction coordinate data')
    parser.add_argument('-o', '--outFile', required=True, help='output file name (including full path)')
    parser.add_argument('-c', '--chrCol', type=int, default=0, help='0-based column containing chromosome info')
    parser.add_argument('-d', '--donorCol', type=int, default=1, help='0-based column containing donor coordinate')
    parser.add_argument('-a', '--acceptorCol', type=int, default=2, help='0-based column containing acceptor coordinate')
    parser.add_argument('-s', '--strandCol', type=int, default=3, help='0-based column containing strand (+/-)')
    parser.add_argument('-v', '--verbose', help='print info about data obtained', action='store_true')
    args = parser.parse_args()
    
    # read in the sequences 
    f_handle = open(args.fastaFile, "rU")
    f_dict = SeqIO.to_dict(SeqIO.parse(f_handle, "fasta"))
    f_handle.close()
    
    o_handle = open(args.outFile, "wb")
    getEquivSeqs(args.txtFile, args.chrCol, args.donorCol, args.acceptorCol, args.strandCol)
    o_handle.close()
