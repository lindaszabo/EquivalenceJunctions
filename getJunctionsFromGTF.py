# Determines equivalent junction sequences at exon-exon boundaries in a transcriptome.
# Sequence of genes on - strand are what would be observed in the transcript.

import sys
import argparse
from Bio import SeqIO
from BCBio import GFF

# usage: python getJunctionsFromGTF.py -f /scratch/PI/horence/linda/EquivJuncPaper/GRCh38.primary_assembly.genome.fa
#                                      -a /scratch/PI/horence/linda/EquivJuncPaper/gencode.v25.annotation.gtf
#                                      -o /scratch/PI/horence/linda/EquivJuncPaper/GRCh38_gencodev25.txt
#                                      -v

# output: tab delimited file with equivJuncSeq, gene name, chr, donor, acceptor  

# TO GENERATE REPORT FROM THIS OUTPUT:

# this removes double-counting cases where the same exon boundaries are used in 2 different transcripts
# and counts number of junctions with a given equivalent junction sequence
# sort -k1 GRCh38_gencodev25.txt | uniq | cut -f1 | sort | uniq -c | sort -k1nr

# limit report to unique junctions (don't double-count if in multiple transcripts)
# sort -k2 GRCm38_gencodevM12.txt | uniq > GRCm38_gencodevM12_uniqueExonsWithEquivJuncs.txt



CHECK_EQUIV_LEN = 30 # number of nt to check on either side of exon boundary for equivalent junction sequences

# tries getting primary gene name, if that doesn't exist tries getting secondary gene name.
# if neither primary or secondary gene name fields exist, an exception will be thrown
# which is caught and handled in the calling function
def getGeneName(use_feature):
    if args.name1 in use_feature.qualifiers:
        use_name = args.name1
    else: 
        use_name = args.name2
        
    return use_feature.qualifiers[use_name][0]
    
# Single-exon genes will be skipped since they do not have sub_features
# for + strand genes location.start is 1 less than is in the gtf file, location.end matches gtf
# param chrSeqRecord: a SeqRecord for 1 chromosome
def parseChrFeatures(chrSeqRecord):
    
    for ftr in chrSeqRecord.features: # this is 1 entire gene
        curDonors = []
        curAcceptors = []
        
        try:
            sftr = None # genes with a single exon are at the top level ftr
            for sftr in ftr.sub_features:
                curStrand = sftr.strand

                if curStrand == 1: # + strand
                    curDonors.append(int(sftr.location.end))  # will have 1 extra for last exon that will just get ignored
                    curAcceptors.append(int(sftr.location.start))  # will have 1 extra at the beginning for start exon that will just get ignored
                else: # -strand
                    curDonors.append(int(sftr.location.start))
                    curAcceptors.append(int(sftr.location.end))
            
            if args.verbose:
                if not sftr:    
                    print "single exon gene:", getGeneName(ftr)
                    
                print "exon donors:", curDonors
                print "exon acceptors:", curAcceptors
                
            for b_ind in xrange(len(curAcceptors)-1):
                if args.verbose:
                    print curDonors[b_ind], curAcceptors[b_ind+1]
                if curStrand == 1: # + strand
                    dseq_prev = str(f_dict[chrSeqRecord.id][curDonors[b_ind]-CHECK_EQUIV_LEN:curDonors[b_ind]].seq).upper()
                    dseq_next = str(f_dict[chrSeqRecord.id][curDonors[b_ind]:curDonors[b_ind]+CHECK_EQUIV_LEN].seq).upper()
                    
                    aseq_prev = str(f_dict[chrSeqRecord.id][curAcceptors[b_ind+1]-CHECK_EQUIV_LEN:curAcceptors[b_ind+1]].seq).upper()
                    aseq_next = str(f_dict[chrSeqRecord.id][curAcceptors[b_ind+1]:curAcceptors[b_ind+1]+CHECK_EQUIV_LEN].seq).upper()
                else:
                    # need reverse complement of sequence for - strand genes, and prev/next are swapped
                    dseq = str(f_dict[chrSeqRecord.id][curDonors[b_ind]-CHECK_EQUIV_LEN:curDonors[b_ind]+CHECK_EQUIV_LEN].seq.reverse_complement()).upper()
                    dseq_prev = dseq[:CHECK_EQUIV_LEN]
                    dseq_next = dseq[CHECK_EQUIV_LEN:]
                    
                    aseq = str(f_dict[chrSeqRecord.id][curAcceptors[b_ind+1]-CHECK_EQUIV_LEN:curAcceptors[b_ind+1]+CHECK_EQUIV_LEN].seq.reverse_complement()).upper()
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
                
                o_handle.write("\t".join([equivJuncSeq, getGeneName(sftr), chrSeqRecord.id, str(curDonors[b_ind]), str(curAcceptors[b_ind+1])])) 
                o_handle.write("\n")
        except Exception as e:
            print "Exception"
            print e
            print "error:", sys.exc_info()[0]
            print "parsing features for", ftr
            
            
if __name__  == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fastaFile', required=True, help='path to fasta file with chromosome sequences')
    parser.add_argument('-a', '--annotationFile', required=True, help='path to gff file with exon annotations')
    parser.add_argument('-o', '--outFile', required=True, help='output file name (including full path)')
    parser.add_argument('-n1', '--name1', help='name of field in gtf to use for gene names', default='gene_name')
    parser.add_argument('-n2', '--name2', help='name of field in gtf to use for gene names if n1 does not exist', default='gene_id')
    parser.add_argument('-v', '--verbose', help='print info about data obtained', action='store_true')
    args = parser.parse_args()
    
    # read in the sequences 
    f_handle = open(args.fastaFile, "rU")
    f_dict = SeqIO.to_dict(SeqIO.parse(f_handle, "fasta"))
    f_handle.close()
    
    # only want exons for now, and only chromosomes that we have sequences for
    limit_info = dict(
        gff_id = f_dict.keys(),
        gff_type = ["exon"])
    
    
    ########### read in the annotations, adding annotations to the sequences 
    a_handle = open(args.annotationFile, "rU")
    if args.verbose:
        print "opening annotation file", args.annotationFile
    
    o_handle = open(args.outFile, "wb")
    
    for rec in GFF.parse(a_handle, base_dict=f_dict, limit_info=limit_info): # each rec is a SeqRecord (for 1 chromosome)
        if args.verbose:
            print "####### starting new chromosome: " + str(rec.id)
        # populate the data structures with genes and exons from this chromosome and pickle the objects
        parseChrFeatures(rec)
    
    o_handle.close()
    a_handle.close()