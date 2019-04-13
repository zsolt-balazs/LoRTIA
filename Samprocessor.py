#!/usr/bin/env python3

###############################################################################
###                            SAMprocessor                                 ###
###                 This script handles the .sam files.                     ###
###                And oh boy, the way it handles them...                   ###
###############################################################################

import pysam
import operator
import os
import subprocess
import pandas as pd
from argparse import ArgumentParser
from Bio.Seq import Seq
from Bio import pairwise2

###############################################################################
###                               ADAPTERS                                  ###
###############################################################################
#TeloPrime 5' TGGATTGATATGTAATACGACTCACTATAG
#PacBio 5' AGAGTACATGGG
#MinION 5' TTTCTGTTGGTGCTGATATTGCTGCCATTAGGCCGGG
#MinION 3' ACTTGCCTGTCGCTCTATCTTC


###############################################################################
###                 Bedtools and Samtools-like functions                    ###
###############################################################################

def input_sorter(in_file, prefix):
    """
    Calls the pysam sort function on the input.
    """
    print("Sorting {}...".format(in_file))
    pysam.sort("-n", "-O" "SAM", "-o", "{}_temp.sam".format(prefix), in_file)

# The code below was meant to run bedtools and get output line-by-line and
# only write out non0 coverage in order to buy space, but it runs incredibly
# slowly.
#
#def runner(bam, tsv, strand):
#    cmd = ["genomeCoverageBed", "-ibam {}".format(bam), "-d", "-split {}".format(strand)]
#    cmd = ["genomeCoverageBed -ibam {} -d -split {}".format(bam, strand)]
#    btools = Popen(cmd, stdout=PIPE, shell = True)
#    while True:
#        line = btools.stdout.readline().decode('utf-8')
#        if not line:
#            break
#        if line[-3:] != "\t0\n":
#            out_appender(line, tsv)
#    btools.stdout.close()

def del0s(bam, strand, tsv):
    """
    Runs bedtools to create the genome coverage and removes the 0 values.
    """
    cmd = "bedtools genomecov -ibam {} -d -split{}".format(bam, strand + tsv)
    subprocess.run(cmd, shell=True)
    with open(tsv.replace(".tsv", "_no0.tsv"), "a") as outfile:
        with open(tsv) as in_tsv:
            for line in in_tsv:
                if line.strip().split("\t")[2] != "0":
                    outfile.write(line)
    os.remove(tsv)
    os.rename(tsv.replace(".tsv", "_no0.tsv"), tsv)

def coverage_counter(outbam, out_stranded_bam):
    """
    Calls del0s to create the coverage files.
    """
    print("Counting coverage...")
    alltsv = outbam.replace("sorted.bam", "allcov.tsv")
    mintsv = outbam.replace("sorted.bam", "minuscov.tsv")
    plustsv = outbam.replace("sorted.bam", "pluscov.tsv")
    del0s(outbam, " > ", alltsv)
    del0s(out_stranded_bam, " -strand - > ", mintsv)
    del0s(out_stranded_bam, " -strand + > ", plustsv)

def output_creator(outsam):
    """
    Uses pysam to sort and index the output and the stranded_only output files.
    """
    print("Sorting and indexing output...")
    outbam = outsam.replace(".sam", "_sorted.bam")
    out_stranded_bam = outsam.replace("out.sam", "stranded_only.bam")
    pysam.sort("-o", outbam, outsam)
    pysam.index(outbam)
    pysam.view("-h",
               "-b",
               "-F 512",
               "-o",
               out_stranded_bam,
               outbam,
               catch_stdout=False)
    pysam.index(out_stranded_bam)
    coverage_counter(outbam, out_stranded_bam)


###############################################################################
###                           Cigar functions                               ###
###############################################################################

def get_left_end(read_seq, cigar, args):
    """
    Gets the left end of the alignment, taking 'check_in_soft' number of 
    nucleotides from the softclip and check_in_match number of nucleotides
    from the mapped part.
    """
    if cigar[0][0] in [4, 5]:
        softleftpos = cigar[0][1]
        if softleftpos > args.check_in_soft - 1:
            checkleft = softleftpos - args.check_in_soft
            cis = args.check_in_soft
        else:
            checkleft = 0
            cis = softleftpos
        leftend_seq = read_seq[checkleft:softleftpos + args.check_in_match]
        left_match = read_seq[softleftpos:softleftpos + args.match_in_first]
    else:
        leftend_seq = "no softclip"
        cis = 0
        left_match = 0
    return leftend_seq, cis, left_match

def get_right_end(read_seq, cigar, args):
    """
    Gets the right end of the alignment, taking 'check_in_soft' number of 
    nucleotides from the softclip and check_in_match number of nucleotides
    from the mapped part.
    """
# here I first calculate the length of the read and I
# substract the length of the closing softclip
    if cigar[-1][0] in [4, 5]:
        softrightpos = (len(read_seq) - cigar[-1][1])
        if args.check_in_soft < cigar[-1][1] - 1:
            checkright = softrightpos + args.check_in_soft
            cis = args.check_in_soft
        else:
            checkright = len(read_seq)
            cis = cigar[-1][1]
        rightend_seq = read_seq[softrightpos - args.check_in_match:checkright]
        right_match = read_seq[softrightpos - args.match_in_first:softrightpos]
    else:
        rightend_seq = "no softclip"
        cis = 0
        right_match = 0
    return rightend_seq, cis, right_match

###############################################################################
###                             Adapter functions                           ###
###############################################################################

def adapter_aligner(sequence, adapter, args):
    """
    Aligns the adapter to a string using the Smith-Waterman algorithm.
    """
    alignments = pairwise2.align.localms(adapter,
                                         sequence,
                                         args.match_score,
                                         args.mismatch_score,
                                         args.gap_open_score,
                                         args.gap_extend_score)
    return alignments
    # alignments is a list of alignment tuples, an alignment has five elements:
    # adapter, query, aln_score, aln_start, aln_end

def pos_wo_gap(alignment, p):
    """
    Returns the position in the alignment, substracting the gaps.
    """
    # the alignment end is reported counting the gaps as well, we change that
    pos_wo_gap = alignment[p] - alignment[1][:alignment[p]].count("-")
    return pos_wo_gap

def in_place_checker(alignments,
                     cigar_info,
                     args,
                     check_in_soft,
                     matchseq,
                     adapter):
    """
    Checks whether the adapter is no further than nts away from the start of
    the alignment.
    """
    ts = check_in_soft + args.shs_for_ts
    soft_match = check_in_soft + args.check_in_match
    adapter_info = (alignments[0][2],
                    pos_wo_gap(alignments[0], 3),
                    pos_wo_gap(alignments[0], 4),
                    "out of place")
    for alignment in alignments:
        if pos_wo_gap(alignment, 4) in range(check_in_soft 
                                              - args.check_from_alignment,
                                              soft_match + 1):
            adapter_info = (alignment[2],
                            pos_wo_gap(alignment, 3),
                            pos_wo_gap(alignment, 4),
                            "correct")
            break
    if adapter_info[3] == "correct":
        exon = 0
        for event in cigar_info:
            event_type, event_count = event
            if event_type == 0:
                if exon + int(event_count) < args.first_exon + 1:
                    exon += int(event_count)
                else:
                    break
            elif event_type == 3:
                match_adapter = adapter_aligner(matchseq, adapter, args)
                if match_adapter:
                    score_limit = args.match_in_first * 0.5 * args.match_score
                    if match_adapter[0][2] >= score_limit:
                        adapter_info = (alignment[2],
                                        pos_wo_gap(alignment, 3),
                                        pos_wo_gap(alignment, 4),
                                        "false exon")
                break
            else:
                for alignment in alignments:
                    if pos_wo_gap(alignment, 4) in range(ts, soft_match + 1):
                        adapter_info = (alignment[2],
                                        pos_wo_gap(alignment, 3),
                                        pos_wo_gap(alignment, 4),
                                        "potential template switching")
                        break
    return adapter_info

def get_adapter_info(sequence,
                     adapter,
                     score_limit,
                     cigar_info,
                     args,
                     check_in_soft,
                     matchseq):
    """
    Sends sequence to the aligner, and gathers adapter information.
    """
    alignments = adapter_aligner(sequence, adapter, args)

    if alignments: # sometimes pairwise2 does not come up with any alignment.
        # Pairwise2 only lists alignments with the best score, so checking the
        # 1st one is enough.

        isadapter = float(alignments[0][2]) > score_limit
        if isadapter:
            adapter_info = in_place_checker(alignments,
                                            cigar_info,
                                            args,
                                            check_in_soft,
                                            matchseq,
                                            adapter)
        else:
            adapter_info = (alignments[0][2],
                            pos_wo_gap(alignments[0], 3),
                            pos_wo_gap(alignments[0], 4),
                            "missing")
    else:
        adapter_info = (0, 0, 0, "missing")
    return adapter_info

def adapter_checker(read, args):
    """
    Retrieves the end sequences and sorts adapter information for each end.
    """
    cigar = read.cigartuples
    seq = read.query_sequence
    leftend, check_in_softl, left_match = get_left_end(seq, cigar, args)
    rightend, check_in_softr, right_match = get_right_end(seq, cigar, args)  
    left_cigar_info = cigar
    right_cigar_info = list(cigar)[::-1]
    if leftend != "no softclip":
        three_left = get_adapter_info(Seq(leftend).complement(),
                                      args.three_adapter,
                                      args.three_score,
                                      left_cigar_info,
                                      args,
                                      check_in_softl,
                                      Seq(left_match).complement())
        five_left = get_adapter_info(leftend,
                                     args.five_adapter,
                                     args.five_score,
                                     left_cigar_info,
                                     args,
                                     check_in_softl,
                                     left_match)
    else:
        three_left = (0, 0, 0, "no softclip")
        five_left = (0, 0, 0, "no softclip")
    if rightend != "no softclip":
        three_right = get_adapter_info(rightend[::-1],
                                       args.three_adapter,
                                       args.three_score,
                                       right_cigar_info,
                                       args,
                                       check_in_softr,
                                       right_match[::-1])
        five_right = get_adapter_info(Seq(rightend).reverse_complement(),
                                      args.five_adapter,
                                      args.five_score,
                                      right_cigar_info,
                                      args,
                                      check_in_softr,
                                      Seq(right_match).reverse_complement())
    else:
        three_right = (0, 0, 0, "no softclip")
        five_right = (0, 0, 0, "no softclip")
    return {"l3": three_left,
            "r3": three_right,
            "l5": five_left,
            "r5": five_right}


###############################################################################
###                              Sam modifiers                              ###
###############################################################################

def out_appender(line, out_file):
    """
    Appends string to a file.
    """
    with open(out_file, "a") as out:
        out.write(line)

def out_writer(dictionary, prefix, suffix):
    """
    Writes a dictionary into a csv file.
    """
    dataframe = pd.DataFrame.from_dict(dictionary, orient="index")
    dataframe.to_csv("{}_{}.tsv".format(prefix, suffix), header=False,
                     sep="\t")

def set_read_strand(summary, old_flag):
    """
    Rewrites flags into 0 for +, 16 for - and 256 for unknown
    """
    is_l3 = summary.get("l3")[3] not in ["missing", "no softclip"]
    is_r3 = summary.get("r3")[3] not in ["missing", "no softclip"]
    is_l5 = summary.get("l5")[3] not in ["missing", "no softclip"]
    is_r5 = summary.get("r5")[3] not in ["missing", "no softclip"]
    if (is_l3 or is_r5) and not (is_r3  or is_l5):
        new_flag = 16
    elif (is_r3 or is_l5) and not (is_l3 or is_r5):
        new_flag = 0
    else:
        new_flag = old_flag + 512
    return new_flag

def intron_finder(read, read_start, read_end, args):
    """
    Finds introns based on the cigar code. Also filters introns based on the 
    first_exon argument and removes introns which contain a long insertion
    preceding the intron (likely triple-chimeric reads).
    """
    introns = ()
    counter = 0
    previous_event = 0, 0
    matches_so_far = 0
    matches_to_come = 0
    is_left_false = ((read.get_tag("l5").split(",")[3] != "correct") or 
                     (read.get_tag("l3").split(",")[3] != "correct"))
    for event_type, event_count in read.cigartuples:
        # event coding: M:0, I:1, D:2, N:3, S:4, H:5
        previous_is_insert = (previous_event[0] == 1 
                              and previous_event[1]
                              > args.insert_before_intron)
        # we set this so that triple-chimeric reads do not appear as intronic
        if event_type == 3:
            intron_start = read_start - 1
            for element_type, element_count in read.cigartuples[:counter]:
                if element_type in [0, 2, 3]:
                    intron_start += element_count
            intron_end = intron_start + event_count + 1
            # we added one, because we actually want the exon start
            introns += (intron_start, intron_end),
            matches_to_come = 0
            if (previous_is_insert 
                or (is_left_false and (matches_so_far < args.first_exon))):
                read.set_tag("ga", ",".join((str(i) for i in introns)), "Z")
                introns = ()
        elif event_type == 0:
            matches_so_far += event_count
            matches_to_come += event_count
        counter += 1
        previous_event = event_type, event_count
    is_right_false = ((read.get_tag("r5").split(",")[3] != "correct") or 
                      (read.get_tag("r3").split(",")[3] != "correct"))
    if matches_to_come < args.first_exon and is_right_false and introns:
        read.set_tag("ga", str(introns[-1]), "Z")
        introns = introns[:-1]
    return introns

def put_in_dict(adapter_sum, end_type, string, dictionary, contig, end):
    """
    Increases adapter value in a dicrionary.
    """
    if adapter_sum.get(end_type)[3] == string:
        if (contig, end) in dictionary:
            dictionary[contig, end] += 1
        else:
            dictionary[contig, end] = 1
    return dictionary

def prepare_new_sam_line(read,
                         args,
                         l3_dict,
                         r3_dict,
                         l5_dict,
                         r5_dict,
                         introns_dict,
                         tsl3_dict,
                         tsr3_dict,
                         tsl5_dict,
                         tsr5_dict):
    """"
    This is the function that summmarizes read info and writes it into a new
    line in the output file
    """
    read_start = read.reference_start + 1
    read_end = read.reference_end
    contig = read.reference_name
    adapter_sum = adapter_checker(read, args)
    read.flag = set_read_strand(adapter_sum, read.flag)
    if read.flag > 511:
        read.mapping_quality = 0
    else:
        read.mapping_quality = 9
    l3_dict = put_in_dict(adapter_sum,
                          "l3",
                          "correct",
                          l3_dict,
                          contig,
                          read_start)
    r3_dict = put_in_dict(adapter_sum,
                          "r3",
                          "correct",
                          r3_dict,
                          contig,
                          read_end)
    l5_dict = put_in_dict(adapter_sum,
                          "l5",
                          "correct",
                          l5_dict,
                          contig,
                          read_start)
    r5_dict = put_in_dict(adapter_sum,
                          "r5",
                          "correct",
                          r5_dict,
                          contig,
                          read_end)
    tsl3_dict = put_in_dict(adapter_sum,
                            "l3",
                            "potential template switching",
                            tsl3_dict,
                            contig,
                            read_start)
    tsr3_dict = put_in_dict(adapter_sum,
                            "r3",
                            "potential template switching",
                            tsr3_dict,
                            contig,
                            read_end)
    tsl5_dict = put_in_dict(adapter_sum,
                            "l5",
                            "potential template switching",
                            tsl5_dict,
                            contig,
                            read_start)
    tsr5_dict = put_in_dict(adapter_sum,
                            "r5",
                            "potential template switching",
                            tsr5_dict,
                            contig,
                            read_end)
    read.set_tag("re", read_end, "i")
    read.set_tag("l3", ",".join((str(i) for i in adapter_sum.get("l3"))), "Z")
    read.set_tag("r3", ",".join((str(i) for i in adapter_sum.get("r3"))), "Z")
    read.set_tag("l5", ",".join((str(i) for i in adapter_sum.get("l5"))), "Z")
    read.set_tag("r5", ",".join((str(i) for i in adapter_sum.get("r5"))), "Z")
    introns = intron_finder(read, read_start, read_end, args)
    if introns:
        for intron in introns:
            if (contig, intron) in introns_dict:
                introns_dict[contig, intron] += 1
            else:
                introns_dict[contig, intron] = 1
    read.set_tag("in", ",".join((str(i) for i in introns)), "Z")
    readline = read.to_string()
    out_appender("".join((readline,"\n")), args.out_file)
    return (l3_dict,
            r3_dict,
            l5_dict,
            r5_dict,
            introns_dict,
            tsl3_dict,
            tsr3_dict,
            tsl5_dict,
            tsr5_dict)

def deal_with_same_name(previous_lines, 
                        args,
                        l3_dict,
                        r3_dict,
                        l5_dict,
                        r5_dict,
                        introns_dict,
                        tsl3_dict,
                        tsr3_dict,
                        tsl5_dict,
                        tsr5_dict):
    """
    Renames non-overlapping reads with the same name and calls 
    prepare_new_sam_line to process reads one-by-one.
    """
    counter = 1
    previous_lines.sort(key=operator.itemgetter(3), reverse=True)
    ranges_by_far = ()
    for read, read_start, read_end, read_span in previous_lines:
        is_in_previous = False
        for ran in ranges_by_far:
            if read_start in ran or read_end in ran:
                is_in_previous = True
                break
        if not is_in_previous:
            if len(previous_lines) > 1:
                read.query_name = "_".join((read.query_name, str(counter)))
                counter += 1
            (l3_dict,
            r3_dict,
            l5_dict,
            r5_dict,
            introns_dict,
            tsl3_dict,
            tsr3_dict,
            tsl5_dict,
            tsr5_dict) = prepare_new_sam_line(read,
                                              args,
                                              l3_dict,
                                              r3_dict,
                                              l5_dict,
                                              r5_dict,
                                              introns_dict,
                                              tsl3_dict,
                                              tsr3_dict,
                                              tsl5_dict,
                                              tsr5_dict)
        ranges_by_far += range(read_start, read_end + 1),
    return (l3_dict,
            r3_dict,
            l5_dict,
            r5_dict,
            introns_dict,
            tsl3_dict,
            tsr3_dict,
            tsl5_dict,
            tsr5_dict)
    
def sam_iterator(args):
    """
    Reads the SAM file line by line, puts the reads with the same name into a
    list and calls deal_with_same_name to process them
    """
    sam = pysam.AlignmentFile("{}_temp.sam".format(args.prefix), "r")
    previous_lines = []
    prog_comment = "@PG	ID:LoRTIA	PN:LoRTIA	VN:0.1	CL:Samprocessor.py {},\
                    \n".format(args)
    outsam = pysam.AlignmentFile(args.out_file, "w", template=sam)
    outsam.close()
    out_appender(prog_comment, args.out_file)
    l3_dict = {}
    r3_dict = {}
    l5_dict = {}
    r5_dict = {}
    introns_dict = {}
    tsl3_dict = {}
    tsr3_dict = {}
    tsl5_dict = {}
    tsr5_dict = {}
    print("Iterating over {}...".format(args.in_file))
    for read in sam:
        if not read.is_unmapped:
            read_start = read.reference_start + 1
            read_end = read.reference_end
            read_span = read_end - read_start
            if (previous_lines and
                    read.query_name != previous_lines[0][0].query_name):
                # this makes sure that the previous reads are processed if the
                # next read does not have the same name
                (l3_dict,
                 r3_dict,
                 l5_dict,
                 r5_dict,
                 introns_dict,
                 tsl3_dict,
                 tsr3_dict,
                 tsl5_dict,
                 tsr5_dict) = deal_with_same_name(previous_lines,
                                                  args,
                                                  l3_dict,
                                                  r3_dict,
                                                  l5_dict,
                                                  r5_dict,
                                                  introns_dict,
                                                  tsl3_dict,
                                                  tsr3_dict,
                                                  tsl5_dict,
                                                  tsr5_dict)
                previous_lines = []
            previous_lines.append((read, read_start, read_end, read_span))
    (l3_dict,
     r3_dict,
     l5_dict,
     r5_dict,
     introns_dict,
     tsl3_dict,
     tsr3_dict,
     tsl5_dict,
     tsr5_dict) = deal_with_same_name(previous_lines,
                                      args,
                                      l3_dict,
                                      r3_dict,
                                      l5_dict,
                                      r5_dict,
                                      introns_dict,
                                      tsl3_dict,
                                      tsr3_dict,
                                      tsl5_dict,
                                      tsr5_dict)
    sam.close()
    print("Generating statistics...")
    out_writer(r3_dict, args.prefix, "r3")
    out_writer(l3_dict, args.prefix, "l3")
    out_writer(l5_dict, args.prefix, "l5")
    out_writer(r5_dict, args.prefix, "r5")
    out_writer(tsr3_dict, args.prefix, "ts_r3")
    out_writer(tsl3_dict, args.prefix, "ts_l3")
    out_writer(tsl5_dict, args.prefix, "ts_l5")
    out_writer(tsr5_dict, args.prefix, "ts_r5")
    out_writer(introns_dict, args.prefix, "in")

def Samprocessor(args):
    """
    Sets argument types and processes the sam file.
    """
    if not os.path.isdir(args.out_path):
        os.mkdir(args.out_path)
    in_prefix = os.path.basename(args.in_file)
    args.prefix = "{}/{}".format(args.out_path,
                                 in_prefix[:len(in_prefix) - 4])
    args.out_file = "{}_out.sam".format(args.prefix)

    input_sorter(args.in_file, args.prefix)
    sam_iterator(args)
    output_creator(args.out_file)
    print("Processed files saved to {}\n".format(args.out_path))

###############################################################################
###                             Main function                               ###
###############################################################################

def main():
    args = parsing()
    Samprocessor(args)

def parsing():
    """
    This part handles the commandline arguments
    """
    parser = ArgumentParser(description="This is the first module of LoRTIA:\
                            a Long-read RNA-Seq Transcript Isofom Annotator")
    parser.add_argument("in_file",
                        help="Input file. Both .sam and .bam files are\
                        accepted.",
                        metavar="input_file")
    parser.add_argument("out_path",
                        help="Output folder. Multiple output files are going\
                        to be created using the input file's prefix (ie. the\
                        part that precedes '.bam' or '.sam')",
                        metavar="output_path")
    parser.add_argument("--match_score", 
                        dest="match_score",
                        help="The alignment scores for each match when \
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: 2",
                        type=float,
                        metavar="[float]", 
                        default=2.0)
    parser.add_argument("--mismatch_score", 
                        dest="mismatch_score",
                        help="The alignment scores for each mismatch when \
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: -3",
                        type=float,
                        metavar="[float]", 
                        default=-3.0)
    parser.add_argument("--gap_open_score", 
                        dest="gap_open_score",
                        help="The alignment scores for each gap opening when\
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: -3",
                        type=float,
                        metavar="[float]", 
                        default=-3.0)
    parser.add_argument("--gap_extend_score", 
                        dest="gap_extend_score",
                        help="The alignment scores for each gap extension \
                        when searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: -3",
                        type=float,
                        metavar="[float]", 
                        default=-3.0)
    parser.add_argument("-3", "--three_adapter", 
                        dest="three_adapter",
                        help="The 3' adapter to look for, the default is a \
                        polyA tail of 30 adenines",
                        metavar="[string]",
                        default=30*"A")
    parser.add_argument("-5", "--five_adapter", 
                        dest="five_adapter",
                        help="The 5' adapter to look for. The default is the\
                        TeloPrime cap adapter: TGGATTGATATGTAATACGACTCACTATAG",
                        metavar="[string]",
                        default="TGGATTGATATGTAATACGACTCACTATAG")
    parser.add_argument("--five_score", 
                        dest="five_score",
                        help="The minimum score that the adapter alignment \
                        should reach to be recognized as a 5' adapter. The \
                        default score is 20.0",
                        type=int,
                        metavar="[float]",
                        default=20.0)
    parser.add_argument("--three_score", 
                        dest="three_score",
                        help="The minimum score that the adapter alignment \
                        should reach to be recognized as a 3' adapter. The \
                        default score is 20.0",
                        type=float,
                        metavar="[float]",
                        default=20.0)
    parser.add_argument("--check_in_soft", 
                        dest="check_in_soft",
                        help="The number of nucleotides from the start of \
                        the soft clip where the adapter sequence is searched \
                        for. It is not advisable to search for adapters too \
                        deep into the soft clip not only because it may \
                        increase running time, but also because it increases \
                        false positive hits. The default is 30.",
                        type=int,
                        metavar="[integer]",
                        default=30)
    parser.add_argument("--check_in_match", 
                        dest="check_in_match",
                        help="The number of nucleotides from the start of \
                        the soft clip where the adapter sequence is searched \
                        for. This lets adapters be discovered even if they \
                        resemble some genomic segments and therefore align to \
                        the genome. The default is 10.",
                        type=int,
                        metavar="[integer]",
                        default=10)
    parser.add_argument("--check_from_alignment", 
                        dest="check_from_alignment",
                        help="The maximum distance of the adapter from the \
                        alignment. Some bases may be inserted or miscalled \
                        between the adapter and the mapping part of the read, \
                        This parameter allows those adapters to be considered \
                        'correct' which do not start exactly at the end of \
                        the mapped part but are a few bases off. The default \
                        value is 3.",
                        type=int,
                        metavar="[integer]",
                        default=3)
    parser.add_argument("--shs_for_ts", 
                        dest="shs_for_ts",
                        help="The minimum length of agreement (Short \
                        Homologous sequence, SHS) between the start of the \
                        match part of the alignment and the adapter that \
                        raises a suspicion of template switching. Putative\
                        template switching artefacts are listed in a \
                        separate file and are excluded form further \
                        statistics. The value has to be lesser than the value\
                        set by --check_in_match. If greater or equal value is\
                        set, the program will not look for signs of template \
                        switching. The default value is 3 nucleotides.",
                        type=int,
                        metavar="[integer]",
                        default=3)
    parser.add_argument("--first_exon", 
                        dest="first_exon",
                        help="Alignment ends are often placed far away from \
                        from the rest of the read if the adapter maps to a \
                        nearby part of the genome. This option sets the length\
                        of the first exon, under which the matching part of\
                        of the alignment should be checked for the presence of\
                        the adapters. The default is 30.",
                        type=int,
                        metavar="[integer]",
                        default=30)
    parser.add_argument("--match_in_first", 
                        dest="match_in_first",
                        help="Alignment ends are often placed far away from \
                        from the rest of the read if the adapter maps to a \
                        nearby part of the genome. With this option, the user\
                        can set how many nucleotides from the matching part \
                        of the alignment should be aligned to the adapter \
                        if at least half of the nucleotides match to the \
                        adapter, the exon will be considered false. The \
                        default is 15.",
                        type=int,
                        metavar="[integer]",
                        default=15)
    parser.add_argument("--insert_before_intron", 
                        dest="insert_before_intron",
                        help="The maximum allowed insert length immediately \
                        before an intron. Triple-chimeric reads are often \
                        as an exon, a long insert and another exon. In these \
                        cases the inserts are usually several hundred nts \
                        long and are unmapped because they stem from another \
                        contig or would be mapped to the complementary strand.\
                        The default value is 20 nts.",
                        type=int,
                        metavar="[integer]",
                        default=20)

    return parser.parse_args()



if __name__== "__main__":
    main()

#Explanation of constants for adapter search: 

#S: softclip      M: Match (alignment)

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#    |____________ADAPTER________|
#    -----------------------------: CHECK_IN_SOFT
#    how many nucleotides are checked from the softclip

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#                                 ADAPTER\
#                  CHECK_IN_MATCH:-----
# how many nucleotides are checked from the beginning of the alignment

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#        |______ADAPTER______|
#        CHECK_FROM_ALIGNMENT:---
# how many bases are allowed between the start of the alignment and the adapter

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#                      FIRST_EXON:--------------------------------------
# if an intron starts in this region, the adapter will be classified as
# 'out of place',

