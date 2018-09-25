
###############################################################################
###                            SAMprocessor                                 ###
###                 This script handles the .sam files.                     ###
###                And oh boy, the way it handles them...                   ###
###############################################################################

from Bio.Seq import Seq
from Bio import pairwise2
import re
from argparse import ArgumentParser


###############################################################################
###                             Set variables                               ###
###############################################################################

parser = ArgumentParser(description = "This is the first module of LoRTIA: a Long-read RNA-Seq Transcript Isofom Annotator")
parser.add_argument("input", help = "Input file", metavar = "FILE")
parser.add_argument("output",
                    help = "Output file", metavar = "FILE")
parser.add_argument("-3", "--threeprime", dest = "three_adapter",
                    help = "The 3' adapter to look for, the default is a polyA tail",
                    metavar = "string", default = 30*"A", required=False)
parser.add_argument("-5", "--fiveprime", dest = "five_adapter",
                    help = "The 5' adapter to look for, the default is the TeloPrime cap adapter",
                    metavar = "string", default = "TGGATTGATATGTAATACGACTCACTATAG",
                    required=False)

args = parser.parse_args()


# Constants of the adapter search
CHECK_IN_SOFT = 30
CHECK_IN_MATCH = 5
CHECK_FROM_ALIGNMENT = 3
FIRST_EXON = 50
FIRST_EXON_QUAL = 0.6

# Explanation of constants for adapter search:

# S: softclip      M: Match (alignment)

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#    \____________ADAPTER________|
#    -----------------------------: CHECK_IN_SOFT
#       how many nucleotides are checked from the softclip

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#                               ADAPTER
#                  CHECK_IN_MATCH:-----
#  how many nucleotides are checked from the beginning of the alignment

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#         \______ADAPTER______|
#         CHECK_FROM_ALIGNMENT:---
# how many bases are allowed between the start of the alignment and the adapter

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
#                      FIRST_EXON:--------------------------------------
# more than FIRST_EXON_QUAL ratio of deletions or introns in this region renders the adapter false

# Scores for adapter alignment
MATCH_SCORE = 2
MISMATCH_SCORE = -3
GAP_OPEN_SCORE = -3
GAP_EXT_SCORE = -3

# Required scores for 5' and 3' adapters
FIVE_SCORE_LIMIT = 20
THREE_SCORE_LIMIT = 20

FIVE_ADAPTER = args.five_adapter
THREE_ADAPTER = args.three_adapter

SAM_FILE = args.input
OUT_FILE = args.output
prog_comment = "@PG	ID:LoRTIA	PN:LoRTIA	VN:0.1	CL:Samprocessor.py{}, {}, {}\n".format(SAM_FILE, OUT_FILE, THREE_ADAPTER, FIVE_ADAPTER)

###############################################################################
###                         Processing the input                            ###
###############################################################################

def cigar_interpreter(cigar_string):
    """
    Creates a list of cigar events and sums them.
    """
# I decided to declare cigar events as global variables,
# because they are needed a lot and I didn't want to pass
# them along many functions (is it bad practice?)
    global M, I, S, D, H, N, cigar_list
    M = 0
    I = 0
    S = 0
    D = 0
    H = 0
    N = 0
    cigar_list = []
    for number, event, in re.findall("(\d+)([ISMDHN])",cigar_string):
        cigar_list.append("_".join((event, number)))
        if event == "M":
            M += int(number)
        elif event == "I":
            I += int(number)
        elif event == "S":
            S += int(number)
        elif event == "D":
            D += int(number)
        elif event == "H":
            H += int(number)
        elif event == "N":
            N += int(number)
    return M, I, S, D, H, N, cigar_list

def get_left_end(read_seq):
    softclipleftpos = int(cigar_list[0][2:])
    if softclipleftpos > CHECK_IN_SOFT - 1:
        checkleft = softclipleftpos-CHECK_IN_SOFT
    else:
        checkleft = 0
    leftend_seq = read_seq[checkleft : softclipleftpos + CHECK_IN_MATCH]
    return leftend_seq

def get_right_end(read_seq):
# here I first calculate the length of the read and I 
# substract the length of the closing softclip
    softcliprightpos = (len(read_seq) - int(cigar_list[len(cigar_list) - 1][2:]) - 1)
    if CHECK_IN_SOFT < int(cigar_list[len(cigar_list)-1][2:]) - 1:
        checkright = softcliprightpos - 1 + CHECK_IN_SOFT
    else:
        checkright = len(read_seq)
    rightend_seq = read_seq[softcliprightpos - CHECK_IN_MATCH : checkright]
    return rightend_seq


###############################################################################
###                             Adapter functions                           ###
###############################################################################

def adapter_aligner(string_sequence, adapter):
    alignment = pairwise2.align.localms(adapter, string_sequence, MATCH_SCORE, MISMATCH_SCORE, GAP_OPEN_SCORE, GAP_EXT_SCORE)
    if alignment:
        return alignment[0][2:]
    # sometimes pairwise2 does not come up with any alignment
    else:
        return 0, -1, -1

def in_place_checker(alignment, cigar_info):
    """
    Checks whether the adapter is no further than nts away from the start of 
    the alignment
    """
    where_it_has_to_be = range(int(alignment[1]), int(alignment[2]) + CHECK_FROM_ALIGNMENT)
    match = 0
    non_match = 0
    for event in cigar_list:
        event_type, event_count = event.split("_")
        if event_type == "M":
            if match + non_match + int(event_count) < FIRST_EXON + 1:
                match += int(event_count)
            else:
                match = FIRST_EXON - non_match
        elif event_type in ("D", "I", "N"):
            if match + non_match + int(event_count) < FIRST_EXON + 1:
                non_match += int(event_count)
            else:
                non_match = FIRST_EXON - match
    return ((CHECK_IN_SOFT in where_it_has_to_be) and (match / (match + non_match) > FIRST_EXON_QUAL))


def adapter_checker(read_sequence, cigar_list):
    """
    Calls the aligner and checks if the alignment has enough score and is in
    the right place for the adapter
    """
    leftend = get_left_end(read_sequence)
    rightend = get_right_end(read_sequence)
    left_cigar_info = cigar_list
    right_cigar_info = cigar_list.reverse
    
    #Well, I guess the following could be written shorter somehow...
    three_left = adapter_aligner(Seq(leftend).reverse_complement(), THREE_ADAPTER)
    isthreeleft = int(three_left[0]) > THREE_SCORE_LIMIT
    three_left += isthreeleft,
    if isthreeleft:
        isthreeleftinplace = in_place_checker(three_left, left_cigar_info)
    else:
        isthreeleftinplace = False
    three_left += isthreeleftinplace,
    
    three_right = adapter_aligner(rightend, THREE_ADAPTER)
    isthreeright = int(three_right[0]) > THREE_SCORE_LIMIT
    three_right += isthreeright,
    if isthreeright:
        isthreerightinplace = in_place_checker(three_right, right_cigar_info)
    else:
        isthreerightinplace = False
    three_right += isthreerightinplace,
    
    five_left = adapter_aligner(leftend, FIVE_ADAPTER)
    isfiveleft = int(five_left[0]) > FIVE_SCORE_LIMIT
    five_left += isfiveleft,
    if isfiveleft:
        isfiveleftinplace = in_place_checker(five_left, left_cigar_info)
    else:
        isfiveleftinplace = False
    five_left += isfiveleftinplace,
        
    five_right = adapter_aligner(Seq(rightend).reverse_complement(), FIVE_ADAPTER)
    isfiveright = int(five_right[0]) > FIVE_SCORE_LIMIT
    five_right += isfiveleft,
    if isfiveright:
        isfiverightinplace = in_place_checker(five_right, right_cigar_info)
    else:
        isfiverightinplace = False
    five_right += isfiverightinplace,
    return {"l3": three_left, "r3": three_right, "l5": five_left, "r5": five_right}
    

###############################################################################
###                          Sam attribute modifiers                        ###
###############################################################################

def set_read_strand(summary):
    """
    Rewrites flags into 0 for +, 16 for - and 256 for unknown
    """
    if (summary.get("l3")[3] or summary.get("r5")[3]) and not (summary.get("r3")[3] or summary.get("l5")[3]):
        new_flag = 16
    elif (summary.get("r3")[3] or summary.get("l5")[3]) and not (summary.get("l3")[3] or summary.get("r5")[3]):
        new_flag = 0
    else:
        new_flag = 256
    return new_flag

def intron_finder(read_start, cigar_list):
    introns = ()
    counter = 0
    for event in cigar_list:
        if event.startswith("N"):
            intron_start = read_start
            for element in cigar_list[:counter]:
                if element.startswith("M") or element.startswith("D") or element.startswith("N"):
                    intron_start += int(element.split("_")[1])
            intron_end = intron_start + int(event.split("_")[1])
            introns += "..".join((str(intron_start), str(intron_end))),
        counter += 1
    return introns

def out_appender(line):
    with open(OUT_FILE, "a") as out:
        out.write(line)
    out.close

###############################################################################
###                             Main function                               ###
###############################################################################

def main():
    """
    Reads the SAM file line by line and writes a new one 
    """
    with open(SAM_FILE) as sam:
        for line in sam:
            line = line.split("\t")
            if line[0].startswith('@'):
                commentline = "\t".join(line)
                out_appender(commentline)
                if line[0] == "@PG":
                    out_appender(prog_comment)
            elif line and not line[0].startswith('@') and not line[1] == "4":
                cigar_interpreter(line[5])
                read_length = len(line[9])
                read_start = int(line[3])
                read_end = int(line[3]) + M + D + N
                adapter_summary = adapter_checker(line[9], cigar_list)
                line[1] = str(set_read_strand(adapter_summary))
                introns = intron_finder(read_start, cigar_list)
                re = ":".join(("re","i",str(read_end)))
                rl = ":".join(("rl","i",str(read_length)))
                l3 = ":".join(("l3","Z",",".join((str(i) for i in adapter_summary.get("l3")))))
                r3 = ":".join(("r3","Z",",".join((str(i) for i in adapter_summary.get("r3")))))
                l5 = ":".join(("l5","Z",",".join((str(i) for i in adapter_summary.get("l5")))))
                r3 = ":".join(("r5","Z",",".join((str(i) for i in adapter_summary.get("r5")))))
                intron = ":".join(("in","Z",",".join((str(i) for i in introns)),"\n"))
                newlyadded = "\t".join((re, rl, l3, r3, l5, r3, intron))
                new_sam_line = "\t".join(("\t".join(line[:12]), newlyadded))
                out_appender(new_sam_line)
    sam.close
    
if __name__== "__main__":
    main()

print("Modified sam file saved as {}".format(OUT_FILE))
#samfile.columns = ["read_name", "flag", "genome", "read_start", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual"]