#!/usr/bin/env python3

import pandas as pd
import os
import sys
from scipy.stats import poisson
from argparse import ArgumentParser
from ast import literal_eval
from Bio import SeqIO


###############################################################################
###                           Common functions                              ###
###############################################################################

def get_cov(position, cov_dict):
    """
    Calculates coverage based on a cov_dict
    """
    if position in cov_dict:
        return cov_dict.get(position)
    else:
        return 0

def coverage(pos_list, args, contig):
    """
    Calculates average coverages from a given distance for a list of positions
    """
    cols = ["contig", "pos", "count"]
    cov = pd.read_csv(args.coverage_file, sep = "\t", names= cols)
    cov = cov.loc[cov.contig == contig].copy()
    cov["pos"] = cov["pos"].astype(int)
    cov["count"] = cov["count"].astype(int)
    cov.set_index("pos", drop=True, inplace=True)
    cov_dict = cov.to_dict()["count"]
    coverages = []
    for pos in pos_list:
        to_avg = []
        if args.distance > args.cov_sample:
            for position in range(pos + args.distance - args.cov_sample,
                                  pos + args.distance):
                to_avg.append(get_cov(position, cov_dict))
        elif args.distance < args.cov_sample:
            for position in range(pos + args.distance,
                                  pos + args.distance - args.cov_sample):
                to_avg.append(get_cov(position, cov_dict))
        else:
            to_avg.append(get_cov(pos + args.distance, cov_dict))
        coverages.append(sum(to_avg)/len(to_avg))
    return coverages

def check_if_qualified(df, minimum, ratio):
    """
    Checks whether the feature position satisfies the minimum count and minimum
    ratio of coverage requirements.
    """
    qual_list = []
    for index, row in df.iterrows():
        is_qual = (row["count"] >= minimum
                   and row["ratio"] >= ratio
                   and row["is_picked"])
        qual_list.append(is_qual)
    return qual_list

def find_features(args):
    """
    Reads in a dataframe from csv, chops it and processes it on contigs.
    """
    df = pd.read_csv(args.feature_file, sep = "\t", names = ["pos", "count"])
    df["pos"] = df["pos"].apply(literal_eval)
    df[["contig", "pos"]] = df["pos"].apply(pd.Series)
    contig_set = list(set(df.contig))
    new_df = ()
    for contig in contig_set:
        current_df = df.loc[df.contig == contig].copy()
        if args.feature == "in":
            current_df = contig_introns(current_df, args, contig)
        else:
            current_df = contig_ends(current_df, args, contig)
        if len(new_df) != 0:
            new_df = pd.concat([new_df, current_df])
        else:
            new_df = current_df
    if args.feature[1] == "3":
        feat = "_tes"
    elif args.feature[1] == "5":
        feat = "_tss"
    else:
        feat = "tron"
    new_df.to_csv(args.feature_file.replace(".tsv", "{}.tsv".format(feat)),
                  index=False,
                  sep="\t")

def Stats(args):
    """
    Sets argument types and runs stat functions for features.
    """
    if not args.feature:
        args.feature = args.feature_file[-6:-4]
    if args.feature == "r5" or args.feature == "r3":
        multiplier = -1
    else:
        multiplier = 1
    print("Calculating {} feature statistics...".format(args.feature))
    args.distance = abs(args.distance) * multiplier
    args.cov_sample = abs(args.cov_sample) * multiplier
    if os.stat(args.feature_file).st_size == 0:
        print("Feature file {} is empty. There is nothing to do here.".format(
              args.feature_file))
    else:
        find_features(args)


###############################################################################
###                             End functions                               ###
###############################################################################

def pick_from_greatests(dictionary, wobble):
    """
    Picks the left- or rightmost positions of the greatests list in a window
    determined by the wobble size. Whether the left or the rightmost positions
    are desired can be set by the user, and the list is ordered accordingly.
    """
    previous = -100
    is_picked_list = []
    for pos, is_greatest in dictionary.items():
        is_picked = False
        if is_greatest:
            if previous not in range(pos - wobble, pos + wobble + 1):
                is_picked = True
            previous = pos
        is_picked_list.append(is_picked)
    return is_picked_list

def check_if_greatest(tuples, wobble):
    """
    Finds the feature position with the highest read support in a window
    determined by the wobble size.
    """
    is_greatest_list = []
    for tup in tuples:
        pos, count = tup
        for position, count2 in tuples:
            is_greatest = True
            if position in range(pos - wobble, pos + wobble + 1):
                is_greatest = count >= count2
                if not is_greatest:
                    break
        is_greatest_list.append(is_greatest)
    return is_greatest_list

def count_average(tuples, window):
    """
    Counts the average values needed for the Poisson and the Pólya-Aeppli
    distributions.
    """
    in_window = []
    for tup in tuples:
        pos, count = tup
        hundred = 0
        for position, count2 in tuples:
            if position in range(pos - window, pos + window + 1):
                hundred += count2
        in_window.append(hundred/(2 * window + 1))
    return in_window

def contig_ends(df, args, contig):
    """
    Processes the dataframe for one contig looking for TSSs or TESs.
    """
    df["pos"] = df["pos"].astype(int)
    if args.feature == "l5" or args.feature == "l3":
    # This makes sure that the leftmost position is taken for each left feature
        df = df.sort_values(by = "pos")
    else:
        df = df.sort_values(by = "pos", ascending = False)
    df["coverage"] = coverage(df["pos"], args, contig)
    pos_count = list(zip(df["pos"], df["count"]))
    df["average"] = count_average(pos_count, args.window)
    df["is_greatest"] = check_if_greatest(pos_count, args.wobble)
    greatests_dict = dict(zip(df["pos"], df["is_greatest"]))
    df["is_picked"] = pick_from_greatests(greatests_dict, args.wobble)
    df["ratio"] = df["count"] / df["coverage"]
    df["is_qualified"] = check_if_qualified(df, args.minimum, args.ratio)
    df["poisp"] = (1 - poisson.cdf(df["count"], df["average"]))
    return df


###############################################################################
###                          Intronic functions                             ###
###############################################################################

def intron_picker(tuples, args):
    """
    Checks whether there is a more frequent intron near the intron and if there
    is, it only picks the intron if it is more abundant than the set 'rare 
    intron' limit.
    """
    is_picked_list = []
    wobble = args.intron_wobble
    for tup in tuples:
        is_picked = True
        left, right, count = tup
        poschecked = []
        for left2, right2, count2 in tuples:
            if left2 in range(left - wobble, left + wobble + 1):
                if right2 in range(right - wobble, right + wobble + 1):
                    poschecked.append((left2, right2, count2))
                if len(poschecked) > 1:
                    counts = [t[2] for t in poschecked]
                    is_picked = count / args.rare_intron >= max(counts)
        is_picked_list.append(is_picked)
    return is_picked_list
        
def check_consensus(left2, right2):
    """
    Finds consensus splice junctions and sets the strand of the intron
    accordingly.
    """
    consensus_list = []
    strand_list = []
    for i in range(len(left2)):
        if left2[i].lower() == "gt" and right2[i].lower() == "ag":
            consensus_list.append("GT/AG")
            strand_list.append("+")
        elif left2[i].lower() == "gc" and right2[i].lower() == "ag":
            consensus_list.append("GC/AG")
            strand_list.append("+")
        elif left2[i].lower() == "at" and right2[i].lower() == "ac":
            consensus_list.append("AT/AC")
            strand_list.append("+")
        elif left2[i].lower() == "ct" and right2[i].lower() == "ac":
            consensus_list.append("GT/AG")
            strand_list.append("-")
        elif left2[i].lower() == "ct" and right2[i].lower() == "gc":
            consensus_list.append("GC/AG")        
            strand_list.append("-")
        elif left2[i].lower() == "gt" and right2[i].lower() == "at":
            consensus_list.append("AT/AC")        
            strand_list.append("-")
        else:
            consensus_list.append("None")        
            strand_list.append(".")
    return consensus_list, strand_list


def get_score(scores, args):
    """
    Gets the alignment scores from the aligner and calculates the maximum
    alignment score between the two exon ends that could have lead to template
    switching.
    """
    pos = args.shs_for_ts + 1
    rscores = []
    score = 0
    while pos <= 2 * args.shs_for_ts:
        score += scores[pos - 1]
        rscores.append(score)
        pos += 1
    rmax = max(rscores)
    if rmax < 0:
        rmax = 0
    pos = args.shs_for_ts
    lscores = []
    score = 0
    while pos >= 1:
        score += scores[pos - 1]
        lscores.append(score)
        pos += -1
    lmax = max(lscores)
    if lmax < 0:
        lmax = 0
    return rmax + lmax

def align(lseq, rseq, args):
    """
    This is a special local aligner that does not allow gaps nor shifts. It
    reports the highest alignment score that borders or involves the centre.
    """
    scores = []
    for i in range(len(lseq)):
        if lseq[i] == rseq[i]:
            scores.append(args.match_score)
        else:
            scores.append(args.mismatch_score)
    return get_score(scores, args)

def intron_seq(df, args, contig):
    """
    Gets exon-intron border sequencences and sends them to the aligner.
    """
    ts = args.shs_for_ts
    for seq_record in SeqIO.parse(args.reference, "fasta"):
        if seq_record.name == contig:
            is_ts_list = []
            lseq_list = []
            rseq_list = []
            l2_list = []
            r2_list = []
            ts_list = []
            for index, row in df.iterrows():
                left = row["left"]
                right = row["right"] - 1
                leftseq = seq_record.seq[left - ts:left + ts]
                rightseq = seq_record.seq[right - ts:right + ts]
                ts_score = align(leftseq, rightseq, args)
                is_ts = ts_score >= args.match_score * ts
                is_ts_list.append(is_ts)
                ts_list.append(ts_score)
                lseq_list.append(leftseq)
                rseq_list.append(rightseq)
                l2_list.append(leftseq[ts:ts + 2])
                r2_list.append(rightseq[ts - 2:ts])
    try:
        df["is_potential_ts"] = is_ts_list
        df["leftseq"] = lseq_list
        df["rightseq"] = rseq_list
        df["left2"] = l2_list
        df["right2"] = r2_list
        df["ts_score"] = ts_list
    except:
        print("ERROR: The specified reference file does not contain contigs \
              from the references used for mapping. Supply the same reference\
              that was used for the mapping and try again.")
        sys.exit()
    return df

def contig_introns(df, args, contig):
    """
    Processes the dataframe for one contig looking for introns.
    """
    df[["left", "right"]] = df["pos"].apply(pd.Series)
    df["left"] = df["left"].astype(int)
    df["right"] = df["right"].astype(int)
    df["rcov"] = coverage(df["right"], args, contig)
    args.distance = abs(args.distance) * -1
    args.cov_sample = abs(args.cov_sample) * -1
    df["lcov"] = coverage(df["left"], args, contig)
    df["coverage"] = (df["lcov"] + df["rcov"]) / 2
    df["ratio"] = df["count"] / df["coverage"]
    df = intron_seq(df, args, contig)
    left_right_count = list(zip(df["left"], df["right"], df["count"]))
    df["is_picked"] = intron_picker(left_right_count, args)
    df["is_qualified"] = check_if_qualified(df, args.minimum, args.ratio)
    df["consensus"], df["strand"] = check_consensus(list(df["left2"]),
                                                    list(df["right2"]))
    args.distance = args.distance * -1
    args.cov_sample = args.cov_sample * -1
    return df


###############################################################################
###                             Main function                               ###
###############################################################################

def main():
    args = parsing()
    Stats(args)


def parsing():
    parser = ArgumentParser(description="This is the second module of \
                            LoRTIA, a Long-read RNA-Seq Transcript Isofom \
                            Annotator. This module calculates the statistics \
                            of transcript features.")
    parser.add_argument("coverage_file",
                        help="The tsv file which contains the coverages.\
                        The tsv file should contain 3 columns: contig, \
                        position and coverage.",
                        metavar="coverage_file")
    parser.add_argument("feature_file",
                        help="A tab-separated values file containing feature\
                        statistics produced by the Samprocessor.",
                        metavar="feature_file")
    parser.add_argument("-w", "--window",
                        dest="window",
                        help="The window that is examined when calculating \
                        the Poisson distribution. Setting low values finds \
                        false positives in a noisy data, while setting high \
                        values leads to false negatives due to the different \
                        transcriptional activity of different genomic regions.\
                        The default value is 50, which translates to a 101 nt\
                        bin (examined nucleotide +/- 50 nucleotides).",
                        type=int,
                        default=50,
                        metavar="[integer]")
    parser.add_argument("-r", "--reference",
                        dest="reference",
                        help="The reference fasta file. Template-switching \
                        in the case of putative introns is going to be checked\
                        according to this file.",
                        default="/mnt/c/Work/LT907985.2/Ref/LT907985.2.fasta",
                        metavar="[reference_fasta]")
    parser.add_argument("-m", "--minimum",
                        dest="minimum", 
                        help="The minimal number of reads for the feature to\
                        be accepted.",
                        type=int,
                        default=2,
                        metavar="[integer]")
    parser.add_argument("-f", "--feature",
                        dest="feature",
                        help="The feature that is examined. Options are \
                        'r5' for reverse strand 5' ends, 'l3' for \
                        reverse strand 3' ends, 'l5' for forward strand 5'\
                        ends, 'r3' for forward strand 3' ends and 'in' for \
                        introns. By default the tsv file's last two characters\
                        before the .tsv extension are considered.",
                        default=False,
                        metavar="[string]")
    parser.add_argument("-b", "--wobble",
                        dest="wobble",
                        help="The window, in which only one of each feature \
                        is expected, and locations with lesser support are \
                        considered to be derivatives of the major. The default\
                        value is 10, which means that only one feature of a \
                        kind can be described in a 21 nt bin (location +/-10 \
                        nt). This only applies to TSSs and TESs.",
                        type=int,
                        default=10,
                        metavar="[integer]")
    parser.add_argument("-i", "--intron_wobble",
                        dest="intron_wobble",
                        help="This option is only important for error-prone \
                        reads. Sequencing errors can disrupt the mapping of \
                        introns. Rare splice juntions can be detected in the \
                        close vicinity of more frequently utilized splice\
                        junctions. The rare splice junctions are likely to \
                        be results of sequencing errors of the more frequent \
                        version. This option regulates the window in which a \
                        rare intron will be considered to have stemmed from a \
                        sequencing error. The default value is 15 nt. That \
                        means that the rare introns which are no further than \
                        15 nt away from more frequent introns, will be \
                        considered to be sequencing errors.",
                        type=int,
                        default=15,
                        metavar="[integer]")
    parser.add_argument("--rare_intron",
                        dest="rare_intron",
                        help="This option is only important for error-prone \
                        reads. Sequencing errors can disrupt the mapping of \
                        introns. Rare splice juntions can be detected in the \
                        close vicinity of more frequently utilized splice\
                        junctions. The rare splice junctions are likely to \
                        be results of sequencing errors of the more frequent \
                        version. This option determines how much rarer an \
                        should be than the most frequent intron in its +/- \
                        'intron_wobble' vicinity, in orderd to be discarded as\
                        a sequencing error. The default value is 0.05.",
                        type=float,
                        default=0.05,
                        metavar="[float]")
    parser.add_argument("--match_score", 
                        dest="match_score",
                        help="The alignment scores for each match when \
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: 2",
                        type=float,
                        metavar="[float]", 
                        default=2.0,
                        required=False)
    parser.add_argument("--mismatch_score", 
                        dest="mismatch_score",
                        help="The alignment scores for each mismatch when \
                        searching for adapters. Penalty scores should be \
                        supplied as negative vaules. The default is: -3",
                        type=float,
                        metavar="[float]", 
                        default=-3.0,
                        required=False)
    parser.add_argument("--shs_for_ts", 
                        dest="shs_for_ts",
                        help="The minimum length of agreement (Short \
                        Homologous sequence = SHS) between the start of the \
                        match part of the alignment and the adapter that \
                        raises a suspicion of template switching. Putative\
                        template switching artefacts are listed in a \
                        separate file and are excluded form further \
                        statistics. The value has to be lesser than the value\
                        set by --check_in_match. If greater or equal value is\
                        set, the program will not look for signs of template \
                        switching. The default value is 3 nucleotides.",
                        type=int,
                        default=3,
                        metavar="[int]")
    parser.add_argument("-t", "--ratio",
                        dest = "ratio",
                        help = "The minimal ratio of the coverage that a \
                        feature has to reach to be accepted. The default value\
                        is 0.001.",
                        type=float,
                        default=0.001,
                        metavar="[float]")
    parser.add_argument("-d", "--distance",
                        dest="distance",
                        help="The distance from the feature position where \
                        coverage should be calculated. The default value is \
                        15. A positive value should be given, the inward \
                        direction is calculated by the program automatically.",
                        type=int,
                        default=15,
                        metavar="[integer]")
    parser.add_argument("-s", "--cov_sample",
                        dest="cov_sample",
                        help="The number of nucleotides where the coverage \
                        should be averaged. This many consecutive nucleotides\
                        will be considered from the 'distance' towards the \
                        feature. Its absolute value has to be smaller than \
                        or equal to the value of 'distance'. The default value\
                        is 5.",
                        type=int,
                        default=5,
                        metavar="[integer]")
    if not parser.parse_args().feature:
        parser.parse_args().feature = parser.parse_args().feature_file[-6:-4]
    return parser.parse_args()
    

if __name__== "__main__":
    main()
