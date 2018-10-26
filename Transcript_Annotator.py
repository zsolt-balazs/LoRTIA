import pysam
import pandas as pd
from argparse import ArgumentParser
from ast import literal_eval

def read_ends(df, read, end, args, feature):
    """
    Checks whether the read has accepted TES or TSS features.
    """
    feature_found = False
    for index, row in df.iterrows():
        if end in range(row["start"] - args.wobble,
                        row["start"] + args.wobble + 1):
            feature_found = True
            position = row["start"]
            break
    if feature_found:
        read.set_tag(feature[:-1], position, "i")
    else:
        read.set_tag(feature[:-1], -1, "i")
    return read, feature_found

def intron_gff_iterator(df, intron, intron_found, real_introns, gap, args):
    """
    Iterates over the intron.gff3 and checks whether the given deletion is an
    intron.
    """
    df = df.loc[df["start"] > intron[0] - args.wobble - 1]
    for index, row in df.iterrows():
        if intron[0] == (row["start"]) and intron[1] == (row["end"] + 2):
            intron_found = True
            real_introns += intron
            break
    if not intron_found and ((intron[1] - intron[0]) > args.gap):
            gap = True
    return real_introns, intron_found, gap

def read_introns(df, read, introns, args):
    """
    Checks whether the read has accepted introns.
    """
    intron_found = False
    gap = False
    real_introns= ()
    if isinstance(introns[0], int):
        real_introns, intron_found, gap = intron_gff_iterator(df,
                                                              introns,
                                                              intron_found,
                                                              real_introns,
                                                              gap,
                                                              args)
    else:
        for intron in introns:
            real_introns, intron_found, gap = intron_gff_iterator(df,
                                                                  intron,
                                                                  intron_found,
                                                                  real_introns,
                                                                  gap,
                                                                  args)
    if intron_found:
        read.set_tag("in", str(real_introns), "Z")
    else:
        read.set_tag("in", None)
    return read, intron_found, gap

def bam_iterator(bam, tr_dict, tss_gff, tes_gff, intron_gff, outbam, args):
    """
    Iterates over reads in a bam file and sets the read tags for each found
    feature.
    """
    for read in bam:
        tss_found = False
        tes_found = False
        contig = read.reference_name
        if read.is_reverse:
            strand = "-"
            tss = "r5"
            tss_pos = read.reference_end
            tes = "l3"
            tes_pos = read.reference_start + 1
        else:
            strand = "+"
            tss = "l5"
            tss_pos = read.reference_start + 1
            tes = "r3"
            tes_pos = read.reference_end
        #GET THE TSS
        if (read.get_tag(tss).split(",")[3] == "correct" 
            or read.get_tag(tss).split(",")[3] =="potential template switching" 
            or read.get_tag(tes).split(",")[3] == "out of place"):
            df = tss_gff.loc[(tss_gff["strand"] == strand)
                             & (tss_gff["contig"] == contig)
                             & (tss_gff["start"] >= tss_pos - args.wobble - 1)]
            read, tss_found = read_ends(df, read, tss_pos, args, "tss")
        #GET THE TES
        if (read.get_tag(tes).split(",")[3] == "correct" 
            or read.get_tag(tes).split(",")[3] =="potential template switching" 
            or read.get_tag(tes).split(",")[3] == "out of place"):
            df = tes_gff.loc[(tes_gff["strand"] == strand) 
                             & (tes_gff["contig"] == contig) 
                             & (tes_gff["start"] >= tes_pos - args.wobble - 1)]
            read, tes_found = read_ends(df, read, tes_pos, args, "tes")
        #GET THE INTRON
        if read.get_tag("in"):
            introns = literal_eval(read.get_tag("in"))
            df = intron_gff.loc[(intron_gff["strand"] == strand)
                                & (intron_gff["contig"] == contig)]
            read, intron_found, gap = read_introns(df, read, introns, args)
        else:
            intron_found = False
        if tss_found and tes_found and not gap:
            if strand == "+":
                leftend = read.get_tag("ts")
                rightend = read.get_tag("te")
            else:
                leftend = read.get_tag("ts")
                rightend = read.get_tag("te")
            if intron_found:
                tr = contig, strand, leftend, read.get_tag("in"), rightend
            else:
                tr = contig, strand, leftend, rightend
            if tr in tr_dict:
                tr_dict[tr] += 1
            else:
                tr_dict[tr] = 1
            read.set_tag("tr", str(tr), "Z")
        outbam.write(read)
    bam.close()
    outbam.close()

def create_gff(tr_dict, tr_gff, args):
    """
    Generates the tr.gff3 file based on tr_dict
    """
    for key, value in tr_dict.items():
        if value >= args.mintr_count:

            if len(key) == 4:
                contig, strand, start, end = key
                l = len(tr_gff)
                if strand == "+":
                    tr_gff.loc[l] = [contig,
                                     "LoRTIA",
                                     "exon",
                                     start,
                                     end,
                                     value,
                                     strand,
                                     ".",
                                     "ID=" + str(key)]
                else:
                    tr_gff.loc[l] = [contig,
                                     "LoRTIA",
                                     "exon",
                                     end,
                                     start,
                                     value,
                                     strand,
                                     ".",
                                     "ID=" + str(key)]
            elif len(key) == 5:
                contig, strand, start, intron, end = key
                l = len(tr_gff)
                intron = literal_eval(intron)
                if strand == "+":
                    tr_gff.loc[l] = [contig,
                                     "LoRTIA",
                                     "exon",
                                     start,
                                     intron[0],
                                     value,
                                     strand,
                                     ".",
                                     "Parent=" + str(key)]
                    if len(intron) > 2:
                        for i in range(2, len(intron))[::2]:
                            l = len(tr_gff)
                            tr_gff.loc[l] = [contig,
                                             "LoRTIA",
                                             "exon",
                                             intron[i],
                                             intron[i+1],
                                             value,
                                             strand,
                                             ".",
                                             "Parent=" + str(key)]
                    l = len(tr_gff)
                    tr_gff.loc[l] = [contig,
                                     "LoRTIA",
                                     "exon",
                                     intron[-1],
                                     end,
                                     value,
                                     strand,
                                     ".",
                                     "Parent=" + str(key)]
                else:
                    tr_gff.loc[l] = [contig,
                                     "LoRTIA",
                                     "exon",
                                     end,
                                     intron[0],
                                     value,
                                     strand,
                                     ".",
                                     "Parent=" + str(key)]
                    if len(intron) > 2:
                        for i in range(2, len(intron))[::2]:
                            l = len(tr_gff)
                            tr_gff.loc[l] = [contig,
                                             "LoRTIA",
                                             "exon",
                                             intron[i],
                                             intron[i+1],
                                             value,
                                             strand,
                                             ".",
                                             "Parent=" + str(key)]
                    l = len(tr_gff)
                    tr_gff.loc[l] = [contig,
                                     "LoRTIA",
                                     "exon",
                                     intron[-1],
                                     start,
                                     value,
                                     strand,
                                     ".",
                                     "Parent=" + str(key)]

    tr_gff.to_csv(args.output_prefix + ".gff3",
                  sep="\t",
                  index = False,
                  header = False)

def annotate_tr(args):
    """
    Sets argument types and runs the bam_iterator
    """
    cols = ["contig",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "info"]
    tss_gff = pd.read_csv(args.gff_prefix + "_tss.gff3", sep="\t", names=cols)
    tes_gff = pd.read_csv(args.gff_prefix + "_tes.gff3", sep="\t", names=cols)
    intron_gff = pd.read_csv(args.gff_prefix + "_intron.gff3", sep="\t",
                             names=cols)
    bam = pysam.AlignmentFile(args.in_bam, "rb")
    tr_dict = {}
    outbam = pysam.AlignmentFile(args.output_prefix + ".bam",
                                 "wb",
                                 template=bam)
    bam_iterator(bam, tr_dict, tss_gff, tes_gff, intron_gff, outbam, args)
    pysam.index(args.output_prefix + ".bam")
    tr_gff = pd.DataFrame(columns=cols)
    print(tr_dict)
    create_gff(tr_dict, tr_gff, args)


###############################################################################
###                             Main function                               ###
###############################################################################

def main():
    args = parsing()
    annotate_tr(args)
    
def parsing():
    """
    This part handles the commandline arguments
    """
    parser = ArgumentParser(description="This is the first module of LoRTIA:\
                            a Long-read RNA-Seq Transcript Isofom Annotator")
    parser.add_argument("in_bam",
                        help="Input file. Both .sam and .bam files are\
                        accepted.",
                        metavar="input_file")
    parser.add_argument("gff_prefix",
                        help="The prefix of the gff files. Endings such as \
                        _tes, _tss and _intron will be assumed.",
                        metavar="gff_prefix")
    parser.add_argument("output_prefix",
                        help="Output prefix used for the .bam and the .gff \
                        files which are created by the program.)",
                        metavar="output_prefix")
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
    parser.add_argument("-g", "--gap",
                        dest="gap",
                        help="The largest allowed gap in the alignment with \
                        which it still constitutes as a transcript. Gaps are \
                        deletions which are not found in the gff of accepted \
                        introns. Such gaps can be present in chimeric reads \
                        which are often artefactual. The default value is \
                        1000.",
                        type=int,
                        default=1000,
                        metavar="[integer]")
    parser.add_argument("--mintr_count",
                        dest="mintr_count",
                        help="The minimum number of reads that has to \
                        support a transcript isoform for the isoform to be \
                        reported. The default is 1.",
                        type=int,
                        default=1,
                        metavar="[integer]")
    return parser.parse_args()


if __name__== "__main__":
    main()