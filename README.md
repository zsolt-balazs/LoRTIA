# LoRTIA - Long-read RNA-Seq Transcript Isoform Annotator toolkit
version 0.9.9

The LoRTIA toolkit annotates transcript features such as transcriptional start sites (TSS), transcriptional end sites (TES), introns and transcript isoforms based on long-read RNA sequencing (direct RNA-Seq and cDNA-Seq) data. The toolkit expects `SAM` or `BAM` files as input and produces `GFF` files as well as other processed files as outputs. The sequencing reads can stem from any long-read sequencing platform (both Nanopore and PacBio reads are accepted). As long as the reads contain adapters which mark the two ends of the full-length transcripts, a complete isoform annotation is possible. The `SAM` or `BAM` files can be produced by any mapper, however when mapping with [minimap2], either run minimap2 with the `-Y` option or filter out secondary alignments before using the LoRTIA toolkit. The toolkit was developed to be run in UNIX environments.

## Contents

- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Credits](#credits)

## <a name="installation"></a>Installation
Clone the repository.
```sh
git clone https://github.com/zsolt-balazs/LoRTIA
```
Make sure you have all the dependencies and you will be able to run the program.

## <a name="dependencies"></a>Dependencies
- [bedtools] (tested v. 2.28) Make sure that `bedtools` is added to your `PATH` and can be run by the `bedtools` command!
- [pandas] (tested v. 0.24.0 and up)
- [Biopython] (tested Release 1.73)
- [pysam] (tested v. 0.15.0)
- [scipy] (tested v. 1.2.2)

If you have your dependencies installed, and also added LoRTIA to your `PATH`, then typing `LoRTIA -h` should show you the help menu.

## <a name="usage"></a>Usage
The following examples will show how to determine isoforms from an `alignments.sam` file that was mapped to `reference.fasta`. The annotations and other processed files will be generated in the `output_folder` with the prefix of the `SAM` file.
When running the LoRTIA pipeline, you have to specify the types of adapters that you used during your library preparation.
The [Lexogen Telo Prime cap selection kit], is the default adapter set of the pipeline. If your sample was prepared with this kit, simply type the following command: 
```sh
LoRTIA -s poisson -f True /path/to/alignments.sam /path/to/output_folder /path/to/reference.fasta
```
If you prepared your library using the [PacBio Isoseq], type the following command: 
```sh
LoRTIA -5 AGAGTACATGGG --five_score 16 --check_in_soft 15 -3 AAAAAAAAAAAAAAA --three_score 18 -s poisson -f True /path/to/alignments.sam /path/to/output_folder /path/to/reference.fasta
```
If you prepared your library using standard [Nanopore cDNA-Seq adapters], type the following command: 
```sh
LoRTIA -5 TGCCATTAGGCCGGG --five_score 16 --check_in_soft 15 -3 AAAAAAAAAAAAAAA --three_score 16 -s poisson -f True /path/to/alignments.sam /path/to/output_folder /path/to/reference.fasta
```
If you prepared your library with a different set of adapters, you will have to specify those when running the program.

For more details, read the [Wiki].

## <a name="credits"></a>Credits
The LoRTIA toolkit is developed by [Zsolt Balázs]. Discussions and earlier work with Attila Szűcs (University of Szeged) contributed a lot of ideas to the implementation the toolkit. The early code was commented on by Tibor Nagy (University of Debrecen). Special thanks to the Department of Medical Biology at the University of Szeged for the early testing feedbacks.

[minimap2]: https://github.com/lh3/minimap2
[bedtools]: https://bedtools.readthedocs.io/en/latest/content/installation.html
[pandas]: https://pandas.pydata.org/pandas-docs/stable/install.html
[Biopython]: http://biopython.org/DIST/docs/install/Installation.html
[pysam]: https://pysam.readthedocs.io/en/latest/installation.html
[scipy]: https://www.scipy.org/install.html
[Lexogen Telo Prime cap selection kit]: https://www.lexogen.com/wp-content/uploads/2015/03/013PF032V0100_TeloPrime.pdf
[PacBio Isoseq]: https://www.pacb.com/blog/introduction-of-the-iso-seq-method-state-of-the-art-for-full-length-transcriptome-sequencing/
[Nanopore cDNA-Seq adapters]: https://nanoporetech.com/resource-centre/guide-cdna-sequencing-oxford-nanopore
[Wiki]: https://github.com/zsolt-balazs/LoRTIA/wiki
[Zsolt Balázs]: https://github.com/zsolt-balazs/
