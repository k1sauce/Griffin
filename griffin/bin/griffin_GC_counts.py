#!/usr/bin/env python
import pysam
import pandas as pd
import numpy as np
import argparse
import sys


# Griffin special sauce, keep if you want to reproduce the
# original results, this function is creating a dict like this:
# {length: {num_GC: count}} for the reads in the bam file
# that overlap the mappable regions provided. Some reads are removed
# based on the specific filtering criteria.
#
# length is the absolute value of the template length
# num_GC is the number of GCs in the fragment
# count is the number of fragments with that length and GC content
#
# This is a dict of dicts, so the outer dict is the length
# and the inner dict is the number of GCs
def collect_reads(sublist):
    # create a dict for holding the frequency of each read length and GC content
    GC_dict = {}
    for length in range(size_range[0], size_range[1] + 1):
        GC_dict[length] = {}
        for num_GC in range(0, length + 1):
            GC_dict[length][num_GC] = 0

    # import the bam file
    # this needs to be done within the loop otherwise it gives a truncated file warning
    bam_file = pysam.AlignmentFile(bam_file_path, "rb", index_filename=index_file_path)

    # this might also need to be in the loop
    # import the ref_seq
    ref_seq = pysam.FastaFile(ref_seq_path)

    for i in range(len(sublist)):
        chrom = sublist.iloc[i][0]
        start = sublist.iloc[i][1]
        end = sublist.iloc[i][2]

        if i % 5000 == 0:
            print(f"Processing interval {i+1} of {len(sublist)}")
            sys.stdout.flush()
        # fetch any read that overlaps the inteterval
        fetched = bam_file.fetch(chrom, start, end)
        for read in fetched:
            # use both fw (positive template length) and rv (negative template length) reads
            if (
                read.is_reverse == False
                and read.template_length >= size_range[0]
                and read.template_length <= size_range[1]
            ) or (
                read.is_reverse == True
                and -read.template_length >= size_range[0]
                and -read.template_length <= size_range[1]
            ):
                # qc filters, some longer fragments are considered 'improper pairs' but I would like to keep these
                if (
                    read.is_paired == True
                    and read.mapping_quality >= map_q
                    and read.is_duplicate == False
                    and read.is_qcfail == False
                ):
                    if read.is_reverse == False:
                        fragment_start = read.reference_start
                        fragment_end = read.reference_start + read.template_length
                    elif read.is_reverse == True:
                        fragment_end = read.reference_start + read.reference_length
                        fragment_start = fragment_end + read.template_length

                    # count the GC content
                    fragment_seq = ref_seq.fetch(
                        read.reference_name, fragment_start, fragment_end
                    )
                    fragment_seq = np.array(list(fragment_seq.upper()))
                    fragment_seq[np.isin(fragment_seq, ["A", "T", "W"])] = 0
                    fragment_seq[np.isin(fragment_seq, ["C", "G", "S"])] = 1
                    rng = np.random.default_rng(fragment_start)
                    fragment_seq[
                        np.isin(
                            fragment_seq, ["N", "R", "Y", "K", "M", "B", "D", "H", "V"]
                        )
                    ] = rng.integers(
                        2,
                        size=len(
                            fragment_seq[
                                np.isin(
                                    fragment_seq,
                                    ["N", "R", "Y", "K", "M", "B", "D", "H", "V"],
                                )
                            ]
                        ),
                    )  # random integer in range(2) (i.e. 0 or 1)
                    fragment_seq = fragment_seq.astype(float)
                    num_GC = int(fragment_seq.sum())
                    GC_dict[abs(read.template_length)][num_GC] += 1

    return GC_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_file", help="sample_bam_file", required=True)
    parser.add_argument(
        "--bam_file_name",
        help="sample name (does not need to match actual file name)",
        required=True,
    )
    parser.add_argument(
        "--index_file_path",
        help="path to the index file for the bam file",
        required=True,
    )
    parser.add_argument(
        "--mappable_regions_path",
        help="highly mappable regions to be used in GC correction, bedGraph or bed foramt",
        required=True,
    )
    parser.add_argument(
        "--ref_seq", help="reference sequence (fasta format)", required=True
    )
    parser.add_argument(
        "--chrom_sizes",
        help="path to chromosome sizes for the reference seq",
        required=True,
    )
    parser.add_argument(
        "--map_q",
        help="minimum mapping quality for reads to be considered",
        type=int,
        required=True,
    )
    parser.add_argument(
        "--size_range",
        help="range of read sizes to be included",
        nargs=2,
        type=int,
        required=True,
    )
    parser.add_argument(
        "--CPU", help="number of CPU for parallelizing", type=int, required=True
    )
    args = parser.parse_args()

    bam_file_path = args.bam_file
    index_file_path = args.index_file_path
    bam_file_name = args.bam_file_name
    mappable_regions_path = args.mappable_regions_path
    ref_seq_path = args.ref_seq
    chrom_sizes_path = args.chrom_sizes
    map_q = args.map_q
    size_range = args.size_range
    CPU = args.CPU

    # import filter
    mappable_intervals = pd.read_csv(mappable_regions_path, sep="\t", header=None)

    # remove non standard chromosomes and X and Y
    chroms = ["chr" + str(m) for m in range(1, 23)]
    mappable_intervals = mappable_intervals[mappable_intervals[0].isin(chroms)]
    GC_dict = collect_reads(mappable_intervals)
    rows = [
        {"length": length, "num_GC": num_gc, "number_of_fragments": count}
        for length, gc_counts in GC_dict.items()
        for num_gc, count in gc_counts.items()
    ]
    GC_df = pd.DataFrame(rows)
    GC_df = GC_df.sort_values(by=["length", "num_GC"])
    GC_df.to_csv(f"{bam_file_name}.GC_counts.txt", sep="\t", index=False)
