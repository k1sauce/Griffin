#!/usr/bin/env python
import sys
import argparse
import pandas as pd
import pysam
import pybedtools
import pyBigWig
import numpy as np


def define_fetch_interval(
    name_to_print,
    sites,
    chroms,
    chrom_sizes_path,
    upstream_bp,
    downstream_bp,
):
    # separate fw and reverse sites
    fw_markers = ["+", 1, "1"]
    rv_markers = ["-", -1, "-1"]
    fw_sites = sites[sites["Strand"].isin(fw_markers)].copy()
    rv_sites = sites[sites["Strand"].isin(rv_markers)].copy()

    undirected_sites = sites[~(sites["Strand"].isin(fw_markers + rv_markers))].copy()

    if len(rv_sites) + len(fw_sites) + len(undirected_sites) == len(sites):
        print(
            name_to_print
            + " (fw/rv/undirected/total): "
            + str(len(fw_sites))
            + "/"
            + str(len(rv_sites))
            + "/"
            + str(len(undirected_sites))
            + "/"
            + str(len(sites))
        )
    else:  # I don't think this should ever happen...
        print("total fw sites:\t\t" + str(len(fw_sites)))
        print("total rv sites:\t\t" + str(len(rv_sites)))
        print("total undirected sites:" + "\t" + str(len(undirected_sites)))
        print("total sites:\t\t" + str(len(sites)))
        sys.exit("Problem with strand column")

    # set up to fetch a window extending across the desired window
    fw_sites["fetch_start"] = fw_sites["position"] + upstream_bp
    fw_sites["fetch_end"] = fw_sites["position"] + downstream_bp

    undirected_sites["fetch_start"] = undirected_sites["position"] + upstream_bp
    undirected_sites["fetch_end"] = undirected_sites["position"] + downstream_bp

    # for reverse sites, flip the window
    rv_sites["fetch_start"] = rv_sites["position"] - downstream_bp
    rv_sites["fetch_end"] = rv_sites["position"] - upstream_bp

    # merge fw and reverse back together and sort them back into the original order
    sites = pd.concat([fw_sites, rv_sites, undirected_sites], ignore_index=True)
    sites = sites.sort_values(by=["Chrom", "position"]).reset_index(drop=True)

    chrom_sizes = pd.read_csv(chrom_sizes_path, sep="\t", header=None)
    chrom_sizes = chrom_sizes[chrom_sizes[0].isin(chroms)]
    chrom_sizes = chrom_sizes.set_index(0)

    adjusted_ends_df = pd.DataFrame()

    for chrom in chroms:
        length = chrom_sizes.loc[chrom][1]
        current = sites[sites["Chrom"] == chrom].copy()
        current["fetch_start"] = np.where(
            current["fetch_start"] < 0, 0, current["fetch_start"]
        )
        current["fetch_end"] = np.where(
            current["fetch_end"] > length, length, current["fetch_end"]
        )
        adjusted_ends_df = pd.concat([adjusted_ends_df, current], ignore_index=True)
    adjusted_ends_df = adjusted_ends_df.sort_values(
        by=["Chrom", "position"]
    ).reset_index(drop=True)
    adjusted_ends_df = adjusted_ends_df.copy()

    return adjusted_ends_df


def collect_fragments(input_list):
    i, chrom, start, end = input_list
    # open the bam file for each pool worker (otherwise individual pool workers can close it)
    bam_file = pysam.AlignmentFile(bam_path, "rb", index_filename=index_file_path)

    # open the ref seq
    ref_seq = pysam.FastaFile(ref_seq_path)

    # make dicts to hold the fetched positions
    columns = np.arange(start, end, 1)
    cov_dict = {m: 0 for m in columns}
    GC_cov_dict = {m: 0 for m in columns}

    # fetch reads
    fetched = bam_file.fetch(
        contig=chrom, start=start, stop=end
    )  # fetch reads that map to the region of interest

    ########################
    # count coverage
    ########################
    for read in fetched:
        # filter out reads
        if (
            abs(read.template_length) >= size_range[0]
            and abs(read.template_length) <= size_range[1]
            and read.is_paired == True
            and read.mapping_quality >= map_q
            and read.is_duplicate == False
            and read.is_qcfail == False
        ):
            # only use fw reads with positive fragment lengths (negative indicates an abnormal pair)
            # all paired end reads have a fw and rv read so we don't need the rv read to find the midpoint.
            if read.is_reverse == False and read.template_length > 0:
                fragment_start = (
                    read.reference_start
                )  # for fw read, read start is fragment start
                fragment_end = read.reference_start + read.template_length
                midpoint = int(np.floor((fragment_start + fragment_end) / 2))

                # count the GC content
                fragment_seq = ref_seq.fetch(
                    read.reference_name, fragment_start, fragment_end
                )
                fragment_seq = np.array(list(fragment_seq.upper()))
                fragment_seq[np.isin(fragment_seq, ["A", "T", "W"])] = 0
                fragment_seq[np.isin(fragment_seq, ["C", "G", "S"])] = 1
                rng = np.random.default_rng(fragment_start)
                fragment_seq[
                    np.isin(fragment_seq, ["N", "R", "Y", "K", "M", "B", "D", "H", "V"])
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

                # check that the site is in the window
                if midpoint >= start and midpoint < end:
                    # count the fragment
                    cov_dict[midpoint] += 1

                    ##get the GC bias
                    read_GC_content = sum(fragment_seq)
                    read_GC_bias = GC_bias[abs(read.template_length)][read_GC_content]

                    # count the fragment weighted by GC bias
                    if not np.isnan(read_GC_bias):
                        GC_cov_dict[midpoint] += 1 / read_GC_bias

                else:  # if fragment doesn't fully overlap
                    continue

                del (read, midpoint, fragment_seq)

            else:
                # print('reverse',read.is_reverse)
                continue

    output = pd.DataFrame(pd.Series(cov_dict, name="uncorrected"))
    output["GC_corrected"] = pd.Series(GC_cov_dict)
    output = output[
        output["uncorrected"] > 0
    ]  # don't waste memory on positions with no coverage
    output["chrom"] = chrom
    output = output.reset_index().rename(columns={"index": "position"})
    output["uncorrected"] = output["uncorrected"].astype(float)
    output["GC_corrected"] = np.round(output["GC_corrected"], 5)

    bam_file.close()
    ref_seq.close()

    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_name", help="name of sample", required=True)
    parser.add_argument("--bam", help="bam file", required=True)
    parser.add_argument(
        "--index_file_path",
        help="path to the index file for the bam file",
        required=True,
    )
    parser.add_argument(
        "--GC_bias", help="GC bias info from griffin_GC_bias", required=True
    )
    parser.add_argument(
        "--reference_genome", help="path to the reference genome", required=True
    )
    parser.add_argument(
        "--chrom_sizes_path", help="path to chrom sizes file", required=True
    )
    parser.add_argument("--site_file", help=".bed file of sites", required=True)
    parser.add_argument(
        "--chroms",
        help="chroms to include when selecting sites",
        nargs="*",
        default=[
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
        ],
    )
    parser.add_argument(
        "--norm_window",
        help="start and end of the window to be used for normalization",
        nargs=2,
        type=int,
        default=(-5000, 5000),
    )
    parser.add_argument(
        "--size_range",
        help="acceptable size range for fragments (to filter out genomic contamination)",
        nargs=2,
        type=int,
        default=(0, 500),
    )
    parser.add_argument(
        "--map_quality", help="minimum mapping quality", type=int, default=60
    )
    args = parser.parse_args()

    # assign arguments to variables
    sample_name = args.sample_name
    bam_path = args.bam
    index_file_path = args.index_file_path
    GC_bias_path = args.GC_bias
    ref_seq_path = args.reference_genome
    chrom_sizes_path = args.chrom_sizes_path
    site_file = args.site_file
    chroms = args.chroms
    norm_window = args.norm_window
    size_range = args.size_range
    map_q = args.map_quality

    ########################################
    # GET GC BIAS
    ########################################
    # 1) open the GC_bias file
    # 2) get rid of extremely low GC bias values, these fragments will now be excluded
    # these fragments are extremely rare so it is difficult to get a good estimate
    # of GC bias
    # 3) convert to a dictionary of {<length>: {<num_GC>: <bias>}}
    GC_bias = pd.read_csv(
        GC_bias_path, sep="\t", usecols=["length", "num_GC", "smoothed_GC_bias"]
    )
    GC_bias["smoothed_GC_bias"] = np.where(
        GC_bias["smoothed_GC_bias"] < 0.05, np.nan, GC_bias["smoothed_GC_bias"]
    )
    GC_bias = {
        length: group.set_index("num_GC")["smoothed_GC_bias"].to_dict()
        for length, group in GC_bias.groupby("length")
    }

    # import the sites from the BED file
    cols_to_use = ['Chrom', 'position', 'Strand']
    all_sites = pd.read_csv(
        site_file,
        sep="\t",
        usecols=lambda c: c in cols_to_use
    )
    # no strand information
    if 'Strand' not in all_sites.columns:
        all_sites['Strand'] = 0
    
    # throw out sites that aren't on the selected chroms
    all_sites = all_sites[all_sites["Chrom"].isin(chroms)]

    # number of bp to fetch upstream and downstream of the site
    upstream_bp = norm_window[0] - size_range[0]  # this should be negative
    downstream_bp = norm_window[1] + size_range[0]  # this should be positive
    all_sites = define_fetch_interval(
        "Total sites",
        all_sites,
        chroms,
        chrom_sizes_path,
        upstream_bp,
        downstream_bp,
    )

    # convert to pybedtools and merge overlapping segments
    all_sites_bed = pybedtools.BedTool.from_dataframe(
        all_sites[["Chrom", "fetch_start", "fetch_end"]]
    )
    all_sites_bed = all_sites_bed.sort()
    all_sites_bed = all_sites_bed.merge()
    print("Intervals to fetch:\t" + str(len(all_sites_bed)))
    print("Total bp to fetch:\t" + str(all_sites_bed.total_coverage()))
    sys.stdout.flush()

    # split the long intervals
    to_fetch = all_sites_bed.to_dataframe()
    to_fetch["length"] = to_fetch["end"] - to_fetch["start"]
    print("Max fetch length: " + str(to_fetch["length"].max()) + " bp")

    to_fetch = to_fetch[["chrom", "start", "end"]]
    to_fetch = to_fetch.sort_values(by=["chrom", "start"]).reset_index(drop=True)
    to_fetch = to_fetch.reset_index()  # add an index column
    split_len = pybedtools.BedTool.from_dataframe(
        to_fetch[["chrom", "start", "end"]]
    ).total_coverage()
    pybedtools.cleanup(verbose=True, remove_all=True)

    results = [collect_fragments(row) for row in to_fetch.values]

    chrom_sizes = pd.read_csv(chrom_sizes_path, sep="\t", header=None)
    chrom_sizes = chrom_sizes[chrom_sizes[0].isin(chroms)]

    uncorrected_bw = pyBigWig.open(f"{sample_name}.uncorrected.bw", "w")
    GC_bw = pyBigWig.open(f"{sample_name}.GC_corrected.bw", "w")

    uncorrected_bw.addHeader([(a, b) for a, b in chrom_sizes.values])
    GC_bw.addHeader([(a, b) for a, b in chrom_sizes.values])

    for i, current in enumerate(results):
        if len(current) > 0 and np.nansum(current["uncorrected"]) > 0:
            uncorrected_bw.addEntries(
                list(current["chrom"]),
                list(current["position"]),
                ends=list(current["position"] + 1),
                values=list(current["uncorrected"]),
            )

        if len(current) > 0 and np.nansum(current["GC_corrected"]) > 0:
            GC_bw.addEntries(
                list(current["chrom"]),
                list(current["position"]),
                ends=list(current["position"] + 1),
                values=list(current["GC_corrected"]),
            )

    uncorrected_bw.close()
    GC_bw.close()
