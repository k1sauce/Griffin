#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse
from matplotlib import pyplot as plt
from pathlib import Path


def median_smoothing(current, fraction):
    bin_size = int(len(current) * fraction)
    if bin_size < 50:
        bin_size = 50
    medians = []

    for i in range(len(current)):
        start = int(i - bin_size / 2)
        end = int(i + bin_size / 2)
        # if the bin starts before the beginning, just take the first bin
        if start < 0:
            start = 0
            end = bin_size
        # if the bin extends beyond the end, take the last bin
        if end >= len(current):
            start = len(current) - bin_size
            end = len(current)
        current_median = np.nanmedian(current["GC_bias"].iloc[start:end])
        medians.append(current_median)
    return medians


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam_file_name",
        help="sample name (does not need to match actual file name)",
        required=True,
    )
    parser.add_argument(
        "--mappable_name",
        help="name of mappable regions file (with .bed removed)",
        required=True,
    )
    parser.add_argument(
        "--genome_GC_frequency",
        help="folder containing GC counts in the reference sequence (made by generate_reference_files.snakemake)",
        required=True,
    )
    parser.add_argument("--out_dir", help="folder for GC bias results", required=True)
    parser.add_argument(
        "--size_range",
        help="range of read sizes to be included",
        nargs=2,
        type=int,
        required=True,
    )
    args = parser.parse_args()

    bam_file_name = args.bam_file_name
    mappable_name = args.mappable_name
    genome_GC_frequency = args.genome_GC_frequency
    out_dir = args.out_dir
    size_range = args.size_range
    out_dir = out_dir.rstrip("/")
    # create output folders if needed
    Path(out_dir + "/GC_plots/").mkdir(parents=True, exist_ok=True)
    Path(out_dir + "/GC_bias/").mkdir(parents=True, exist_ok=True)
    # For now I'm going to keep the smoothing bin size as a set variable
    GC_smoothing_step = 20
    # input is the out_file from the previous step
    in_file = out_dir + "/GC_counts/" + bam_file_name + ".GC_counts.txt"
    # output is smoothed version
    smoothed_out_file = out_dir + "/GC_bias/" + bam_file_name + ".GC_bias.txt"
    # plot files
    plot_file1 = out_dir + "/GC_plots/" + bam_file_name + ".GC_bias.pdf"
    plot_file2 = out_dir + "/GC_plots/" + bam_file_name + ".GC_bias.summary.pdf"
    plot_file3 = out_dir + "/GC_plots/" + bam_file_name + ".GC_bias.key_lengths.pdf"
    # import the GC info from the genome
    frequency_prefix = genome_GC_frequency + "/" + mappable_name + "."

    # load the GC frequency data
    GC_freq = pd.DataFrame()
    for i in range(size_range[0], size_range[1] + 1):
        current_path = frequency_prefix + str(i) + "bp.GC_frequency.txt"
        current_data = pd.read_csv(current_path, sep="\t")
        GC_freq = pd.concat([GC_freq, current_data], ignore_index=True)
    GC_freq["GC_content"] = GC_freq["num_GC"] / GC_freq["length"]
    GC_freq = GC_freq.sort_values(by=["GC_content", "length"]).reset_index(drop=True)

    # Import GC counts from the sample
    GC_df = pd.read_csv(in_file, sep="\t")

    # Calculate the GC content of the sample
    GC_df["GC_content"] = GC_df["num_GC"] / GC_df["length"]
    GC_df = GC_df.sort_values(by=["GC_content", "length"]).reset_index(drop=True)
        
    # calculate the GC bias without a loop
    new_df = pd.merge(
        GC_df,
        GC_freq,
        how="left",
        left_on=["length", "GC_content"],
        right_on=["length", "GC_content"],
        suffixes=("", "_freq"),
    )
    new_df["GC_bias"] = (new_df["number_of_fragments"] / new_df["number_of_fragments_freq"])
    # calculate the GC bias for each length using groupby
    new_df["GC_bias"] = new_df.groupby("length")["GC_bias"].transform(
        lambda x: x / np.nanmean(x)
    )
    new_df = new_df.rename(columns={"number_of_fragments_freq": "number_of_positions"})
    # remove the num_GC_freq column
    new_df = new_df.drop(columns=["num_GC_freq"])
    new_df.sort_values(by=["GC_content", "length"], inplace=True)
    new_df.reset_index(drop=True, inplace=True)

    # smooth GC bias by size bin
    new_df2 = pd.DataFrame()
    for length in new_df["length"].unique():
        # get a bin of similar sized fragments
        min_len = int(length - (GC_smoothing_step / 2))
        max_len = int(length + (GC_smoothing_step / 2))
        current = new_df[
            (new_df["length"] >= min_len) & (new_df["length"] <= max_len)
        ].copy()

        # perform smoothing
        fit = median_smoothing(current, 0.05)
        current["smoothed_GC_bias"] = fit

        # only keep smoothed values for the selected length
        current = current[current["length"] == length]

        # get rid of values for GC contents that are never observed
        current["smoothed_GC_bias"] = np.where(
            current["number_of_positions"] == 0, np.nan, current["smoothed_GC_bias"]
        )

        # normalize to a mean of 1
        current["smoothed_GC_bias"] = current["smoothed_GC_bias"] / np.nanmean(
            current["smoothed_GC_bias"]
        )
        new_df2 = pd.concat([new_df2, current], ignore_index=True)

    new_df = new_df2

    # export results
    new_df2.to_csv(smoothed_out_file, sep="\t", index=False)

    # generate one plot per size bin
    # set up a figure for plotting
    plot_indexes = np.arange(
        size_range[0] + GC_smoothing_step,
        size_range[1] + GC_smoothing_step,
        GC_smoothing_step,
    )
    lengths_to_plot = plot_indexes
    x_dim = 6
    y_dim = int(np.ceil(len(plot_indexes) / 6))
    empty_plots = int(x_dim * y_dim - len(plot_indexes))
    plot_indexes = np.append(plot_indexes, [np.nan for m in range(empty_plots)])
    plot_indexes = np.reshape(plot_indexes, (y_dim, x_dim))
    fig, axes = plt.subplots(
        y_dim, x_dim, figsize=(5 * x_dim, 3.5 * y_dim), sharex=True, sharey=True
    )
    axes = axes.reshape(
        y_dim, x_dim
    )  # make sure the axes array is two dimension (just in case it has less than 7 value)

    # do the plotting
    min_len = 0
    for max_len in lengths_to_plot:
        if max_len % 20 == 0:
            print(max_len)

        # pick the axis
        current_index = np.where(plot_indexes == max_len)
        current_index = (current_index[0][0], current_index[1][0])
        current_ax = axes[current_index]

        # pick the data
        current1 = new_df2[
            (new_df2["length"] > min_len) & (new_df2["length"] <= max_len)
        ].copy()

        # plot the smoothed data over top
        for length2 in current1["length"].unique():
            current2 = current1[current1["length"] == length2]
            current_ax.plot(
                current2["GC_content"],
                current2["smoothed_GC_bias"],
                label=str(length2) + "bp",
            )

        current_ax.set_title(str(min_len) + "bp to " + str(max_len) + "bp")
        current_ax.legend(ncol=2)
        min_len = max_len

    for i in range(x_dim):
        axes[y_dim - 1, i].set_xlabel("GC content")

    for i in range(y_dim):
        axes[i, 0].set_ylabel("coverage bias")

    ylim = axes[0, 0].get_ylim()
    old_title = axes[0, 0].get_title()
    axes[0, 0].set_title(bam_file_name + "\n" + mappable_name + "\n" + old_title)
    fig.tight_layout()
    plt.savefig(plot_file1)
    plt.close("all")

    # key lengths
    selected_lengths = np.arange(100, 201, GC_smoothing_step)
    fig, ax = plt.subplots(1)
    for i, length in enumerate(selected_lengths):
        current = new_df2[new_df2["length"] == length]
        ax.plot(
            current["GC_content"], current["smoothed_GC_bias"], label=str(length) + "bp"
        )

    ax.legend(ncol=2, bbox_to_anchor=[1, 1], loc="upper left")
    ax.set_xlabel("GC content")
    ax.set_ylabel("coverage bias")
    ax.set_title(bam_file_name + "\n" + mappable_name)
    fig.tight_layout()
    fig.savefig(plot_file3)
    plt.close("all")

    # summary figure
    selected_lengths = np.arange(size_range[0], size_range[1], GC_smoothing_step)
    fig, ax = plt.subplots(1)
    for length in selected_lengths:
        current = new_df2[new_df2["length"] == length]
        ax.plot(
            current["GC_content"], current["smoothed_GC_bias"], label=str(length) + "bp"
        )
    ax.legend(ncol=2, bbox_to_anchor=[1, 1], loc="upper left")
    ax.set_xlabel("GC content")
    ax.set_ylabel("coverage bias")
    ax.set_title(bam_file_name + "\n" + mappable_name)
    fig.tight_layout()
    fig.savefig(plot_file2)
    plt.close("all")
