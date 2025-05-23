#!/usr/bin/env python
import os
import sys
import argparse
import pandas as pd
import pyBigWig
import numpy as np
import time
from scipy.signal import savgol_filter
from scipy.stats import zscore


def fetch_bw_values(bw_path, current_sites, site_name, name):
    bw = pyBigWig.open(bw_path)
    fw_markers = ["+", 1, "1"]
    rv_markers = ["-", -1, "-1"]
    results = pd.DataFrame(
        np.zeros([len(current_sites), norm_window[1] - norm_window[0]])
    )
    start_time = time.time()
    for i in range(len(current_sites)):
        chrom, start, end, strand = current_sites.iloc[i][
            ["Chrom", "fetch_start", "fetch_end", "Strand"]
        ]

        values = bw.values(chrom, start, end, numpy=True)
        values = np.nan_to_num(
            values
        )  # turn nan into zero because bw doesn't store zero

        # if the window extends beyond the end of the chromosome add np.nan to fill the gap
        if len(values) < (norm_window[1] - norm_window[0]):
            ###################
            position = current_sites.iloc[i]["position"]
            array_start = position + norm_window[0]
            array_end = position + norm_window[1]

            temp_series = pd.Series(
                np.full(norm_window[1] - norm_window[0], np.nan),
                index=np.arange(array_start, array_end),
            )
            temp_series[np.arange(start, end)] = values

            values = temp_series.values
            del (temp_series, position, array_start, array_end)
            ###################
        if strand in rv_markers:
            values = values[::-1]
        results.iloc[i] = pd.Series(values)

    results.columns = np.arange(norm_window[0], norm_window[1])
    return results


def sum_bins(results, name):
    # If a bin has an np.nan value, the whole bin will become np.nan
    summed = np.sum(
        results.values.reshape(len(results), int(len(results.columns) / step), step),
        axis=2,
    )
    summed = pd.DataFrame(summed)
    summed.columns = norm_columns
    return summed


def exclude_regions(results, excluded_regions, name):
    results = np.where(
        excluded_regions > 0, np.nan, results
    )  # if any bp in the bin were excluded (1) exclude the bin
    results = pd.DataFrame(results)
    results.columns = norm_columns
    return results


def exclude_zero_mappability(results, mappability_values, name):
    results = np.where(
        mappability_values > 0, results, np.nan
    )  # only retain positions where the mappability is >0
    results = pd.DataFrame(results)
    results.columns = norm_columns
    return results


def make_outlier_mask(results, site_name):
    max_value = results.max().max()
    min_cutoff = 2  # minimum coverage that must be retained even if it is an outlier
    print(site_name, "max_bin_coverage is", max_value, "midpoints")

    scores = pd.DataFrame(zscore(results.values, axis=None, nan_policy="omit"))
    outlier_mask = pd.DataFrame(np.where(scores < 10, 1, np.nan))
    outlier_mask.columns = norm_columns

    if (results * outlier_mask).max().max() < min_cutoff:
        print(
            site_name,
            "low coverage, resetting the outlier cutoff to " + str(min_cutoff),
        )
        outlier_mask = pd.DataFrame(np.where(results <= min_cutoff, 1, np.nan))
        outlier_mask.columns = norm_columns

    outlier_cutoff = (results * outlier_mask).max().max()
    print(site_name, "masking sites with  >", outlier_cutoff, "midpoints")
    return (outlier_mask, outlier_cutoff)


def normalize_and_smooth(results,site_name,name, norm_columns, save_columns):
    #get the mean midpoints per valid position in each site
    mean_reads_per_bp_in_normalization_window = np.nanmean(results[norm_columns],axis = 1)/step
    mean_reads_per_bp_in_saved_window = np.nanmean(results[save_columns],axis = 1)/step
        
    #normalize individual sites to 1 to remove CNA
    if CNA_normalization.lower() == 'true':
        print(site_name,name,'normalizing CNAs')
        mean_data = np.nanmean(results.values,axis = 1, keepdims=True)
        #replace zero with nan so there aren't any infinities in the output
        mean_data = np.where(mean_data==0,np.nan,mean_data)
        results[norm_columns] = results[norm_columns]/mean_data  
        
    #take the mean of all sites
    if not individual.lower()=='true':
        print(site_name,name,'averaging sites')
        results = pd.DataFrame(pd.Series(np.nanmean(results[norm_columns], axis = 0), index=norm_columns)).T
        results.columns = norm_columns
        mean_reads_per_bp_in_normalization_window = np.nanmean(mean_reads_per_bp_in_normalization_window)
        mean_reads_per_bp_in_saved_window = np.nanmean(mean_reads_per_bp_in_saved_window)
        
    #smooth the sites
    if smoothing.lower()=='true':
        print(site_name,name,'smoothing')
        #savgol window should be approx one fragment length but it must be odd
        savgol_window=np.floor(smoothing_length/step)
        if savgol_window%2==0:
            savgol_window=savgol_window+1
        savgol_window=int(savgol_window)

        results[norm_columns] = savgol_filter(results[norm_columns], savgol_window, 3)
    
    #normalize the average site to 1
    print(site_name,name,'correcting for read depth')
    mean_value = np.nanmean(results[norm_columns])
    results[norm_columns] = results[norm_columns]/mean_value
    
    #save only plot columns 
    results = results[save_columns].copy()
    
    results['mean_reads_per_bp_in_normalization_window'] = mean_reads_per_bp_in_normalization_window
    results['mean_reads_per_bp_in_saved_window'] = mean_reads_per_bp_in_saved_window
    
    return(results)


def calculate_features(results):
    results["mean_coverage"] = results[save_columns].mean(
        axis=1, skipna=False
    )  # pandas skips na by default
    results["central_coverage"] = results[center_columns].mean(axis=1, skipna=False)
    fft_res = np.fft.fft(results[fft_columns])
    results["amplitude"] = np.abs(fft_res[:, fft_index])
    return results


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


def merge_sites(site_file, chroms):
    # get the site lists and define the fetch interval
    # import the sites from the BED file
    cols_to_use = ["Chrom", "position", "Strand"]
    all_sites = pd.read_csv(site_file, sep="\t", usecols=lambda c: c in cols_to_use)
    # no strand information
    if "Strand" not in all_sites.columns:
        all_sites["Strand"] = 0

    # throw out sites that aren't on the selected chroms
    all_sites = all_sites[all_sites["Chrom"].isin(chroms)]
    # throw out sites that aren't on the selected chroms
    all_sites = all_sites[all_sites["Chrom"].isin(chroms)]

    all_sites = define_fetch_interval(
        "Total sites",
        all_sites,
        chroms,
        chrom_sizes_path,
        norm_window[0],
        norm_window[1],
    )
    # dict to hold results
    results_dict = results_dict_template.copy()

    # fetch coverage and sum into bins of length step
    for key in results_dict.keys():
        results_dict[key]["coverage"] = fetch_bw_values(
            results_dict[key]["input_path"], all_sites, site_name, key
        )
        results_dict[key]["coverage"] = sum_bins(results_dict[key]["coverage"], key)

    # exclude specified regions
    if exclude_path:
        regions_to_exclude = fetch_bw_values(
            "excluded_regions.bw", 
            all_sites, 
            site_name, 
            "to_exclude"
        )
        regions_to_exclude = sum_bins(regions_to_exclude, "to_exclude")
        for key in results_dict.keys():
            results_dict[key]["coverage"] = exclude_regions(
                results_dict[key]["coverage"], regions_to_exclude, key
            )
        del regions_to_exclude

    # exclude zero mappability
    if exclude_zero_mappability_parameter.lower() == "true":
        print(site_name + " - excluding zero mappability.")
        # fetch excluded_regions
        mappability_values = fetch_bw_values(
            mappability_bw, all_sites, site_name, "mappabilty"
        )
        # replace zero with np.nan for mappability
        # when summing bins, any bin with one or more zeros will now be np.nan
        mappability_values[all_positions] = np.where(
            mappability_values[all_positions] == 0,
            np.nan,
            mappability_values[all_positions],
        )
        mappability_values = sum_bins(mappability_values, "mappability")

        for key in results_dict.keys():
            results_dict[key]["coverage"] = exclude_zero_mappability(
                results_dict[key]["coverage"], mappability_values, key
            )
        del mappability_values

    # mask out bins with coverage >10 SD above the mean
    if exclude_outliers_parameter.lower() == "true":
        print(site_name + " - excluding outliers.")
        outlier_mask, outlier_cutoff = make_outlier_mask(
            results_dict["uncorrected"]["coverage"], site_name
        )
        sys.stdout.flush()

        for key in results_dict.keys():
            results_dict[key]["coverage"] = results_dict[key]["coverage"] * outlier_mask
    else:
        outlier_cutoff = "NA"

    
    # normalize to a mean of 1 and smooth
    #normalize to a mean of 1 and smooth
    for key in results_dict.keys():
        results_dict[key]['coverage'] = normalize_and_smooth(results_dict[key]['coverage'],site_name,key,norm_columns,save_columns)

    # get features
    for key in results_dict.keys():
        results_dict[key]["coverage"] = calculate_features(
            results_dict[key]["coverage"]
        )

    # get metadata
    for key in results_dict.keys():
        results_dict[key]["coverage"]["outlier_cutoff"] = outlier_cutoff
        results_dict[key]["coverage"][
            "exclude_zero_mappability"
        ] = exclude_zero_mappability_parameter
        results_dict[key]["coverage"]["correction"] = key
        results_dict[key]["coverage"]["number_of_sites"] = len(all_sites)
        results_dict[key]["coverage"]["site_name"] = site_name
        results_dict[key]["coverage"]["smoothing"] = smoothing
        results_dict[key]["coverage"]["CNA_normalization"] = CNA_normalization
        results_dict[key]["coverage"]["sample"] = sample_name
        results_dict[key]["coverage"] = results_dict[key]["coverage"].copy()

    # if saving individual sites, keep the locations
    if individual.lower() == "true":
        current_sites = current_sites.drop(columns=["site_name"])
        for key in results_dict.keys():
            results_dict[key]["coverage"] = results_dict[key]["coverage"].merge(
                current_sites, left_index=True, right_index=True, validate="one_to_one"
            )

    return results_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_name", help="name of sample", required=True)
    parser.add_argument(
        "--uncorrected_bw_path",
        help="uncorrected bigWig from griffin_coverage",
        required=True,
    )
    parser.add_argument(
        "--GC_corrected_bw_path",
        help="GC_corrected bigWig from griffin_coverage",
        required=True,
    )
    parser.add_argument(
        "--mappability_bw",
        help="bigWig file of genome wide mappability scores",
        required=True,
    )
    parser.add_argument(
        "--chrom_sizes_path", help="path to chrom sizes file", required=True
    )
    parser.add_argument("--site_file", help=".bed file of sites", required=True)
    parser.add_argument("--site_name", help="name of the site file", required=True)
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
        "--save_window",
        help="start and end of the window to be saved in the outputs",
        nargs=2,
        type=int,
        default=(-1000, 1000),
    )
    parser.add_argument(
        "--center_window",
        help="start and end of the window to be used for calculating the central coverage feature",
        nargs=2,
        type=int,
        default=(-1000, 1000),
    )
    parser.add_argument(
        "--fft_window",
        help="start and end of the window to be used for calculating the amplitude feature (default is for bin size 15)",
        nargs=2,
        type=int,
        default=(-960, 960),
    )
    parser.add_argument(
        "--fft_index",
        help="index of the fft component to be saved as the amplitude feature (default is for 15bp bins and -960 to 960 fft_window)",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--smoothing_length",
        help="window length for the Savitzky-Golay smoothing function (should be approximately the mean fragment length)",
        type=int,
        default=165,
    )
    parser.add_argument(
        "--exclude_path",
        help='path to bed files of regions to filter out, excluded regions, centromeres, gaps, patches, alternative haplotypes',
    )
    parser.add_argument(
        "--step", help="step size when calculating coverage", type=int, default=5
    )
    parser.add_argument(
        "--CNA_normalization",
        help="whether to normalize each site individually to the copy number within the normalization window",
        default="False",
        required=True,
    )
    parser.add_argument(
        "--individual",
        help="save individual site coverage. TRUE WILL RESULT IN HUGE OUTPUT FILES. (True/False)",
        default="False",
        required=True,
    )
    parser.add_argument(
        "--smoothing",
        help="whether to use a savgol filter to smooth sites (True/False)",
        default="True",
        required=True,
    )
    parser.add_argument(
        "--exclude_outliers",
        help="whether to exclude bins with extreme outlier coverage >10SD above the mean (True/False)",
        default="True",
        required=True,
    )
    parser.add_argument(
        "--exclude_zero_mappability",
        help="whether to exclude bins with zero mappability (True/False)",
        default="True",
        required=True,
    )
    parser.add_argument(
        "--number_of_sites", help="number of sites to analyze", default="none"
    )
    args = parser.parse_args()

    sample_name = args.sample_name
    uncorrected_bw_path = args.uncorrected_bw_path
    GC_corrected_bw_path = args.GC_corrected_bw_path
    mappability_bw = args.mappability_bw
    chrom_sizes_path = args.chrom_sizes_path

    chroms = args.chroms
    norm_window = args.norm_window
    save_window = args.save_window
    center_window = args.center_window
    fft_window = args.fft_window
    fft_index = args.fft_index
    smoothing_length = args.smoothing_length
    exclude_path = args.exclude_path
    step = args.step
    CNA_normalization = args.CNA_normalization
    individual = args.individual
    smoothing = args.smoothing
    exclude_outliers_parameter = args.exclude_outliers
    exclude_zero_mappability_parameter = args.exclude_zero_mappability
    number_of_sites = args.number_of_sites
    site_file = args.site_file
    site_name = args.site_name
    
    results_dict_template = {
        "uncorrected": {"input_path": uncorrected_bw_path},
        "GC_corrected": {"input_path": GC_corrected_bw_path},
    }
    norm_window = [
        int(np.ceil(norm_window[0] / step) * step),
        int(np.floor(norm_window[1] / step) * step),
    ]  # round to the nearest step inside the window
    save_window = [
        int(np.ceil(save_window[0] / step) * step),
        int(np.floor(save_window[1] / step) * step),
    ]  # round to the nearest step inside the window
    center_window = [
        int(np.ceil(center_window[0] / step) * step),
        int(np.floor(center_window[1] / step) * step),
    ]  # round to the nearest step inside the window
    fft_window = [
        int(np.ceil(fft_window[0] / step) * step),
        int(np.floor(fft_window[1] / step) * step),
    ]  # round to the nearest step inside the window
    all_positions = np.arange(norm_window[0], norm_window[1])
    norm_columns = np.arange(norm_window[0], norm_window[1], step)
    save_columns = np.arange(save_window[0], save_window[1], step)
    center_columns = np.arange(center_window[0], center_window[1], step)
    fft_columns = np.arange(fft_window[0], fft_window[1], step)
    smoothing_length = int(
        np.round(smoothing_length / step) * step
    )  # round fragment length to the nearest step

    if exclude_path:
        # read in excluded regions, only keep the ones on the selected chromosomes
        merged_exclude_regions = pd.read_csv(exclude_path, sep='\t', header=0)
        merged_exclude_regions = merged_exclude_regions[
            merged_exclude_regions["chrom"].isin(chroms)
        ]
        # convert to bw, add header and add entries
        excluded_regions_bw = pyBigWig.open("excluded_regions.bw", "w")
        chrom_sizes = pd.read_csv(chrom_sizes_path, sep="\t", header=None)
        chrom_sizes = chrom_sizes[chrom_sizes[0].isin(chroms)]
        excluded_regions_bw.addHeader([(a, b) for a, b in chrom_sizes.values])
        excluded_regions_bw.addEntries(
            list(merged_exclude_regions["chrom"]),
            list(merged_exclude_regions["start"]),
            ends=list(merged_exclude_regions["end"]),
            values=[1.0 for m in range(len(merged_exclude_regions))],
        )
        excluded_regions_bw.close()

    results = merge_sites(site_file, chroms)

    for key in results_dict_template.keys():
        current_results = pd.DataFrame()
        current_results = pd.concat(
            [current_results, results[key]["coverage"]], ignore_index=True
        )
        current_results.to_csv(
            f'{sample_name}.{key}.coverage.tsv',
            sep="\t", 
            index=False, 
            float_format="%.5f"
        )

    if exclude_path:
        os.remove("excluded_regions.bw")
