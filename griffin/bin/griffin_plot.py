#!/usr/bin/env python
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sample", help="sample name", required=True)
parser.add_argument(
    "--uncorrected_coverage_tsv",
    help="path to uncorrected coverage tsv file",
    required=True,
)
parser.add_argument(
    "--gc_corrected_coverage_tsv",
    help="path to gc corrected coverage tsv file",
    required=True,
)
parser.add_argument(
    "--save_window",
    help="start and end of window to be plotted",
    nargs=2,
    type=int,
    default=(-1000, 1000),
)
parser.add_argument(
    "--step", help="step size when calculating coverage", type=int, default=5
)
parser.add_argument(
    "--individual",
    help="if individual sites were saved in previous steps. (True/False)",
    action="store_true",
)
args = parser.parse_args()

sample = args.sample
uncorrected_coverage_tsv = args.uncorrected_coverage_tsv
gc_corrected_coverage_tsv = args.gc_corrected_coverage_tsv
save_window = args.save_window
step = args.step
individual = args.individual

save_window = [
    int(np.ceil(save_window[0] / step) * step),
    int(np.floor(save_window[1] / step) * step),
]  # round to the nearest step inside the window
save_columns = np.arange(save_window[0], save_window[1], step)
str_save_columns = [str(m) for m in save_columns]
print("save_window rounded to step:", save_window)


# dict to hold results grouped by correction type
print("Importing data")
results_dict = {"uncorrected": pd.DataFrame(), "GC_corrected": pd.DataFrame()}

current = pd.read_csv(uncorrected_coverage_tsv, sep="\t")
if individual:
    current = current.groupby("site_name")[str_save_columns].mean()
    current = current.reset_index()  # retain site_name
    current["sample"] = sample

results_dict["uncorrected"] = pd.concat(
    [results_dict["uncorrected"], current], ignore_index=True
)

current = pd.read_csv(gc_corrected_coverage_tsv, sep="\t")
if individual:
    current = current.groupby("site_name")[str_save_columns].mean()
    current = current.reset_index()  # retain site_name
    current["sample"] = sample

results_dict["GC_corrected"] = pd.concat(
    [results_dict["GC_corrected"], current], ignore_index=True
)

site_names = results_dict["uncorrected"]["site_name"].unique()

# if info about individual sites was kept, the averaging process can take quite a while. Save for later use.
if individual:
    for i, key in enumerate(results_dict.keys()):
        data = results_dict[key].copy()
        data.to_csv(f"{key}.mean_data.txt", sep="\t", index=False)


# generate plots
num_plots = len(results_dict.keys())
for j, site_name in enumerate(site_names):
    fig, axes = plt.subplots(1, num_plots, figsize=(4 * num_plots, 3.5), sharey="row")
    for i, key in enumerate(results_dict.keys()):
        data = results_dict[key].copy()
        ax = axes[i]
        for sample in data["sample"].unique():
            current = data[
                (data["sample"] == sample) & (data["site_name"] == site_name)
            ]
            ax.plot(save_columns, current[str_save_columns].T, label=sample)
            ax.tick_params(labelleft=True)
        ax.set_title(site_name + " " + key)
        ax.set_xlabel("distance from site")

    axes[0].set_ylabel("normalized coverage")

    if len(data["sample"].unique()) < 15:
        axes[-1].legend(bbox_to_anchor=[1, 1], loc="upper left")
    else:
        axes[-1].legend(bbox_to_anchor=[1, 1], loc="upper left", ncol=2)

    fig.tight_layout()
    plt.savefig(f"{site_name}.pdf")
    plt.close("all")
