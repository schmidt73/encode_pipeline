import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

def parse_args():
    p = argparse.ArgumentParser(
        description="Analyze and design gRNAs for cCRE regions."
    )
    
    p.add_argument(
        "ccres", 
        help="BED file containing cCREs."
    )

    p.add_argument(
        "candidate_ccre_grnas", 
        help="List of candidate gRNAs that overlap cCRE regions."
    )

    p.add_argument(
        "--output",
        help="Output directory.",
        default="candidate_regions",
    )

    return p.parse_args()

def has_monopolymer(sequence):
    for N in list("ATCG"):
        if (N * 4) in sequence:
            return True
    return False

if __name__ == "__main__":
    args = parse_args()

    ccres = {}
    with open(args.ccres, 'r') as f:
        for line in f:
            chrm, start, end, rdhs_id, ccre_id, ccre_type = line.split()
            start, end = int(start), int(end)
            center = start + (end - start) // 2
            ccres[rdhs_id] = (chrm, int(start), int(end), ccre_type, center)
            
    candidate_grnas = pd.read_csv(args.candidate_ccre_grnas)
    candidate_grnas['cutting_position'] = candidate_grnas.apply(
        lambda r:  r.grna_start + 17 if r.grna_strand == '+' else r.grna_start - 6,
        axis=1
    )

    candidate_grnas['ccre_center'] = candidate_grnas.ccre_id.apply(
        lambda ccre_id: ccres[ccre_id][4]
    )

    candidate_grnas['ccre_distance'] = np.abs(candidate_grnas.ccre_center - candidate_grnas.cutting_position)
    candidate_grnas = candidate_grnas.sort_values(["ccre_id", "ccre_distance"], ascending=True).groupby("ccre_id").apply(
        lambda df: df.head(20).drop(columns=["ccre_id"])
    ).reset_index().drop(columns=["level_1"])

    print(candidate_grnas)

    fig, ax = plt.subplots()
    grnas_per_region = candidate_grnas.groupby("ccre_id").ccre_type.count()
    sns.histplot(data=grnas_per_region, ax=ax)
    ax.set_xlabel("Number of gRNAs per cCRE")

    try:
        os.mkdir(args.output)
    except FileExistsError:
        pass

    candidate_grnas.to_csv(f"{args.output}/unfiltered_grnas.csv")

    candidate_grnas = candidate_grnas[candidate_grnas.specificity >= 0.2]
    candidate_grnas = candidate_grnas[~candidate_grnas.grna_sequence.apply(has_monopolymer)]

    candidate_grnas.to_csv(f"{args.output}/filtered_grnas.csv")

    fig, ax = plt.subplots()
    grnas_per_region = candidate_grnas.groupby("ccre_id").ccre_type.count()
    sns.histplot(data=grnas_per_region, ax=ax)
    ax.set_xlabel("Number of gRNAs per cCRE after filtering")

    plt.show()
