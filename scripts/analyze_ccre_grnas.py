import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from loguru import logger

from tqdm import tqdm

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

    ccres_to_center = {}
    with open(args.ccres, 'r') as f:
        for line in f:
            chrm, start, end, rdhs_id, ccre_id, ccre_type = line.split()
            start, end = int(start), int(end)
            center = start + (end - start) // 2
            # ccres[rdhs_id] = (chrm, int(start), int(end), ccre_type, center)
            ccres_to_center[rdhs_id] = center

    logger.info("Loaded cCREs")
            
    candidate_grnas = pd.read_csv(args.candidate_ccre_grnas)#.set_index("ccre_id")
    logger.info("Loaded candidate gRNAs")

    forward_strand_rows = candidate_grnas['grna_strand'] == '+'

    candidate_grnas['cutting_position'] = candidate_grnas['grna_start']
    candidate_grnas.loc[forward_strand_rows, 'cutting_position'] += 17
    candidate_grnas.loc[~forward_strand_rows, 'cutting_position'] -= 6
    logger.info("Computed cutting position")

    candidate_grnas['ccre_center'] = candidate_grnas.ccre_id.map(ccres_to_center)
    candidate_grnas['ccre_distance'] = np.abs(candidate_grnas.ccre_center - candidate_grnas.cutting_position)
    logger.info("Computed distance to ccre center")

    candidate_grnas = candidate_grnas[candidate_grnas.specificity >= 0.2]
    candidate_grnas = candidate_grnas[~candidate_grnas.grna_sequence.map(has_monopolymer)]

    logger.info("Processing cCRE groups...")
    ccre_dfs = []
    for ccre_id, df in tqdm(candidate_grnas.groupby("ccre_id")):
        sorted_ccre_df = df.sort_values("ccre_distance", ascending=True).head(20)
        ccre_dfs.append(sorted_ccre_df)

    candidate_grnas = pd.concat(ccre_dfs)

    # fig, ax = plt.subplots()
    # grnas_per_region = candidate_grnas.groupby("ccre_id").ccre_type.count()
    # sns.histplot(data=grnas_per_region, ax=ax)
    # ax.set_xlabel("Number of gRNAs per cCRE")

    try:
        os.mkdir(args.output)
    except FileExistsError:
        pass

    # candidate_grnas.to_csv(f"{args.output}/unfiltered_grnas.csv")

    candidate_grnas.to_csv(f"{args.output}/filtered_grnas.csv")

    # fig, ax = plt.subplots()
    # grnas_per_region = candidate_grnas.groupby("ccre_id").ccre_type.count()
    # sns.histplot(data=grnas_per_region, ax=ax)
    # ax.set_xlabel("Number of gRNAs per cCRE after filtering")

    # plt.show()
