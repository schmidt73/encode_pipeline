import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def parse_args():
    p = argparse.ArgumentParser(
        description="Design gRNAs for cCRE regions."
    )
    
    p.add_argument(
        "candidate_ccre_grnas", 
        help="List of candidate gRNAs that overlap cCRE regions."
    )

    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()

    candidate_grnas = pd.read_csv(args.candidate_ccre_grnas)
    grnas_per_region = candidate_grnas.groupby("ccre_id").ccre_type.count()

    print(candidate_grnas.specificity)
    candidate_grnas_good_specificity = candidate_grnas[candidate_grnas.specificity >= 0.2]
    print(candidate_grnas_good_specificity)

    fig, ax = plt.subplots()
    sns.histplot(data=grnas_per_region, ax=ax)
    ax.set_xlabel("Number of gRNAs per cCRE")

    fig, ax = plt.subplots()
    grnas_good_specificity_per_region = candidate_grnas_good_specificity.groupby("ccre_id").ccre_type.count()
    sns.histplot(data=grnas_good_specificity_per_region, ax=ax)
    ax.set_xlabel("Number of gRNAs per cCRE\nwith >= 0.2 Specificity")

    plt.show()
