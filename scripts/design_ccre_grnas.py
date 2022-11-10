import argparse
import pandas as pd

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
    pass
