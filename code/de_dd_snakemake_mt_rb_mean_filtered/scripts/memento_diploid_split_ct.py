import pandas as pd

var_path  = snakemake.input["var_combined"]
mean_path = snakemake.input["mean_combined"]
out_var   = snakemake.output["var"]
out_mean  = snakemake.output["mean"]
celltype  = snakemake.wildcards["celltype"]

def subset_ct(path, ct):
    df = pd.read_csv(path, index_col=0)
    # column format: sg^species^donor_id^celltype
    cols = [c for c in df.columns if c.split("^")[-1] == ct]
    return df[cols]

subset_ct(var_path,  celltype).to_csv(out_var)
subset_ct(mean_path, celltype).to_csv(out_mean)
