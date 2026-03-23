import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import memento
from statsmodels.stats.multitest import fdrcorrection



# read in metadata with indv and rep info
data_path = r"/project/gilad/awchen55/differentialDispersion/data/crossfilt_clustered_data/cell_annotations/"
# metadata = pd.read_csv(data_path +'indv_rep_meta_data.csv')
# metadata = metadata.set_index('cell_id')
# 
# celltypes = ['ENS_glia', 'Glioblast',  'Neuroblast', 'Neuronal_IPC', 'Radial_glia']
# 
# for ct in celltypes:

# read data
adata = sc.read(data_path + 'parent_crossfilt_celltypes_of_interest.h5ad')

# add indv and rep data
#adata.obs = adata.obs.join(metadata, how = 'left')

# setup memento
adata.obs['capture_rate'] = 0.2
memento.setup_memento(adata, q_column='capture_rate')
memento.create_groups(adata, label_columns=['species','individual','replicate','celltype'])

# compute moments
memento.compute_1d_moments(adata,
min_perc_group=0.95) # percentage of groups that satisfy the condition for a gene to be considered. 

# get var
mean, var, counts = memento.get_1d_moments(adata)
# var = var.set_index('gene')
# var.columns = [col.split('^')[1] + '_' + col.split('^')[2] for col in var.columns]


var.to_csv(data_path + 'parent_crossfilt_memento_var_estimate_cot.csv')
mean.to_csv(data_path + 'parent_crossfilt_memento_mean_estimate_cot.csv')
# counts.to_csv(data_path + 'parent_crossfilt_memento_counts_estimate_cot.csv')


    




