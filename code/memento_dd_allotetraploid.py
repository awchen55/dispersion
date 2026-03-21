import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import memento



# Read data and subset for hybrid lines
data_path = r"/project/gilad/awchen55/differentialDispersion/data/crossfilt_clustered_data/cell_annotations/"
data_path_out = r"/project/gilad/awchen55/differentialDispersion/data/crossfilt_clustered_data/memento_results/"

cell_types = ['ENS_glia','Glioblast','Neuroblast','Neuronal_IPC','Radial_glia']
q = 0.2

for ct in cell_types:
    # load data
    adata = sc.read(data_path + 'hybrid_crossfilt_' + ct + '.h5ad')
    # set chimp to 0 and human to 1
    adata.obs['species'] = adata.obs['species'].apply(lambda x: 0 if x == 'chimp' else 1)

    # run memento
    result_1d = memento.binary_test_1d(
            adata=adata, 
            capture_rate=q, 
            treatment_col='species', 
            num_cpus=12,
            num_boot=5000)
    
    # FDR correction
    result_1d["de_pval_adj"] = memento.util._fdrcorrect(result_1d["de_pval"])
    result_1d["dv_pval_adj"] = memento.util._fdrcorrect(result_1d["dv_pval"])
    
    # save data
    result_1d.to_pickle(data_path_out + 'hybrid_crossfilt_memento_results_' + ct + '.pkl')


