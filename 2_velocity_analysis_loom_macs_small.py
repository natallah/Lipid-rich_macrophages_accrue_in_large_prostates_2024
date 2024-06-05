import scvelo as scv
import scanpy as sc
#import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
#cr.settings.verbosity = 2

adata = sc.read_h5ad('macs_10_10.h5ad')

#read in cell IDs
sample_obs = pd.read_csv("macs_cellID_obs.csv")

#get cell IDs for each sample
cellID_obs_003 = sample_obs[sample_obs["x"].str.contains("_003")]
cellID_obs_004 = sample_obs[sample_obs["x"].str.contains("_004")]
cellID_obs_006 = sample_obs[sample_obs["x"].str.contains("_006")]
cellID_obs_007 = sample_obs[sample_obs["x"].str.contains("_007")]
cellID_obs_008 = sample_obs[sample_obs["x"].str.contains("_006")]
cellID_obs_009 = sample_obs[sample_obs["x"].str.contains("_090")]
cellID_obs_010 = sample_obs[sample_obs["x"].str.contains("_010")]
cellID_obs_013 = sample_obs[sample_obs["x"].str.contains("_013")]
cellID_obs_1144 = sample_obs[sample_obs["x"].str.contains("_1144")]
cellID_obs_1196 = sample_obs[sample_obs["x"].str.contains("_1196")]

# load files for spliced/unspliced matrices for each sample:
ldata003 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/003_small_vel/003_small.loom')
ldata004 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/004_small_vel/004_small.loom')
ldata006 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/006_small_vel/006_small.loom')
ldata007 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/007_small_vel/007_small.loom')
ldata008 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/008_small_vel/008_small.loom')
ldata009 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/009_small_vel/009_small.loom')
ldata010 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/010_small_vel/010_small.loom')
ldata013 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/013_small_vel/013_small.loom')
ldata1144 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1144_small_vel/1144_small.loom')
ldata1196 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1196_small_vel/1196_small.loom')


#add gene id

# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata003.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_003' for bc in barcodes]
ldata003.obs.index = barcodes
ldata003 = ldata003[np.isin(ldata003.obs.index, cellID_obs_003["x"])]

barcodes = [bc.split(':')[1] for bc in ldata004.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_004' for bc in barcodes]
ldata004.obs.index = barcodes
ldata004 = ldata004[np.isin(ldata004.obs.index, cellID_obs_004["x"])]

barcodes = [bc.split(':')[1] for bc in ldata006.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_006' for bc in barcodes]
ldata006.obs.index = barcodes
ldata006 = ldata006[np.isin(ldata006.obs.index, cellID_obs_006["x"])]

barcodes = [bc.split(':')[1] for bc in ldata007.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_007' for bc in barcodes]
ldata007.obs.index = barcodes
ldata007 = ldata007[np.isin(ldata007.obs.index, cellID_obs_007["x"])]

barcodes = [bc.split(':')[1] for bc in ldata008.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_008' for bc in barcodes]
ldata008.obs.index = barcodes
ldata008 = ldata008[np.isin(ldata008.obs.index, cellID_obs_008["x"])]

barcodes = [bc.split(':')[1] for bc in ldata009.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_090' for bc in barcodes]
ldata009.obs.index = barcodes
ldata009 = ldata009[np.isin(ldata009.obs.index, cellID_obs_009["x"])]

barcodes = [bc.split(':')[1] for bc in ldata010.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_010' for bc in barcodes]
ldata010.obs.index = barcodes
ldata010 = ldata010[np.isin(ldata010.obs.index, cellID_obs_010["x"])]

barcodes = [bc.split(':')[1] for bc in ldata013.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_013' for bc in barcodes]
ldata013.obs.index = barcodes
ldata013 = ldata013[np.isin(ldata013.obs.index, cellID_obs_013["x"])]

barcodes = [bc.split(':')[1] for bc in ldata1144.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1144' for bc in barcodes]
ldata1144.obs.index = barcodes
ldata1144 = ldata1144[np.isin(ldata1144.obs.index, cellID_obs_1144["x"])]

barcodes = [bc.split(':')[1] for bc in ldata1196.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1196' for bc in barcodes]
ldata1196.obs.index = barcodes
ldata1196 = ldata1196[np.isin(ldata1196.obs.index, cellID_obs_1196["x"])]

# make variable names unique
ldata003.var_names_make_unique()
ldata004.var_names_make_unique()
ldata006.var_names_make_unique()
ldata007.var_names_make_unique()
ldata008.var_names_make_unique()
ldata009.var_names_make_unique()
ldata010.var_names_make_unique()
ldata013.var_names_make_unique()
ldata1144.var_names_make_unique()
ldata1196.var_names_make_unique()

# concatenate the three loom
ldata = ldata003.concatenate([ldata004, ldata006,ldata007,ldata008,ldata009,ldata010,ldata013,ldata1144,ldata1196])

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)
adata.write('macs_small_velocity1.h5ad', compression='gzip')

# pre-processing
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
adata.write('macs_small_velocity2.h5ad', compression='gzip')
# compute velocity
#scv.tl.velocity(adata, mode='stochastic')
scv.tl.recover_dynamics(adata, var_names='all')
adata.write('macs_small_velocity3.h5ad', compression='gzip')
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

adata.write('macs_small_velocity.h5ad', compression='gzip')
