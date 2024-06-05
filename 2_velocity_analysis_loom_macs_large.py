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
cellID_obs_012 = sample_obs[sample_obs["x"].str.contains("_012")]
cellID_obs_1157 = sample_obs[sample_obs["x"].str.contains("_1157")]
cellID_obs_1195 = sample_obs[sample_obs["x"].str.contains("_1195")]
cellID_obs_766 = sample_obs[sample_obs["x"].str.contains("_766")]
cellID_obs_1652 = sample_obs[sample_obs["x"].str.contains("_1652")]
cellID_obs_1595 = sample_obs[sample_obs["x"].str.contains("_1595")]
cellID_obs_1579 = sample_obs[sample_obs["x"].str.contains("_1579")]
cellID_obs_1344 = sample_obs[sample_obs["x"].str.contains("_1344")]
cellID_obs_1338 = sample_obs[sample_obs["x"].str.contains("_1338")]
cellID_obs_0118 = sample_obs[sample_obs["x"].str.contains("_0118")]


gyr = ['#88CCEE','#CC6677', '#DDCC77','#117733','#332288','#AA4499','#44AA99','#999933','#882255','#661100','#6699CC','#888888','#CD1076','#2F4F4F']
# load files for spliced/unspliced matrices for each sample:
ldata012 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/012_large_vel/012_large.loom')
ldata1157 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1157_large_vel/1157_large.loom')
ldata1195 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1195_large_vel/1195_large.loom')
ldata766= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/766_large_vel/766_large.loom')
ldata1652= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1652_large_vel/1652_large.loom')
ldata1595= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1595_large_vel/1595_large.loom')
ldata1579= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1344_actually_1579_large_vel/1344_actually_1579_large.loom')
ldata1344= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1579_actually_1344_large_vel/1579_actually_1344_large.loom')
ldata1338= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1338_large_vel/1338_large.loom')
ldata0118= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/0118_large_vel/0118_large.loom')


# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata766.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_766' for bc in barcodes]
ldata766.obs.index = barcodes
ldata766 = ldata766[np.isin(ldata766.obs.index, cellID_obs_766["x"])]

barcodes = [bc.split(':')[1] for bc in ldata1652.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1652' for bc in barcodes]
ldata1652.obs.index = barcodes
ldata1652 = ldata1652[np.isin(ldata1652.obs.index, cellID_obs_1652["x"])]

barcodes = [bc.split(':')[1] for bc in ldata1579.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1579' for bc in barcodes]
ldata1579.obs.index = barcodes
ldata1579 = ldata1579[np.isin(ldata1579.obs.index, cellID_obs_1579["x"])]

barcodes = [bc.split(':')[1] for bc in ldata1595.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1595' for bc in barcodes]
ldata1595.obs.index = barcodes
ldata1595 = ldata1595[np.isin(ldata1595.obs.index, cellID_obs_1595["x"])]

barcodes = [bc.split(':')[1] for bc in ldata1344.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1344' for bc in barcodes]
ldata1344.obs.index = barcodes
ldata1344 = ldata1344[np.isin(ldata1344.obs.index, cellID_obs_1344["x"])]

barcodes = [bc.split(':')[1] for bc in ldata1338.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1338' for bc in barcodes]
ldata1338.obs.index = barcodes
ldata1338 = ldata1338[np.isin(ldata1338.obs.index, cellID_obs_1338["x"])]

barcodes = [bc.split(':')[1] for bc in ldata0118.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_0118' for bc in barcodes]
ldata0118.obs.index = barcodes
ldata0118 = ldata0118[np.isin(ldata0118.obs.index, cellID_obs_0118["x"])]

barcodes = [bc.split(':')[1] for bc in ldata012.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_012' for bc in barcodes]
ldata012.obs.index = barcodes
ldata012 = ldata012[np.isin(ldata012.obs.index, cellID_obs_012["x"])]

barcodes = [bc.split(':')[1] for bc in ldata1157.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1157' for bc in barcodes]
ldata1157.obs.index = barcodes
ldata1157 = ldata1157[np.isin(ldata1157.obs.index, cellID_obs_1157["x"])]

barcodes = [bc.split(':')[1] for bc in ldata1195.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_1195' for bc in barcodes]
ldata1195.obs.index = barcodes
ldata1195 = ldata1195[np.isin(ldata1195.obs.index, cellID_obs_1195["x"])]

# make variable names unique
ldata012.var_names_make_unique()
ldata1157.var_names_make_unique()
ldata1195.var_names_make_unique()
ldata766.var_names_make_unique()
ldata1652.var_names_make_unique()
ldata1579.var_names_make_unique()
ldata1595.var_names_make_unique()
ldata1344.var_names_make_unique()
ldata1338.var_names_make_unique()
ldata0118.var_names_make_unique()

ldatalarge = ldata012.concatenate([ldata1157, ldata1195,ldata766,ldata1652,ldata1579,ldata1595,ldata1344,ldata1338,ldata0118])

# merge matrices into the original adata object
adatalarge = scv.utils.merge(adata, ldatalarge)

#plot a umap
sc.pl.umap(adatalarge, color='seurat_clusters', frameon=False, legend_loc='on data', palette= gyr, title='', save='_macs_large_celltypes.pdf')

#inspect spliced and unspliced reads
scv.pl.proportions(adatalarge, groupby='seurat_clusters', save='_macs_large_celltypes.pdf')

# pre-processing
scv.pp.filter_and_normalize(adatalarge)
scv.pp.moments(adatalarge)

# compute velocity
#scv.tl.velocity(adata, mode='stochastic')
scv.tl.recover_dynamics(adatalarge, var_names='all')
scv.tl.velocity(adatalarge, mode='dynamical')
scv.tl.velocity_graph(adatalarge)

adatalarge.write('macs_large_velocity.h5ad', compression='gzip')
#adata = scv.read('macs_large_velocity.h5ad')

scv.pl.velocity_embedding(adatalarge, basis='umap', frameon=False, save='_macs_large_embedding.pdf')
scv.pl.velocity_embedding_grid(adatalarge, basis='umap', color='seurat_clusters', palette=gyr, save='_macs_large_embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adatalarge, basis='umap',color='seurat_clusters', legend_loc='none', palette=gyr,save='_macs_large_embedding_stream_May2024.pdf', title='')

