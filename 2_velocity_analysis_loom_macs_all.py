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

#766
#1652
#1579
#1595
#1344
#1338
#0118

# load files for spliced/unspliced matrices for each sample:
ldata766= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/766_large_vel/766_large.loom')
ldata1652= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1652_large_vel/1652_large.loom')
ldata1595= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1595_large_vel/1595_large.loom')
ldata1579= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1344_actually_1579_large_vel/1344_actually_1579_large.loom')
ldata1344= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1579_actually_1344_large_vel/1579_actually_1344_large.loom')
ldata1338= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1338_large_vel/1338_large.loom')
ldata0118= scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/0118_large_vel/0118_large.loom')
ldata003 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/003_small_vel/003_small.loom')
ldata004 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/004_small_vel/004_small.loom')
ldata006 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/006_small_vel/006_small.loom')
ldata007 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/007_small_vel//007_small.loom')
ldata008 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/008_small_vel/008_small.loom')
ldata009 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/009_small_vel/009_small.loom')
ldata010 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/010_small_vel/010_small.loom')
ldata013 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/013_small_vel/013_small.loom')
ldata1144 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1144_small_vel/1144_small.loom')
ldata1196 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1196_small_vel/1196_small.loom')
ldata012 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/012_large_vel/012_large.loom')
ldata1157 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1157_large_vel/1157_large.loom')
ldata1195 = scv.read('/depot/tlratlif/data/P20_2019/scRNA_seq/step6_velocityAnalysis/quantification/velocytoResults/1195_large_vel/1195_large.loom')

#rename barcodes in order to merge:
#766
#1652
#1579
#1595
#1344
#1338
#0118
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

ldata766.var_names_make_unique()
ldata1652.var_names_make_unique()
ldata1579.var_names_make_unique()
ldata1595.var_names_make_unique()
ldata1344.var_names_make_unique()
ldata1338.var_names_make_unique()
ldata0118.var_names_make_unique()
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
ldata012.var_names_make_unique()
ldata1157.var_names_make_unique()
ldata1195.var_names_make_unique()



ldata = ldata003.concatenate([ldata004, ldata006,ldata007,ldata008,ldata009,ldata010,ldata013,ldata1144,ldata1196,ldata012,ldata1157,ldata1195,ldata766,ldata1652,ldata1579,ldata1595,ldata1344,ldata1338,ldata0118])

# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

#plot a umap
gyr = ['#88CCEE','#CC6677', '#DDCC77','#117733','#332288','#AA4499','#44AA99','#999933','#882255','#661100','#6699CC','#888888','#CD1076','#2F4F4F']

sc.pl.umap(adata, color='seurat_clusters', frameon=False, legend_loc='on data', palette= gyr, title='', save='large_small_macs_celltypes.pdf')

#inspect spliced and unspliced reads
scv.pl.proportions(adata, groupby='seurat_clusters', save='largesmall_macs_celltypes.pdf')

# pre-processing
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# compute velocity
#scv.tl.velocity(adata, mode='stochastic')
scv.tl.recover_dynamics(adata, var_names='all')
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

adata.write('largesmall_macs_velocity.h5ad', compression='gzip')
#adata = scv.read('largesmall_macs_velocity.h5ad')

scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='largesmall_macs_embedding.pdf')
scv.pl.velocity_embedding_grid(adata, basis='umap', color='seurat_clusters', palette=gyr, save='largesmall_macs_embedding_grid.pdf', title='', scale=0.25)
scv.pl.velocity_embedding_stream(adata, legend_loc=none, basis='umap',color='seurat_clusters', palette=gyr,save='_macs_largesmall_embedding_streami_May2024.pdf', title='')



