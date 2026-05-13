#!~/miniconda3/envs/scvi-env/bin/python

import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor="white")

results_file = "/public/home/chenjiaminggroup/wufan/20250207_scRNAseq_IBD_PBMC_WF/analysis_res/write/"

sampleName = pd.read_table("/public/home/chenjiaminggroup/wufan/20250207_scRNAseq_IBD_PBMC_WF/analysis_res/sampleName.txt",header=None)
samples = sampleName[0].to_list()
adatas = {}

## Read samples and merge all samples
for sample_id in samples:
    path = "/public/home/chenjiaminggroup/wufan/20250207_scRNAseq_IBD_PBMC_WF/"+sample_id+"_scRNA/outs/filtered_feature_bc_matrix.h5"
    sample_adata = sc.read_10x_h5(path)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
# print(adata.obs["sample"].value_counts())
print(adata)
# AnnData object with n_obs × n_vars = 346505 × 38606
#     obs: 'sample', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'log1p_total_counts_hb', 'pct_counts_hb', 'n_genes'
#     var: 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts', 'n_cells'


## Filter cells and genes
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
# percent_top=[20]

mito_filter = 10
n_counts_filter = 5000
fig, axs = plt.subplots(ncols = 2, figsize = (8,4))
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',ax = axs[0], show=False)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',ax = axs[1], show = False)
#draw horizontal red lines indicating thresholds.
axs[0].hlines(y = mito_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed',linewidth=0.5)
axs[1].hlines(y = n_counts_filter, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed',linewidth=0.5)
axs[1].hlines(y = 400, xmin = 0, xmax = max(adata.obs['total_counts']), color = 'red', ls = 'dashed',linewidth=0.5)
axs[1].vlines(x = 25000, ymin = 0, ymax = max(adata.obs['n_genes_by_counts']), color = 'red', ls = 'dashed',linewidth=0.5)
fig.tight_layout()
plt.savefig('QC.scatter_plots.png', dpi=150)
plt.close(fig)

sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    stripplot = False, # jitter=0.4,
    multi_panel=True, show=False, 
    save="QC.violin.png")

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
    show=False,
    save="QC.scatter.png")


n0 = adata.shape[0]
print(f'Original cell number: {n0}')
print("\n")


# for i in [25, 50, 75]:
#     percentile_value = np.percentile(adata.obs['pct_counts_mt'], i)
#     print(f'{i}th percentile of pct_counts_mt: {percentile_value}')

# for i in [25, 50, 75]:
#     percentile_value = np.percentile(adata.obs['n_genes_by_counts'], i)
#     print(f'{i}th percentile of n_genes_by_counts: {percentile_value}')

# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
sc.pp.filter_cells(adata, max_genes=5000)
n1 = adata.shape[0]
print(f'Higher treshold, n_genes_by_counts: 5000; filtered-out-cells: {n0-n1}, remain {n1} cells')

sc.pp.filter_cells(adata, min_genes=400)
n2 = adata.shape[0]
print(f'Lower treshold, n_genes_by_counts: 200; filtered-out-cells: {n1-n2}, remain {n2} cells')

adata = adata[adata.obs['pct_counts_mt']<15]
n3 = adata.shape[0]
print(f'Higher treshold, pct_counts_mt: 10%; filtered-out-cells: {n2-n3}, remain {n3} cells')

# adata = adata[~((adata.obs['total_counts']>3000) & (adata.obs['n_genes_by_counts']<400))]
adata = adata[adata.obs['total_counts']<=25000]
n4 = adata.shape[0]
print(f'Removing the outlier cells in scatter plot: {n3-n4}, remain {n4} cells, last cells before rm doublet')
print("\n")


g0 = adata.shape[1]
sc.pp.filter_genes(adata, min_cells=3)
print(f'Gene treshold, min_cells: 3; filtered-out-genes: {g0-adata.shape[1]}, remain {adata.shape[1]} genes')


sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    stripplot = False, # jitter=0.4,
    multi_panel=True, show=False,
    save="QC.violin.afterFilter.png")

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
    show=False,
    save="QC.scatter.afterFilter.png")


# adata = adata[~(adata.obs['sample'].isin(['B118CP','B232BP','B272AP']))]


## Doublet detection
sc.pp.scrublet(adata, batch_key="sample")

print(adata.obs.groupby("sample")["predicted_doublet"].value_counts().unstack(fill_value=0))
print('\n')


adata = adata[~adata.obs['predicted_doublet']].copy()
n5 = adata.shape[0]
print(f'Removing the doublet cells by scrublet: {n4-n5}, remain {n5} cells, last cells after rm doublet')
print("\n")

sc.pl.scrublet_score_distribution(adata, scale_hist_obs = 'log', scale_hist_sim = 'log', show=False, save="QC.doublets.png")



adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata) #log(normalize+1)
adata.raw = adata

# adda metadata information
sampleInfor = pd.read_csv('/public/home/chenjiaminggroup/wufan/20250207_scRNAseq_IBD_PBMC_WF/analysis_res/sampleInformation_clear.txt', sep='\t')
merge_df = adata.obs.merge(sampleInfor, on='sample', how='left').copy()
merge_df.index = adata.obs.index
adata.obs = merge_df
del(merge_df)
del(sampleInfor)

adata.obs.drop(columns=['Description','predicted_doublet','Lib_time','Tissue','Type','IBD'], inplace = True)
adata = adata[adata.obs['Group'] != '_']
# 291825 cells

mt_genes = adata.var[adata.var['mt']].index.to_list()
ribo_genes = adata.var[adata.var['ribo']].index.to_list()
hb_genes = adata.var[adata.var['hb']].index.to_list()
adata = adata[:,~adata.var.index.isin(mt_genes + ribo_genes + hb_genes)].copy()

### remove batch effect by scVI
sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, subset=True,
    layer="counts", flavor="seurat_v3", batch_key="sample",)

scvi.model.SCVI.setup_anndata(adata, layer="counts", 
    categorical_covariate_keys=["sample","Group"],
    continuous_covariate_keys=["pct_counts_mt"])

num_layers = 2
num_latent = 30
disp="gene-batch"
g_like = "zinb"
model = scvi.model.SCVI(adata, n_layers=num_layers, n_latent=num_latent, gene_likelihood=g_like,
        dispersion=disp)

print("scVI model infor: ", model,'\n')

model.train(max_epochs=1000, early_stopping=True)

model_dir = os.path.join(results_file, "scvi_model_plusGroup_rmTHMTinHVG"+'_'+str(num_layers)+'_'+str(num_latent)+'_'+disp+'_'+g_like)
model.save(model_dir, overwrite=True)

train_elbo = model.history['elbo_train'][1:]
test_elbo = model.history['elbo_validation']
ax = train_elbo.plot()
test_elbo.plot(ax=ax).figure.savefig("model_scvi_plusGroup_rmTHMTinHVG_loss.png")


### Obtaining model outputs

SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
print(latent.shape)

## Nearest neighbor graph constuction and visualization
# adata_test = adata.copy()
# for n_n in [5,10,15,20,25,40,50]: # original is 15
#     sc.pp.neighbors(adata_test, use_rep='X_scVI', n_neighbors = n_n, knn=True)
#     sc.tl.umap(adata_test, min_dist = 0.3) 
#     sc.pl.umap(adata_test, color=["sample", "Group", "Sex", "Drug","pct_counts_mt","pct_counts_ribo"],
#         size=2, show=False, wspace=0.5,hspace=0.5,ncols=3,
#         save='UMAP.NNG.plot.scVI_RmBEs_plusGroup'+f"_neighbors{n_n}"+'.png')
# del(adata_test)

sc.pp.neighbors(adata, use_rep='X_scVI', n_neighbors = 30, knn=True)
sc.tl.umap(adata, min_dist = 0.5) 
sc.pl.umap(adata, color=["sample", "Group", "Sex", "Drug","pct_counts_mt","pct_counts_ribo"],
    size=2, show=False, wspace=0.5,hspace=0.5,ncols=3,
    save="UMAP.NNG.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_neighbors30.png")

for res in [0.15, 0.5, 0.18, 0.2, 0.22]:
    sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, n_iterations = 2, flavor="igraph")
    sc.pl.umap(adata, color=[f"leiden_res_{res:4.2f}"], show=False, legend_loc = 'on data',
        save=f"leiden_res_{res:4.2f}"+"_clustering.RmBE_plusGroup_rmTHMTinHVG.png")

# sc.tl.leiden(adata, key_added="leiden_res_0.22", resolution=0.22, n_iterations = 2, flavor="igraph")
# sc.pl.umap(adata, color=["leiden_res_0.22"], show=False, legend_loc = 'on data',
#     save="leiden_res_0.22_clustering.RmBE_plusGroup_rmTHMTinHVG.png")

adata.write(results_file + 'pbmcMerge.afterNNG.scVI_RmBE_plusGroup_rmTHMTinHVG.20250214.h5ad')




for res in [0.15, 0.18, 0.22, 0.2]:
    # sc.tl.leiden(adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph", n_iterations=2)
    sc.tl.rank_genes_groups(adata, groupby=f"leiden_res_{res:4.2f}", method="wilcoxon", use_raw = True)
    sc.pl.rank_genes_groups_dotplot(adata, groupby=f"leiden_res_{res:4.2f}", standard_scale="var", n_genes=5,
        show=False, save=f"dotplot_markers.leiden_res_{res:4.2f}.RmBE_plusGroup_rmTHMTinHVG.png")
    # sc.pl.rank_genes_groups_stacked_violin(adata, groupby=f"leiden_res_{res:4.2f}", n_genes=5,
    #     show=False, save=f"violin_markers.leiden_res_{res:4.2f}.RmBE_plusGroup.png")
    # sc.pl.rank_genes_groups_heatmap(adata, groupby=f"leiden_res_{res:4.2f}", n_genes=5,
    #     show=False, save=f"heatmap_markers.leiden_res_{res:4.2f}.RmBE_plusGroup.png")
    print("Remove batch's marker genes:")
    print(pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).head(5))
    print("\n")
    markersDF1 = sc.get.rank_genes_groups_df(adata, group = None,pval_cutoff=0.001,log2fc_min=1)
    markersDF1.to_csv(f"rank_genes_groups_{res:4.2f}.marker_genes_RmBE_plusGroup_rmTHMTinHVG.csv", index=False)


adata.obs.groupby('sample')['leiden_res_0.22'].value_counts().unstack(fill_value=0)


marker_genes = {
    "Neu": ['CSF3R',"FCGR3B","NAMPT", "NEAT1","AQP9","G0S2"],
    "Tcell": ["CD3E","CD3D"],
    'CD4+T': ["CD4", "IL7R", "CD40LG"],
    "NKT": ["NKG7", "GZMA", "GNLY", "SYNE1","GZMB","GZMH","GZMK"],
    "Precursor": ["CD34", "KIT",'CDK6'],
    "Mega": ["PPBP", "PF4","NRGN","GP1BB"],
    "T+Neu mix":["CD3E","CD3D",'CSF3R',"FCGR3B"],
    "CD8T": ["CD8A", "CD8B","NELL2"],
    "B": ["CD79A", "BANK1", "MS4A1","IGHD","FCRL1","IGHM","PAX5"],
    "CD14+Mono": ["LYZ", "CD14", "FCN1","CST3","DMXL2"],
    "CD16+Mono":['FCGR3A'],
    'pDC':['IRF8','TCF4','UGCG','PTPRS'],
    'Basophil':['ENPP3','CCR3','IL3RA'],
    "plasma": ["JCHAIN",'IRF4',"XBP1","CD74","LILRA4",],
    "DC": ["CD1C", "CD83"],
}

# 'CD1C' in adata.raw.var_names

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.22", standard_scale="var", use_raw = True, 
show=False, save="dotplot_markers.leiden_res_0.22.RmBE_plusGroup_rmTHMTinHVG_markerGenes.png")
# sc.pl.stacked_violin(adata, marker_genes, groupby='leiden_res_0.20', show=False, save="stacked_violin_violin_markers.leiden_res_0.20.RmBE_plusGroup_markerGenes.png")

# sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.30", standard_scale="var", use_raw = True, 
# show=False, save="dotplot_markers.leiden_res_0.30.RmBE_plusGroup_rmTHMTinHVG_markerGenes.cycle3.png")


marker_genes2 = {
    "CD14+ Mono": ["FCN1", "CD14"],
    "CD16+ Mono": ["TCF7L2", "FCGR3A", "LYN"],
    # Note: DMXL2 should be negative
    "cDC2": ["CST3", "COTL1", "LYZ", "DMXL2", "CLEC10A", "FCER1A"],
    "Erythroblast": ["MKI67", "HBA1", "HBB"],
    # Note HBM and GYPA are negative markers
    "Proerythroblast": ["CDK6", "SYNGR1", "HBM", "GYPA"],
    "NK": ["GNLY", "NKG7", "CD247", "FCER1G", "TYROBP", "KLRG1", "FCGR3A"],
    "ILC": ["ID2", "PLCG2", "GNLY", "SYNE1"],
    "Naive CD20+ B": ["MS4A1", "IL4R", "IGHD", "FCRL1", "IGHM"],
    # Note IGHD and IGHM are negative markers
    "B cells": [
        "MS4A1",
        "ITGB1",
        "COL4A4",
        "PRDM1",
        "IRF4",
        "PAX5",
        "BCL11A",
        "BLK",
        "IGHD",
        "IGHM",
    ],
    "Plasma cells": ["MZB1", "HSP90B1", "FNDC3B", "PRDM1", "IGKC", "JCHAIN"],
    # Note PAX5 is a negative marker
    "Plasmablast": ["XBP1", "PRDM1", "PAX5"],
    "CD4+ T": ["CD4", "IL7R", "TRBC2"],
    "CD8+ T": ["CD8A", "CD8B", "GZMK", "GZMA", "CCL5", "GZMB", "GZMH", "GZMA"],
    "T naive": ["LEF1", "CCR7", "TCF7"],
    "pDC": ["GZMB", "IL3RA", "COBLL1", "TCF4"],
}

# sc.pl.dotplot(adata, marker_genes2, groupby="leiden_res_0.20", standard_scale="var",show=False, save="dotplot_markers.leiden_res_0.20.RmBE_plusGroup_markerGenes2.png")
# sc.pl.stacked_violin(adata, marker_genes2, groupby='leiden_res_0.20', show=False, save="stacked_violin_violin_markers.leiden_res_0.20.RmBE_plusGroup_markerGenes2.png")
sc.pl.dotplot(adata, marker_genes2, groupby="leiden_res_0.22", standard_scale="var",show=False, save="dotplot_markers.leiden_res_0.22.RmBE_plusGroup_rmTHMTinHVG_markerGenes2.png")

# using celltypist to test; step 5.1


# final cell type
adata.obs["cell_type1"] = adata.obs["leiden_res_0.22"].map(
    {
        "0": "Neutrophil",
        "1": "Neutrophil",
        "2": "Neu-del",
        "3": "T(DNT+CD4T)",
        "4": "NKT",
        "5": "Precursor",
        "6": "Mega+NKT-del",
        "7": "Mix-del",
        "8": "CD8T",
        "9": "B",
        "10": "Mono",
        "11": "pDC",
        "12": "Basophil",
        "13": "plasmaB",
        "14": "DC"
    }
)

sc.pl.umap(adata,
           color=['CSF3R','FCGR3B','CD3D','CD4','CD8B','NKG7','CD34','PPBP','BANK1',
                'CD14','LYZ','JCHAIN','TCF4','ENPP3','CD83','leiden_res_0.22','cell_type1'],
           frameon=False,
           ncols=3,
           show=False,save="UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_celltype1.png")

adata.write(results_file+'pbmcMerge.afterFindMarkers_RmBE_plusGroup_rmTHMTinHVG.20250214.celltype.h5ad')

sc.pl.umap(adata,
           color=['leiden_res_0.22','cell_type1'],
           frameon=False,
           ncols=3,legend_loc = 'on data',
           show=False,save="UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_celltype1.1.png")


# rmCells = adata.obs[adata.obs['doublet_score']>0.1].index.to_list()
# adata = adata[~adata.obs.index.isin(rmCells)].copy()
# sc.pl.umap(adata,
#            color=['CSF3R','FCGR3B','CD3D','CD4','CD8B','NKG7','CD34','PPBP','BANK1',
#                 'CD14','LYZ','JCHAIN','TCF4','ENPP3','CD83','leiden_res_0.22','cell_type1'],
#            frameon=False,
#            ncols=3, 
#            show=False,save="UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_celltype1.rmDS0.1.png")
# sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.22", standard_scale="var", use_raw = True, 
# show=False, save="dotplot_markers.leiden_res_0.22.RmBE_plusGroup_rmTHMTinHVG_markerGenes.rmDS0.1.png")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

my_list = adata.obs[adata.obs['leiden_res_0.50'].isin(['5','9','16']) |
                    adata.obs['leiden_res_0.22'].isin(['2','7'])].index.to_list()
rmCells = pd.Series(my_list).unique().tolist()
# >>> len(rmCells)
# 12144

adata = adata[~adata.obs.index.isin(rmCells)].copy()
adata.write(results_file+'pbmcMerge.afterFindMarkers_RmBE_plusGroup_rmTHMTinHVG.20250214.celltype.rm259.h5ad')

adata.obs["cell_type2"] = adata.obs["leiden_res_0.22"].map(
    {
        "0": "Neutrophil",
        "1": "Neutrophil",
        "3": "T(DNT+CD4T)",
        "4": "NKT",
        "5": "Precursor",
        "6": "Mega(Platelets)",
        "8": "CD8T",
        "9": "B",
        "10": "Mono",
        "11": "pDC",
        "12": "Basophil",
        "13": "plasmaB",
        "14": "DC"
    }
)

adata.obs["leiden_res_0.22.2"] = adata.obs["leiden_res_0.22"].map(
    {
        "0": "10",
        "1": "10",
        "3": "0",
        "4": "1",
        "5": "8",
        "6": "9",
        "8": "2",
        "9": "3",
        "10": "5",
        "11": "7",
        "12": "11",
        "13": "4",
        "14": "6"
    }
)

desired_order = ["T(DNT+CD4T)","NKT","CD8T","B","plasmaB","Mono","DC","pDC","Precursor","Mega(Platelets)","Neutrophil","Basophil"]
adata.obs['cell_type2'] = pd.Categorical(
    adata.obs['cell_type2'],
    categories=desired_order,
    ordered=True
)

adata.uns['leiden_res_0.22.2_colors'] = ['#1f77b4', '#17becf', '#279e68', '#e377c2', '#aa40fc', '#8c564b', '#ff7f0e', '#b5bd61', '#d62728', '#aec7e8', '#ffbb78', '#98df8a']
adata.uns['cell_type2_colors'] = ['#1f77b4', '#17becf', '#279e68', '#e377c2', '#aa40fc', '#8c564b', '#ff7f0e', '#b5bd61', '#d62728', '#aec7e8', '#ffbb78', '#98df8a']

sc.pl.umap(adata,
           color=['leiden_res_0.22.2','cell_type2'],
           frameon = False,# legend_fontsize=7.5,
           legend_loc = 'on data',
           show=False, save="UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_celltype1.2.png")

sc.pl.umap(adata,
           color=['leiden_res_0.22.2','cell_type2'],
           frameon = False,  # legend_fontsize='x-small',
           legend_loc = 'right margin', size=2, add_outline=True, # outline_color='black',
           show=False, save="UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_celltype1.2.1.png")

desired_order = ["T(DNT+CD4T)","NKT","CD8T",
"B","plasmaB",
"Mono","DC","pDC",
"Precursor",
"Mega(Platelets)",
"Neutrophil","Basophil"]

marker_genes = {
    "Tcell": ["CD3E","CD3D","IL7R"],
    'CD4T': ["CD4", "CD40LG"],
    "NKT": ["NKG7", "GZMA", "GNLY", "SYNE1","GZMH","GZMK"],
    "CD8T": ["CD8A", "CD8B","NELL2"],
    "B": ["CD79A", "BANK1", "MS4A1","IGHD","FCRL1","IGHM","PAX5"],
    "plasmaB": ["JCHAIN",'IRF4',"XBP1","CD74","LILRA4",],
    "CD14+Mono": ["LYZ", "CD14", "FCN1","CST3","DMXL2"],
    "CD16+Mono":['FCGR3A'],
    "DC": ["CD1C", "CD83"],
    'pDC':['IRF8','TCF4','UGCG','PTPRS'],
    "Precursor": ["CD34", "KIT",'CDK6'],
    "Mega(platelets)": ["PPBP", "PF4","NRGN","GP1BB"],
    "Neutrophil": ['CSF3R',"FCGR3B","NAMPT", "NEAT1","AQP9","G0S2"],
    'Basophil':['ENPP3','CCR3','IL3RA'],
}

fig = sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.22.2", standard_scale="var", use_raw = True, 
return_fig=False, show = False)
ax = fig["mainplot_ax"]

for l in ax.get_xticklabels():
    l.set_style("italic")
    g = l.get_text()
    # Change settings (e.g. color) of certain ticklabels based on their text (here gene name)
    if g == "MS4A1":
        l.set_color("#A97F03")
plt.savefig("dotplot_markers.leiden_res_0.22.2.RmBE_plusGroup_rmTHMTinHVG_markerGenes.2.png", 
            dpi=150, bbox_inches="tight")
plt.close()

# adata.var.drop(columns=['italic_name'], inplace = True)
adata.write(results_file+'pbmcMerge.afterFindMarkers_RmBE_plusGroup_rmTHMTinHVG.20250214.celltype.rm259.h5ad')



for i in adata.obs['leiden_res_0.20'].unique():
    tempAdata = adata[adata.obs['leiden_res_0.20'].isin([str(i)]),:].copy()
    sc.pp.neighbors(tempAdata, use_rep='X_scVI', n_neighbors = 30, knn=True)
    sc.tl.umap(tempAdata, min_dist = 0.5) 
    sc.tl.leiden(tempAdata, key_added=f"C{i}_leiden_res_0.22", resolution=0.22, n_iterations = 2, flavor="igraph")
    sc.pl.umap(tempAdata, color=[f"C{i}_leiden_res_0.22"], show=False, legend_loc = 'on data',
        save=f"UMAP_C{i}_leiden_res_0.22_clustering.RmBE_plusGroup_rmTHMTinHVG.cycle3.png")
    sc.tl.rank_genes_groups(tempAdata, groupby=f"C{i}_leiden_res_0.22", method="wilcoxon", use_raw = True)
    sc.pl.rank_genes_groups_dotplot(tempAdata, groupby=f"C{i}_leiden_res_0.22", standard_scale="var", n_genes=5,
        show=False, save=f"dotplot_markers.C{i}_leiden_res_0.22.RmBE_plusGroup_rmTHMTinHVG.cycle3.png")
    markersDF1 = sc.get.rank_genes_groups_df(tempAdata, group = None, pval_cutoff=0.001, log2fc_min=1)
    markersDF1.to_csv(f"rank_genes_groups_C{i}_0.22.marker_genes_RmBE_plusGroup_rmTHMTinHVG.cycle3.csv", index=False)
    tempAdata.raw = None
    tempAdata.write(results_file+f'C{i}_pbmcMerge.afterFindMarkers_RmBE_plusGroup_rmTHMTinHVG.20250214.cycle3.h5ad')



