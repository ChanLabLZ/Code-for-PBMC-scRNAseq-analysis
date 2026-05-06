#!~/miniconda3/bin/python

import scanpy as sc
import os
import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/public/home/chenjiaminggroup/wufan/20250207_scRNAseq_IBD_PBMC_WF/analysis_res/")

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor="white")

results_file = "/public/home/chenjiaminggroup/wufan/20250207_scRNAseq_IBD_PBMC_WF/analysis_res/write/"
adata = sc.read_h5ad(results_file+'pbmcMerge.afterFindMarkers_RmBE_plusGroup_rmTHMTinHVG.20250401.h5ad')

adata.obs.groupby('sample')['leiden_res_0.20'].value_counts().unstack(fill_value=0)

marker_genes = {
    "Tcell": ["CD3E","CD3D","IL7R"],
    'CD4T': ["CD4", "CD40LG"],
    "NKT": ["NKG7", "GZMA", "GNLY", "SYNE1","GZMH","GZMK"],
    "CD8T": ["CD8A", "CD8B","NELL2"],
    "B": ["CD79A", "BANK1", "MS4A1","IGHD","FCRL1","IGHM","PAX5"],
    "plasmaB": ["JCHAIN",'IRF4',"XBP1","LILRA4",],
    "CD14+Mono": ["LYZ", "CD14", "FCN1","CST3","DMXL2"],
    "CD16+Mono":['FCGR3A'],
    "DC": ["CD1C", "CD83",'HLA-DPA1','HLA-DPB1','HLA-DRA'],
    'pDC':['IRF8','TCF4','UGCG','PTPRS'],
    "Precursor": ["CD34", "KIT",'CDK6'],
    "Mega(platelets)": ["PPBP", "PF4","NRGN","GP1BB"],
    "Neutrophil": ['CSF3R',"FCGR3B","NAMPT", "NEAT1","AQP9","G0S2"],
    'Basophil':['ENPP3','CCR3','IL3RA'],
}

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.20", standard_scale="var", use_raw = True, 
show=False, save="250401.dotplot_markers.leiden_res_0.20.RmBE_plusGroup_rmTHMTinHVG_markerGenes.png")

adata.obs["cell_type2"] = adata.obs["leiden_res_0.20"].map(
    {
        "0": "Neutrophil",
        "1": "plasmaB",
        "2": "NKT",
        "3": "Precursor",
        "4": "T(DNT+CD4T)",
        "5": "CD8T",
        "6": "Mega(Platelets)",
        "7": "Mono",
        "8": "DC",
        "9": "B",
        "10": "Basophil",
        "11": "pDC",
    }
)


adata.obs["leiden_res_0.22_new"] = adata.obs["leiden_res_0.20"].map(
    {
        "0": "10",
        "1": "4",
        "2": "1",
        "3": "8",
        "4": "0",
        "5": "2",
        "6": "9",
        "7": "5",
        "8": "6",
        "9": "3",
        "10": "11",
        "11": "7",
    }
)

desired_order = ["T(DNT+CD4T)","NKT","CD8T","B","plasmaB","Mono","DC","pDC","Precursor","Mega(Platelets)","Neutrophil","Basophil"]
adata.obs['cell_type2'] = pd.Categorical(
    adata.obs['cell_type2'],
    categories=desired_order,
    ordered=True
)
desired_order = ["0","1","2","3","4","5","6","7","8","9","10","11"]
adata.obs['leiden_res_0.22_new'] = pd.Categorical(
    adata.obs['leiden_res_0.22_new'],
    categories=desired_order,
    ordered=True
)


adata.uns['leiden_res_0.22_new_colors'] = ['#1f77b4', '#17becf', '#279e68', '#e377c2', '#aa40fc', '#8c564b', '#ff7f0e', '#b5bd61', '#d62728', '#aec7e8', '#ffbb78', '#98df8a']
adata.uns['cell_type2_colors'] = ['#1f77b4', '#17becf', '#279e68', '#e377c2', '#aa40fc', '#8c564b', '#ff7f0e', '#b5bd61', '#d62728', '#aec7e8', '#ffbb78', '#98df8a']

sc.pl.umap(adata,
           color=['leiden_res_0.22_new','cell_type2'],
           frameon = False,# legend_fontsize=7.5,
           legend_loc = 'right margin',
           show=False, save="250401.UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_celltype1.2.png")
sc.pl.umap(adata,
           color=['leiden_res_0.22_new','cell_type2'],
           frameon = False,# legend_fontsize=7.5,
           legend_loc = 'on data',
           show=False, save="250401.UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_celltype1.2.2.png")
sc.pl.umap(adata,
           color=['leiden_res_0.22_new','cell_type2'],
           frameon = False,  # legend_fontsize='x-small',
           legend_loc = 'right margin', size=2, add_outline=True, # outline_color='black',
           show=False, save="250401.UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_celltype1.2.1.png")

marker_genes = {
    "Tcell": ["CD3E","CD3D","IL7R"],
    'CD4T': ["CD4", "CD40LG"],
    "NKT": ["NKG7", "GZMA", "GNLY", "SYNE1","GZMH","GZMK"],
    "CD8T": ["CD8A", "CD8B","NELL2"],
    "B": ["CD79A", "BANK1", "MS4A1","IGHD","FCRL1","IGHM","PAX5"],
    "plasmaB": ["JCHAIN",'IRF4',"XBP1","LILRA4",],
    "Mono": ["LYZ", "CD14", "FCN1","CST3","DMXL2"],
    "DC": ["CD83",'HLA-DPA1','HLA-DPB1','HLA-DRA'],
    'pDC':['IRF8','TCF4','UGCG','PTPRS'],
    "Precursor": ["CD34", "KIT",'CDK6'],
    "Mega(platelets)": ["PPBP", "PF4","NRGN","GP1BB"],
    "Neutrophil": ['CSF3R',"FCGR3B","NAMPT", "NEAT1","AQP9","G0S2"],
    'Basophil':['ENPP3','CCR3','IL3RA'],
}

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.22_new", standard_scale="var", use_raw = True, 
show=False, save="250401.dotplot_markers.leiden_res_0.20.RmBE_plusGroup_rmTHMTinHVG_markerGenes.png")

fig = sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.22_new", standard_scale="var", use_raw = True, 
return_fig=False, show = False)
ax = fig["mainplot_ax"]
for l in ax.get_xticklabels():
    l.set_style("italic")
    # g = l.get_text()
    # # Change settings (e.g. color) of certain ticklabels based on their text (here gene name)
    # if g == "MS4A1":
    #     l.set_color("#A97F03")
plt.savefig("250401.dotplot_markers.leiden_res_0.20.RmBE_plusGroup_rmTHMTinHVG_markerGenes.1.png", 
            dpi=150, bbox_inches="tight")
plt.close()


# read again

# the UMAP in different groups
ncols = 4
nrows = 1
figsize = 4
wspace = 0.5
fig, axs = plt.subplots(nrows=nrows, ncols=ncols,
    figsize=(ncols * figsize + figsize * 0.1 * (ncols-1), nrows * figsize),)
plt.subplots_adjust(wspace=0.1)
# print("axes:", axs)
sc.pl.umap(adata[adata.obs.Group.isin(["UNT"]),:], color='cell_type2', ax=axs[0], title='UNT', frameon=False, legend_loc = None, show=False,)
sc.pl.umap(adata[adata.obs.Group.isin(["NR"]), :], color='cell_type2', ax=axs[1], title='NR', frameon=False, legend_loc = None, show=False,)
sc.pl.umap(adata[adata.obs.Group.isin(["SOR"]), :], color='cell_type2', ax=axs[2], title='SOR', frameon=False, legend_loc = None, show=False,)
sc.pl.umap(adata[adata.obs.Group.isin(["R"]), :], color='cell_type2', ax=axs[3], title='R', frameon=False, legend_loc = 'right margin', show=False,)
plt.savefig("UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_celltype2.splitGroup.png", dpi=150, bbox_inches="tight")
plt.close()

adata.obs["lineage"] = adata.obs["leiden_res_0.22_new"].map(
    {
        "0": "T",
        "1": "T",
        "2": "T",
        "3": "B",
        "4": "plasmaB",
        "5": "Mono",
        "6": "Mono",
        "7": "Mono",
        "8": "Precursor",
        "9": "Mega",
        "10": "Granu",
        "11": "Granu"
    }
)

ncols = 4
nrows = 1
figsize = 4
wspace = 0.5
fig, axs = plt.subplots(nrows=nrows, ncols=ncols,
    figsize=(ncols * figsize + figsize * 0.1 * (ncols-1), nrows * figsize),)
plt.subplots_adjust(wspace=0.1)
sc.pl.umap(adata[adata.obs.Group.isin(["UNT"]),:], color='lineage', ax=axs[0], title='UNT', frameon=False, legend_loc = None, show=False,)
sc.pl.umap(adata[adata.obs.Group.isin(["NR"]), :], color='lineage', ax=axs[1], title='NR', frameon=False, legend_loc = None, show=False,)
sc.pl.umap(adata[adata.obs.Group.isin(["SOR"]), :], color='lineage', ax=axs[2], title='SOR', frameon=False, legend_loc = None, show=False,)
sc.pl.umap(adata[adata.obs.Group.isin(["R"]), :], color='lineage', ax=axs[3], title='R', frameon=False, legend_loc = 'right margin', show=False,)
plt.savefig("UMAP.plot.scVI_RmBEs_plusGroup_rmTHMTinHVG_lineage.splitGroup.png", dpi=150, bbox_inches="tight")
plt.close()



for i in adata.obs['leiden_res_0.22_new'].unique():
    tempAdata = adata[adata.obs['leiden_res_0.22_new'].isin([str(i)]),:].copy()
    sc.pp.neighbors(tempAdata, use_rep='X_scVI', n_neighbors = 30, knn=True)
    sc.tl.umap(tempAdata, min_dist = 0.5) 
    sc.tl.leiden(tempAdata, key_added=f"C{i}_leiden_res_0.22", resolution=0.22, n_iterations = 2, flavor="igraph")
    sc.pl.umap(tempAdata, color=[f"C{i}_leiden_res_0.22"], show=False, legend_loc = 'on data',
        save=f"250401.UMAP_C{i}_leiden_res_0.22_clustering.RmBE_plusGroup_rmTHMTinHVG.cycle4.png")
    sc.tl.rank_genes_groups(tempAdata, groupby=f"C{i}_leiden_res_0.22", method="wilcoxon", use_raw = True)
    sc.pl.rank_genes_groups_dotplot(tempAdata, groupby=f"C{i}_leiden_res_0.22", standard_scale="var", n_genes=5,
        show=False, save=f"250401.dotplot_markers.C{i}_leiden_res_0.22.RmBE_plusGroup_rmTHMTinHVG.cycle4.png")
    markersDF1 = sc.get.rank_genes_groups_df(tempAdata, group = None, pval_cutoff=0.001, log2fc_min=1)
    markersDF1.to_csv(f"250401.rank_genes_groups_C{i}_0.22.marker_genes_RmBE_plusGroup_rmTHMTinHVG.cycle4.csv", index=False)
    tempAdata.raw = None
    tempAdata.write(results_file+f'250401.C{i}_pbmcMerge.afterFindMarkers_RmBE_plusGroup_rmTHMTinHVG.20250401.cycle4.h5ad')



# identify the DEG for each cell type;
adata.obs['Group2'] = adata.obs['Group']
adata.obs['Group2'] = adata.obs['Group2'].cat.add_categories(['NR+SOR'])
adata.obs.loc[adata.obs['Group'].isin(['NR','SOR']),'Group2'] = 'NR+SOR'
adata.obs['Group2'] = pd.Categorical(
    adata.obs['Group2'],
    categories=["UNT", "NR+SOR", "R"],
    ordered=True)
adata.obs['Group3'] = adata.obs.apply(lambda row: f"{row['Group2']}_{row['cell_type2']}", axis=1)
adata.raw = adata.copy()
for i in adata.obs['cell_type2'].unique().tolist():
    # for j in adata.obs['Group2'].unique().to_list():
    a = 'UNT_'+i
    b = 'NR+SOR_'+i
    c = 'R_'+i

    sc.tl.rank_genes_groups(adata, groupby="Group3", groups= [a], reference = b, method="wilcoxon", use_raw = True)
    markersDF1 = sc.get.rank_genes_groups_df(adata, group = None, pval_cutoff=0.01)
    markersDF1.to_csv(f"250401.rank_genes_groups_{a}_vs_{b}_Group3_RmBE_plusGroup_rmTHMTinHVG.csv", index=False)

    sc.tl.rank_genes_groups(adata, groupby="Group3", groups= [a], reference = c, method="wilcoxon", use_raw = True)
    markersDF1 = sc.get.rank_genes_groups_df(adata, group = None, pval_cutoff=0.01)
    markersDF1.to_csv(f"250401.rank_genes_groups_{a}_vs_{c}_Group3_RmBE_plusGroup_rmTHMTinHVG.csv", index=False)

    sc.tl.rank_genes_groups(adata, groupby="Group3", groups= [b], reference = c, method="wilcoxon", use_raw = True)
    markersDF1 = sc.get.rank_genes_groups_df(adata, group = None, pval_cutoff=0.01)
    markersDF1.to_csv(f"250401.rank_genes_groups_{b}_vs_{c}_Group3_RmBE_plusGroup_rmTHMTinHVG.csv", index=False)


tt = adata.obs

tt['group_old'] = tt['Group']
tt['Group'] = tt['Group'].replace('SOR', 'NR')

tt['Group'] = pd.Categorical(tt['Group'], categories=['UNT', 'NR', 'R'], ordered=True)
adata.obs = tt


tt['cell_type2'] = tt['cell_type2'].cat.add_categories(['DNT'])
tt.loc[tt['cell_type3'] == '0_DNT', 'cell_type2'] = 'DNT'
tt['cell_type2'] = tt['cell_type2'].replace('T(DNT+CD4T)', 'CD4T')
tt['cell_type2'] = pd.Categorical(tt['cell_type2'], 
        categories=['DNT', 'CD4T', 'CD8T', 'NKT','B', 'plasmaB', 'Mono',
        'DC', 'pDC', 'Precursor', 'Mega(Platelets)', 'Neutrophil','Basophil'], 
        ordered=True)
adata.obs = tt

adata_sub = adata[(adata.obs['cell_type3'] != '1_Neu-del') & (adata.obs['cell_type3'] != '7_doublet-del')].copy()
adata_sub


# ### test marker genes
# test_genes = {
#     'TNFγ':['IFNG'],
#     "Receptor": ["IFNGR1","IFNGR2"],
#     'Downstream':['JAK1','JAK2','STAT4','MAPK3','MAPK1','PTK2B','CRKL']
# }

# sc.pl.dotplot(adata_sub, test_genes, groupby='cell_type3', standard_scale="var", use_raw = True, swap_axes=True,
# show=False, save="250501.dotplot_markers.celltype3.RmBE_plusGroup_rmTHMTinHVG_markerGenes.IFN.png")



# ### Just test the counts of cell polirazed cell by R.
# #

adata_sub.raw = None
adata_sub.write(results_file + 'pbmcMerge.afterFindMarkers_RmBE_plusGroup_rmTHMTinHVG.20250801.h5ad')


