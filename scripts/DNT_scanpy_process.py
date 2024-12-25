import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.pct_counts_mt < 10, :]
adata = adata[adata.obs.n_genes_by_counts < 4000, :]

gene = pd.read_csv("/work/haoq/Pancancer_T/gene_symbol/HGNC-20231122.txt",sep = "\t")
gene.columns = gene.columns.str.replace(' ', '_')

gene = gene[(gene.Locus_group.isin(["protein-coding gene","non-coding RNA","other"])) & 
     (gene.Status == "Approved")]

selected_genes = adata.var.HGNCsymbol.isin(gene.Approved_symbol)
adata = adata[:, selected_genes]

# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(('RPL', 'RPS'))

adata = adata[:, adata.var['mt'] == False]
adata = adata[:, adata.var['ribo'] == False]

adata.obs['batch'] = adata.obs['ResearchID'].astype(str) + "_" + adata.obs['Patient_ID'].astype(str)
patient_name =list(adata.obs.groupby("batch").count()[list(adata.obs.groupby("batch").count().iloc[:,0]>5)].index)
adata =adata[[i in patient_name for i in adata.obs['batch']]]

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = "diseaseTCGA")
sc.pl.highly_variable_genes(adata)

sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='ResearchID')

sc.tl.paga(adata)
sc.pl.paga(adata) 
sc.tl.umap(adata, init_pos='paga')
sc.tl.leiden(adata)

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon',use_raw=False)

adata.write("/work/haoq/Pancancer_T/h5ad/pancancerT_symbol_leiden.h5ad", compression="gzip")
