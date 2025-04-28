import scanpy as sc
from muon import prot as pt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from IPython.display import display, HTML
from itertools import chain
from adjustText import adjust_text

from pySPaRTAN import pySPaRTAN



adata=sc.read_10x_h5("../data/CytAssist_FFPE_Protein_Expression_Human_Tonsil_AddOns_filtered_feature_bc_matrix.h5", gex_only=False)
# making the variable name unique
adata.var_names_make_unique()
adata

set(adata.var["feature_types"])



# The ADT and RNA data are both stored in this adata object, differentiated by the "feature_types" variable. We will separate them into two distinct objects for ease of processing:
# extracting gene expression raw counts
RNA = adata[:, adata.var["feature_types"] == "Gene Expression"]



# extracting ADT raw counts
ADT = adata[:,adata.var["feature_types"]=="Antibody Capture"].copy()

# removing the isotype control for antibody
ADT=ADT[:,[x for x in ADT.var_names if "control" not in x]]


# ================ Quality control ================

# We calculate QC matrix and inspect the following features:
# * the number of genes expressed in the count matrix
# * the total counts per cell
# * the percentage of counts in mitochondrial genes

RNA.var['mt'] = RNA.var_names.str.startswith('MT-')
# sc.pp.calculate_qc_metrics(RNA, qc_vars=['mt'], percent_top=None, log1p=False,inplace=True)

sc.pl.violin(
    RNA,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)


# Based on QC metric plots, we filter out cells with fewer than 1,000 or more than 5,000 unique gene sequences, as well as cells with more than 30% mitochondrial counts. Additionally, we exclude genes that are not present in at least 3% of cells. we also remove all mitochondrial genes from the dataset.
# filtering mRNA cells
sc.pp.filter_cells(RNA, min_genes=1000)
RNA=RNA[RNA.obs.query("n_genes_by_counts < 5000 and pct_counts_mt<30").index]

# filtering mRNA genes
sc.pp.filter_genes(RNA, min_cells=0.03*RNA.n_obs)

# removing all mitochondrial genes from the dataset
RNA=RNA[:, RNA.var['mt']==False]

# filtering cells of ADT according to RNA
ADT = ADT[RNA.obs_names, :]

    