# Dependencies
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse

# Import h5ad file
adata = sc.read_h5ad("raw/snRNA-seq-submission.h5ad")

# extract sparse matrix
mat = adata.raw.X

# Convert to dataframe
df = pd.DataFrame.sparse.from_spmatrix(mat)
df.to_csv('file_name.csv', index=False)

