# how to update an old version of an h5ad (annData) file
# this step has already been done for the objects so this is just for reference

import scanpy as sc

#please note the file extension will be different on your machine 
adata = sc.read_h5ad("tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad")

adata.write("droplet-marrow.h5ad")


