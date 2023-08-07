import pickle
from collections import ChainMap #hopefully this works, if not we can merge with comprehension
import pandas as pd
import anndata
import numpy as np
import scanpy as sc

#Get all the entries written
data = []
with open('shared_dict', 'rb') as fr:
    try:
        while True:
            data.append(pickle.load(fr))
    except EOFError:
        pass

#Combine all the indivudal written entries into big dictionary
total_dict = dict(ChainMap(*data))

#Load in PTR csv. It has a good header (sample names), but the first column should just be row attributes instead of a real collumn
df = pd.read_csv('out.csv')
df.fillna(0,inplace=True)
#change the row index values to the first column
df.index = [i for i in df.iloc[:,0]]
#remove that first column and put the df into anndata
df = df.iloc[:,1:]
adata_PTR = anndata.AnnData(df)
#Now we transpose it (already a AnnData object)
adata_PTR_trans = adata_PTR.transpose()

#Get the obs and var atrributes from it
obs_len = adata_PTR_trans.n_obs 
var_len = adata_PTR_trans.n_vars
obs = adata_PTR_trans.obs_names
var = adata_PTR_trans.var_names

#Create a new zeros df for abundances
abund_df = pd.DataFrame(np.zeros((obs_len, var_len)))
abund_df.index = list(obs)
abund_df.columns = list(var)

#sort df for faster access
abund_df.sort_index(inplace=True)

#Edit this new datafram with abundances
#May want to look into optimizations
for tup, count in total_dict.items():
    abund_df.loc[tup[1],tup[0]] = count

#convert this df to a adata object and make it the PTR object a layer of it
adata_abund = anndata.AnnData(abund_df)
adata_abund.layers["abundance"] = adata_abund.X #Need to refer to the abundance as a layer, so make it an alias for the base layer X
adata_abund.layers["PTR"] = adata_PTR_trans.X

print(adata_abund.to_df(layer="abundance"))
print(adata_abund.to_df(layer="PTR"))

#ADD CLUSTER LABELS
ct = np.random.choice(["Healthy","Sick"], size=(adata_abund.n_obs,))
adata_abund.obs["clusters"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency

print(adata_abund.obs["clusters"])

#TEST FILTRATION AND NORMALIZATION
#sc.pp.filter_genes(adata_abund, min_cells = 5) #do filtration once you have a larger dataset, can tweak parameters to decide species to exclude based on # counts
sc.pp.normalize_total(adata_abund, layer ="abundance") #normalization does this for the base layer X as well

print(adata_abund.to_df(layer="abundance"))
print(adata_abund.to_df(layer="PTR"))

sc.tl.pca(adata_abund, random_state=0)
sc.pp.neighbors(adata_abund, random_state=0)

#Getting VELOCITY Kernel
print("KERNEL TIME")
from cellrank.kernels import VelocityKernel
vk = VelocityKernel(adata_abund, gene_subset=adata_abund.var, vkey="PTR",xkey="abundance")
vk.compute_transition_matrix()

#Getting CONNECTIVITY Kernel
from cellrank.kernels import ConnectivityKernel
ck = ConnectivityKernel(adata_abund).compute_transition_matrix()

#Combined kernel (linear combination of matrices)
kernel = 0.8*vk + 0.2*ck

#ESTIMATOR
#from cellrank.tl.estimators import GPCCA #cellrank.tl depreciated recently???
import cellrank as cr

g = cr.estimators.GPCCA(kernel)

g.compute_schur() #it was 20 components in the example. This is the part where you'd want to use KRYLOV instead for very large sample counts
g.plot_spectrum(real_only=True)
g.compute_macrostates(n_states=2, n_cells=2, cluster_key="clusters") #if you have a "clusters" annotation in adata_abund obs (classified as health and sick, for instance), put it here as cluster_key="clusters"
