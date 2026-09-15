import Statescope.Statescope as scope
import pandas as pd
import anndata
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pickle 

Bulk = pd.read_csv(
    "/DIT/WITH/GENECOUNT/MATRIX/gene_counts_matrix.csv",
    sep=',',
    index_col='symbol')

#Drop dupicates, I had this problem with rRNAs
Bulk = Bulk[~Bulk.index.duplicated(keep="first")]

Tumor_purities = pd.read_csv('/DIT/WITH/TUMORPURITY/MATRIX/tumor_purity_absolute.csv', sep = ',', index_col = 'Sample')


##Match the tumorpurity samples with the bulk seq sample
Tumor_purities.index = ['-'.join(x.split('-')[0:3]) for x in Tumor_purities.index]
Bulk = Bulk[Tumor_purities.index]

ss = scope.Initialize_Statescope(
    Bulk,
    TumorType="PDAC",
    Ncores=16
)
    
Expectation_single = pd.DataFrame(np.nan, index=ss.Samples, columns=ss.Celltypes)

Expectation_single.loc[:,'Epithelial']  = Tumor_purities.loc[ss.Samples,'ABSOLUTE Purity']
#Purity of 1 is not alowed, to clip to 0.99
Expectation_single = Expectation_single.clip(upper=0.99)

##This took about 3 hours running on a RTX5080 and 32 CPU threads. RAM usage was suprisingly low
deconv_res = ss.Deconvolution(Expectation=Expectation_single) 

#Save the ss object:

#Refine the results (don't know why, manual says so). Need to figure it out also ran it twice
refined_res = ss.Refinement()

#saved with
with open("/home/bram/Documents/Deconvolution_PDAC_nicole/Deconvalution results/Statescope_PDAC_refined.pkl", "wb") as f: pickle.dump(ss, f) 
#opened with
with open("/home/bram/Documents/Deconvolution_PDAC_nicole/Deconvalution results/Statescope_PDAC_refined.pkl", "rb") as f: ss = pickle.load(f) 

##folowed by running state discovery I also don't know exactly how this works.
ss.StateDiscovery()
