#/usr/bin/env python3
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
import pandas as pd
from io import  StringIO
import os

# Set matplotlib to interactive
plt.ion()

## Biomart
def get_biomart(atts=None):
    """
    Get biomart id mappings

    >>> bmrt = get_biomart()
    """
    import biomart
    server = biomart.BiomartServer("http://www.ensembl.org/biomart")
    # default atts
    if not atts:
        atts = ['external_gene_name','external_gene_source','ensembl_gene_id',
                'ensembl_transcript_id','ensembl_peptide_id']
    hge = server.datasets['hsapiens_gene_ensembl']
    s = hge.search({'attributes': atts}, header=1)
    data = pd.read_table(StringIO(s.content.decode()))
    return data

b = get_biomart()
b=b[b['Source of gene name']=='HGNC Symbol']
g=b['Gene name']
b=b[~g.duplicated()]
#ask server to not send duplicate ids
b['network_id'] = b['Gene stable ID'].apply(lambda x: int(x.replace('ENSG','')))

# Load network
# als complex netwerk -> complex als entiteit meenemen?
# voor brute analyze best niet
# TODO hoe gaan we daarmee om
networktable = pd.read_table(os.path.expanduser('~/Dropbiz/Personal/Tools/ninklings/data/reactome_FI.txt'))
networktable = networktable[networktable.Gene1.isin(b['Gene name'])]
networktable = networktable[networktable.Gene2.isin(b['Gene name'])]
networktable = networktable[networktable.Score >= .75]
print('unique genes:',len(set(networktable.Gene1)|set(networktable.Gene2)))
networkgenes = networktable[['Gene1','Gene2']]
networkgenes.to_csv(os.path.expanduser('~/Dropbiz/Personal/Tools/ninklings/data/reactome_FI_filteredEdges.txt'),
                    sep='\t')

cosmicgenes = set(pd.read_table(os.path.expanduser('/Users/cvneste/Dropbiz/Personal/Tools/ninklings/data/cosmic_20180125.tsv'))['Gene Symbol'])

reactome_FI = nx.Graph()
networkgenes.T.apply(lambda x: reactome_FI.add_edge(x['Gene1'],x['Gene2']))

admatrix = nx.adjacency_matrix(reactome_FI)

reactome_FI_cosmic = reactome_FI.subgraph(cosmicgenes)

for comp in nx.connected_component_subgraphs(reactome_FI_cosmic):
    print(comp.number_of_nodes())
    if comp.number_of_nodes() == 534:
        reactome_FI_cosmic = comp
        break

admatrix = nx.adjacency_matrix(reactome_FI_cosmic)
admatrix = admatrix.todense()
admatrix = admatrix.astype(np.float32)
plt.matshow(admatrix)

#laplacian exponential diffusion kernel
degreematrix = np.diag(np.array(admatrix.sum(axis=1)).flatten())
laplacian = degreematrix - admatrix
alpha = 1
expm(-alpha*np.array(laplacian)) # => wel een kernel
admatrix.tofile('Downloads/admatrix.txt', sep='\t')

#restart random walk
restartProb = 0.1
inversematrix = np.linalg.inv(degreematrix - ((1-restartProb)*admatrix))
restartmatrix = inversematrix*degreematrix # => op factor na probabiliteiten => geen kernel
plt.matshow(restartmatrix)
#vertrekken van gen dat e.g. gemuteerd is (restartmatrix is niet symmetrisch)
# restartmatrix -> converteren naar afstandsmatrix (1 - ) en diagonaal op 0 zetten
# daarmee kan je clusteren -> groepen van genen in elkaars buurt
# goed ter controle, enrichment ook controleren -> biologische functie

#genes of interest impact verspreiden
# vector bouwen even groot als netwerk (1 voor gene of interest, 0 anders),
# elementwise vermenigvuldigen met restartmatrix of expm
# sum -> definieert subnetwerk
