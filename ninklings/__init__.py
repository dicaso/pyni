#/usr/bin/env python3
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
import pandas as pd
from io import  StringIO
import os

# Set matplotlib to interactive mode when executed from interactive shell
if 'ps1' in vars(os.sys): plt.ion()

# Package data
pkgdatadir = '{}/../data/%s'.format(os.path.dirname(__file__)) #TODO change to subdir ninklings when making public

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

class Netwink():
    """
    Ninkling's starting class.
    Does the setup for all down stream network analyses
    """

    def __init__(self,savedir,name=None,cosmicOnly=True,fullInit=True):
        self.location = savedir
        self.name = name if name else os.path.basename(savedir)
        if fullInit:
            self.check_location()
            self.set_annotation()
            self.load_network(cosmicOnly)

    def check_location(self):
        if not os.path.exists(self.location):
            os.mkdir(self.location)
    
    def set_annotation(self):
        self.annotationFile = os.path.join(self.location, 'ensembl_annotation.tsv')
        if os.path.exists(self.annotationFile):
            self.annotation = pd.read_table(self.annotationFile)
        else:
            self.annotation = get_biomart()
            self.annotation = self.annotation[
                self.annotation['Source of gene name']=='HGNC Symbol'
            ]
            self.annotation = self.annotation[
                ~self.annotation['Gene name'].duplicated()
            ] #TODO ask server to not send duplicate ids
            self.annotation['network_id'] = self.annotation['Gene stable ID'].apply(lambda x: int(x.replace('ENSG','')))
            self.annotation.to_csv(self.annotationFile,sep='\t')

    def load_network(self,cosmicOnly=False):
        # als complex netwerk -> complex als entiteit meenemen?
        # voor brute analyze best niet
        # TODO hoe gaan we daarmee om
        self.networkFile = os.path.join(self.location, 'reactome_FI_filteredEdges.tsv')
        if os.path.exists(self.networkFile):
            self.networkgenes = pd.read_table(self.networkFile)
        else:
            networktable = pd.read_table(os.path.expanduser(pkgdatadir % 'reactome_FI.txt'))
            networktable = networktable[networktable.Gene1.isin(self.annotation['Gene name'])]
            networktable = networktable[networktable.Gene2.isin(self.annotation['Gene name'])]
            networktable = networktable[networktable.Score >= .75]
            print('unique genes:',len(set(networktable.Gene1)|set(networktable.Gene2)))
            self.networkgenes = networktable[['Gene1','Gene2']]
            self.networkgenes.to_csv(self.networkFile, sep='\t')
        # Construct network
        self.graph = nx.Graph()
        self.networkgenes.T.apply(lambda x: self.graph.add_edge(x['Gene1'],x['Gene2']))
        if cosmicOnly:
            cosmicgenes = set(pd.read_table(
                os.path.expanduser(pkgdatadir % 'cosmic_20180125.tsv')
            )['Gene Symbol'])
            self.graph = self.graph.subgraph(cosmicgenes)
            components = {
                comp.number_of_nodes():comp
                for comp in nx.connected_component_subgraphs(self.graph)
            }
            self.graph = components[max(components)]
        self.admatrix = nx.adjacency_matrix(self.graph).todense()
        self.admatrix = self.admatrix.astype(np.float32)

    def export_admatrix(self):
        self.admatrix.tofile(os.path.join(self.location,'admatrix.tsv'), sep='\t')

    def plot_admatrix(self):
        plt.matshow(self.admatrix)

# Kernels
class Kernel():
    def __init__(self,netwink):
        self.netwink = netwink

    def __repr__(self):
        if 'computedMatrix' in vars(k):
            return self.computedMatrix.__repr__()

    def compute(self):
        raise NotImplementedError('implement in inheriting classes')

    def visualize(self):
        plt.matshow(self.computedMatrix)
        
## laplacian exponential diffusion kernel
class ExponentialDiffusionKernel(Kernel):
    def compute(self,alpha = 1):
        from scipy.linalg import expm
        self.degreematrix = np.diag(np.array(self.netwink.admatrix.sum(axis=1)).flatten())
        self.laplacian = self.degreematrix - self.netwink.admatrix
        self.computedMatrix = expm(-alpha*np.array(self.laplacian)) # => wel een kernel
        return self

## restart random walk
class RestartRandomWalk(Kernel):
    def compute(self,restartProb = 0.1):
        self.degreematrix = np.diag(np.array(self.netwink.admatrix.sum(axis=1)).flatten())
        self.inversematrix = np.linalg.inv(self.degreematrix - ((1-restartProb)*self.netwink.admatrix))
        self.restartmatrix = self.inversematrix*self.degreematrix # => op factor na probabiliteiten => geen kernel
        self.computedMatrix = self.restartmatrix
        return self
        
#vertrekken van gen dat e.g. gemuteerd is (restartmatrix is niet symmetrisch)
# restartmatrix -> converteren naar afstandsmatrix (1 - ) en diagonaal op 0 zetten
# daarmee kan je clusteren -> groepen van genen in elkaars buurt
# goed ter controle, enrichment ook controleren -> biologische functie

#genes of interest impact verspreiden
# vector bouwen even groot als netwerk (1 voor gene of interest, 0 anders),
# elementwise vermenigvuldigen met restartmatrix of expm
# sum -> definieert subnetwerk
