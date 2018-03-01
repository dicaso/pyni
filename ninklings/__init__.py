#/usr/bin/env python3
import networkx as nx
import numpy as np
from scipy.linalg import expm
import pandas as pd
from io import  StringIO
import os, logging, random, shelve
import matplotlib.pyplot as plt
import seaborn as sns
from bidali.LSD.dealer.external.ensembl import get_biomart

# Set matplotlib to interactive mode when executed from interactive shell
if 'ps1' in vars(os.sys): plt.ion()

# Configuration
from .config import config

class Netwink():
    """
    Ninkling's starting class.
    Does the setup for all down stream network analyses
    """

    def __init__(self,savedir,name=None,loglevel='DEBUG',cosmicOnly=True,fullInit=True):
        self.location = savedir
        self.name = name if name else os.path.basename(savedir)
        self.logger = logging.Logger(self.name,level=loglevel)
        if fullInit:
            self.check_location()
            self.setup_logging(loglevel)
            self.shelve = shelve.open(os.path.join(self.location,'distrostore'))
            self.set_annotation()
            self.load_network(cosmicOnly)
            self.logger.debug('%s fully created',self)
        else: self.logger.warning('%s created but no network loaded',self)

    def __del__(self):
        self.shelve.close()

    def __repr__(self):
        return '<Netwink object: {}>'.format(self.name)

    def check_location(self):
        if not os.path.exists(self.location):
            os.mkdir(self.location)

    def setup_logging(self,level='DEBUG'):
        stdout = logging.StreamHandler(stream=os.sys.stdout)
        stdout.setLevel(level)
        logfile = logging.FileHandler(os.path.join(self.location,self.name+'.log'))
        logfile.setLevel(level)
        self.logger.addHandler(stdout)
        self.logger.addHandler(logfile)
        
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
        # if complexes in network -> condir complex as its own entity?
        # for brute force analysis best not
        # TODO how are we going to deal with it?
        self.networkFile = os.path.join(self.location, 'reactome_FI_filteredEdges.tsv')
        if os.path.exists(self.networkFile):
            self.networkgenes = pd.read_table(self.networkFile)
        else:
            networktable = pd.read_table(os.path.join(config['ninklings']['datadir'],'reactome_FI.txt'))
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
                os.path.join(config['ninklings']['datadir'],'cosmic_20180125.tsv')
            )['Gene Symbol'])
            self.graph = self.graph.subgraph(cosmicgenes)
            components = {
                comp.number_of_nodes():comp
                for comp in nx.connected_component_subgraphs(self.graph)
            }
            self.graph = components[max(components)]
        self.admatrix = nx.adjacency_matrix(self.graph).todense()
        self.admatrix = self.admatrix.astype(np.float32)
        # Set weights to default of 1
        self.weights = np.ones(self.admatrix.shape)

    def load_weights(self,weights):
        """
        Load weights such as context relevant correlation values

        TODO what to do with NA weights?
        """
        self.weights = weights

    def load_correlations(self,nodesCorrelationData,method='pearson',expSmoothening=1,**kwargs):
        """
        Calculate correlations based on the data in nodesCorrelationData pd.DataFrame
        and load them as weights

        if nodesCorrelationData is not pd.DataFrame, it is assumed to be a filename,
        that will be used to load the DataFrame with pd.read_csv and kwargs provided.

        correlation method => any method accepted by DataFrame.corr (pearson, kendall, spearman)

        expSmoothening => all corr values are transformed to the power of expSmoothening
          values furhter away from 1 will go exponentially towards 0 with higher expSmoothening value
        """
        if type(nodesCorrelationData) != pd.DataFrame:
            nodesCorrelationData = pd.read_csv(nodesCorrelationData,**kwargs)
        # Filter non-node data
        nodesCorrelationData = nodesCorrelationData.reindex(self.get_nodes_series())
        correlations = nodesCorrelationData.T.corr(method=method)
        self.correlationSigns = np.sign(correlations)
        weights = correlations.abs() ** expSmoothening
        self.load_weights(weights)
        

    def get_nodes_series(self):
        return pd.Series(self.graph.nodes)

    def export_admatrix(self):
        self.admatrix.tofile(os.path.join(self.location,'admatrix.tsv'), sep='\t')

    def plot_admatrix(self):
        plt.matshow(self.admatrix)

    def apply_gene_scores(self,genescores,fillna=0):
        """
        Apply gene scores to the adjacency matrix. 
        Genescores should be dict with as keys the node names.
        When a node/gene is not in genescores it gets a score of `fillna`

        In the returned matrix, each gene connected to the row's gene
        gets the score of that gene. The returned matrix is not symmetrical.
        """
        scoresDiag = np.matrix(np.diag(
            self.get_nodes_series().map(genescores).fillna(fillna).astype(np.float32).as_matrix()
        ))
        self.scoresMatrix = scoresDiag * np.matrix(self.admatrix)
        
    
    def subset_adjmatrix(self,geneset,externalEdges=False):
        """
        For a given set of genes, returns the original adjancency matrix
        with 1 indicating edges within the geneset, or if externalEdges is
        True, edges between 1 of the genes of interest and any other gene.
        Genes that are not in the geneset or not connected to a geneset gene
        if externalEdges is True are 0.
        """
        return np.multiply(
            self.admatrix,
            self.get_geneset_selection_matrix(geneset,externalEdges)
        )

    def get_geneset_selection_matrix(self,geneset,externalEdges=False):
        """
        For a given set of genes, returns a matrix with row and columns
        1 for genes of the set, and 0 for all others if externalEdges is True.
        If externalEdges is False only 1 for combinations of genes within
        the geneset, and 0 for all others.
        """
        genesetVector = self.get_nodes_series().isin(geneset).astype(np.float32).as_matrix()
        if externalEdges:
            return np.logical_or(
                np.outer(genesetVector,np.ones(len(genesetVector))),
                np.outer(np.ones(len(genesetVector)),genesetVector)
            ).astype(np.float32)
        else:
            return np.outer(genesetVector,genesetVector)

# Kernels
class Kernel():
    def __init__(self,netwink):
        self.netwink = netwink

    def __repr__(self):
        if 'computedMatrix' in vars(self):
            return self.computedMatrix.__repr__()

    def compute(self):
        """
        sets self.computedMatrix
        returns self
        """
        raise NotImplementedError('implement in inheriting classes')

    def score_geneset(self,geneset,ax=None):
        """
        Calculate the sum of a geneset.
        netwink.apply_gene_scores(genescores) has to been run before
        
        The function spreads the impact of the genes in the geneset of interest
        by constructing a vector as big as the network (1 for gene of interest, 0 otherwise),
        elementwise multpication with kernel computed matrix
        sum defines the subnetwork, but permutation test is needed to test how relevant
        its sum is, compared to scores randomly same-sized genesets get.
        """
        score = np.multiply(
            self.computedMatrix,
            np.multiply(self.netwink.scoresMatrix,self.netwink.subset_adjmatrix(geneset))
        ).sum()
        if ax: ax.axvline(score,c='r')
        return score

    def permutate_geneset_scores(self,geneset,numberOfPermutations=1000,seed=None,ax=None):
        """
        Takes the length of the passed geneset, then makes `numberOfPermutations` genesets
        of the same length and calculates their scores. Returns a list with the scores.

        seed can be used to (re)set the random state
        """
        genesetLength = len(geneset)
        allGenes = set(self.netwink.get_nodes_series())
        if seed: random.seed(seed)
        shelveKey = 'Geneset length {} - permutations {} - seed {}'.format(
            genesetLength,numberOfPermutations,seed
            #optionally also include a hash of the netwink object to be sure it did not change
        )
        try: scores = self.netwink.shelve[shelveKey]
        except KeyError:
            scores = [
                self.score_geneset(random.sample(allGenes,genesetLength))
                for i in range(numberOfPermutations)
            ]
            self.netwink.shelve[shelveKey] = scores
            #logger.
        if ax:
            sns.distplot(scores,ax=ax)
        return scores
        
    def visualize(self):
        plt.matshow(self.computedMatrix)

    def convert_to_distance(self):
        """
        convert a kernel computed matrix to a distance matrix by subtrating from 1
        and equaling the diagonale positions to 0 distance
        
        the distance matrix can be used for clustering related genes
        this is good to use in a check, e.g. checking enrichment of biological functions
        """
        # Take average values in case computedMatrix is not symmetrical
        averageMatrix = (self.computedMatrix + self.computedMatrix.T)/2
        # Substract from 1 to get distance instead of connectedness and set diagonal distances to 0
        distanceMatrix = np.multiply(
            1 - averageMatrix,
            np.logical_not(np.matrix(np.diag(np.ones(len(averageMatrix)))))
        )
        return distanceMatrix
        
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
        """
        computed restartmatrix is not symmetrical
        """
        self.degreematrix = np.diag(np.array(self.netwink.admatrix.sum(axis=1)).flatten())
        self.inversematrix = np.linalg.inv(self.degreematrix - ((1-restartProb)*self.netwink.admatrix))
        self.restartmatrix = self.inversematrix*self.degreematrix # => op factor na probabiliteiten => geen kernel
        self.computedMatrix = self.restartmatrix
        return self
