
# -*- coding: utf-8 -*-
"""Kernels.

Kernel base (interface) class and implementing classes.
"""

from .config import config
import random, numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import multiprocessing as mp
import itertools as it
import networkx as nx

# Kernels
class Kernel():
    def __init__(self,netwink,kernelspec={}):
        """
        netwink: the network object with which the kernel is associated
        kernelspec: dict with parameters for kernel computation
        """
        self.netwink = netwink
        self.spec = kernelspec

    def __repr__(self):
        if 'computedMatrix' in vars(self):
            return self.computedMatrix.__repr__()

    def __str__(self):
        return '{}_{}_{}'.format(
            self.netwink.name,
            self.__class__.__name__,
            '_'.join(['{}_{}'.format(k,self.spec[k]) for k in sorted(self.spec)]) if self.spec else 'default'
        )

    def __hash__(self):
        return hash(str(self))

    def compute(self):
        """
        sets self.computedMatrix
        returns self
        """
        raise NotImplementedError('implement in inheriting classes')

    def score_geneset(self,geneset,method='vector',ax=None,c='r'):
        """
        Calculate the sum of a geneset.
        netwink.apply_gene_scores(genescores) has to been run before
        
        The function spreads the impact of the genes in the geneset of interest
        by constructing a vector as big as the network (1 for gene of interest, 0 otherwise),
        elementwise multpication with kernel computed matrix
        sum defines the subnetwork, but permutation test is needed to test how relevant
        its sum is, compared to scores randomly same-sized genesets get.
        """
        if method == 'vector':
            diffusion =  self.netwink.scoresVector @ self.computedMatrix
            #score matrix repeat tot size computed matrix, nadien col sums resultaat (diffusie matrix moet wel row sums 1 hebben)
            #genormaliseerde correlatie waarden als gewichten
            score = diffusion.T[self.netwink.get_nodes_series().isin(geneset)].sum()
        elif method == 'matrix':
            score = np.multiply(
                self.computedMatrix,
                np.multiply(self.netwink.scoresMatrix,self.netwink.subset_adjmatrix(geneset))
            ).sum()
        if ax: ax.axvline(score,c=c)
        return score

    @staticmethod
    def mp_score_geneset(geneset,diffusion):
        """Static method for calculating geneset score

        For multiprocessing use

        Args:
            geneset (pd.Series): boolean pd.Series indicating which scoresVector genes are in the geneset
            scoresVector (matrix): vector with scores
            computedMatrix (matrix): diffusion kernel

        Returns:
            score (float)
        """
        score = diffusion.T[geneset].sum()
        return score

    def plot_scores_diffusion_scatter(self):
        """Scatter plot of original scores and diffused scores
        """
        fig, ax = plt.subplots()
        ax.scatter(
            np.array(self.netwink.scoresVector[:,None]),
            np.array((self.netwink.scoresVector @ self.computedMatrix).T)
        )
        return ax

    def plot_geneset_scores(self,geneset,filename=None,cmap='hot_r',border_cmap='hot_r',border_scores=None,penwidth=10):
        """Plot a geneset network with Graphiz

        Args:
            geneset (set): Set of genes that will be plotted.
            filename (str): filename path, should end in '.svg'.
            cmap (str): matplotlib colormap name.
            border_cmap (str): colormap for the border.
            border_scores (pd.Series): Instead of the original scores, 
                put in a different metric for border color.
            penwidth (int): the size of the border line around the gene.

        Todo:
            * add gradient legend -> dot fillcolor="orange:yellow"
              or generating color bar legend with matplotlib and then referencing that image
        """
        # Collect data
        nodes =  pd.Index(self.netwink.get_nodes_series())
        originalScores = self.netwink.scoresVector
        diffusedScores = (originalScores @ self.computedMatrix).T
        # Setup colors
        import pydot
        import matplotlib as mpl
        from matplotlib.colors import rgb2hex
        from bidali.visualizations import labelcolor_matching_backgroundcolor
        vmin = originalScores.min()
        vmax = originalScores.max()
        norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
        cmap = plt.get_cmap(cmap)
        # Setup graph
        subgraph = self.netwink.graph.subgraph(geneset)
        dotgraph = nx.drawing.nx_pydot.to_pydot(subgraph)
        for n in dotgraph.get_nodes():
            gene = n.get_name()
            gene_loc = nodes.get_loc(gene)
            n.set_penwidth(10)
            n.set_style('filled')
            fillcolortuple = cmap(norm(diffusedScores[gene_loc,0]))
            n.set_fillcolor(rgb2hex(fillcolortuple))
            n.set_color(rgb2hex(cmap(norm(originalScores[gene_loc]))))
            n.set_fontcolor(rgb2hex(labelcolor_matching_backgroundcolor(fillcolortuple)))
            
        # Create color legend
        graphlegend = pydot.Cluster(
            graph_name="legend", label="Color legend", fontsize="15", color="blue",
            style="filled", fillcolor="lightgrey", rankdir="TB"
        )
        minnode = pydot.Node('min', label="Min: {:.2f}".format(vmin), style="filled",
                             fontcolor=rgb2hex(labelcolor_matching_backgroundcolor(cmap(norm(vmin)))),
                             fillcolor=fillcolor=rgb2hex(cmap(norm(vmin))),
                             shape="Mrecord", rank="same"
        )
        graphlegend.add_node(minnode)
        maxnode = pydot.Node('max', label="Max: {:.2f}".format(vmax), style="filled",
                             fontcolor=rgb2hex(labelcolor_matching_backgroundcolor(cmap(norm(vmax)))),
                             fillcolor=rgb2hex(cmap(norm(vmax))), shape="Mrecord", rank="same"
        )
        graphlegend.add_node(maxnode)
        graphlegend.add_edge(pydot.Edge(minnode, maxnode, style="invis"))
        dotgraph.add_subgraph(graphlegend)
        
        # Graph output
        if filename:
            dotgraph.write_svg(filename)
        else:
            import io
            import matplotlib.image as mpimg
            from PIL import Image
            png = io.BytesIO(dotgraph.create_png())
            img = Image.open(png)
            #img.thumbnail(resize, Image.ANTIALIAS)
            img.show()
            #imgplot = plt.imshow(mpimg.imread(png))
            #imgplot.get_figure().axes[0].axis('off')
            #return imgplot

    def plot_gene_neighbourhood(self,gene,depth=2,**kwargs):
        """Plot network centered around a gene

        Args:
            gene (str): Gene name.
            depth (int): Connection distance with gene.
        """
        geneset = {gene} if isinstance(gene, str) else set(gene)
        for level in range(depth):
            for g in geneset.copy():
                geneset.update(set(self.netwink.graph.neighbors(g)))
        self.plot_geneset_scores(geneset,**kwargs)
        
    def permutate_geneset_scores(self,geneset,numberOfPermutations=1000,seed=None,ax=None,multiprocessing=True):
        """
        Takes the length of the passed geneset, then makes `numberOfPermutations` genesets
        of the same length and calculates their scores. Returns a list with the scores.

        seed can be used to (re)set the random state

        Args:
            geneset (set): Set of genes.
            numberOfPermutations (int): Number of similarly sized random geneset calculations.
            seed (float): Random seed.
            ax (plt.Axes): Ax to plot results on.
            multiprocessing (bool): Use a pool of processes to calculate random scores.
                Number of processes configurable in pyni config.
        """
        genesetLength = len(geneset)
        allGenes = set(self.netwink.get_nodes_series())
        if seed: random.seed(seed)
        shelveKey = 'Geneset length {} - permutations {} - seed {} - kernel hash {}'.format(
            genesetLength,numberOfPermutations,seed,hash(self)
            #optionally also include a hash of the netwink nodes+edges to be sure it did not change
        )
        try: scores = self.netwink.shelve[shelveKey]
        except KeyError:
            if multiprocessing:
                with mp.Pool(processes=int(config['pyni']['threads'])) as pool:
                    scores = pool.starmap(
                        Kernel.mp_score_geneset,
                        zip(
                            (self.netwink.get_nodes_series().isin(random.sample(allGenes,genesetLength)) for i in range(numberOfPermutations)),
                            it.repeat(self.netwink.scoresVector @ self.computedMatrix)
                        )
                    )
            else:
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
    def compute(self):
        """
        Spec parameters:
        alpha (default: 1)
        """
        # Set parameters
        alpha = 1 if 'alpha' not in self.spec else self.spec['alpha']
        # Compute
        from scipy.linalg import expm
        self.degreematrix = np.diag(np.array(self.netwink.admatrix.sum(axis=1)).flatten())
        self.laplacian = self.degreematrix - self.netwink.admatrix
        self.computedMatrix = expm(-alpha*np.array(self.laplacian)) # => wel een kernel
        return self

## restart random walk
class RestartRandomWalk(Kernel):
    def compute(self):
        """
        computed restartmatrix is not symmetrical

        Spec parameters:
        restartProb (default: 0.1)
        """
        # Set parameters
        restartProb = 0.1 if 'restartProb' not in self.spec else self.spec['restartProb']
        # Compute
        self.degreematrix = np.diag(np.array(self.netwink.admatrix.sum(axis=1)).flatten())
        #self.inversematrix = np.linalg.inv(self.degreematrix - ((1-restartProb)*self.netwink.admatrix))
        self.inversematrix = np.linalg.pinv(self.degreematrix - ((1-restartProb)*self.netwink.admatrix)) #pseudo inverse
        self.restartmatrix = self.inversematrix*self.degreematrix # => op factor na probabiliteiten => geen kernel
        # Normalize row sums to 1 => this way each row shows how much a gene keeps from its score and gets from other genes' scores
        self.restartmatrix = self.restartmatrix / self.restartmatrix.sum(axis=1) 
        self.computedMatrix = self.restartmatrix
        return self
