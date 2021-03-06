
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

    def score_geneset(self,geneset,method='vector',signedDiffusion=False,ax=None,c='r'):
        """
        Calculate the sum of a geneset.
        netwink.apply_gene_scores(genescores) has to been run before
        
        The function spreads the impact of the genes in the geneset of interest
        by constructing a vector as big as the network (1 for gene of interest, 0 otherwise),
        elementwise multpication with kernel computed matrix
        sum defines the subnetwork, but permutation test is needed to test how relevant
        its sum is, compared to scores randomly same-sized genesets get.

        Todo:
            * use correlation signs for multiplying when scores have relevant signs
            * check passing signedDiffusion upstream
        """
        if method == 'vector':
            if signedDiffusion:
                diffusion =  self.netwink.scoresVector @ np.multiply(self.computedMatrix,self.netwink.correlationSigns)
            else: diffusion =  self.netwink.scoresVector @ self.computedMatrix
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
    def mp_diffuse_score_geneset(geneset,scores,computedMatrix,nperm=1000):
        """Static method for calculating diffused geneset score.
        Opposing to mp_score_geneset this function expects the original computedMatrix
        and then calculates diffused scores.

        Ready for multiprocessing use.

        Args:
            geneset (pd.Series): boolean pd.Series indicating which scoresVector genes are in the geneset
            scoresVector (matrix): vector with scores
            computedMatrix (matrix): diffusion kernel
            nperm (int): If nperm, calculate random permutations of scoresVector

        Returns:
            score (float)
        """
        diffusion = scores @ computedMatrix
        score = diffusion.T[geneset].sum()
        if nperm:
            randomScoresVector = scores.copy()
            random_scores = [
                (randomScoresVector @ computedMatrix).T[geneset].sum()
                for i in range(nperm)
                # random shuffling for every iteration as a side effect:
                if np.random.shuffle(randomScoresVector) or True
            ]
            return score, random_scores
        else: return score

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

    def plot_scores_diffusion_scatter(self,diffusionDiff=False,returnDiff=False,annotate=0.05,arrowColor='red',figsize=(7,7)):
        """Scatter plot of original scores and diffused scores

        Args:
            diffusionDiff (bool): Show the difference with original score
                after diffusion, instead of the diffused score.
            annotate (float): Float between 0-1 indicating upper and lower 
                quantile of diffusion difference that will be annotated.
            arrowColor (color): Color of the annotation arrows.
            figsize (float,float): Figure size passed to plt.subplots.
        """
        fig, ax = plt.subplots(figsize=figsize)
        x = np.array(self.netwink.scoresVector[:,None])
        y = np.array((self.netwink.scoresVector @ self.computedMatrix).T)
        if diffusionDiff:
            y = y - x
        ax.scatter(x, y)
        if diffusionDiff and annotate:
            annotdata = pd.Series(y.flatten(), index=self.netwink.get_nodes_series())
            annotbool = (annotdata > annotdata.quantile(1-annotate)) | (annotdata < annotdata.quantile(annotate))
            textBoxes = [
                ax.text(self.netwink.scoresVector[i],annotdata[i], annotdata.index[i], ha='center', va='center')
                #ax.annotate(annotdata.index[i],(self.netwink.scoresVector[i],annotdata[i]))
                for i,a in enumerate(annotbool) if a
            ]
            from adjustText import adjustText #adjust_text needs to be last action on fig before visualization
            adjustText.adjust_text(
                textBoxes,x=self.netwink.scoresVector[annotbool],y=annotdata[annotbool],
                ax=ax, arrowprops=dict(arrowstyle='->', color=arrowColor)
            )
        return ax if not returnDiff else (ax,y)

    def plot_geneset_scores(
            self, geneset, filename=None, cmap='Greens', vmin=None, vmax=None,
            border_cmap=None, border_scores=None, penwidth=10, legend_title=None,
            color_edges = True, edge_cmap='RdGy', edge_filter=.2, makeColorbars=False,
            colorbar_title=None, border_colorbar_title=None, dotprog = 'dot'
    ):
        """Plot a geneset network with Graphiz

        Args:
            geneset (set): Set of genes that will be plotted.
            filename (str): filename path, should end in '.svg'.
            cmap (str): matplotlib colormap name.
            vmin (float): minimum value for cmap normalization.
                lower scores will not be distinguishable.
            vmax (float): maximum value for cmap normalization, 
                higher scores will not be distinguishable.
            border_cmap (str): colormap for the border.
            border_scores (pd.Series): Instead of the original scores, 
                put in a different metric for border color.
            penwidth (int): the size of the border line around the gene.
            legend_title (str): Title for legend.
            color_edges (bool): Use netwink correlation data for coloring edges.
            edge_cmap (str): If color_edges, which cmap to use.
            edge_filter (float): If color_edges, edges with absolute correlation
                smaller than edge_filter will be filtered. Allows to make more
                comprehensible graphs.
            makeColorbars (bool): If true, make matplotlib figures with the colorbars used.

        Todo:
            * add gradient legend -> dot fillcolor="orange:yellow"
              or generating color bar legend with matplotlib and then referencing that image
            * look into other (python) graph visualizers such as igraph, graph-tool
        """
        # Collect data
        nodes =  pd.Index(self.netwink.get_nodes_series())
        originalScores = self.netwink.scoresVector
        diffusedScores = (originalScores @ self.computedMatrix).T
        # Setup colors
        import pydot
        import matplotlib as mpl, tempfile
        from matplotlib.colors import rgb2hex
        from bidali.visualizations import labelcolor_matching_backgroundcolor
        vmin = vmin if vmin else originalScores.min()
        vmax = vmax if vmax else originalScores.max()
        norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
        cmap = plt.get_cmap(cmap)
        fill_color = lambda x: cmap(norm(x))
        if border_scores is not None:
            border_vmin = border_scores.min()
            border_vmax = border_scores.max()
            border_norm = mpl.colors.Normalize(vmin=border_vmin,vmax=border_vmax)
            border_cmap = plt.get_cmap(border_cmap) if border_cmap else cmap
            border_color = lambda x: border_cmap(border_norm(x))
        else:
            border_color = fill_color
        if makeColorbars:
            # Reference: https://matplotlib.org/examples/api/colorbar_only.html
            cbarfig, cbax = plt.subplots(
                nrows=1 if border_scores is None else 2, ncols=1, figsize=(5,2)
            )
            cbarfig.tight_layout()
            if border_scores is not None:
                cbax, cboax = cbax
                mpl.colorbar.ColorbarBase(
                    cboax, cmap=border_cmap, norm=border_norm, orientation='horizontal'
                )
                cboax.set_xlabel(border_colorbar_title or 'Border score colorbar')
            mpl.colorbar.ColorbarBase(cbax, cmap=cmap, norm=norm, orientation='horizontal')
            cbax.set_xlabel(colorbar_title or 'Center score colorbar')
            # Write out to tempfile
            tf = tempfile.NamedTemporaryFile(suffix='.png')
            print(tf.name)
            tf.close()
            cbarfig.savefig(tf.name)
            plt.close(cbarfig) #closing figure as downstream only file ref needed
            graphlegend = pydot.Cluster(
                graph_name="legend", label=legend_title if legend_title else "Color legend",
                fontsize="15", color="blue", style="filled", fillcolor="lightgrey", rankdir="LR"
            )
            colorbarnode = pydot.Node('colorbars', label='', shape="box", image=tf.name)
            graphlegend.add_node(colorbarnode)
        # Setup graph
        subgraph = self.netwink.graph.subgraph(geneset)
        dotgraph = nx.drawing.nx_pydot.to_pydot(subgraph)
        for n in dotgraph.get_nodes():
            gene = n.get_name()
            gene_loc = nodes.get_loc(gene)
            n.set_penwidth(10)
            n.set_style('filled')
            fillcolortuple = fill_color(diffusedScores[gene_loc,0])
            n.set_fillcolor(rgb2hex(fillcolortuple))
            n.set_color(rgb2hex(
                border_color(
                    (originalScores if border_scores is None else border_scores)[gene_loc]
                )
            ))
            n.set_fontcolor(rgb2hex(labelcolor_matching_backgroundcolor(fillcolortuple)))
            # If there are negative scores, set different shape
            if diffusedScores[gene_loc,0] < 0:
                n.set_shape('diamond')
        if color_edges:
            for e in dotgraph.get_edges():
                edge_norm = mpl.colors.Normalize(vmin=-1,vmax=1)
                edge_cmap = plt.get_cmap(edge_cmap)
                edge_color = lambda x: edge_cmap(edge_norm(x))
                try:
                    corr = self.netwink.correlations[e.get_source()][e.get_destination()]
                    if abs(corr) < edge_filter:
                        dotgraph.del_edge(e.get_source(),e.get_destination())
                    else:
                        e.set_color(rgb2hex(edge_color(corr)))
                        if corr < 0: e.set_style('dashed')
                except KeyError:
                    # No correlation info as source or target not in expression data
                    e.set_style('dotted')
            
        # Create color legend
        #graphlegend = pydot.Cluster(
        #    graph_name="legend", label=legend_title if legend_title else "Color legend",
        #    fontsize="15", color="blue", style="filled", fillcolor="lightgrey", rankdir="LR"
        #)
        #minnode = pydot.Node('min', label="Min: {:.2f}".format(vmin), style="filled",
        #                     fontcolor=rgb2hex(labelcolor_matching_backgroundcolor(cmap(norm(vmin)))),
        #                     fillcolor=rgb2hex(cmap(norm(vmin))),
        #                     shape="Mrecord", rank="same"
        #)
        #graphlegend.add_node(minnode)
        #maxnode = pydot.Node('max', label="Max: {:.2f}".format(vmax), style="filled",
        #                     fontcolor=rgb2hex(labelcolor_matching_backgroundcolor(cmap(norm(vmax)))),
        #                     fillcolor=rgb2hex(cmap(norm(vmax))), shape="Mrecord", rank="same"
        #)
        #graphlegend.add_node(maxnode)
        #graphlegend.add_edge(pydot.Edge(minnode, maxnode, style="invis"))
        dotgraph.add_subgraph(graphlegend)
        
        # Graph output
        if filename:
            if filename.endswith('.svg'):
                dotgraph.write_svg(filename, prog = dotprog)
            elif filename.endswith('.png'):
                dotgraph.write_png(filename, prog = dotprog)
            else:
                raise TypeError('Only svg or png file types can be written.',filename)
        else:
            import io
            import matplotlib.image as mpimg
            from PIL import Image
            png = io.BytesIO(dotgraph.create_png(prog = dotprog))
            img = Image.open(png)
            #img.thumbnail(resize, Image.ANTIALIAS)
            img.show()
            #imgplot = plt.imshow(mpimg.imread(png))
            #imgplot.get_figure().axes[0].axis('off')
            #return imgplot
            
        return dotgraph

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
        return self.plot_geneset_scores(geneset,**kwargs)
        
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

    def permutate_unique_geneset_scores(self,geneset,numberOfPermutations=1000,seed=None,ax=None,multiprocessing=False):
        """
        Randomizes the scores for `numberOfPermutations` and calculates the null distribution of 
        scores for the specific geneset.

        seed can be used to (re)set the random state

        Args:
            geneset (set): Set of genes.
            numberOfPermutations (int): Number of similarly sized random geneset calculations.
            seed (float): Random seed.
            ax (plt.Axes): Ax to plot results on.
            multiprocessing (bool): Use a pool of processes to calculate random scores.
                Number of processes configurable in pyni config.
        """
        genesetHash = hash(str(sorted(geneset)))
        allGenes = set(self.netwink.get_nodes_series())
        if seed: random.seed(seed)
        shelveKey = 'Geneset {} - permutations {} - seed {} - kernel hash {}'.format(
            genesetHash,numberOfPermutations,seed,hash(self)
            #optionally also include a hash of the netwink nodes+edges to be sure it did not change
        )
        try: scores = self.netwink.shelve[shelveKey]
        except KeyError:
            if multiprocessing:
                raise NotImplementedError #TODO refactor code for more optimal use of multiprocessing
            else:
                score, scores = Kernel.mp_diffuse_score_geneset(
                    self.netwink.get_nodes_series().isin(geneset),
                    self.netwink.scoresVector,
                    self.computedMatrix,
                    nperm=numberOfPermutations
                )
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
