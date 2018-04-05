
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
        self.computedMatrix = self.restartmatrix #TODO possibly normalise so that row sums are 1
        return self
