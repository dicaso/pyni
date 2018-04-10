# -*- coding: utf-8 -*-
"""Benchmarking.

Module for comparing pyni netwerk enrichment results
with other genset enrichment methods
"""
from . import Netwink
from .simulations import SimScores
from .config import config

class GenesetBenchmark:
    """Geneset Benchmark.

    Compares results of two methods.
    """
    def __init__(self,savedir,genesets,trueGenesets,netwink_kwargs={},simul_kwargs={}):
        """
        genesets and trueGenesets should be dict like,
        with geneset names as keys and geneset lists as values.
        """
        # Set base benchmark attributes
        self.genesets = genesets
        self.truth = trueGenesets
        self.trueGenes = {g for t in self.truth for g in self.truth[t]}
        # Load netwink
        netwink_defaults = {}
        netwink_defaults.update(netwink_kwargs)
        self.netwink = Netwink(savedir,**netwink_defaults)
        # Simulate scores
        self.sims = SimScores(
            self.netwink.get_nodes_series(),
            self.trueGenes,
            **simul_kwargs
        )
        # Prepare R data
        from bidali.retro import genesets2indices_r
        self.geneLabels = self.netwink.get_nodes_series().as_matrix()
        self.genesets_r = genesets2indices_r(self.genesets,self.geneLabels)

    def run_ni(self, kernel = 'random_walk', kernel_compute_kwargs = {}):
        """
        Run pyni algorithm
        """
        self.niresults = self.netwink.gsea(
            self.genesets, self.sims.data.set_index('gene').score, kernel = kernel, minSize = 10,
            kernel_compute_kwargs = kernel_compute_kwargs
        )
    
    def run_opponent(self):
        """
        Runs the competing method
        """
        from bidali.retro import fgsea
        self.opresults = fgsea(self.genesets_r, self.sims.data.set_index('gene').score, minSize = 10)

    def compare_methods(self, permutations = 30, kernel = 'random_walk', simul_kwargs = {}, kernel_compute_kwargs = {}):
        """
        Runs run_ni and run_opponent x times (permutations)
        and compares results
        """
        import pandas as pd
        self.niresultsList = []
        self.opresultsList = []
        for i in range(permutations):
            self.sims = SimScores(
                self.netwink.get_nodes_series(),
                self.trueGenes,
                **simul_kwargs
            )
            print('Running ni {} ...'.format(i))
            self.run_ni(kernel = kernel, kernel_compute_kwargs = kernel_compute_kwargs)
            print('Running op {} ...'.format(i))
            self.run_opponent()
            self.niresultsList.append(self.niresults)
            self.opresultsList.append(self.opresults)
        # summarize comparison
        rankpositions = {}
        probabilities = {}
        for true_geneset in self.truth:
            rankpositions[(true_geneset,'ni')] = []
            rankpositions[(true_geneset,'op')] = []
            probabilities[(true_geneset,'ni')] = []
            probabilities[(true_geneset,'op')] = []
            for nir,opr in zip(self.niresultsList, self.opresultsList):
                rankpositions[(true_geneset,'ni')].append(nir.index.get_loc(true_geneset))
                rankpositions[(true_geneset,'op')].append(opr.set_index('pathway').index.get_loc(true_geneset))
                probabilities[(true_geneset,'ni')].append(nir.loc[true_geneset].prob)
                probabilities[(true_geneset,'op')].append(opr.set_index('pathway').loc[true_geneset].pval)
        self.rankcomparison = pd.DataFrame(rankpositions)
        self.probcomparison = pd.DataFrame(probabilities)

def main(netwink_kwargs={}):
    """
    Runs benchmarking test for pyni' algorithms
    over a set of random 'true' genesets
    """
    import os, random
    from bidali.LSD import get_msigdb6
    import numpy as np, pandas as pd
    import matplotlib.pyplot as plt
    netwink_kwargs.setdefault('cosmicOnly', True)
    mdb = get_msigdb6()

    # Select genesets that are big enough
    nw = Netwink('/tmp/nw_checksize',**netwink_kwargs)
    print('Original number of genesets',len(mdb['H']))
    mdb['H'] = {
        gs:mdb['H'][gs] for gs in mdb['H']
                if len(set(mdb['H'][gs])&set(nw.get_nodes_series())) >= 15
    }
    print('Filtered number of genesets',len(mdb['H']))
    del nw

    # Select true geneset combinations
    trueGenesetCombinations = {('HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT')}
    trueGenesetCombinations.update({tuple(random.sample(mdb['H'].keys(),2)) for i in range(5)})
    trueGenesetResults = {}
    for ti,trueGenesetsTuple in enumerate(trueGenesetCombinations):
        print('Processing simulated true genesets',trueGenesetsTuple)
        trueGenesets = set(trueGenesetsTuple)
        bm = GenesetBenchmark(
            os.path.join(config['pyni']['projectdir'],'ni_benchmark_{}'.format(ti)),
            genesets = mdb['H'],
            trueGenesets = {k:mdb['H'][k] for k in mdb['H'].keys() & trueGenesets},
            netwink_kwargs = netwink_kwargs
        )
        # Set restart probabilities
        restartProbs = (np.logspace(0,1,5)/-10)+1.1
        kernel_kwargs_list = [
          {'restartProb':i} for i in restartProbs
        ]
        probcomparisons = []
        rankcomparisons = []
        for rrwkernel in kernel_kwargs_list:
            print('Restart probability:',rrwkernel['restartProb'])
            bm.compare_methods(kernel_compute_kwargs = rrwkernel, permutations = 5)
            probcomparisons.append(bm.probcomparison)
            rankcomparisons.append(bm.rankcomparison)
        # Plot summaries
        means_probcomparisons = pd.DataFrame(pd.concat([df.mean() for df in probcomparisons]),columns=['probMean'])
        means_probcomparisons['restartProb'] = np.repeat(restartProbs, 4)
        means_probcomparisons['probVar'] = pd.concat([df.var() for df in probcomparisons])
        fig,ax1=plt.subplots()
        for gs in trueGenesets:
            for method in ('ni','op'):
                ax1.errorbar(
                    means_probcomparisons.loc[gs,method].restartProb,
                    means_probcomparisons.loc[gs,method].probMean,
                    yerr = means_probcomparisons.loc[gs,method].probVar ** .5,
                    fmt='o',
                    label = '{} - {}'.format(gs,method)
                )
        ax1.legend()
        ax1.invert_xaxis()
        ax1.set_title('Probability comparison - {}'.format(trueGenesetsTuple))
        ax1.set_xlabel('Random restart probability')
        means_rankcomparisons = pd.DataFrame(pd.concat([df.mean() for df in rankcomparisons]),columns=['rankMean'])
        means_rankcomparisons['restartProb'] = np.repeat(restartProbs, 4)
        means_rankcomparisons['rankVar'] = pd.concat([df.var() for df in rankcomparisons])
        fig,ax2=plt.subplots()
        for gs in trueGenesets:
            for method in ('ni','op'):
                ax2.errorbar(
                    means_rankcomparisons.loc[gs,method].restartProb,
                    means_rankcomparisons.loc[gs,method].rankMean,
                    yerr = means_rankcomparisons.loc[gs,method].rankVar ** .5,
                    fmt='o',
                    label = '{} - {}'.format(gs,method)
                )
        ax2.legend()
        ax2.invert_xaxis()
        ax2.set_title('Rank comparison - {}'.format(trueGenesetsTuple))
        ax2.set_xlabel('Random restart probability')
        bm.rankcomparisons = means_rankcomparisons
        bm.probcomparisons = means_probcomparisons
        trueGenesetResults[trueGenesetsTuple] = (bm, ax1, ax2)

    return trueGenesetResults

    
