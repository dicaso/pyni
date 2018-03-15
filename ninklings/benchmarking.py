# -*- coding: utf-8 -*-
"""Benchmarking.

Module for comparing ninklings netwerk enrichment results
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

    def run_ni(self, kernel = 'random_walk'):
        """
        Run ninklings algorithm
        """
        self.niresults = self.netwink.gsea(
            self.genesets, self.sims.data.set_index('gene').score, kernel = kernel, minSize = 10
        )
    
    def run_opponent(self):
        """
        Runs the competing method
        """
        from bidali.retro import fgsea
        self.opresults = fgsea(self.genesets_r, self.sims.data.set_index('gene').score, minSize = 10)

    def compare_methods(self, permutations = 30, kernel = 'random_walk', simul_kwargs = {}):
        """
        Runs run_ni and run_opponent x times (permutations)
        and compares results
        """
        self.niresultsList = []
        self.opresultsList = []
        for i in range(permutations):
            self.sims = SimScores(
                self.netwink.get_nodes_series(),
                self.trueGenes,
                **simul_kwargs
            )
            print('Running ni {} ...'.format(i))
            self.run_ni(kernel = kernel)
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

def main(args=None):
    """
    Runs benchmarking test for ninklings' algorithms
    """
    import os
    from bidali.LSD import get_msigdb6
    mdb = get_msigdb6()
    trueGenesets = {
        'HALLMARK_E2F_TARGETS',
        'HALLMARK_G2M_CHECKPOINT'
    }
    bm = GenesetBenchmark(
        os.path.join(config['ninklings']['projectdir'],'ni_benchmark'),
        genesets = mdb['H'],
        trueGenesets = {k:mdb['H'][k] for k in mdb['H'].keys() & trueGenesets},
        netwink_kwargs = {'cosmicOnly': True}
    )
    bm.compare_methods()

    return bm
