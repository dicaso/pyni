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

    def run_ni(self):
        """
        Run ninklings algorithm
        """
        pass
    
    def run_opponent(self):
        """
        Runs the competing method
        """
        from bidali.retro import fgsea
        self.opresults = fgsea(self.genesets_r, self.sims.data.set_index('gene').score)

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
    bm.run_opponent()

    return bm
