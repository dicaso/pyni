# -*- coding: utf-8 -*-
"""Simulations.

Random simulations for testing.
"""

import random, numpy as np, pandas as pd
import scipy.special as sps

class SimScores:
    """
    Make a similuted scores object.
    Produces the scores to use, and
    methods for calculating false positives/negatives
    """
    def __init__(
            self, allGenes, trueGenes, distro='gamma', onlyPositiveScores = True,
            gamma_shape_t = 2, gamma_scale_t = 2, gamma_shape_f = 1, gamma_scale_f = .5
    ):
        """
        Sets up simulated scores
        distro choices are: gamma, uniform
        """
        self.genes = pd.Series(allGenes)
        self.truth = self.genes.isin(trueGenes)
        self.shape = (sum(self.truth),len(self.genes)-sum(self.truth))
        self.data = pd.DataFrame({'gene': self.genes, 'truth': self.truth})
        self.distro = distro
        if self.distro == 'gamma':
            self.distro_params = {
                't': (gamma_shape_t, gamma_scale_t),
                'f': (gamma_shape_f, gamma_scale_f)
            }
            self.data.loc[self.truth,'score'] = np.random.gamma(gamma_shape_t, gamma_scale_t, self.shape[0])
            self.data.loc[~self.truth,'score'] = np.random.gamma(gamma_shape_f, gamma_scale_f, self.shape[1])
        elif self.distro == 'uniform':
            raise NotImplementedError
        else:
            raise Exception('%s distro not supported' % distro)
        if not onlyPositiveScores: raise NotImplementedError

    def plot_scores(self):
        """Plot scores
        
        Visually screen the similated scores.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        if self.distro == 'gamma':
            for torf in ('f','t'):
                gamma_shape, gamma_scale = self.distro_params[torf]
                color = 'g' if torf == 't' else 'r'
                count, bins, ignored = ax.hist(
                    self.data.score[self.truth if torf=='t' else ~self.truth], 'auto',
                    normed=True, color=color, alpha=.5 if torf=='t' else .8
                )
                y = bins**(gamma_shape-1)*(np.exp(-bins/gamma_scale) /
                                           (sps.gamma(gamma_shape)*gamma_scale**gamma_shape))
                ax.plot(bins, y, linewidth=2, color=color)
            
        else:
            raise NotImplementedError

        return ax

    def empirical_param_calc(self,scores,trueGenes=None,distro='gamma'):
        """Calculate distro param values

        Parameter value calculation based on real data.
        If no trueGenes are given, only calculates the average distribution
        """
        raise NotImplementedError
