# -*- coding: utf-8 -*-
"""Simulations.

Random simulations for testing.
"""

import random, numpy as np, pandas as pd
import scipy.special as sps
from scipy import stats

class SimScores:
    """
    TODO refactor -> SimScores should be interface class
    Make a similuted scores object.
    Produces the scores to use, and
    methods for calculating false positives/negatives

    >>> SimScores.empirical_param_calc(testdata.logFC,testdata['adj.P.Val']<=.05)
    """
    def __init__(
            self, allGenes, trueGenes, distro='gamma', onlyPositiveScores = True,
            gamma_shape_t = 1.4, gamma_scale_t = .5, gamma_loc_t = .15, gamma_shape_f = 1, gamma_scale_f = .5,
            unif_low_t = 0, unif_high_t = .05, unif_low_f = 0, unif_high_f = 1
    ):
        """
        Sets up simulated scores
        distro choices are: gamma, uniform

        TODO ?include gaussian for t-stat score?
        """
        self.genes = pd.Series(allGenes)
        self.truth = self.genes.isin(trueGenes)
        self.shape = (sum(self.truth),len(self.genes)-sum(self.truth))
        self.data = pd.DataFrame({'gene': self.genes, 'truth': self.truth})
        self.distro = distro
        if self.distro == 'gamma':
            self.distro_params = {
                't': (gamma_shape_t, gamma_scale_t, gamma_loc_t),
                'f': (gamma_shape_f, gamma_scale_f, 0)
            }
            self.data.loc[self.truth,'score'] = np.random.gamma(gamma_shape_t, gamma_scale_t, self.shape[0]) + gamma_loc_t
            self.data.loc[~self.truth,'score'] = np.random.gamma(gamma_shape_f, gamma_scale_f, self.shape[1])
        elif self.distro == 'uniform':
            self.distro_params = {
                't': (unif_low_t, unif_high_t),
                'f': (unif_low_f, unif_high_f)
            }
            self.data.loc[self.truth,'score'] = np.random.uniform(self.distro_params['t'][0], self.distro_params['t'][1], self.shape[0])
            self.data.loc[~self.truth,'score'] = np.random.uniform(self.distro_params['f'][0], self.distro_params['f'][1], self.shape[1])
        else:
            raise Exception('%s distro not supported' % distro)
        if not onlyPositiveScores: raise NotImplementedError

    def plot_scores(self,density=False):
        """Plot scores
        
        Visually screen the similated scores.
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax_y2 = ax.twinx()
        if self.distro == 'gamma':
            for torf in ('f','t'):
                gamma_shape, gamma_scale, gamma_loc = self.distro_params[torf]
                color = 'g' if torf == 't' else 'r'
                count, bins, ignored = ax.hist(
                    self.data.score[self.truth if torf=='t' else ~self.truth], 'auto',
                    density=density, color=color, alpha=.5 if torf=='t' else .8
                )
                bins = bins - gamma_loc
                y = bins**(gamma_shape-1)*(np.exp(-bins/gamma_scale) /
                                           (sps.gamma(gamma_shape)*gamma_scale**gamma_shape))
                ax_y2.plot(bins, y, linewidth=2, color=color)
            
        elif self.distro == 'uniform':
            for torf in ('f','t'):
                unif_low, unif_high = self.distro_params[torf]
                color = 'g' if torf == 't' else 'r'
                count, bins, ignored = ax.hist(
                    self.data.score[self.truth if torf=='t' else ~self.truth], 'auto',
                    density=density, color=color, alpha=.5 if torf=='t' else .8
                )
                ax_y2.plot(bins, np.ones_like(bins), linewidth=2, color=color)

        return ax

    @staticmethod
    def empirical_param_calc(scores,trueGenes=None,distro='gamma',onlyPositiveScores=True):
        """Calculate distro param values

        Parameter value calculation based on real data.
        If no trueGenes are given, only calculates the average distribution.
        trueGenes should be bool pd.Series with same index as a pd.Series scores
        """
        if distro == 'gamma':
            if onlyPositiveScores: scores = scores.abs()
            if trueGenes is not None:
                fit_alpha_t, fit_loc_t, fit_beta_t = stats.gamma.fit(scores[trueGenes])
                fit_alpha_f, fit_loc_f, fit_beta_f = stats.gamma.fit(scores[~trueGenes])
                print('false gamma location was {} and should be close to 0!'.format(fit_loc_f))
                return {
                    'gamma_shape_t': fit_alpha_t,
                    'gamma_scale_t': fit_beta_t,
                    'gamma_loc_t': fit_loc_t,
                    'gamma_shape_f': fit_alpha_f,
                    'gamma_scale_f': fit_beta_f,                    
                }
            else:
                fit_alpha, fit_loc, fit_beta = stats.gamma.fit(scores)
                print('gamma location was {}. Should be close to 0!'.format(fit_loc))
                return {
                    'gamma_shape': fit_alpha,
                    'gamma_scale': fit_beta
                }
