Tutorials
=========

Simulation tutorial
-------------------

Setting up the network:

    import pyni as ni
    from pyni.simulations import SimScores
    from bidali.LSD import get_msigdb6
    netwink = ni.Netwink('/tmp/pyni_cosmic',cosmicOnly=True)
    netwink.plot_admatrix()
    netkernel = ni.RestartRandomWalk(netwink,kernelspec={'restartProb':0.9}).compute()
    mdb = get_msigdb6()
    sims = SimScores(
      netwink.get_nodes_series(),
      set(mdb['H']['HALLMARK_E2F_TARGETS']) &
      set(mdb['H']['HALLMARK_G2M_CHECKPOINT'])
    )
    
Working with genesets:

    #niresults = netwink.gsea(
    #        mdb['H'], sims.data.set_index('gene').score, kernel = 'random_walk', minSize = 10,
    #)
    netwink.apply_gene_scores(sims.data.set_index('gene').score)
    ax = netkernel.plot_scores_diffusion_scatter()
    netkernel.plot_geneset_scores(set(mdb['H']['HALLMARK_E2F_TARGETS']))
    #genesetscore_p53 = netkernel.score_geneset(gs_p53,ax=ax)#,c='g')
    #nulldistro = netkernel.permutate_geneset_scores(randomgeneset,ax=ax)

    
Loading weights:

    from bidali.LSD.dealer.external.cohorts import get_NRC
    nrc = get_NRC()
    netwink.load_correlations(nrc.exprdata)
    netkernel = ni.RestartRandomWalk(netwink,kernelspec={'restartProb':0.9}).compute()
    netkernel.plot_geneset_scores(set(mdb['H']['HALLMARK_E2F_TARGETS']))
