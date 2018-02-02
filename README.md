# ninklings
Network ink enrichment stats

## Get started

    git clone https://github.ugent.be/cvneste/ninklings.git && cd ninklings
    git checkout -b $USER
    mkvirtualenv ninklings
    pip install -e .     # installs package linked to the git repo
    python setup.py test # runs all the tests

## Working on the project
`workon ninklings`

    import ninklings as ni, random
    import matplotlib.pyplot as plt
    netwink = ni.Netwink('/tmp/testnet')
    #netkernel = ni.ExponentialDiffusionKernel(netwink).compute()
    netkernel = ni.RestartRandomWalk(netwink).compute()
    genesetscores = {g:random.randint(10,100) for g in random.sample(set(netwink.get_nodes_series()),50)}
    randomgeneset = random.sample(set(netwink.get_nodes_series()),20)
    netwink.apply_gene_scores(genesetscores)
    fig, ax = plt.subplots()
    genesetscore = netkernel.score_geneset(randomgeneset,ax=ax)
    nulldistro = netkernel.permutate_geneset_scores(randomgeneset,ax=ax)

## Working on the documentation

   pip install ninklings[documentation]
   #sphinx-quickstart #was used for initiating the documentation
   