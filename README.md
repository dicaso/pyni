# ninklings -Network ink enrichment stats
<img title="ninklings logo" src="ni_logo.svg" width="300">

## Requirements

 - https://virtualenvwrapper.readthedocs.io/en/latest/

## Get started

    git clone https://github.ugent.be/cvneste/ninklings.git && cd ninklings
    git checkout -b $USER
    mkvirtualenv ninklings
    pip install -e .     # installs package linked to the git repo
    python setup.py test # runs all the tests

### Configuration file

`ninklings.cfg` in your working directory, `~/ninklings.cfg`, or `/etc/ninklings.cfg` (in this order of precedence)

    [ninklings]
    datadir = ...path-to-repository.../ninklings/data

## Working on the project
`workon ninklings`

### Setting up network

    import ninklings as ni, random
    import matplotlib.pyplot as plt
    netwink = ni.Netwink('/tmp/testnet',cosmicOnly=True)
    netwink.plot_admatrix()
    netkernel = ni.RestartRandomWalk(netwink).compute()    #ni.ExponentialDiffusionKernel(netwink).compute()

### Working with genesets

    from bidali.LSD import get_msigdb6
    mdb = get_msigdb6()
    gs_p53 = {g for g in mdb['H']['HALLMARK_P53_PATHWAY'] if g in set(netwink.get_nodes_series())} #only 26 out of 200 in cosmic
    genescores = {g:random.randint(10,100) for g in random.sample(set(netwink.get_nodes_series()),len(gs_p53))} #TODO good random model
    randomgeneset = random.sample(set(netwink.get_nodes_series()),len(gs_p53))
    netwink.apply_gene_scores(genescores)
    fig, ax = plt.subplots()
    genesetscore_random = netkernel.score_geneset(randomgeneset,ax=ax)
    genesetscore_p53 = netkernel.score_geneset(gs_p53,ax=ax)#,c='g')
    nulldistro = netkernel.permutate_geneset_scores(randomgeneset,ax=ax)

### Loading weights

    from bidali.LSD.dealer.external.cohorts import get_FischerData
    fd = get_FischerData()
    netwink.load_correlations(fd.exprdata)

### Hierarchical clustering

    from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
    Z=linkage(netkernel.convert_to_distance(),'ward')
    fig = plt.figure(figsize=(25, 10))
    dn = dendrogram(Z)
    root,nodeslist = to_tree(Z,rd=True)
    netwink.get_nodes_series()[netwink.get_nodes_series().index.isin(root.left.pre_order())]

## Working on the documentation

    pip install ninklings[documentation]
    #sphinx-quickstart #was used for initiating the documentation
    cd docs && make html

Before publication, follow http://dont-be-afraid-to-commit.readthedocs.io/en/latest/documentation.html
to make it public


# Meeting 19-02-2018

disconneted components, full lines or columns that are 0
check up front for those rows and cols
remove them, for calculations, then reinserting

=> in the base class

for inversion matrices very important

* normalization
cosine normalization and eucledian transformation
=> hard to explain why normalization would be needed

eucledian distance
meta distance

hotnet paper

ward only when eucledian distance


weight matrix for overlay of adjacency matrix
whenever you have weights you should use them

normal expression in that tissue, or independent
  gene expression network
filter the network for our condition

in undirectional network use absolute values of the correlation

sanity check: check clustering in
visualized network (cytoscape e.g.)

with random walk weights are in any case normolized
with weight they become proportional

PLANNING
--------
- weight matrices (mogelijks verschillende manieren om dat te doen) (dissapear node -> sink node)

- ranks toepassen -> monotoon, continue
 - subranks afzonderlijk toepassen en dan integreren

- opmerkelijke stijgers en dalers in je rankings

- analyze laten shiften tussen volledig prior en volledig colleratie netwerk

CAMDA challenge
---------------
Maarten gaat erop werken