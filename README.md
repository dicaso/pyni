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
    netwink = ni.Netwink('/tmp/testnet',cosmicOnly=True)
    #netkernel = ni.ExponentialDiffusionKernel(netwink).compute()
    netkernel = ni.RestartRandomWalk(netwink).compute()
    genesetscores = {g:random.randint(10,100) for g in random.sample(set(netwink.get_nodes_series()),50)}
    randomgeneset = random.sample(set(netwink.get_nodes_series()),20)
    netwink.apply_gene_scores(genesetscores)
    fig, ax = plt.subplots()
    genesetscore = netkernel.score_geneset(randomgeneset,ax=ax)
    nulldistro = netkernel.permutate_geneset_scores(randomgeneset,ax=ax)
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