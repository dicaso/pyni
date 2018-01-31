# ninklings
Network ink enrichment stats

## Get started

   git clone https://github.ugent.be/cvneste/ninklings.git && cd ninklings
   mkvirtualenv ninklings
   pip install -e .     # installs package linked to the git repo
   python setup.py test # runs all the tests

## Working on the project
`workon ninklings`

	import ninklings as ni
	netwink = ni.Netwink('/tmp/testnet')
	ni.ExponentialDiffusionKernel(netwink).compute()
   