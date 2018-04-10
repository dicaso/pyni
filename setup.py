from setuptools import setup

package = "pyni"
version = "0.0.7"

setup(name = package,
      version = version,
      description="Python network inference package",
      url='https://github.com/dicaso/pyni',
      author = 'Christophe Van Neste',
      author_email = 'christophe.vanneste@ugent.be',
      license = 'GNU GENERAL PUBLIC LICENSE',
      packages = ['pyni'],
      python_requires='>3.6',
      install_requires = [
          'matplotlib',
          'pillow',
          'pandas',
          'scipy',
          'requests',
          'networkx',
          'pydot',
          'bidali'
      ],
      extras_require = {
          'development': ['twine','Sphinx']
          #twine for uploading to pypi, e.g.: twine upload --repository pypi dist/pyni-0.0.5.tar.gz
      },
      package_data = {
          'pyni': [
              'data/reactome_FI.txt',
              'data/reactome_FI_filteredEdges.txt',
              'data/cosmic_20180125.tsv'
          ]
      },
      include_package_data = True,
      zip_safe = False,
      entry_points = {
          'console_scripts': ['ni=pyni.__main__:main'],
      },
      test_suite = 'nose.collector',
      tests_require = ['nose']
)

#To install with symlink, so that changes are immediately available:
#pip install -e .
