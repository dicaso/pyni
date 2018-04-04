from setuptools import setup

package = "pyni"
version = "0.0.3"

setup(name = package,
      version = version,
      description="Python network inference package",
      url='https://github.com/beukueb/pyni',
      author = 'Christophe Van Neste',
      author_email = 'christophe.vanneste@ugent.be',
      license = 'GNU GENERAL PUBLIC LICENSE',
      packages = ['pyni'],
      install_requires = [
          'matplotlib',
          'pandas',
          'scipy',
          'requests',
          'networkx',
          'pydot',
          'bidali'
      ],
      extras_require = {
          'documentation': ['Sphinx']
      },
      package_data = {
          'pyni': [
          ]
      },
      include_package_data = False,
      zip_safe = False,
      entry_points = {
          'console_scripts': ['ni=pyni.__main__:main'],
      },
      test_suite = 'nose.collector',
      tests_require = ['nose']
)

#To install with symlink, so that changes are immediately available:
#pip install -e .
