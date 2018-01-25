from setuptools import setup

package = "ninklings"
version = "0.1.1"

setup(name = package,
      version = version,
      description="Network ink enrichment stats",
      url='https://github.ugent.be/cvneste/ninklings',
      author = 'Christophe Van Neste',
      author_email = 'christophe.vanneste@ugent.be',
      license = 'GNU GENERAL PUBLIC LICENSE',
      packages = ['ninklings'],
      install_requires = [
          'matplotlib',
          'pandas',
          'requests',
          'networkx'
      ],
      package_data = {
          'ninklings': [
          ]
      },
      include_package_data = False,
      zip_safe = False,
      entry_points = {
          'console_scripts': ['ni=ninklings.__main__:main'],
      },
      test_suite = 'nose.collector',
      tests_require = ['nose']
)

#To install with symlink, so that changes are immediately available:
#pip install -e .
