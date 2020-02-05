import os
from distutils.core import setup


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    version_file = open(os.path.join(setup_dir, 'phylorank', 'VERSION'))
    return version_file.readline().strip()


setup(
    name='phylorank',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['phylorank', 'phylorank.plot'],
    scripts=['bin/phylorank'],
    package_data={'phylorank': ['VERSION']},
    url='http://pypi.python.org/pypi/phylorank/',
    license='GPL3',
    description='Assigns taxonomic ranks based on evolutionary divergence.',
    install_requires=["biolib >= 0.0.46", 'dendropy', 'numpy', 'mpld3', 'scipy'],
)
