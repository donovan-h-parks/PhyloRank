from distutils.core import setup

import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'phylorank', 'VERSION'))
    return versionFile.readline().strip()

setup(
    name='phylorank',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['phylorank', 'phylorank.plot'],
    scripts=['bin/phylorank'],
    package_data={'phylorank' : ['VERSION']},
    url='http://pypi.python.org/pypi/phylorank/',
    license='GPL3',
    description='Assigns taxonomic ranks based on evolutionary divergence.',
    install_requires=[
        "biolib >= 0.0.46"],
)
