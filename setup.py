from distutils.core import setup

import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'autorank', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='autorank',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['autorank'],
    scripts=['bin/autorank'],
    package_data={'autorank' : ['VERSION']},
    url='http://pypi.python.org/pypi/autorank/',
    license='GPL3',
    description='Assigns taxonomic ranks based on evolutionary divergence.',
    long_description=open('README.md').read(),
    install_requires=[
        "scikit-bio >= 0.4.0"],
)
