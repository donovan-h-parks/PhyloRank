from distutils.core import setup

import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'phylorank', 'VERSION'))
    return versionFile.readline().strip()

data_files = []
js_files = []
js_dir = './phylorank/mpld3/mpld3/js'
for f in os.listdir(js_dir):
    js_files.append(os.path.join(js_dir, f))
data_files.append((js_dir, js_files))

setup(
    name='phylorank',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['phylorank', 
                'phylorank.plot', 
                'phylorank.mpld3.mpld3', 
                'phylorank.mpld3.mpld3',
                'phylorank.mpld3.mpld3.mplexporter',
                'phylorank.mpld3.mpld3.mplexporter.renderers'],
    scripts=['bin/phylorank'],
    package_data={'phylorank' : ['VERSION']},
    data_files=data_files,
    url='http://pypi.python.org/pypi/phylorank/',
    license='GPL3',
    description='Assigns taxonomic ranks based on evolutionary divergence.',
    install_requires=[
        "biolib >= 0.1.0"],
)
