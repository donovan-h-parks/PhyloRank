import os
import re
from setuptools import setup, find_packages


def read_meta():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'phylorank/__init__.py')
    with open(path) as fh:
        hits = re.findall(r'__(\w+)__ ?= ?["\'](.+)["\']\n', fh.read())
    return {k: v for k, v in hits}


def readme():
    with open('README.md') as fh:
        return fh.read()


meta = read_meta()

setup(
    author=meta['author'],
    author_email=meta['author_email'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    data_files=[("", ["LICENSE"])],
    description=meta['description'],
    entry_points={
        'console_scripts': [
            'phylorank = phylorank.__main__:main'
        ]
    },
    install_requires=['biolib >= 0.1.0', 'numpy', 'matplotlib',
                      'dendropy>=4.1.0', 'scipy', 'mpld3>=0.5.2'],
    keywords='phylogenetics taxonomy relative evolutionary divergence tree '
             'decorate decoration',
    license=meta['license'],
    long_description=readme(),
    long_description_content_type='text/markdown',
    name=meta['title'],
    packages=find_packages(),
    project_urls={
        "Bug Tracker": "https://github.com/dparks1134/PhyloRank/issues",
        "Documentation": "https://github.com/dparks1134/PhyloRank",
        "Source Code": "https://github.com/dparks1134/PhyloRank",
    },
    python_requires='>=3.6',
    url=meta['url'],
    version=meta['version']
)
