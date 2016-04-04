#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
import sys, os
import enrichme


def publish():
    """Publish to PyPi"""
    os.system("python setup.py bdist_wheel sdist upload")
    
if sys.argv[-1] == "publish":
    publish()
    sys.exit()

setup(
    name='enrichme',
    version=enrichme.__version__,
    description='Test enrichment of genome-wide statistics in gene (or other) categories while correcting for gene-length and clustering.',
    long_description=open('README.rst').read() + '\n\n' +
                     open('HISTORY.rst').read(),
    url='https://github.com/feilchenfeldt/enrichme',
    author='Hannes Svardal',
    author_email='hannes.svardal@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
	'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    keywords='Enrichment',
    py_modules=['enrichme','pandas_util'],
    #packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=[
        "pandas"
    ],
    tests_require=[
            'pytest',
    ],
    entry_points={
        'console_scripts': [
            'enrichme=enrichme:enrichme'
        ],
    },
)
