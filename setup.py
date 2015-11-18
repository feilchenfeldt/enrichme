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
    name='Enrichme',
    version=enrichme.__version__,
    description='An enrichment library',
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
    py_modules=['enrichment','pandas_util'],
    #packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=[
        "numpy >=1.6.1",
        "scipy >=0.13.0",
        "matplotlib >= 1.4.3",
        "pandas >= 0.17.0"
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
