#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools


setuptools.setup(
    name='urban_metagenomics_tutorial',
    version='0.1.0',
    description="Metagenomics utilities",
    author="David C. Danko",
    author_email='dcdanko@gmail.com',
    url='https://github.com/dcdanko/urban_metagenomics_tutorial',
    packages=setuptools.find_packages(),
    package_dir={
        'urban_metagenomics_tutorial': 'urban_metagenomics_tutorial',
        'squid': 'squid',
    },
    install_requires=[
        'click',
        'biopython',
        'capalyzer',
        'gimmebio.kmers',
    ],
    entry_points={
        'console_scripts': [
            'umt=urban_metagenomics_tutorial.cli:main',
            'squid=squid.cli:main',
        ]
    },
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
)
