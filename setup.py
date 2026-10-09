import os
from setuptools import setup, find_packages
import sys

version_py = os.path.join(os.path.dirname(__file__), 'Pyntegrate', 'version.py')
version = open(version_py).read().split('=')[-1].strip().replace('"','')

requirements = open(os.path.join(os.path.dirname(__file__), 'requirements.txt')).readlines()
long_description = open('README.md').read()

setup(
        name='Pyntegrate',
        version=version,
        description="Integrative analysis of high-thoughput sequencing data",
        long_description=long_description,
        long_description_content_type='text/markdown',
        license="MIT",
        install_requires=requirements,
        packages=find_packages(),
        python_requires='>=3.9,<3.13',
        classifiers=[
            'Intended Audience :: Science/Research',
            'Intended Audience :: Developers',
            'Intended Audience :: System Administrators',
            'Operating System :: POSIX',
            'Operating System :: MacOS :: MacOS X',
            'Environment :: Console',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Medical Science Apps.',
        ]
    )
