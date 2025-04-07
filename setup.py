#!/usr/bin/env python
from setuptools import setup


__version__ = "2022.03"


classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 3.6
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""

with open('README.md') as f:
    long_description = f.read()

classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='sequence-processing-pipeline',
      version=__version__,
      long_description=long_description,
      license="BSD",
      description='Sequence processing pipeline',
      author="Qiita development team",
      author_email="qiita.help@gmail.com",
      url='https://github.com/biocore/mg-scripts',
      setup_requires=['numpy', 'cython'],
      include_package_data=True,
      packages=['sequence_processing_pipeline'],
      package_data={
          'sequence_processing_pipeline': [
              'templates/*',
              'scripts/*'
              ]
          },
      install_requires=[
        'click', 'requests', 'pandas', 'flake8', 'nose', 'coverage',
        'pgzip', 'jinja2', 'metapool @ https://github.com/biocore/'
        'metagenomics_pooling_notebook/archive/master.zip'
        ],
      entry_points={
          'console_scripts': ['demux=sequence_processing_pipeline.scripts.cli'
                              ':demux'],
      })
