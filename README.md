# Sequence Processing Pipeline

A Jupyter notebook to assist wet lab shotgun pipeline.
A packaged Python-based implementation of Knight Lab's sequencing process pipeline.

## Installation

To install this package, first clone this repository from GitHub:

```bash
git clone https://github.com/biocore/mg-scripts.git
```

Create a Python3 Conda environment in which to run the notebook:

```bash
conda create --yes -n spp python='python=3.9' scikit-learn pandas numpy nose pep8 flake8 matplotlib jupyter notebook 'seaborn>=0.7.1' pip openpyxl 'seqtk>=1.4' click scipy fastq-pair
```

Activate the Conda environment:

```bash
source activate sp_pipeline
```

Change directory to the cloned repository folder and install:

```bash
cd mg-scripts
pip install -e .
```

This will automatically install https://github.com/biocore/metagenomics_pooling_notebook.git, a dependency of mg-scripts and the sequence_processing_pipeline.

## Running Unittests

Change directory to the downloaded repository folder:

```bash
cd mg-scripts
nosetests --with-coverage --cover-inclusive --cover-package sequence_processing_pipeline
```

## Getting Started

Review Pipeline.py and main.py to learn how to import and access package functionality:

```bash
cd mg-scripts/sequence_processing_pipeline
more Pipeline.py
more main.py
```

Adjust configuration settings as needed:

```bash
cd mg-scripts/sequence_processing_pipeline
vi configuration.json
```

Please note that the setting 'minimap2_databases' is expected to be a list of paths to individual .mmi files for QCJob.
For NuQCJob, minimap2_databases is expected to be the path to a directory containing two subdirectories: 'metagenomic'
and 'metatranscriptomic'. Each directory should contain or symlink to the appropriate .mmi files needed for that Assay
type.

Additional TellSeq-related notes:
'spades-cloudspades-0.1', 'tellread-release-novaseqX' or similar directories must be placed in a location available to SPP.
Their paths should be made known to SPP in the configuration files. (See examples for details).
Additional scripts found in sequence_processing_pipeline/contrib were contributed by Daniel and Omar and can be similarly located and configured.

