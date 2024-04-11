# SPERS

SPatialomics Enhanced Research Suite (we can change the name later)

## Install

```shell
# create a conda env for the install
conda create -n spers python=3.12
conda activate spers

# clone this repo
git clone https://github.com/beardymcjohnface/SPERS.git
cd SPERS/

# install
pip install -e .

# check
spers --help
spers ficture --help

# to update - just use git
git pull
```

## Ficture wrapper

Wrapper for the ficture package, including file conversion.

```shell
# for xenium
spers ficture --input transcripts.csv.gz --platform xenium

# for cosmx
spers ficture --input transcripts.csv.gz --platform cosmx
```

__Current steps covered:__
- csv_to_tsv
- make_spatial_minibatch
