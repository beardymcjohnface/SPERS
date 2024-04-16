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

__Run the test dataset__
```shell
spers ficture-test
```

__Current steps covered:__
- Convert CosMX or Xenium CSV to TSV
- Hex binning
- Minibatch binning for model training
- Model generation
- Trained model plot
- Transcript rescoring (with overlapping hex bins)
- Final plot