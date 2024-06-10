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


__Get bins of clustered transcripts__

This will let you analyse with single cell tools. TODO: add this as default output.

```python
import pandas as pd
from spers.ficture.workflow.scripts.hex_bin import transcript_to_hex_bins

transcripts_df = pd.read_csv(
    "spers.out/ficture/results/transcripts.tsv.gz",
    sep="\t",
    compression="gzip",
    index_col="transcript_id")
scored_df = pd.read_csv(
    "spers.out/ficture/results/scored_transcripts/transcripts.rescored.tsv.gz",
    sep="\t",
    compression="gzip",
    index_col="transcript_id")

transcripts_df = pd.concat((transcripts_df, scored_df), axis=1, join="inner")
transcripts_df = transcript_to_hex_bins(transcripts_df, hex_width=16) # I used 24 but 16 um is probably better IDK
transcripts_df = transcripts_df.groupby(["hex_id", "topK", "gene"]).size().reset_index()
transcripts_df["bin"] = transcripts_df.groupby(["hex_id", "topK"]).ngroup()
transcripts_df = transcripts_df.loc[:, ["bin", "gene", 0]]
transcripts_df = pd.pivot_table(transcripts_df, index="bin", values=0, columns="gene", fill_value=0)
transcripts_df.to_csv("spers.out/ficture/results/transcripts.scored.binned.csv")
```

You will probably want to filter low count bins when analysing.

```r
spers = read.csv(
  "spers.out/ficture/results/transcripts.scored.binned.csv", 
  row.names = 1)
spers = t(as.matrix(spers))
spers = CreateSeuratObject(
  counts=spers,
  min.cells = 3,
  min.features = 16
)
# ...
```

Have a look at our preliminary comparison [HERE](https://bioinf.cc/docs/2024_WEHI.png)
