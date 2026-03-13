# moshpit-notebooks
Analysis notebooks accompanying the MOSHPIT manuscript.

## Contents
1. [TARA Oceans Expedition analysis](tara/tara.ipynb)
2. [Cocoa fermentation time course analysis](cocoa/cocoa.ipynb)
3. [Mock community analysis](mag-mock/MAGMock.ipynb)
4. [CAMI II Toy Human Microbiome Project dataset analysis](cami2-thmp/thmp.ipynb)

## Setup
To run the notebooks included in this repository you will need a working QIIME 2 environment created 
using the **moshpit** distribution. Please follow the steps below to create and update one:

1. Create a fresh 2026.1 environment as described in the [QIIME 2 library](https://library.qiime2.org/quickstart/moshpit)
2. Activate the environment:
   ```bash
   conda activate qiime2-moshpit-2026.1
   ```
3. Install the other required dependencies:
   ```bash
   conda install -c bioconda -c conda-forge plotly
   ```
   ```bash
   pip install kaleido
   ```
4. You should be good to go!
