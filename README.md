# moshpit-notebooks
Analysis notebooks accompanying the MOSHPIT manuscript.

## Contents
1. [TARA Oceans Expedition analysis](tara/tara.ipynb)
2. [Cocoa fermentation time course analysis](cocoa/cocoa.ipynb)
3. [Mock community analysis](mag-mock/MAGMock.ipynb)

## Setup
To run the notebooks included in this repository you will need a working QIIME 2 environment created 
using the **metagenome** distribution. Please follow the steps below to create and update one:

1. Create a fresh 2024.10 environment as described in the [QIIME 2 documentation](https://docs.qiime2.org/2024.10/install/native/#qiime-2-metagenome-distribution)
2. Activate the environment:
   ```bash
   conda activate qiime2-metagenome-2024.10
   ```
3. Install the other required dependencies:
   ```bash
   conda install -c bioconda -c conda-forge mason plotly
   ```
   ```bash
   pip install kaleido
   ```
4. You should be good to go!
