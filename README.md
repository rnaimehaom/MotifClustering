# Motif Clustering

A simple tool to cluster protein motifs using protein sequence based features

## Installation 

We assume here conda is already installed on your machine. If not, install conda (). Then, when you are ready type the following command to install this script :

First, clone the git repository to the desired location on your machine:
```
cd "path/to/the/desired/location"
git clone git@github.com:DjampaKozlowski/MotifClustering.git 
cd MotifClustering/
```

Create a new conda environment named 'moclu_env' and install all the dependecies :
```
conda env create -f moclu_env.yml
```

Activate the envrionment :
```
conda activate moclu_env
```

Install the 'motifclustering' package:
```
pip install -e .
``` 
