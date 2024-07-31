# sync_osc_MSTMs
Code used for the article "Synchrony, oscillations, and phase relationships in collective neuronal activity: a highly-comparative overview of methods" by F. Baroni and B.D. Fulcher (2024).
The code is licensed under the [GNU GPL v3 license](http://www.gnu.org/licenses/gpl-3.0.html) (or later).

Most code is in Matlab. Some functions use `tightPosition` for creating figures and hence require a modern version of Matlab (R2022b or newer), but most of the code will run with an older version.
Some of the scripts and functions require the Wavelet Toolbox, the Image Processing Toolbox, or the Statistics and Machine Learning Toolbox, but most of the code will run with a basic Matlab installation.
Some of the scripts and functions require a Python environment recognized by Matlab. In particular, these are the functions that perform [FOOOF](https://github.com/fooof-tools/fooof) model estimation, and those that estimate intrinsic dimensionality using [DADApy](https://github.com/sissa-data-science/DADApy).

This repository includes code for:
- the generation of synthetic spike trains, inside the directory `synth_train_generation`;
- the analysis of real and synthetic spike trains, including the generation of all paper figures included in Baroni & Fulcher 2024, inside the directory `data_analysis`.
The README file includes a description of all high-level scripts and functions. 

# Synthetic Spike Train Generation

# Synthetic Spike Train Analysis
## Applying the library of Multineuron Spike Train Measures (MSTMs)

## Correlation of each MSTM with each generative parameter

## Organizing the MSTM library using inter-MSTM correlations and clustering

# Real Spike Train Analysis
## Applying the library of Multineuron Spike Train Measures (MSTMs)
## Organizing data segments using clustering








