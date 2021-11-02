# Joint estimation of neonatal mortality and vital registration completeness across Mexico

This repository contains code used to produce the results described in Chapter 2 of the doctoral thesis, "Assessing local health outcomes using spatially-resolved health surveillance data."


## Downloading data

All input data for this project was downloaded directly from the [public data portal](http://en.www.inegi.org.mx/datos/) of the Mexican Institute for Geography and Statistics, INEGI.


## Code structure

This repository is organized into the following folders:
- `data_analysis/`: One-off and exploratory scripts used to assess data availability and quality
- `data_prep/`: Scripts used to clean input data and prepare it for modeling
- `model/`: Functions and scripts used to run the joint birth-history and CRVS model. The main execution script, `run_mex_model.R`, calls portable analysis and spatial modeling functions defined in `model_functions.R`. The core statistical model, writted using the Template Model Builder library, is defined in `vr_bh_tmb_model.cpp`.
- `viz/`: Scripts to analyze the raw data and model results, including all figures included in the chapters.


## Running this code

This code is designed to be executed in R version 3.6 or above. Paths to top-level folders should be set in `config.yaml`; these paths will be loaded in executed scripts and used to locate input data objects. The path to this repository on the user's own computer should be set at the top of each execution script.
