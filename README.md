# Synthetic Difference-in-Differences Analysis in R

This repository contains R code developed for the analysis of
land-use dynamics using the Synthetic Difference-in-Differences
(SDID) methodology.

The code demonstrates programming and data-analysis skills including:

- panel-data preparation and transformation;
- implementation of Synthetic Difference-in-Differences;
- extraction and management of synthetic-control weights;
- construction of treated and synthetic trajectories;
- leave-one-out robustness analysis;
- spatial placebo tests;
- estimation and visualisation of treatment effects;
- publication-ready graphics using ggplot2.

## Methodology

The empirical analysis is implemented using the `synthdid` R package.

The workflow constructs synthetic counterfactual trajectories and
compares them with observed outcomes for treated areas. Robustness
checks include leave-one-out donor analyses and spatial placebo tests.

## Example output

The figure below reports the main estimation results and robustness
checks for one study area.

![Estimation results](figures/crater_area_results.png)

## Software and packages

- R
- synthdid
- dplyr
- ggplot2
- gridExtra
- scales

## Data availability

The underlying dataset is not included in this repository because
redistribution is not permitted.

The repository is intended to demonstrate the code used for data
management, econometric analysis and visualisation.
