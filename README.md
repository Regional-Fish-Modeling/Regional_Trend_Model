# Regional Trend Model

## Objectives

## Run

1. `run_gompertz.R`
2. `2_Random walk panel plot.R`

## Model Structure

*Density Dependent*

Abundance allows for density-dependent change over time following Gompertz function where `r` represents the intrinsic population growth rate (lumps survival and recruitment) and `K` represents the semi-equalibrium abundance (Hostetler and Chandler 2015 Ecology).

Population growth `r` is modeled as a function of seasonal covariates and a random site effect (intercept). The seasonal effects are currently fixed but we could explore random effects by site or HUC. Ultimately in the future it would be great to add spatial dependence to these.

`K` is modeled as a random effect by site but I'd like to add drainage area as a covariate.

*Density Independent*

The density independent effects allow for immigration to rescue populations following local extirpation. They `rho` are currently modeled as random by site-year, but I'm not sure if this will converge. It represents the average number of immigrants to the population independent of the current population size. It could be a fixed value or random by site and or year rather than site-year. It could also be random at the HUC (8, 10, 12) level, but it gets messy to read too much into it or spend too much time on covariates since it can be a function of so many things related to the surrounding populations and network connectivity.
