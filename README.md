# VonBertalanffy
R script for the estimation of the VonBertalanffy growth curve

The function uses the age-length data by species to estimate the growth parameters by mean of the Von Bertalanffy growth curve.
The results of the analysis are automatically saved in the user defined working directory.

The analysis is conducted by sex. Combined sexes ('C') are composed aggregating all the available sexes.
The analysis by sex is conducted aggregating undetermined sex ('I') with each other (e.g. M = M + I; F = F + I).
