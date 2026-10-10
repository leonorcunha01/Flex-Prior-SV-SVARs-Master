# Flex-Prior-SV-SVARs-Master
Flexible Prior Distributions for Stochastic Volatility in Bayesian Structural Vector Autoregressions

## Description
This is the repository for the Master Thesis 'Flexible Prior Distributions for Stochastic Volatility in Bayesian Structural Vector Autoregressions'. 

The folder SV_SVARs contains the relevant code, Thesis folder contains the thesis pdf.

The author is Leonor Cunha (Msc. Econometrics), supervised by Dr. Jan Hendryk Prüser & Prof. Dr. Christoph Hanck, at TU Dortmund in 2025. 

### Project Description
In Bayesian estimation of SV-identified Structural VARs, the variance of log-volatility processes tends to be time-constant. This thesis investigates whether
such models can be extended with time-variant priors for this parameter, namely
a Horseshoe prior, to better model periods of strong, sudden changes in variables’ behavior such as during the COVID-19 pandemic.

The flexible VAR model is compared to a more standard specification using
a tight Inverse-Gamma prior for the variance. MCMC samplers for each variant
are presented. Challenges arising with this framework are discussed and dealt
with, namely a scaling issue with the impact matrix and the corresponding
implications for VAR tools. 

The variants are compared first in a Monte Carlo
study and then in an empirical application to global oil markets in the pandemic
period, using instrumental variables and narrative measures for shock labeling.
IRFs, FEVDs and HDs are computed and results compared to findings from
prominent literature on the subject.

The flexible prior variant is found to perform slightly better in the simulation
study, but in the empirical application the evidence is not strong enough to
determine a clear winner.

### Code Attribution
This repository uses or adapts code from publicly available implementations associated with the following papers:
- Bertsche, Braun (2020). *Identification of Structural Vector Autoregressions by Stochastic Volatility*
- Chan, Koop and Yu (2024), *Large Order-Invariant Bayesian VARs with Stochastic Volatility*
- Prüser, J. (2021). *The horseshoe prior for time-varying parameter VARs and Monetary Policy*

Where applicable, the original authors' code has been modified to fit the purposes and structure of this repository. Please see the respective papers and original code repositories for the underlying methodology and implementations.

## Repository structure
The repository structure is the following:

```
Flex-Prior-SV-SVARs-Master/
│
├── SV_SVARs/
│   ├── Empirical/
│   │   ├── Braun_oil_app/       # Adapted code from Bertsche & Braun (2020)
│   │   ├── data_raw/            # Raw data and seasonal adjustment script
│   │   ├── estimation/          # Estimation results
│   │   ├── figures/             # MATLAB figures
│   │   ├── tables/              # MATLAB tables
│   │   └── ...                  # Scripts for individual analysis steps
│   │
│   ├── functions/               # Custom MATLAB functions
│   │   ├── func_SVAR7.m         # Main model estimation / sampler
│   │   └── ...
│   │
│   ├── simulations/
│   │   ├── results/             # Saved .mat files (parameters and variables)
│   │   ├── results_aux/         # Auxiliary saved .mat files
│   │   ├── MC_01_study.m        # Simulation study 1
│   │   ├── MC_02_study.m        # Simulation study 2
│   │   └── ...
│   │
│   ├── figures/                 # MATLAB figures from simulations
│   ├── tables/                  # Summary tables
│   └── Main.m                   # Main script running the complete analysis
│
├── README.md
│
└── Thesis/
    └── sv_svars_thesis.pdf
```

#### Notes
Matlab data files (.mat files) with the simulation results and the empirical datasets are included, so the analysis can be run directly and briefly; the simulations can be run from scratch, but this requires Matlab parallel computing toolbox and takes several hours
