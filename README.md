# PoresIdentifiability

Code to perform inference and identifiability analysis for pore closing experiments. Repository will form supplementary material for a preprint.


## Module

### Getting started

Ensure `Julia` is installed and download the repository, in its entirity, to your machine. You should then run `InstallRequiredPackages.jl` from the `Module` folder.

Use the following commands to add the module to your current search path, and load the module:
```
  push!(LOAD_PATH,"path/to/module/folder")  # Add to load path
  using PoreIdentifiability
```
You may wish to add the first line above to your `startup.jl`.

### Module

The module `PoreIdentifiability` provides access to the following functions, each (to be) thoroughly documented.
  - `ExperimentalData()` to load experimental data
  - `PorePDE()` to solve the PDE
  - `SummaryStatistic()` to compute summary statistics from the PDE solution
  - `MLE()` to estimate the maximum likelihood estimate
  - `ProfileLogLikelihood()` to profile the log-likelihood function
