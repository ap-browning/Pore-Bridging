# Pores Bridging

Code to perform the model-based data analysis with the Porous-Fisher to explore pore bridging by osteoblastic cells in 3D printed scaffolds. This repository is supplementary material for the preprint "Model-based data analysis of osteoblastic tissue growth in shallow 3D printed scaffolds" available on [bioRxiv](https://www.biorxiv.org).

The majority of the code contains the `Julia` module `Pores` (`Module/Pores.jl`) that exports the experimental data and functions used to perform the analysis. The `Results` folder contains code that performs the analysis displayed in the main document (the `Results/Supplementary` folder does the same for the supplementary material document). The `Figures` folder contains scripts to produce and export a `.pdf` for each figure in the document. 


## Installation

Ensure `Julia` is installed (see *Required software*) and download the repository, in its entirety, to your machine. You should then run `Install_Required_Packages.jl` from the `Module` folder.

Use the following commands to add the module to your current search path, and load the module:
```
  push!(LOAD_PATH,"/path/to/module/folder/")  # Add to load path
  using Pores                                 # Load module
```
If using Windows, ensure to escape the backslashes in the path: `C:\\path\\to\\module\\folder`, or use Unix style forward slashes.


## Module
The module `Pores` provides access to the following functions, each thoroughly documented.
  - `SolvePDE()` to solve the Porous-Fisher model numerically using the method of lines. Output is a callable function, i.e. `u(t)` that can be called to give the full, spatial, solution at any time `t` as a 2D array.
  - `SummaryStatistics()` to summarise the numerical, spatial, solution `u(t)` with four summary statistics (see main and supplementary documents).
  - `Maximise()` perform maximisation on a function (i.e., maximum likelihood estimation).
  - `FitModel()` to fit a simple model to data (i.e., the quadratic noise function).
  - `Profile()` to profile a log-likelihood function.
  - `LogLikelihood()` to calculate the log-likelihood for a given set of parameters `θ`, data `Data`, and pore size `L`.
  - `GetΘ()` to convert partitioned parameter space (θ,β) into full parameter space Θ. Here, θ might be unknown parameters (i.e., θ = (D,λ,K,u₀)), and β known, fixed, parameters.
  - `ConfidenceInterval()` to compute an approximate confidence interval from a profile likelihood.

The module also exports the following variables:
  - `Data` a `DataFrame` containing the full experimental data set (this is also available at `Module/Data.csv`).
  - `P` an array of all pore sizes.
  - `T` an array of all observation times.
  - `S` a set of all summary statistics.
  

## Results

All results in the main (and supporting material) document can be obtained by running the corresponding script in the `Results` folder (or `Results/Supplementary` folder). The output of each script is stored in a `JLD2` file in the `/Saved` directories (which are called when producing figures).

In general, each `JLD2` file of results contains two variables. For the `MLE_Individual` results:
  - `MLE_Individual` contains the result.
  - `MLE_Individual_Settings` contains information about the parameters used to produce the result, and the time used to produce the result.

The `Options.jl` file contains global options used to produce all results (i.e., parameter bounds, maximisation tolerances, etc).

Before running scripts, ensure the working directory of the `Julia` session, `pwd()` is set to the results folder:
```
  cd("/path/to/module/folder")
```
Approximate runtimes for the full computation of the results for each model (using a 3.7GHz Quad-Code i7 desktop running Windows 10), and figures produced, are given below.

| Script        | Description                                         | Runtime    |
| ------------- | --------------------------------------------------- | ---------- |
| `MLE.jl`      | Compute individual and combined MLEs                | ~24 hours  |
| `Profiles.jl` | Compute individual and combined profile likelihoods | ~9.5 hours |
| `Profiles.jl` | Compute bivariate profiles                          | ~16 hours  |


Note that the code uses the `.Threads` module to run the analysis simultaneously on CPU threads. Use the `nthreads()` command (run `using .Threads` first) to verify the number of threads in the `JULIA_NUM_THREADS` environment variable. For more information on setting the number of threads in Julia visit [julialang.org](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_NUM_THREADS-1).


## Figures

Each figure in the main document (except for purely experimental or graphical figures) can be produced from the corresponding script in the `Figures` folder (and `Figures/Supplementary`), which uses the `PyPlot` package (which calls `matplotlib`). The `Preferences.jl` file contains graphical tweaks to obtain the figure style used in the main document.


## Required software

  - `Julia` can be downloaded from [julialang.org](https://julialang.org/downloads/) or on macOS using `homebrew`: just run `brew cask install julia` in terminal.
  - All `Julia` packages used are available from the standard package installed. Run `Module/Install_Required_Packages.jl` in `Julia` to ensure all required packages are installed.
