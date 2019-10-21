# README:  GEOCARB-calibration

## Purpose

To improve estimates of Earth-system sensitivity by assimilating deep-time paleocliamte CO2 proxy data with the GEOCARBSULFvolc long-term carbon cycle model. This is done via Bayesian model calibration (Markov chain Monte Carlo).

## Installation and compilation

Markov chain Monte Carlo requires a great many model executions. So, the model is written in Fortran and called from R. After pulling the codes, you will need to navigate to the `fortran` directory and compile the core Fortran model:
* In a terminal window, navigate to the `fortran` directory (`cd fortran`)
* If the `obj` subdirectory is not present, then create it (`mkdir obj`)
* Make sure the file named `Makefile` has the proper directory path to your Fortran-90 compiler (line 5). Two common ones are given, with one commented out (line 6).
* Run the Makefile to compile the model (in the terminal window, enter the command: `make`).

## Workflow

Within the `R` directory, the following routines are needed.

1. Install the relevant packages (`install_packages.R`)
1. Fit skew-normal (and normal) distributions to each data point, including for the supplemental experiment using the data set of Park and Royer (2011) (processData_skewNormal.R and `processData_normal.R`, and `processData_skewNormal_PR2011.R`)
1. Run the sensitivity analysis (`sensitivity_driver.R)`
1. Run the main result calibration (`calib_driver_general.R`)
  1. This requires to run it once from the default parameters (`DO_PARAM_INIT = FALSE`),
  1. then create initial conditions files for the model parameters and the transition covariance matrix (`make_initialization_files_sn-mixture.R`),
  1. and start a new chain from those initialization files, setting filename.`covarinit = "covar_init_sn-mix_12Aug2019.rds"` and `filename.paraminit = "param_init_sn-mix_12Aug2019.rds"`, the initialization files generated in the last step. These initial conditions were found by running a chain with all 69 parameters for 12 million iterations.
  1. Continue that process until the chain is sufficiently spun up.
  1. The batch scripts `rgeocarb_mcmc_aci.pbs` and `rgeocarb_sens_aci.pbs` was used to perform the experiments in this study on the Penn State ACI machines. Feel free to use and modify these scripts for working on otehr high-performance computing systems.
  1. Once the chains are sufficiently converged (by eye and/or by Gelman and Rubin diagnostic (1992)), run a set of parallel chains (changing `nnode_mcmc` in `calib_driver_general.R`). After the initial 12 million iterations (that is, starting from the 12Aug2019 files above), the chains are all okay to not worry about burn-in at that point. The first 15 million iterations from the warm-up chain are ignored for all experiments.
1. The supplemental experiments are all run from those same initial conditions, with appropriate settings in `calib_driver_general.R`.
1. Chop off for burn-in and thin the chains to maintain independent MCMC samples (`process_results.R`).
  1. For the chains that use all 69 parameters (instead of the 7 for Park and Royer (2011)), an interactive routine called `pick_chains.R` will plot the 5 Markov chains' history plots and give the user the option of choosing which ones to keep for analysis. This is because some chains get stuck periodically in local posterior modes. You could keep all chains and just end up with a very high lag for thinning (to maintain independent samples), or could trim out the stuck sections, but we have so many samples already that just discarding the stuck chains is the simplest and fastest option.
1. Run and model ensemble and compute relevant statistics (`analysis.R`)
1. Make plots (`plotting.R` and `plotting_sobol.R`)

## Redistribution and use

Please redistribute, please modify and please use! This code is part of GEOCARB-calibration.
GEOCARB-calibration is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GEOCARB-calibration is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with GEOCARB-calibration.  If not, see <http://www.gnu.org/licenses/>.

## Questions?

There will always be questions. Whether you are trying to get the model to run, or wondering what a line of code does, we would love to talk to you. Do not hesitate to reach out to us at:
* Tony Wong (aewsma@rit.edu)
* Ying Cui (cuiy@montclair.edu)
