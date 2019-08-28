# README

## Workflow

1. Install the relevant packages (install_packages.R)
1. Fit skew-normal (and normal) distributions to each data point, including for the supplemental experiment using the data set of Park and Royer (2011) (processData_skewNormal.R and processData_normal.R, and processData_skewNormal_PR2011.R)
1. Run the sensitivity analysis (sensitivity_driver.R)
1. Run the main result calibration (calib_driver_sn-mixture.R)
  1. This requires to run it once from the default parameters (DO_PARAM_INIT = FALSE),
  1. then create initial conditions files for the model parameters and the transition covariance matrix (make_initialization_files_sn-mixture.R),
  1. and start a new chain from those initialization files (calib_driver_sn-mixture.R, with DO_PARAM_INIT = FALSE and set filename.covarinit and filename.paraminit to the initialization files generated in the last step).
  1. Be careful not to delete your old results! You will want them for plotting the full history plot.
  1. Continue that process until the chain is sufficiently spun up.
  1. The batch scripts `rgeocarb_mcmc_aci.pbs` and `rgeocarb_mcmc_all.pbs` are good examples of submission scripts for working on high-performance computing systems (experiments in this study were performed on the Penn State ACI machines)
  1. Once the chains are sufficiently converged (by eye and/or by Gelman and Rubin diagnostic (1992)), run a set of parallel chains (changing nnode_mcmc in calib_driver_sn-mixture.R)
1. Chop off for burn-in and thin the chains to maintain independent MCMC samples (burnin_thinning.R)
1. Run and model ensemble and compute relevant statistics (analysis.R)
1. Make plots (plotting.R and plotting_sobol.R)

## Questions?

Tony Wong (aewsma@rit.edu) and Ying Cui (cuiy@montclair.edu)
