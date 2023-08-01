This file relates to the research done in the following paper:

S. Borgmans, S. M. J. Rogge, L. Vanduyfhuys, V. Van Speybroeck, xxx (2023)

This file is part of the longterm storage (publically available under the [CC BY-SA license](https://creativecommons.org/licenses/by-sa/4.0/)) of the input files relevant in this work. In the remainder of this file, we outline in detail the workflow used in this work including references to the relevant input and output files. In case the provided workflow would appear incomprehensible, and lacking to reproduce our results, please contact Sander.Borgmans@ugent.be

The paper reports on a new python software package, named OGRe, to construct an optimal sampling grid of umbrella sampling simulations. This optimization procedure significantly increases the accuracy, and provides an easy-to-use protocol based on three clear metrics that originate from the underlying weighted histogram analysis method (WHAM) equations to calculate the free energy surface (FES).

This protocol was benchmarked on four analytic potentials:
* two 1D potentials: skewed bimodal, trimodal
* two 2D potentials: double well, Ackley

and subsequently applied on two physical systems:
* proton hopping in a zeolite, from https://doi.org/10.1038/s41467-023-36666-y, comparing with the published results
* COF-5, to establish its efficacy to characterize layer stacking in 2D COFs


# Software
The following software packages are used to perform all relevant calculations.

- Yaff (version 1.6.0), with adapted logger functionality, see https://github.com/SanderBorgmans/yaff/tree/logger_overhaul
- OGRe (version 0.0.1, yaff branch), https://github.com/SanderBorgmans/OGRe/tree/yaff
- ThermoLIB (version 1.4.1), with updated wham_2d_scf functionality, see https://github.ugent.be/lvduyfhu/ThermoLIB/tree/boltmanfactor_sb/thermolib
- Schnet (see conda yml file, for MLP simulations)


# Workflow
The workflow for each system can be divided into three parts:
* Input generation
* Simulations
* Post-processing

where the last two steps are repeated until convergence.


Since the input and simulation scripts differ slightly for each material, as they define the initial parameters and define the energy evaluation, they are stored for each system separately. In contrast, the post-processing scripts are identical, and are thus stored in a common script folder. Additionally, as both OGRe and Schnet were used in a conda environment, separate `environment.yml` files are provided in the scripts folder.

data/longterm/
	scripts/
	benchmarking/
		1D/
			skewed_bimodal/
			trimodal/
			scripts/
		2D/
			ackley/
			double_well/
			scripts/
	applications/
		proton_hopping/
		layer_stacking/


In contrast to normal OGRe use, where a single parameter combination is postulated, the benchmarking required many different parameter combination to be considered. As such, instead of storing and replacing trajectories in-place, a database was used, where all trajectories pertaining to the same grid spacing and kappa growth factors are combined. 


The following subsections give the details of each of the parts (and subparts) in the workflow.

## STEP 1 - Benchmarking
### Step 1a - setup
For each potential a `potential.py` file needs to be constructed, which defines the analytic potential as a class, with a function `internal_compute` computing the potential and its partial derivatives towards its arguments. This is illustrated by the `potential.py` scripts for each of the benchmarking potentials. Second, a `database_init.py` is defined, which stores all simulation parameters for this system, and considers several user arguments to set the hyperparameters.

### Step 1b - parameter screening
Then the `scan.sh` script is used to screen several parameter combinations, creating the setup through the `database_init.py` script. Afterward, this bash script loops over a `database_simulate_load.py` script which attempt to load the relevant simulations from the database, complemented by `database_simulate.py` for the entries that are missing (storing these new simulations into the database for later use). Finally, every iteration ends with a post-processing step to determine whether additional simulations are required through the `database_post.py` script. An `adapt_cores.py` script is also provided which monitors the number of simulations required, and adapts the number of cores in the run scripts, as the simulations are run under the hood with the worker module to minimize the number of submitted jobs.

### Step 1c - analysis
Finally, figures documenting the performance of each parameter combination can be created through the error_plots.py script. 

## STEP 2 - Applications
### Step 2a - setup
Similar to the benchmarking a certain setup is required for the individual applications, which is specific for each application. The OGRe software package is intended to be used through a `custom_cv.py` script, which defines a `yaff.pes.colvar.CollectiveVariable` object as the collective variable (cv), in combination with a force field to evaluate the energy at each instance (defined through `pars_*.txt` files in yaff format). However, to illustrate that OGRe can also interface with other simulation engines, the proton hopping system required interfacing with plumed and schnetpack for the energy evaluation through a machine learned potential. To this end, the internal `OGRe_Simulation` class is extended, as illustrated in the `database_simulate_load.py` script for the proton hopping application. 

### Step 2b - simulations
In correspondence to the benchmarking simulations, a particular parameter combination was fixed in the `database_init.py` scripts, and the OGRe protocol was applied until converged. Again, the worker module was used to run each series of simulations, through the `run.pbs` script, and the post-processing was performed through the `database_post.sh` script. The remaining simulations are always stored as a list in the `run.txt` file, which is translated to a `grid_restart.txt` file when the `database_simulate_load.py` script is applied. 

The final `layerxx.txt` files are also provided, to allow for a reproduction. Evidently, as there might be some numerical variations through the yaff engine, it is advised to rerun the post-processing on the final data to establish whether the obtained collection of simulations is representative.


### Step 2c - analysis
The analysis was mainly facilitated through the final FES, which is either calculated for every iteration when a CONSISTENCY_THR is defined, or manually through the `database_fes.sh` script. 

## Step 3 - Compression
Several ways exist to compress the trajectory data created by the simulations, as the required disk space for all possible combination can accumulate to untractable amounts. Either `compress_database_trajs.py` can be used to extract the `cv_values` array from each trajectory (with a chunked hdf5 storage for optimal compression), or the OGRe provided `ogre_compress_iteration.py` can be used to create a single `.h5` file with the cv trajectory data of a full layer.