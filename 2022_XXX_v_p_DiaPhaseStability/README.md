

﻿This file relates to the research done in the following paper:

S. Borgmans, S. M. J. Rogge, J. S. De Vos, V. Van Speybroeck, _xxxx_, 2022, XX, XX, XXXX
https://doi.org/

This file is part of the longterm storage (publically available under the [CC BY-SA license](https://creativecommons.org/licenses/by-sa/4.0/)) of the input files relevant in this work. In the remainder of this file, we outline in detail the workflow used in this work including references to the relevant input and output files. In case the provided workflow would appear incomprehensible, and lacking to reproduce our results, please contact Sander.Borgmans@ugent.be

The paper deals with the characterization of the relative phase stability in diamondoid covalent organic frameworks, as a function of interpenetration and temperature, considering four different materials:

* COF-300, COF-320, NPN-1, and NPN-3

on which the same workflow is applied. Below, we address the workflow that is applied on each COF individually.

# Software
The following software packages are used to perform all relevant calculations.

- Gaussian (Gaussian 16, Revision A.03) # check with Juul
- HORTON (version 2.1.0) # check with Juul
- QuickFF (version 2.2.0) # check with Juul
- TAMkin (version 1.2.6) # check with Juul
- Yaff (version 1.6.0), with adapted logger functionality, see https://github.com/SanderBorgmans/yaff/tree/logger_overhaul
- Ndfsampler (version 0.0.1), https://github.com/SanderBorgmans/ndfsampler
- Mepsa (version 1.4), http://bioweb.cbm.uam.es/software/MEPSA/
- ThermoLIB (version 1.4.1), with updated wham_2d_scf functionality, see https://github.ugent.be/lvduyfhu/ThermoLIB/tree/boltmanfactor_sb/thermolib
- Dask (version 2021.03.0)

# Workflow

The workflow can be divided into four main parts:
* Structure generation
* Force field derivation and validation
	* ai calculations
	* charge partitioning
	* quickff
	* validation
* Enhanced sampling simulations
* Minimal free energy path analysis

Since the script files are identical for each material, they are consistently stored in a separate scripts folder. Then, for each material, a distinction is made between input and output files, with a subdirectory for each degree of interpenetration and temperature combination, as shown in the directory structure below:


	data/longterm/
		FES_sampling/
            scripts/
            results/
                COF-300/
                    50K/
                        2fold/
                        3fold/
                        ...
                    100K/
                        ...
                    ...
                COF-320/
                    ...
                        ...
                NPN-1/
                    ...
                        ...
                NPN-3/
                    ...
                        ...
			input/
                COF-300/
                    ...
                        ...
                COF-320/
                    ...
                        ...
                NPN-1/
                    ...
                        ...
                NPN-3/
                    ...
                        ...
        FF_development/

		initial_structures/
            COF-300/
            COF-320/
            NPN-1/
            NPN-3/

		MFEP/
            scripts/
            results/
                COF-300/
                    ...
                        ...
                COF-320/
                    ...
                        ...
                NPN-1/
                    ...
                        ...
                NPN-3/
                    ...
                        ...

The following subsections give the details of each of the parts (and subparts) in the workflow.

## STEP 1 - Initial structure generation
In order to consistently compare all possible degrees of interpenetration with respect to each other, all the structures of the different degrees of interpenetration have to be created in a consistent manner. Each structure in `initial_structures/material` is labelled as follows: `dia-cX_block1_block2.cif`, where `X` indicated the degree of interpenetration This is done through:

*SEE JUUL*


## STEP 2 - Force field development (OUTDATED)

This should be updated in accordance with Juul's protocol. (This is the old manual protocol)

### Step 2a - Cluster force fields
#### (i) *Ab initio* input
A geometry optimization followed by a frequency calculation was performed using Gaussian (Gaussian 16, Revision A.03). The molecular structure was generated using MarvinSketch, which was exported to an XYZ file. Then, a Gaussian route section and molecule specification header (with charge and multiplicity) was appended to the XYZ file to arrive at the required Gaussian COM file. Finally, a corresponding job script was created using the `js-g16.sh` script. In summary,  the input, program command and generated output for this step are:

**input**
`gaussian.com`, `gaussian.sh`, `js-g16.sh`, `g09-restart.py`

**command line**
`qsub gaussian.sh`

**output**
`gaussian.fchk`, `gaussian.log`

If the optimization did not converge (as seen in the log file without `Normal termination`), this procedure was repeated, using the `g09-restart.py` script to create a new Gaussian COM file.


#### (ii) Charge partitioning
The fixed atomic charges were estimated with the Minimal Basis Iterative Stockholder partitioning scheme using HORTON (version 2.1.0). The only required input is the Gaussian FCHK file generated in **Step 2a _(i)_**. HORTON then generates the output HDF5 file containing the raw charges for each individual atom in the specified array (array named *mbis* in this case). In summary, the input, program command and generated output for this step are:

**input**
`gaussian.fchk`

**command line**
`horton-wpart.py --grid=veryfine gaussian.fchk horton_out.h5 mbis > horton.log`

**output**
`horton_out.h5`, `horton.log`

#### (iii) QuickFF
- **Atom type assignation**

Instead of the automatic QuickFF atom type assignation protocol, a custom protocol was used, using the `define_atomtypes.py` script. Starting from the `gaussian.fchk` file, CON3F (`c3f.py`) was used to create a chk file, in which the atom types are subsequently defined through the previously mentioned script.

**input**
`gaussian.fchk`

**command line**
`c3f.py gaussian.fchk init.chk; python define_atomtypes.py`

**output**
`sample.chk`

- **Conversion of atomic charges to Yaff parameter file for electrostatics**

Then, the script `qff-ei-input.py` (part of the QuickFF (version 2.2.0) package) was applied to generate a Yaff (version 1.6.0) parameter file containing the electrostatic contribution to the force field in terms of charges of atom types.  Next to the Horton HDF5 output (`horton_out.h5`), this script also requires some user-specified keyword arguments:

* `--bci`
To convert averaged atomic charges to bond charge increments
* `--gaussian`
To characterize the atomic charges as gaussian distributed charges with radii taken from Chen and Martinez.

In summary, the input, program command and generated output for this step are:

**input**
`sample.chk`, `horton_out.h5`

**command line**
`qff-input-ei.py --bci --gaussian sample.chk horton_out.h5:/charges`

**output**
`pars_ei.txt`

- **Derivation of covalent force field terms**

Finally, a covalent force field was estimated using QuickFF using `sample.chk` (containing the equilibrium geometry, Hessian and atom types) and `pars_ei.txt` (containing the electrostatic contribution) as input.

To this end, the script `qff.py` (part of the QuickFF package) was executed using the QuickFF configuration defined in `config.txt`). More information on this configuration files can be found in the [online documentation of QuickFF](http://molmod.github.io/QuickFF/ug.html). QuickFF then generates as output the covalent force field (`pars_yaff.txt`), the system file (system.chk that can be used by Yaff) and a log file (`quickff.log`).

As specified in the `config.txt` files, QuickFF will also generate the trajectories.pp file. This is a intermediary output file (containing the perturbation trajectories in a pickled (binary) file). However, this file has no longterm storage value and is therefore omitted in here.

**input**
`config.txt`, `sample.chk`, `pars_ei.txt`

**command line**
`qff.py -c config.txt sample.chk`

**output**
`pars_cov.txt`, `quickff.log`, `system.chk`


#### (iv) Validation
In order to validate each cluster force field, a basic frequency analysis is performed, comparing the force field and *ab initio* normal mode frequencies. To ensure positive frequencies for the force field, first an optimization is performed using `yopt.py` outputted in `traj.h5`, afterwards the normal modes are calculated through the `freq_analysis.py` script using TAMkin (version 1.2.6). The normal mode frequencies from the Gaussian calculation are calculated from the *ab initio* Hessian in `sample.chk`. The comparison is finally made through `plot_spectrum.py`.

**input**
`sample.chk`, `pars_cov.txt`, `pars_ei.txt`

**command line**
`python yopt.py > yopt.log; python freq_analysis.py > freq_analysis.log; python plot_spectrum.py`

**output**
`yopt.log`, `freq_analysis.log`, `opt.chk`, `comparison_AI_FF.pdf`, `traj.h5`


### Step 2b - Combining cluster force fields into periodic force field
Starting from the building blocks, the diamondoid framework structure constructed, together with a fully flexible *ab-initio* based force field. By combining the cluster force field terms, taking a weighted average based on the number of atoms inside each of the core regions of the building blocks (see supporting information for more details), a periodic force field is created. This is executed through the `combine_pars.py` script. Aside from the parameter files of each of the building blocks, a periodic structure file (chk format) is required, where the atom types have been defined through `define_atomtypes.py`, which identifies how the building blocks are linked to each other.

**input**
`periodic.chk`, `pars_cov_block1.chk`, `pars_ei_block1.txt`, `pars_cov_block2.chk`, `pars_ei_block2.txt`

**command line**
`python define_atomtypes.py; python combine_pars.py`

**output**
`pars_cov.txt`, `pars_ei.txt`, `sample.chk`


Subsequently, the van der Waals contribution to the force field was taken from the MM3 force field of Allinger et al. and added to the covalent and electrostatic contributions *a posteriori*. To automate this process, the cluster structures are loaded into Molden which can extract the MM3 atom types. The so-generated `.xyz` files are then automatically converted into a `Yaff` format through the `mm3_from_tinker.py` script. Finally, the parameter files are combined into a single `pars_mm3.txt` file through manually removing the redundant lines (easily automated).


### Step 2c - Additional force field terms
The cluster force fields, and in turn the derived periodic force field, can give rise to significant deviations between the ab initio cluster model and the optimal force field geometry. When considering the rotation of the imine linked building block in COF-300 and COF-320, it was deemed appropriate to add an additional term. To this end, an *ab inito* rotation scan was performed to fit his additional term, replacing the original torsion term. This procedure is outlined in the `longterm/FF_development/rotational_barrier_diaphasestability.ipynb` file, using a pyiron workflow (https://pyiron.org/).


## Step 3 - Enhanced sampling simulations
After obtaining consistent structural models (in **Step 1**) and their corresponding periodic force fields (in **Step 2**), enhanced molecular dynamics simulations were performed to map the 2D free energy surface describing the relative phase stability of diamondoid COFs as a function of the size and shape of the 1D channels. To this end, umbrella sampling (US) was used, where the ndfsampler software was used to create grid files (containing the umbrella positions and their umbrella strengths) that were subsequently refined based on the sampling overlap of neighbouring simulations.

### Step 3a - Initial grid file
The initial grid files were generated for each material through the `sampling_grid_init.py` script. This script requires the specification of several parameters, some of which can be provided through the command line, while others should be adapted in the script file. An overview of the parameters:

**command line arguments**
- CONFINEMENT_THR: (float) between 0 and 1, indicating the minimum fraction of the phase space that remains in the neighbourhood of the umbrella position (neighbourhood is defined as position +- grid_step/2 in each direction)
- OVERLAP_THR: (float) between 0 and 1, indicating the minimum relative overlap between the sampled phase space of two neighbouring simulations to be considered OK
- KAPPA_GROWTH_FACTOR: (float) between 0 and +inf, scale factor to increase the umbrella strength for each grid point that did not meet the CONFINEMENT_THR
- MAX_GRID_REFINEMENT_IT: (int) maximum number of grid refinements, you should start at 1

**script arguments**
- edges: (dict) specifying the rectangular boundaries of the phase space in units
- kappas: (list) specifying the initial umbrella strengths in kJ/mol/unit
- spacings: (list) the initial grid spacings in units
- units: the units in which all previous parameters are defined, corresponds to the natural unit in which the collective variables are defined
- MD parameters (trivial)

**command line**
`python sampling_grid_init.py CONFINEMENT_THR OVERLAP_THR KAPPA_GROWTH_FACTOR MAX_GRID_REFINEMENT_IT`

**output**
`grid00.txt`, `data.yml`, `run.txt`

The `run.txt` initially is a copy of the `grid00.txt` file, as this file contains all the simulations which have to be run. After each iteration of the post processing step, this `run.txt` file will be updated accordingly.


### Step 3b - Running the enhanced sampling simulations
Aside from the input defined in **Step 3a** (grid files and MD parameters), and the structure definition (structure and the force field parameters), the collective variables and their corresponding potential energy contributions have to be defined as well. This is specified in the `custom_cv.py` script, which has been provided for each material/interpenetration combination. Then, with the provided job script `US_run.sh`, the `US_simulate.py` script is executed on a full node, where Dask distributes the individual US simulations over the available resources that have been assigned by the job script. While the individual simulations can be executed separately as well (for which the standard Yaff version is OK), Dask allows for very easy HPC computation.


**input**
`init.chk`, `run.txt`, `pars.txt`, `data.yml`

**command line**
`python US_simulate.py`

**output**
`trajs/traj_*.h5`, `logs/log_*.txt`

For each individual grid point a trajectory file and log file is generated, as in standard Yaff practice. A collection of all trajectories, defined by the grid points, can then be used in the post processing step.


### Step 3c - Post processing
To update the grid files, and possibly refine the sampling grid, the `US_post.py` is used. This script generates maps of the grid locations and whether there are confinement or overlap issues between certain grid locations.

**command line**
`python US_post.py`

**output**
`grid_points_*.pdf`, `overlap_*.pdf`, `run.txt`

Aside from a new `run.txt` for the next iteration, the original `gridXX.txt` files are updated accordingly.

### Step 3d - FES calculation
After iterating **Steps 3a-3c** several times, the post processing script will reach convergence. At that point you can calculate the free energy surface (FES) using either the provided adaptive WHAM script `WHAM_calc_FES.py` (based on the `wham.yml` input file) or the ThermoLIB script `ThermoLIB_calc_FES.py`. A comparison between the two can be found in the Supporting Information (SI).

#### (i) WHAM code

The former uses the WHAM code from Grossfield *et al.* (http://membrane.urmc.rochester.edu/?page_id=126), and requires you to specify the required parameters for the adaptive algorithm ,as described in the SI, in `wham.yml`. This describes the limits of the rectangular regions taken into account in the WHAM analysis. It can be run through provided job script `WHAM_adaptive_algorithm.sh` or through the command line:


**input**
`trajs/`, `data.yml`, `wham.yml`, `grid*.txt`

**command line**
`python WHAM_calc_FES.py`

**output**
`combined_fes/`

The resulting `combined_fes` folder contains the information on all the individual WHAM regions that were calculated, and the final combined FES. To write the resulting free energy to a readable text file, the `WHAM_write_FES.py` should be executed in the `combined_fes` folder.

#### (ii) ThermoLIB code

The latter uses our in-house ThermoLIB software package, that was adapted for this manuscript to accommodate regions with prohibitively large Boltzmann factors, as elaborated in the SI.

**input**
`trajs/`, `data.yml`, `grid*.txt`

**command line**
`python ThermoLIB_calc_FES.py`

**output**
`colvars/`, `fes.pdf`, `fes.txt`, `metadata`

#### (iii) Grid refinement
Additionally, the MAX_GRID_REFINEMENT_IT can be incremented to introduce finer grid points (with half the original step size) to increase the resolution of the free energy surface. After increasing the MAX_GRID_REFINEMENT_IT parameter in the `data.yml` file, re-execute **Step 3c** to generate a new `run.txt` file and repeat **Steps 3a-3c** until convergence. Moreover, when the rectangular region was not sufficiently large, the grid files can be expanded through the `sampling_grid_expand.py` script.


#### (iv) Data compression
As a lot of individual trajectories will be generated, a compression script `US_compress_data.py` has been provided, that removes all redundant information, and only stores the CV data as a function of time per grid iteration. While this no longer allows for a reconstruction of the original trajectories, it retains all the relevant information for the FES calculation at a fraction of the original storage requirement.

**input**
`trajs/traj_ITERATION_*.h5`, `data.yml`, `grid*.txt`

**command line**
`python US_compress_data.py ITERATION`

**output**
`trajs/compressed_ITERATION.h5`

To avoid data corruption, the individual trajectories are not automatically deleted, but should be by the user.

## Step 4 - Minimal free energy path analysis
Finally, after obtaining the 2D free energy profile, a minimal free energy path (MFEP) can be derived using the MEPSA software package (http://bioweb.cbm.uam.es/software/MEPSA/). While this can be used as a GUI, loading each free energy map individually, a script `derive_MFEP.py` was created to automate the derivation of the minimal free energy path the different (meta)stable states for all structures. Similar to the GUI, this takes the free energy map as input, and, based on a user-prompted identification of the (meta)stable states the MFEPs between these states is calculated and the relevant forward and backward transition barriers are calculated.

**input**
`fes.txt`

**command line**
`python derive_MFEP.py`

**output**
`MFEP.pdf`, `path_transitions.yml`, `path_forward_transition_barriers.yml`, `path_backward_transition_barriers.yml`, `phase_idx.yml`

All the `.yml` files contain user readable information of the different (meta)stable phases, how they were identified by the user, the transition barriers, and their relevant CV paths.

The user-prompted phase identification consists of an image of the minima identified by the MEPSA code superimposed on the FES image, where the minima have been clustered, and the used has to identify which of the clusters corresponds to the 'SQ', 'sq', and 'rect' phases.
