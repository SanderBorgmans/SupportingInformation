

﻿This file relates to the research done in the following paper:

S. Borgmans, S. M. J. Rogge, J. S. De Vos, P. Van Der Voort, V. Van Speybroeck, _xxxx_, 2022, XX, XX, XXXX
https://doi.org/

This file is part of the longterm storage (publically available under the [CC BY-SA license](https://creativecommons.org/licenses/by-sa/4.0/)) of the input files relevant in this work. In the remainder of this file, we outline in detail the workflow used in this work including references to the relevant input and output files. In case the provided workflow would appear incomprehensible, and lacking to reproduce our results, please contact Sander.Borgmans@ugent.be

The paper deals with the characterization of the relative phase stability in diamondoid covalent organic frameworks, as a function of interpenetration and temperature, considering four different materials:

* COF-300, COF-320, NPN-1, and NPN-3

on which the same workflow is applied. Below, we address the workflow that is applied on each COF individually.

# Software
The following software packages are used to perform all relevant calculations.

- Gaussian (Gaussian 16, Revision C.01)
- HORTON (version 2.0.0)
- Molden (version 6.8)
- QuickFF (version 2.2.4)
- TAMkin (version 1.2.6)
- Yaff (version 1.6.0), with adapted logger functionality, see https://github.com/SanderBorgmans/yaff/tree/logger_overhaul
- Ndfsampler (version 0.0.1, legacy branch), https://github.com/SanderBorgmans/ndfsampler/tree/legacy_diaphase
- Mepsa (version 1.4), http://bioweb.cbm.uam.es/software/MEPSA/
- ThermoLIB (version 1.4.1), with updated wham_2d_scf functionality, see https://github.ugent.be/lvduyfhu/ThermoLIB/tree/boltmanfactor_sb/thermolib


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

complemented by two additional parts:
* CV measurement - for the calculation and connection of the collective variables to system parameters
* Water loading simulations - for the calculations of the FES for water filled frameworks

Since the script files are identical for each material, they are consistently stored in a separate scripts folder. Then, for each material, a distinction is made between input and output files, with a subdirectory for each degree of interpenetration and temperature combination, as shown in the directory structure below:


	data/longterm/
	
		FES_sampling/
			input/
				COF-300/
					1fold/
					2fold/
					...
				COF-320/
					...
				NPN-1/
					...
				NPN-3/
					...
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
					
        FF_development/
			abinitio_input/
				scripts/
				... # folder per building block
			cluster_forcefield/
				scripts/
				... # folder per building block
			rotational_barriers/
			
		initial_structures/
			input/
				sbus/
					... # folder per building block
				topologies/
				
			output/
				COF-300/
					... # folder per degree of interpenetration
				COF-320/
					... # folder per degree of interpenetration
				NPN-1/
					... # folder per degree of interpenetration
				NPN-3/
					... # folder per degree of interpenetration
			scripts/

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

		CV_measurement
			link_CV_pore/
			generate_CV/
			AuthorYear/

		water_loading/
			input/
			scripts/
			results/
				loading_XX/
				...

The following subsections give the details of each of the parts (and subparts) in the workflow.

## STEP 1 - Cluster force field development

### Step 1a - QuickFF cluster force fields
#### (i) *Ab initio* input

The *ab initio* hessian for each SBU is calculated with Gaussian (Gaussian 16, Revision C.01). For each molecular building block, an XYZ file is generated with Avogadro and converted to a Gaussian COM file using the `g09-xyztocom.py` script. An initial fast optimization with the PBE functional is followed by an optimization with the B3LYP functional with additional Grimme D3 dispersion correction. If the optimization does not converge (as seen in the log file without `Normal termination`), this procedure is repeated, using the `g09-restart.py` script to create a new Gaussian COM file. Afterwards, the hessian is computed by a frequency calculation. For each COM file, a corresponding job script can be created using the `js-g16.py` script.

**input**
`SBU.xyz`

**command line**
`python g09-xyztocom.py SBU.xyz`
`python js-g16.py SBU.com`
`qsub SBU.sh`
`python g09-restart.py SBU.log`
`python g09-restart.py -l B3LYP -d GD3 -o SBU_b3lyp.com SBU.log`
`python g09-restart.py -j 'freq(noraman) -o SBU_freq.com SBU_b3lyp.com'`

**output**
`SBU.com`, `SBU.sh`, `SBU.fchk`, `SBU.log`
`SBU_b3lyp.com`, `SBU_b3lyp.sh`, `SBU_b3lyp.fchk`, `SBU_b3lyp.log`
`SBU_freq.com`, `SBU_freq.sh`, `SBU_freq.fchk`, `SBU_freq.log`


#### (ii) Charge partitioning
The fixed atomic charges are estimated with the Minimal Basis Iterative Stockholder partitioning scheme using HORTON (version 2.0.0). The only required input is the Gaussian FCHK file generated in **Step 1a _(i)_**. HORTON then generates the output HDF5 file containing the raw charges for each individual atom.

**input**
`SBU_freq.fchk`

**command line**
`horton-wpart.py --grid=veryfine SBU_freq.fchk horton_out.h5 mbis > horton.log`

**output**
`horton_out.h5`, `horton.log`

#### (iii) QuickFF
- **Atom type assignation**

Starting from the `SBU_freq.fchk` file, CON3F (`c3f.py`) was used to create an CHK file, from which the bonds are detected using the detect_bonds() method in the `molmod` package. From the molecular graph, all sets of equivalent atoms are detect using the `detect_ffatypes.py` script and each set was assigned a different atom type. The cluster termination has to be manually defined by means of the indices of the atoms that are bonded to the core of the building block and one index of an atom in that core. The atom types of the atoms in the termination end with `_term`, while the atoms in the core of the SBU are given a label identifying the SBU itself (`_BDC`, `_BPDC`, `_TAM`, `_TNM` or `_TNA`). The final atom types are written to the `SBU_ffatypes.txt` file, which is finally used by the `read_atomtypes.py` script to generate the CHK file `system_ffatype.chk` that contains the atom types.

**input**
`SBU_freq.fchk`

**command line**
`c3f.py convert SBU_freq.fchk system.chk`
`python detect_ffatypes -i CORE -o TERM1 -o TERM2 ... --output SBU_ffatypes.txt system.chk`
`python read_atomtypes.py system.chk SBU_ffatypes.txt > ffatypes.log`

**output**
`system.chk`, `SBU_ffatypes.txt`, `system_ffatype.chk`

- **Conversion of atomic charges to Yaff parameter file for electrostatics**

Then, the script `qff-ei-input.py` (part of the QuickFF (version 2.2.4) package) is applied to generate a Yaff (version 1.6.0) parameter file containing the electrostatic contribution to the force field in terms of charges of atom types.  Next to the Horton HDF5 output (`horton_out.h5`), this script also requires some user-specified keyword arguments:

* `--bci`
To convert averaged atomic charges to bond charge increments
* `--gaussian`
To characterize the atomic charges as gaussian distributed charges with radii taken from Chen and Martinez.

**input**
`system_ffatype.chk`, `horton_out.h5`

**command line**
`qff-input-ei.py --bci --gaussian system_ffatype.chk horton_out.h5:/charges`

**output**
`pars_ei.txt`

- **Description of the Van der Waals interactions**

The Van der Waals interactions are described using a Buckingham potential with the MM3 parameters from Allinger *et al.*. To obtain these parameters, the MM3 atom types are identified using the Tinker package as implemented in Molden (version 6.8) and written to `mm3.xyz`. Once these atom types are identified, the parameters are read and converted to the Yaff parameter file `pars_mm3.txt`using the `mm3_from_tinker.py` script.

**input**
`system_ffatype.chk`, `mm3.xyz`

**command line**
`python mm3_from_tinker.py system_ffatype.chk mm3.xyz`

**output**
`pars_mm3.txt`


- **Derivation of covalent force field terms**

Finally, a covalent force field is estimated using QuickFF using `system_ffatype.chk` (containing the equilibrium geometry and atom types) , `SBU_freq.fchk` (containing the *ab initio* hessian) and the reference force field parameter files `pars_ei.txt` and `pars_mm3.txt` (containing the electrostatic and Van der Waals contribution respectively) as input.

To this end, the script `qff-derive-cov.py` is executed using the QuickFF configuration defined in `config_quickff.txt`). More information on this configuration files can be found in the [online documentation of QuickFF](http://molmod.github.io/QuickFF/ug.html). QuickFF then generates as output the covalent force field (`pars_yaff.txt`), the system file (`system.chk` that can be used by Yaff) and a log file (`quickff.log`).

As specified in the `config_quickff.txt` file, QuickFF will also generate the `trajectories.pp` file. This is an intermediary output file (containing the perturbation trajectories in a pickled (binary) file). However, this file has no longterm storage value and is therefore omitted in here.

**input**
`config.txt`, `system_ffatype.chk`, `SBU_freq.fchk`, `pars_ei.txt`, `pars_mm3.txt`

**command line**
`qff-derive-cov.py -c config_quickff.txt -e pars_ei.txt -m pars_mm3.txt -p system_ffatype.chk SBU_freq.fchk > quickff.log`

**output**
`pars_yaff.txt`, `quickff.log`, `system.chk`


#### (iv) Validation

In order to validate each cluster force field, the rest values of the force field and *ab initio* optimized systems are compared. Furthermore, a basic frequency analysis is performed. First, the *ab initio* optimized system is relaxed using the derived QuickFF force field with the `optimize.py`script, ensuring that only positive frequencies are obtained.

**input**
`system.chk`, `pars_yaff.txt`, `pars_ei.txt`, `pars_mm3.txt`

**command line**
`python optimize.py -c pars_yaff.txt -e pars_ei.txt -m pars_mm3.txt -o system_opt.chk system.chk > optimization.log`

**output**
`system_opt.chk`, `optimization.log`

- **Rest value comparison**

The rest values (bonds, bends, dihedral angles and out-of-plane distances) for both the *ab initio* and force field optimized clusters are determined and compared with each other using the `validate_rest_values.py` script. All values are printed in the TXT files and a summary is given in the `rv_ai_ff.txt`file. For all types of internal coordinates, a figure is made to show the agreement between both structures.

**input**
`system.chk`, `system_opt.chk`

**command line**
`python validate_rest_values.py -o rv_ai_ff.txt --ff_bond bond_ff.txt --ff_bend bend_ff.txt --ff_dihed dihed_ff.txt --ff_oop oop_ff.txt --ai_bond bond_ai.txt --ai_bend bend_ai.txt --ai_dihed dihed_ai.txt --ai_oop oop_ai.txt --fig_bond bond_ai_ff.pdf --fig_bend bend_ai_ff.pdf --fig_dihed dihed_ai_ff.pdf --fig_oop oop_ai_ff.pdf system.chk system_opt.chk > comparison.log`

**output**
`rv_ai_ff.txt`, `bond_ff.txt`, `bend_ff.txt`, `dihed_ff.txt`, `oop_ff.txt`, `bond_ai.txt`, `bend_ai.txt`, `dihed_ai.txt`, `oop_ai.txt`, `bond_ai_ff.pdf`, `bend_ai_ff.pdf`, `dihed_ai_ff.pdf`, `oop_ai_ff.pdf`, `comparison.log`

- **Frequency comparison**

A basic frequency analysis is performed, comparing the force field and *ab initio* normal mode frequencies. The *ab initio* normal mode frequencies are obtained from the *ab initio* hessian that is stored in `system.chk` and formatted using the `nma.py` script. For the force field frequencies, a normal mode analysis is performed through the `nma.py` script using TAMkin (version 1.2.6). The comparison is finally made through the `validate_frequencies.py` script.

**input**
`system.chk`, `system_opt.chk`, `pars_yaff.txt`, `pars_ei.txt`, `pars_mm3.txt`

**command line**
`python nma.py --freq freqs_ai.txt system.chk > nma_ai.log`
`python nma.py -c pars_yaff.txt -e pars_ei.txt -m pars_mm3.txt --freq freqs_ff.txt system_opt.chk > nma_ff.log`
`python validate_frequencies.py -o freqs_ai_ff.txt -f freqs_ai_ff.pdf freqs_ai.txt freqs_ff.txt`

**output**
`freqs_ai.txt`,  `freqs_ff.txt`, `freqs_ai_ff.txt`, `freqs_ai_ff.pdf`, `nma_ff.log`, `nma_ai.log`

### Practical implementation
The whole procedure from Step 1a _(ii)_ - 1a _(iv)_ is automatically performed by a single bash script, `derive_FF.sh` for convenience.

### Step 1b - Additional force field terms
The cluster force fields, and in turn the derived periodic force field, can give rise to significant deviations between the ab initio cluster model and the optimal force field geometry. When considering the rotation of the imine linked building block in COF-300 and COF-320, it was deemed appropriate to add an additional term. To this end, an *ab inito* rotation scan was performed to fit his additional term, replacing the original torsion term. This procedure is outlined in the `longterm/FF_development/rotational_barrier_diaphasestability.ipynb` file, using a pyiron workflow (https://pyiron.org/).

### Step 1c - UFF cluster force fields

Besides the system-specific QuickFF cluster force fields, the fully transferable UFF force field is used to derive an additional force field for the clusters. This is used for an initial optimization of the periodic structure to avoid collapse of closely placed atoms, which would be troublesome to describe with the Buckingham potential in the Van der Waals part of the QuickFF force field. As UFF uses a Lennard-Jones potential, nearly overlapping atoms can be pulled apart during an optimization.

This UFF force field is generated using the `create_uff.py` and `uff.py` scripts, which identify the UFF atom types and assign the correct parameters. The detected atom types and bond orders, from which the UFF parameters are derived, are written out to the `ffatypes_uff.txt` and `bonds.txt`files. If necessary, the atom types and bond orders can be redefined in these files, which are read in in subsequent runs of the `create_uff.py` script. The parameters are written in the Yaff parameter files `pars_cov_uff.txt` and `pars_lj_uff.txt`.

**input**
`system.chk`

**command line**
`python create_uff.py -f ffatypes_uff.txt -b bonds.txt -c pars_cov_uff.txt -l pars_lj_uff.txt system.chk`

**output**
`pars_cov_uff.txt`, `pars_lj_uff.txt`, `ffatypes_uff.txt`, `bonds.txt`

## STEP 2 - Initial structure generation

### STEP 2a - Structure assembly

The initial structures are generated using our in-house structure assembly software. A first version of this software will be published soon. For the `SBU`, `Topology`, `ParametersCombination` and `Construct` modules, which are necessary to run the `run.py` and `dia_cN.py` scripts, we refer to the upcoming publication.

The program requires two `input` types: SBUs and topologies. The atomic representation of each SBU is given in the CHK files, whereas the points of extension are defined in the SBU files. The cluster force fields of each SBU, which are derived in **STEP 1**, are used to generate a force field for the periodic structure. In order to consistently compare all possible degrees of interpenetration with respect to each other, the interpenetrated **dia** nets should be defined in a consistent manner. This is done in the `dia_cN.py`script.

Each structure is labelled as follows: `dia-cX_block1_block2`, where `X` indicates the degree of interpenetration and `block1`and `block2`are the SBUs that are used to generate the structure. For NPN-1 and NPN-2, a single building block is used and `block2`is given the label `None`. Besides the initial structures, the program also derives a force field for the periodic structure from the cluster force fields of the individual building blocks. Both a QuickFF force field `pars_cluster.txt`and a UFF force field `pars_uff.txt` are generated.

**input**
`blockN.chk`, `blockN.sbu`, `pars_yaff.txt`, `pars_ei.txt`, `pars_mm3.txt`, `pars_cov_uff.txt`, `pars_lj_uff.txt`, `dia_cN.py`

**command line**
`python run.py`

**output**
`dia-cX_block1_block2.chk`, `pars_cluster.txt`, `pars_uff.txt`

### STEP 2b - Structure optimization

Once each structure is generated, it is relaxed using the following Yaff protocol. First, the system is optimized using the UFF force field to push atoms that are placed closely together apart (e.g. in high degrees of interpenetration) and prevent collapse. Secondly, the structure is optimized with the QuickFF force field. When appropriate (for COF-300 and COF-320) the `pars_polysix.txt` were added to the `pars_cluster.txt` force field for the optimization and subsequent simulations.

**input**
`dia-cX_block1_block2.chk`, `pars_cluster.txt`, `pars_uff.txt`(, `pars_polysix.txt`)

**command line**
`qsub opt.sh` or `python opt.py`

**output**
`dia-cX_block1_block2_opt.chk`


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


# Complementary 1 - CV measurement
An additional folder `CV_measurement` was added in accordance with the reported collective variables of Table 1 in the Supporting Information, with a subfolder per paper. These folders contain all the reference structures from literature and the scripts to calculate the 'experimental' values for our collective variables.

**input**
`example.chk`

**command line**
`python measure_CV.py example.chk`

As output, the relevant information is printed to the terminal. The conversion of the cif files to chk files was performed with CON3F. As the collective variable requires atom types to be specified, `define_atomtypes_nohess.py` scripts are present per material which determine and store the atom types in the chk file. For the NPN COFs an additional script `adapt_chk.py` is provided, as there certain atoms needed to be replaced by their periodic equivalent for the CV function to work. The lattice vectors of each of these materials can be readily accessed from the cif or chk files.

Aside from the calculation of the CV for each reference structure, scripts were provided to perform a specific umbrella sampling simulation to drive to system towards a specific CV combination

**input**
`custom_cv.py`, `init.chk`, `pars.txt`

**command line**
python simulate.py cv1 cv2 kappa1 kappa2

Finally, a script was provided to calculate the system parameters (cell parameters, pore volume, CVs) for a given trajectory file, to provide an overview of how these parameters relate to each other.

# Complementary 2 - water loading
To replicate the experimental phase behaviour of COF-300(7), additional simulations were performed with a TIP-4P based force field. To this end, the relevant force field files were taken from the supporting information of doi:10.1038/s41563-021-00977-6. To couple the existing force field files of the framework to the water force field, the `pars_framework_water_coupling.txt` was created, with LJCROSS for all water-framework atom pairs, and LJ or MM3 parameters with a zero contribution for the water and framework atoms repsectively. The US simulations are performed similar to the other US simulations, but the force field is replaced by a GhostForceField, takinjg ghost atoms into account for an adequate representation of the water molecules. A script `insert_molecules.py` was provided to generate initial water filled structures.

