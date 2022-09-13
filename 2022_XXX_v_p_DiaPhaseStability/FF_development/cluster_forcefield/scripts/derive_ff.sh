#!/bin/sh
#
#PBS -N _derive_ff
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=1
#PBS -m n

### Script to create force field and perform some validation
### This script will run all necessary scripts to create a force field for your system if the necessary files are provided
### Subsequently, a validation of the force field is done by comparing the ab initio and force field frequencies and rest values

### Necessary files
###      Ab initio input: either Gaussian FCHK File or Vasp XML file from frequency job
###      Definition of force field atom types (.py script with an afilters object)
###      Tinker MM3 xyz file with MM3 force field atom type numbers, can be created with Molden
###
### These files should be given as input via the options
### 		-d  --directory: new working directory, where all necessary files are provided
### 		-a  --abinitio: Gaussian FCHK File or VASP XML file
###          	-f  --ffatypes: Definition of FFatypes
###          	-t  --tinker: Tinker XYZ File
###
###     Example
###     	local : bash derive_ff.sh -d path/to/BDBA_boronate -a BDBA_boronate_freq.fchk -f get_ffatype_filters.py -t BDBA_boronate_mm3.xyz
###	        hpc   : qsub -v directory="path/to/BDBA_boronate",abinitio="BDBA_boronate_freq.fchk",ffatypes="get_ffatype_filters.py",tinker="BDBA_boronate_mm3.xyz" derive_ff.sh
###
### 	If a file with the name config.txt is present in the folder, this will be used as quickff config file


### Dependencies: Con3F, horton, QuickFF, TAMkin
### In October 2019 these packages are supported on the HPC Clusters golett, phanpy and swalot

date

scripts=.

###################################    INPUT    ############################

POSITIONAL=()
while [[ $# -gt 0 ]] # While number of arguments is greater than 0
do
key="$1"

case $key in
	-d|--directory)
	directory="$2"
	shift
	shift
	;;
	-a|--abinitio)
	abinitio="$2"
	shift
	shift
	;;
	-f|--ffatypes)
	ffatypes="$2"
	shift
	shift
	;;
	-t|--tinker)
	tinker="$2"
	shift
	shift
	;;
	*)
	POSITIONAL+=("$1")
	shift
	;;
esac
done

echo "Creating force field from following files:"
echo ""
echo "Directory              : ${directory}"
echo ""
echo "FFatype definitions    : ${ffatypes}"
echo "Tinker MM3 definitions : ${tinker}"
echo "Ab Initio input        : ${abinitio}"
echo ""
echo ""

################################    PREPARATION    #########################

echo "Switching to directory ${directory}"
cd "$directory"

if [ ! -f ${ffatypes} ]; then
	echo "File ${ffatypes} not in directory ${directory}"
	exit 1
fi
if [ ! -f ${tinker} ]; then
        echo "File ${tinker} not in directory ${directory}"
        exit 1
fi
if [ ! -f ${abinitio} ]; then
        echo "File ${abinitio} not in directory ${directory}"
        exit 1
fi

echo "Storing backup"
cd ../
dirname="${directory##*/}"
cp -r "${dirname}" "${dirname}"_BACKUP
cd "${dirname}"

echo "Creating directory structure"
[[ -d AI ]] || mkdir AI
[[ -d Derivation ]] || mkdir Derivation
[[ -d Validation ]] || mkdir Validation

cd Derivation
[[ -d programs ]] || mkdir programs
[[ -d trajectories ]] || mkdir trajectories
cd ../

cd Validation
[[ -d programs ]] || mkdir programs
[[ -d frequencies ]] || mkdir frequencies
[[ -d rest_values ]] || mkdir rest_values
cd ../

# Put everything in the AI folder (assumed these are AI calculations), except for the necessary files
mv *.* AI
cd AI
mv "${ffatypes}" ../
mv "${tinker}" ../
cp "${abinitio}" ../
if [ -f config.txt ]; then
	mv config.txt ../
fi
cd ../

echo ""

############################    DERIVE FORCE FIELD    ######################

echo "FORCE FIELD DERIVATION"

if [[ "${machine}" == node* ]] || [[ "${machine}" == gligar*.gastly.os ]]; then
	module load yaff/1.6.0-intel-2020a-Python-3.8.2
	module load Con3F/1.0-20190329-intel-2019a-Python-3.7.2
fi
c3f.py convert ${abinitio} system.chk

if [[ "${machine}" == node* ]] || [[ "${machine}" == gligar*.gastly.os ]]; then
	module purge
        module load yaff/1.6.0-intel-2020a-Python-3.8.2
fi

### Atom Types
echo "    Defining atom types"
python "${scripts}"/read_atomtypes.py system.chk "${ffatypes}" > Derivation/programs/ffatypes.log # FFatypes in system_ffatype.chk

mv "${ffatypes}" Derivation
rm "${ffatypes}"c
mv system.chk system_ai.chk # All information stored in system_ffatype.chk

### Van der Waals Interaction
echo "    Creating Van der Waals part of Force Field using tinker MM3"
python "${scripts}"/mm3_from_tinker.py system_ffatype.chk "${tinker}" # pars_mm3.txt
mv "${tinker}" Derivation

### Electrostatic Interaction
if [[ "${machine}" == node* ]] || [[ "${machine}" == gligar*.gastly.os ]]; then
        module purge
	module load horton/2.0.0-intel-2015b-Python-2.7.10
fi

echo "    Partitioning charge distribution"
horton-wpart.py --grid veryfine ${abinitio} horton_out.h5 mbis > Derivation/programs/horton.log # Charges in horton_out.h5

if [[ "${machine}" == node* ]] || [[ "${machine}" == gligar*.gastly.os ]]; then
        module purge
        module load QuickFF/2.2.4-intel-2020a-Python-3.8.2
fi

echo "    Creating electrostatic part of Force Field using Horton"
qff-input-ei.py --bci --gaussian system_ffatype.chk horton_out.h5:/charges # pars_ei.txt
mv horton_out.h5 Derivation/programs

### Covalent interaction

echo "    Creating covalent part of Force Field using QuickFF"
if [ -f config.txt ]; then
	echo "QuickFF config file detected: config.txt"
	config="config.txt"
else
	config="${scripts}"/config_quickff.txt
fi

python3 "${scripts}"/qff-derive-cov.py -c "${config}" -e pars_ei.txt -m pars_mm3.txt -p system_ffatype.chk "${abinitio}" > Derivation/programs/quickff.log # pars_yaff.txt

# Cleaning Up
rm "${abinitio}" # Already backup in AI folder
mv system_ffatype.chk Derivation/
mv trajector*.png Derivation/trajectories
mv trajector*.pp Derivation/trajectories

echo ""

###########################    VALIDATE FORCE FIELD    #####################

echo "FORCE FIELD VALIDATION"

echo "    Optimizing structure using Force Field"
python "${scripts}"/optimize.py -c pars_yaff.txt -e pars_ei.txt -m pars_mm3.txt -o system_opt.chk system.chk > Validation/programs/optimization.log

### Frequencies

echo "Frequency comparison"

# Force field frequencies are written to frequencies_FF.txt
if [[ "${machine}" == node* ]] || [[ "${machine}" == gligar*.gastly.os ]]; then
        module purge
        module load TAMkin/1.2.6-intel-2020a-Python-3.8.2
	module load yaff/1.6.0-intel-2020a-Python-3.8.2
fi

echo "    Performing FF NMA"
python "${scripts}"/nma.py -c pars_yaff.txt -e pars_ei.txt -m pars_mm3.txt --freq freqs_ff.txt system_opt.chk > Validation/programs/nma_ff.log
echo "    Performing AI NMA"
python "${scripts}"/nma.py --freq freqs_ai.txt system_ai.chk > Validation/programs/nma_ai.log

# Plot spectrum
python "${scripts}"/validate_frequencies.py -o freqs_ai_ff.txt -f freqs_ai_ff.pdf freqs_ai.txt freqs_ff.txt

# Cleaning up
mv freqs* Validation/frequencies
mv system_ai.chk AI
mv system_opt_nma.chk Validation
mv system_ai_nma.chk Validation

### Rest values

echo "Rest value comparison"

python "${scripts}"/validate_rest_values.py -o rv_ai_ff.txt --ff_bond bond_ff.txt --ff_bend bend_ff.txt --ff_dihed dihed_ff.txt --ff_oop oop_ff.txt --ai_bond bond_ai.txt --ai_bend bend_ai.txt --ai_dihed dihed_ai.txt --ai_oop oop_ai.txt --fig_bond bond_ai_ff.pdf --fig_bend bend_ai_ff.pdf --fig_dihed dihed_ai_ff.pdf --fig_oop oop_ai_ff.pdf system.chk system_opt.chk > Validation/programs/comparison.log

mv rv_ai_ff.txt Validation/rest_values
mv bond* Validation/rest_values
mv bend* Validation/rest_values
mv dihed* Validation/rest_values
mv oop* Validation/rest_values

date

echo "####################################################################################"
echo "####################################################################################"

