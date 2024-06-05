# GROMACS Molecular Dynamics Introduction Tutorial

This GROMACS introductory project is based on the [official GROMACS introduction to molecular dynamics tutorial](https://tutorials.gromacs.org/md-intro-tutorial.html), adapted slightly in directory structure, molecular simulation trajectory analysis and Python/Jupyter Labs data analysis techniques.

The project involves setting up the topology, energy minimization, NVT & NPT equilibration and finally molecular dynamics simulation of the small protein **Coagulation Factor Xa** in ionized water.

## Introduction

Coagulation Factor X is a serine protease that sits at a pivotal point in the coagulation cascade. It has a role in the three major pathways of the cascade; the intrinsic, extrinsic and common pathways.<sup>1</sup> 

Factor X is a proenzyme and is activated into Factor Xa by hydrolysis, which can occur via the extrinsic or intrinsic pathways<sup>1</sup>.

The 3D structure of Factor Xa is available from the protein data bank (PDB), with the PDB code **1FJS**. The [crystal structure](https://www.rcsb.org/3d-view/1FJS/1) can be downloaded as a PDB file.

## Simulation Commands

### Cleaning the Input Structure
1. Navigate to `1fjs-protein/protein/data/topology` directory
2. Delete non-protein heteroatoms belonging to ligands, cofactors, ions etc. (labelled "HETATM" in the PDB file) by extracting all non-HETATM lines into temp file
	`grep -v HETATM ../input/protein/1fjs.pdb > protein/1fjs_protein_tmp.pdb`
3. Delete non-protein atoms connectivity/bonds (labelled "CONECT" in the PDB file by extracting all non-CONECT lines from temp file into final input structure file
	`grep -v CONECT protein/1fjs_protein_tmp.pdb > protein/1fjs_protein.pdb`

### Generating a Topology
1. The first GROMACS tool we use is `gmx pdb2gmx`, which is used to convert a PDB file to a GROMACS coordinate file, and generate other necessary file like the topology file (`.top`) and a position restraint file (`itp`)
2. We have also specified a water model for GROMACS to use (the TIP3P model, `tip3p`) and a force field (the CHARMM27 all-atom force field, `"charmm27"`)
	`gmx pdb2gmx -f 1fjs_protein.pdb -o ../processed/1fjs_processed.gro -p topol.top -i ../position-restraints/porse.itp -water tip3p -ff "charmm27"`

### Defining the Simulation Box
1. The next GROMACS tool we use is `gmx editconf`, which is used to edit the configuration of molecular structures e.g. defining the simulation box and centering the molecule
2. The command defines the simulation box as a rhombic dodecahedron (`-bt dodecahedron`) and centers the protein in the box (`-c`)
3. A rhombic dodecahedron is chosen as it is more space-efficient than a cubic box, requiring fewer solvent molecules and so reducing computational overhead
4. The command also defines a solute-box edge distance of 1.0 nm (`-d 1.0`) ensures there is at least a 2.0 nm buffer of solvent between any two periodic images of the protein. This is important as if the protein interacts with its periodic image, the calculated forces will be invalid
	`gmx editconf -f 1fjs_processed.gro -o ../simulation-parameters/1fjs_newbox.gro -c -d 1.0 -bt dodecahedron`

### Solvating the Simulation Box with Water
1. The next GROMACS tool we use is `gmx solvate` which is used for adding solvent molecules to a simulation box
2. The protein & simulation box configuration (`-cp`) are contained within the previously generated `1fjs_newbox.gro` file, and the solvent configuration (`-cs`) is part of the standard GROMACS installation
3. Solvent configuration `spc216.gro` is a generic equilibrated 3-point solvent model box, and because out water model (TIP3P) is a three-point water model, we can use it here
4. The solvated system will be output to `1fjs_solv.gro`, and the `topol.top` topology file will also be updated
	`gmx solvate -cp ../simulation-parameters/1fjs_newbox.gro -cs spc216.gro -o ../solvent/1fjs_solv.gro -p topol.top`

Adding Ions
	touch ../solvent/ions.mdp
	gmx grompp -f ../solvent/ions.mdp -c ../solvent/1fjs_solv.gro -p topol.top -o ../solvent/ions.tpr
	printf "SOL\n" | gmx genion -s ../solvent/ions.tpr -o ../solvent/1fjs_solv_ions.gro -conc 0.15 -p topol.top -pname NA -nname CL -neutral

Navigate to 1fjs-protein/energy-minimization/data/processed

Energy minimisation
	gmx grompp -f  ../input/energy-minimization/emin-charmm.mdp -c ../input/topology/solvent/1fjs_solv_ions.gro -p ../input/topology/protein/topol.top -o em.tpr

	gmx mdrun -v -deffnm em

Navigate to 1fjs-protein/energy-minimization/data-analysis

Energy minimisation Data Analysis
	printf "Potential\n0\n" | gmx energy -f ../data/processed/em.edr -o potential.xvg -xvg none

Navigate to 1fjs-protein/protein/topology/protein

Equilibration Run - Temperature (NVT Ensemble Equilibration)
	gmx grompp -f ../../../../nvt-equilibration/data/input/nvt-charmm.mdp -c ../../../../energy-minimization/data/processed/em.gro -r ../../../../energy-minimization/data/processed/em.gro -p topol.top -o nvt.tpr

	mv nvt.tpr ../../../../nvt-equilibration/data/processed
	
	mv mdout.mdp ../../../../nvt-equilibration/data/processed


Navigate to 1fjs-protein/nvt-equilibration/data/processed
	gmx mdrun -ntmpi 1 -v -deffnm nvt
	
Navigate to 1fjs-protein/nvt-equilibration/data-analysis
	echo "Temperature" | gmx energy -f ../data/processed/nvt.edr -o temperature.xvg -xvg none -b 20
	
Navigate to 1fjs-protein/data/input/topology/protein

Equilibration Run - Pressure (NPT Ensemble Equilibration)
	gmx grompp -f ../../../../npt-equilibration/data/input/npt-charmm.mdp -c ../../../../energy-minimization/data/processed/em.gro -r ../../../../energy-minimization/data/processed/em.gro -p topol.top -o npt.tpr

	mv npt.tpr ../../../../npt-equilibration/data/processed
	
	mv mdout.mdp ../../../../npt-equilibration/data/processed

Navigate to 1fjs-protein/nvt-equilibration/data/processed
	gmx mdrun -ntmpi 1 -v -deffnm npt
	
Navigate to 1fjs-protein/nvt-equilibration/data-analysis
	echo "Pressure" | gmx energy -f ../data/processed/npt.edr -o pressure.xvg -xvg none
	
	echo "Density" | gmx energy -f ../data/processed/npt.edr -o density.xvg -xvg none

Navigate to 1fjs-protein/data/input/topology/protein
	gmx grompp -f ../../../../molecular-dynamics/data/input/md-charmm.mdp -c ../../../../npt-equilibration/data/processed/npt.gro -t ../../../../npt-equilibration/data/processed/npt.cpt -p topol.top -o md.tpr
	
	mv md.tpr ../../../../molecular-dynamics/data/processed
	
	mv mdout.mdp ../../../../molecular-dynamics/data/processed

Navigate to 1fjs-protein/molecular-dynamics/data/processed
	gmx mdrun -ntmpi 1 -v -deffnm md

## References
[1] Camire, R.M. (2021) ‘Blood coagulation factor X: Molecular biology, inherited disease, and engineered therapeutics’, *Journal of Thrombosis and Thrombolysis*, 52(2), pp. 383–390. doi:10.1007/s11239-021-02456-w.<br>
