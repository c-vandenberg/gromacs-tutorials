### GROMACS Molecular Dynamics Introduction Tutorial

This GROMACS introductory project is based on the [official GROMACS introduction to molecular dynamics tutorial](https://tutorials.gromacs.org/md-intro-tutorial.html), adapted slightly in directory structure, molecular simulation trajectory analysis and Python/Jupyter Labs data analysis techniques.

The project involves setting up the topology, energy minimization, NVT & NPT equilibration and finally molecular dynamics simulation of a small protein (1FJS) in ionized water.

## Simulation Commands
Navigate to 1fjs-protein/protein/data/topology directory

Cleaning the input structure
	grep -v HETATM ../input/protein/1fjs.pdb > protein/1fjs_protein_tmp.pdb
	grep -v CONECT protein/1fjs_protein_tmp.pdb > protein/1fjs_protein.pdb

Navigate to 1fjs-protein/protein/data/topology/protein directory
Generating a topology
	gmx pdb2gmx -f 1fjs_protein.pdb -o ../processed/1fjs_processed.gro -p topol.top -i ../position-restraints/porse.itp -water tip3p -ff "charmm27"

Solvating the simulation system
	gmx editconf -f 1fjs_processed.gro -o ../simulation-parameters/1fjs_newbox.gro -c -d 1.0 -bt dodecahedron

Filling the box with water
	gmx solvate -cp ../simulation-parameters/1fjs_newbox.gro -cs spc216.gro -o ../solvent/1fjs_solv.gro -p topol.top

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

