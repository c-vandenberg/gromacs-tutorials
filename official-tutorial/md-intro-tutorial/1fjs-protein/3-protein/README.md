# 3. Generating Topologies, Defining Simulation Box & Solvating System


<div align="center">
  <img src="https://github.com/c-vandenberg/gromacs-tutorials/assets/60201356/cd2d1faf-87e3-4078-acf9-683c3b465e66" alt="1fjs-protein-crystal-structure" width="">
</div>

**Fig 1** Coagulation Factor Xa Crystal Structure

## 3.1 Introduction

Prior to any simulation run, we must define our system by:
1. Cleaning the input protein crystal structure of non-protein atoms
2. Generating the protein topology by converting the cleaned PDB file to a GROMACS coordinate file (`.gro`)
3. Defining the simulation box
4. Solvating the simulation box with water
5. Adding ions to the solvated system

## 3.2 Simulation Commands

### 3.2.1 Cleaning the Input Structure
1. Navigate to `1fjs-protein/protein/data/topology` directory
2. Delete non-protein heteroatoms belonging to ligands, cofactors, ions etc. (labelled "HETATM" in the PDB file) by extracting all non-HETATM lines into temp file
	`grep -v HETATM ../input/protein/1fjs.pdb > protein/1fjs_protein_tmp.pdb`
3. Delete non-protein atoms connectivity/bonds (labelled "CONECT" in the PDB file by extracting all non-CONECT lines from temp file into final input structure file
	`grep -v CONECT protein/1fjs_protein_tmp.pdb > protein/1fjs_protein.pdb`

### 3.2.2 Generating a Topology
1. The first GROMACS tool we use is `gmx pdb2gmx`, which is used to convert a PDB file to a GROMACS coordinate file (`.gro`), and generate other necessary file like the topology file (`.top`) and a position restraint file (`itp`)
2. We have also specified a water model for GROMACS to use (the TIP3P model, `tip3p`) and a force field (the CHARMM27 all-atom force field, `"charmm27"`)
	`gmx pdb2gmx -f 1fjs_protein.pdb -o ../processed/1fjs_processed.gro -p topol.top -i ../position-restraints/porse.itp -water tip3p -ff "charmm27"`

### 3.2.3 Defining the Simulation Box
1. The next GROMACS tool we use is `gmx editconf`, which is used to edit the configuration of molecular structures e.g. defining the simulation box and centering the molecule
2. The command defines the simulation box as a rhombic dodecahedron (`-bt dodecahedron`) and centers the protein in the box (`-c`)
3. A rhombic dodecahedron is chosen as it is more space-efficient than a cubic box, requiring fewer solvent molecules and so reducing computational overhead
4. The command also defines a solute-box edge distance of 1.0 nm (`-d 1.0`) ensures there is at least a 2.0 nm buffer of solvent between any two periodic images of the protein. This is important as if the protein interacts with its periodic image, the calculated forces will be invalid
	`gmx editconf -f 1fjs_processed.gro -o ../simulation-parameters/1fjs_newbox.gro -c -d 1.0 -bt dodecahedron`

### 3.2.4 Solvating the Simulation Box with Water
1. The next GROMACS tool we use is `gmx solvate` which is used for adding solvent molecules to a simulation box
2. The protein & simulation box configuration (`-cp`) are contained within the previously generated `1fjs_newbox.gro` file, and the solvent configuration (`-cs`) is part of the standard GROMACS installation
3. Solvent configuration `spc216.gro` is a generic equilibrated 3-point solvent model box, and because out water model (TIP3P) is a three-point water model, we can use it here
4. The solvated system will be output to `1fjs_solv.gro`, and the `topol.top` topology file will also be updated
	`gmx solvate -cp ../simulation-parameters/1fjs_newbox.gro -cs spc216.gro -o ../solvent/1fjs_solv.gro -p topol.top`

### 3.2.5 Adding Ions to the Solvated System
1. From the `topol.top` file, we can see that the protein has a net charge of -2e (+1e for chain A, and -3e for chain L). Since life does not exist at a net charge, we must add ions to our system and approximate physiological conditions using 0.15 M NaCl
2. The next GROMACS tool we use is `gmx grompp` which is the GROMACS preprocessor command. It prepares an input file for MD simulations by processing it and generating an atomic-level binary file (`.tpr`) that GROMACS can use to run simulations
3. A `.mdp` file is normally used to prepare a `.tpr` binary file for energy minimization or an MD simulation using `gmx grompp`, here we just need an empty `.mdp` file to generate the `ions.tpr` file
4. We also pass in the solvated system GROMACS coordinate file (`1fjs_solv.gro`), and the topology file (`topol.top`)
	`touch ../solvent/ions.mdp`
	`gmx grompp -f ../solvent/ions.mdp -c ../solvent/1fjs_solv.gro -p topol.top -o ../solvent/ions.tpr`
5. The next GROMACS tool we use is `gmx genion` which is used to add ions to the system
6. We pipe in the `SOL` string into the `gmx genion` to specify that only solvent molecules are to be replaced, ensuring no protein atoms are replaced. This is not necessarily needed as GROMACS will prompt you for this if you run the command without it
7. We pass in the previously generated `ions.tpr` file, along with the solvated system GROMACS coordinate file (`1fjs_solv.gro`), and the topology file (`topol.top`)
8. We specify that the cations are Na<sup>+</sup> (`-pname NA`), the anions are Cl<sup>-</sup> (`-nname CL`), and approximate physiological conditions using 0.15 M concentration (`-conc 0.15`)
9. Finally, we ensure GROMACS adds enoughb ions to neutralize the total charge of the system (`-neutral`)
	`printf "SOL\n" | gmx genion -s ../solvent/ions.tpr -o ../solvent/1fjs_solv_ions.gro -conc 0.15 -p topol.top -pname NA -nname CL -neutral`
