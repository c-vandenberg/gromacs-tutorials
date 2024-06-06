# GROMACS Molecular Dynamics Introduction Tutorial

This GROMACS introductory project is based on the [official GROMACS introduction to molecular dynamics tutorial](https://tutorials.gromacs.org/md-intro-tutorial.html), adapted slightly in directory structure, molecular simulation trajectory analysis and Python/JupyterLab data analysis techniques.

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
1. The first GROMACS tool we use is `gmx pdb2gmx`, which is used to convert a PDB file to a GROMACS coordinate file (`.gro`), and generate other necessary file like the topology file (`.top`) and a position restraint file (`itp`)
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

### Adding Ions to the Solvated System
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

### Energy minimisation
1. Before we can begin molecular dynamics, we must ensure that the solvated, electronically neutral system has no steric clashes or inappropriate geometry via energy minimization
2. Navigate to `1fjs-protein/energy-minimization/data/processed`
3. We again use `gmx grompp` to generate the atomic-level binary file (`em.tpr`) that GROMACS can use to run energy minimization
	`gmx grompp -f  ../input/energy-minimization/emin-charmm.mdp -c ../input/topology/solvent/1fjs_solv_ions.gro -p ../input/topology/protein/topol.top -o em.tpr`
4. The next GROMACS tool we use is `gmx mdrun` which is used for running simulations or computations from an atomic-level input binary file (`.tpr` file
5. We specify GROMACS to use verbose terminal output (`-v`), and we specify that the input and all output files have the base name `em` (`-deffnm em`)
	gmx mdrun -v -deffnm em
6. To ensure the energy minimization run was successful, in the end output, the potential energy (`epot`) must be negative & on the order of 1e5, and the maximum force `Fmax` < `emtol` (defined in `emin-charmm.mdp`)

### Energy minimisation Data Analysis
1. Navigate to `1fjs-protein/energy-minimization/data-analysis`
2. The GROMACS tool `gmx energy` is used for extracting and analysing energy terms from an energy (`.edr`) file 
3. Pipe in `Potential` string into `gmx energy` command. Again this isn't needed as GROMACS will prompt you for it if you run it without it
4. Analyse potential energy & time data via Python & JupyterLab
	`printf "Potential\n0\n" | gmx energy -f ../data/processed/em.edr -o potential.xvg -xvg none`

Navigate to 1fjs-protein/protein/topology/protein

### Equilibration Run - Temperature (NVT Ensemble Equilibration)
1. Energy minimization has given us a reasonable starting structure in terms of protein geometry and solvent orientation
2. However, before the MD simulation we must equilibrate the solvent and ions around the protein. If we were to attempt an MD simulation without this, the system may collapse as the solvent is mostly optimised within itself, not necessarily with the solute, and ions have been randomly placed
3. Within the `topol.top` file, there are `porse.itp` position restraint files that apply a position restraining force on the heavy atoms of the protein (i.e. anything that is not H). Movement is permitted, but only after overcoming a substantial energy penalty
4. These position restraints allow us to relax the solvent and ions around the protein, without the added variable of structural changes in the protein
5. To use position restraints, we have added the following line to the NVT equilibration simulation parameter file (`nvt-charmm.mdp`)
	`define                  = -DPOSRES  ; position restrain the protein`
6. Equilibration of temperature is conduted under a NVT ensemble. This ensemble is also called the 'isothermal-isochoric' ensemble
7. The NVT equilibration simulation parameter file (`nvt-charmm.mdp`) uses a modified Berendsen thermostat to control temperature, and specifies a 100 ps NVT equilibration
8. Navigate to `1fjs-protein/protein/topology/protein`
	`gmx grompp -f ../../../../nvt-equilibration/data/input/nvt-charmm.mdp -c ../../../../energy-minimization/data/processed/em.gro -r ../../../../energy-minimization/data/processed/em.gro -p topol.top -o nvt.tpr`
	`mv nvt.tpr ../../../../nvt-equilibration/data/processed`
	`mv mdout.mdp ../../../../nvt-equilibration/data/processed`
9. Navigate to `1fjs-protein/nvt-equilibration/data/processed`
	`gmx mdrun -ntmpi 1 -v -deffnm nvt`
10. The `-ntmpi 1` flag specifies the number of Message Passing Interface (MPI) ranks to use. Here we are using only one, so we are not taking advantage of parallel execution. If we were to specify 4 MPI ranks for example (`-ntmpi 4`), if on a single computer (as is the case here), the simulation workload could be distributed between separate CPU cores

### NVT Ensemble Equilibration Data Analysis
1. Navigate to `1fjs-protein/nvt-equilibration/data-analysis`
2. Analyse temperature & time data via Python & JupyterLab
	`echo "Temperature" | gmx energy -f ../data/processed/nvt.edr -o temperature.xvg -xvg none -b 20`

### Equilibration Run - Pressure (NPT Ensemble Equilibration)
1. The NVT ensemble equilibration stabilises the temperature of the system 
2. Prior to the MD simulation, we must also stabilise the pressure (and thus the density) of the system
3. Equilibration of pressure is conducted under an NPT ensemble. This ensemble is also called the 'isothermal-isobaric' ensemble, and most closely resembles experimental conditions
4. The NPT equilibration simulation parameter file (`npt-charmm.mdp`) uses a Berendsen barostat to control pressure, and specifies a 100 ps NPT equilibration
5. Navigate to `1fjs-protein/data/input/topology/protein`
	`gmx grompp -f ../../../../npt-equilibration/data/input/npt-charmm.mdp -c ../../../../energy-minimization/data/processed/em.gro -r ../../../../energy-minimization/data/processed/em.gro -p topol.top -o npt.tpr`
	`mv npt.tpr ../../../../npt-equilibration/data/processed`
	`mv mdout.mdp ../../../../npt-equilibration/data/processed`
6. Navigate to `1fjs-protein/nvt-equilibration/data/processed`
	`gmx mdrun -ntmpi 1 -v -deffnm npt`

### NPT Ensemble Equilibration Data Analysis
1. Navigate to `1fjs-protein/nvt-equilibration/data-analysis`
2. Extract pressure & time data
	`echo "Pressure" | gmx energy -f ../data/processed/npt.edr -o pressure.xvg -xvg none`
3. Extract density & time data
	`echo "Density" | gmx energy -f ../data/processed/npt.edr -o density.xvg -xvg none`
4. The average system density is close to the experimental value of 100 kg/m<sup>3</sup> and the expected density of the TIP3P model of 1001 kg/m<sup>3</sup>. The density values are also very stable over time, indicating that the system is now well-equilibrated with respect to pressure & density

### Molecular Dynamics (MD) Simulation
1. After the two equilibration phases, the system is now well-equilibrated at the desired temperature (NVT equilibration) and pressure (NPT equilibration)
2. We are now ready to release the position restraints and run the MD simulation
3. The MD simulation parameter file (`md-charmm.mdp`) uses velocity-rescaling temperature coupling as the thermostat and stochastic cell rescaling as the barostat. It also specifies a 1 ns MD simulation
4. Navigate to `1fjs-protein/data/input/topology/protein`
	`gmx grompp -f ../../../../molecular-dynamics/data/input/md-charmm.mdp -c ../../../../npt-equilibration/data/processed/npt.gro -t ../../../../npt-equilibration/data/processed/npt.cpt -p topol.top -o md.tpr`
	`mv md.tpr ../../../../molecular-dynamics/data/processed`
	`mv mdout.mdp ../../../../molecular-dynamics/data/processed`
5. Navigate to `1fjs-protein/molecular-dynamics/data/processed`
	`gmx mdrun -ntmpi 1 -v -deffnm md`

### Molecular Dynamics (MD) Simulation Post-Processing and VMD Processing
1. For post-processing MD simulation analysis, the first GROMACS tool we will use is `gmx trjconv` which is used to manipulate trajectory files e.g. strip out coordinates, correct for periodicity (i.e. fix any parts of the protein that have extended beyond a simulation box with PBC on one side, and re-entered on the opposite side) or manually alter the trajectory (time units, frame frequency etc.)
2. The GROMACS documentation has a [suggested workflow](https://manual.gromacs.org/2021/user-guide/terminology.html?highlight=periodic%20boundary) for `gmx trjconv` to fix periodicity effects
3. We will pipe in the `1\n1\n"` string to select atom 'Group 1' ('Protein') to be centered and atom 'Group 1' ('Protein') to be output to the `center.xtc` file. If we omit this piping GROMACS will prompt us for centering and output groups
4. We use the `-center` flag to center the chosen group of atoms and the `-pbc mol` flag corrects for periodic boundary conditions (PBC). The `mol` option means that whole molecules are kept together when correcting for PBC
5. Navigate to `1fjs-protein/molecular-dynamics/data/processed`
	`printf "1\n1\n" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol`
6. Consolidate files into VMD directory
	`mkdir ../vmd && cp md_center.xtc ../vmd && cd../vmd`
	`cp ../../../protein/data/topology/simulation-parameters/1fjs_newbox.gro .`
7. In directory `fjs-protein/molecular-dynamics/data/vmd` open MD simulation in VMD and modify graphical representation once open
	`vmd 1fjs_newbox.gro md_center.xtc`
8. Render image via File > Renderer > Render the current scene using > Tachyon > Start Rendering
9. Render MP4 via Extensions > Visualizations > Movie Maker then Renderer > Tachyon and Movie Settings > Trajectory (Un-tick 'Delete image files')
10. Once `.ppm` files are generated, navigate to `md-intro-tutorial/scripts` and run `./vmd_movie_processor.sh <vmd_directory_absolute_path>`

### Molecular Dynamics (MD) Simulation Data Analysis
**Protein-Periodic Boundary Minimum Distance**
1. Another GROMACS MD simulation post-processing tool is `gmx mindist`. This calculates the minimum distance between a given group of atoms and its periodic image (i.e. the distance to the periodic boundary). Here we pipe in `1` for the 'Protein' group
	`printf "1\n" | gmx mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg`
2. The distance between the protein and its periodic image should not be smaller than the cut-off distance used to describe non-bonded interactions (in this case, 1.2 nm as per `md-charmm.mdp`)
3. From the data analysis in Python & JupyterLabs we can see that the distance doesn't drop below 1.4 nm

**Root Mean Square Deviation (RMSD)**
1. The post-processing tool `gmx rms` calculates the root mean square deviation (RMSD) of a group of atoms from a reference structure over time
2. RMSD is a measure of the average distance between the atoms of superimposed proteins (usually the backbone of atoms). It is commonly used to assess the similarity between protein structures or to monitor structural changes over time in MD simulations
3. Here the reference structure is the backbone ('Group 4') of the energy minimized topology (`em.tpr`)
    `printf '4\n1\n' | gmx rms -s ../../../energy-minimization/data/processed/em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns -xvg none`

**Radius of Gyration (Rg)**
1. The post-processing tool `gmx gyrate` calculates the radius of gyration (Rg) of a specified group of atoms over time
2. The Rg of a protein is a measure of its compactness. If a protein is stably folded, it will likely maintain a relatively steady value of Rg. If a protein unfolds, its Rg will change over time
3. Again we pipe in `1` for the 'Protein' group
	`echo "1" | gmx gyrate -f md_center.xtc -s md.tpr -o gyrate.xvg -xvg none`

**Index File**
1. The post-processing tool `gmx make_ndx` is used to create and manipulate index groups. Index groups are sets of atoms that are used for various analyses and operations within GROMACS
2. The string `splitch 1` that is piped in splits chain 1 into separate index groups based on criteria such as residues
	`printf "splitch 1\nq\n" | gmx make_ndx -f md.tpr -o`
3. We can now calculate the hydrogen bonds between the two protein chains using `gmx hbond`, with the `-num` flag specifying GROMACS to output the number of hydrogen bonds as a function of time
	
**Report Methods**
1. Once we have run the simulation, it is good practice to report what type of simulation we have performed, as well as the basic system information
2. This can be achieved with the `gmx report-methods` GROMACS tool
	`gmx report-methods -s md.tpr`

## References
[1] Camire, R.M. (2021) 'Blood coagulation factor X: Molecular biology, inherited disease, and engineered therapeutics', *Journal of Thrombosis and Thrombolysis*, 52(2), pp. 383–390.<br>
