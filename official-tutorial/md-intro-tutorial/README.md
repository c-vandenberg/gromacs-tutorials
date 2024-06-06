# GROMACS Molecular Dynamics Introduction Tutorial

This GROMACS introductory project is based on the [official GROMACS introduction to molecular dynamics tutorial](https://tutorials.gromacs.org/md-intro-tutorial.html), adapted slightly in directory structure, molecular simulation trajectory analysis and Python/JupyterLab data analysis techniques.

The project involves setting up the topology, energy minimization, NVT & NPT equilibration and finally molecular dynamics simulation of the small protein **Coagulation Factor Xa** in ionized water.

## Introduction

Coagulation Factor X is a serine protease that sits at a pivotal point in the coagulation cascade. It has a role in the three major pathways of the cascade; the intrinsic, extrinsic and common pathways.<sup>1</sup> 

Factor X is a proenzyme and is activated into Factor Xa by hydrolysis, which can occur via the extrinsic or intrinsic pathways<sup>1</sup>.

The 3D structure of Factor Xa is available from the protein data bank (PDB), with the PDB code **1FJS**. The [crystal structure](https://www.rcsb.org/3d-view/1FJS/1) can be downloaded as a PDB file.

## Simulation Commands

### Equilibration Run - Pressure (NPT Ensemble Equilibration)
1. The NVT ensemble equilibration stabilises the temperature of the system 
2. Prior to the MD simulation, we must also stabilise the pressure (and thus the density) of the system
3. Equilibration of pressure is conducted under an NPT ensemble. This ensemble is also called the 'isothermal-isobaric' ensemble, and most closely resembles experimental conditions
4. The NPT equilibration simulation parameter file (`npt-charmm.mdp`) uses the Berendsen barostat to control pressure, and specifies a 100 ps NPT equilibration
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
