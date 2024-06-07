# 5. Molecular Dynamics Simulation

## 5.1. Introduction

After the NVT & NPT equilibration runs, the system is now well-equilibrated at the desired temperature (NVT equilibration) and pressure/density (NPT equilibration). We are now ready to release the position restraints and run the production molecular dynamics simulation for data collection.

A molecular dynamics (MD) simulation predicts the motion of atoms/molecules at a given temperature for a given simulation time length. 
* The basic idea is that by heating the system to a certain temperature, molecules in the system have sufficient energy to cross potential energy barriers and move across the PES to adopt new energy minimum structures/conformations. 
* This is similar to **Monte Carlo Simulation** energy minimization, with the main difference being we are generating structures/conformations sequentially in time (i.e. dynamics) in MD.

At the heart of MD simulations are **Newton's equations of motions**, which describe how the positions & velocities of particles evolve over time. MD simulations make use of **integration algorithms** such as the **Velocity Verlet algorithm** to integrate and solve Newton's equations of motion.

The general steps of the Velocity Verlet algorithm **<sup>1</sup>** are:
1. **Initial Position & Velocity**
    * Each particle/atom in the molecule is given an initial position & velocity
    * Using a force field (e.g. CHARMM 27), the forces acting on each atom at this initial position are calculated (sum of the bonding & non-bonding interactions)
    * From these force field calculations, the potential energy at this initial position is calculated, similar to energy minimization
2. **Update Positions**
    * The Velocity Verlet algorithm updates the positions of each particle/atom in the molecule using the current positions, velocities & accelerations (accelerations are derived from the forces):

    <br>
    <div align="center">
      <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Br%7D%28t%20%2B%20%5CDelta%20t%29%20%3D%20%5Cmathbf%7Br%7D%28t%29%20%2B%20%5Cmathbf%7Bv%7D%28t%29%5CDelta%20t%20%2B%20%5Cfrac%7B1%7D%7B2%7D%20%5Cmathbf%7Ba%7D%28t%29%20%5CDelta%20t%5E2", alt="velocity-verlet-algorithm-position-update-equation"/>
    </div>

    where:
    * ![new_position](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Br%7D%28t%20%2B%20%5CDelta%20t%29) is the new position at time ![new_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t%20%2B%20%5CDelta%20t)
    * ![current_position](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Br%7D%28t%29) is the current position at time *t*
    * ![current_velocity](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%28t%29) is the velocity at time *t*
    * ![timestep](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5CDelta%20t) is the timestep
    * ![acceleration](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Ba%7D%28t%29%20%3D%20%5Cfrac%7B%5Cmathbf%7BF%7D%28t%29%7D%7Bm%7D) is the acceleration at time *t*, ![force_at_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%29) is the force acting on the particle/atom at time *t*, and *m* is the mass of the particle/atom

      <br>

3. **Velocity Update (Half-step)**
    * The Velocity Verlet algorithm calculates the velocities of each particle/atom in the molecule at the half-time step
    * The half-step velocity update is performed after the positions are updated but before the new forces (and hence new accelerations) are calculated. It uses the current accelerations to calculate & update the velocities of each particle/atom

    <br>
    <div align="center">
      <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%5Cleft%28t%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7B2%7D%5Cright%29%20%3D%20%5Cmathbf%7Bv%7D%28t%29%20%2B%20%5Cfrac%7B1%7D%7B2%7D%20%5Cmathbf%7Ba%7D%28t%29%20%5CDelta%20t", alt="velocity-verlet-algorithm-half-step-velocity-equation"/>
    </div>


    where:
    * ![half_step_velocity](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%5Cleft%28t%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7B2%7D%5Cright%29) is the velocity at the half-step ![half_step_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7B2%7D)
    * ![current_velocity](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%28t%29) is the velocity at time *t*
    * ![half_step_vacceleration](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Ba%7D%28t%29%20%3D%20%5Cfrac%7B%5Cmathbf%7BF%7D%28t%29%7D%7Bm%7D) is the acceleration at time *t*, where ![force_at_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%29) is the force acting on the particle/atom at time *t*, and *m* is the mass of the particle/atom

      <br>

4. **Force Calculation**
   * The Velocity Verlet algorithm then calculates the forces ![force_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%20%2B%20%5CDelta%20t%29) at the new position ![position_at_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Br%7D%28t%20%2B%20%5CDelta%20t%29):

      <br>
      <div align="center">
         <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%20%2B%20%5CDelta%20t%29%20%3D%20-%5Cfrac%7B%5Cpartial%20U%28%5Cmathbf%7Br%7D%29%7D%7B%5Cpartial%20%5Cmathbf%7Br%7D%28t%20%2B%20%5CDelta%20t%29%7D", alt="velocity-verlet-algorithm-force-equation"/>
      </div>

    where:
    * ![force_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%20%2B%20%5CDelta%20t%29) is the force at time ![current_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bt%7D%28t%20%2B%20%5CDelta%20t%29)
    * ![negative_pot_energy_gradient](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20-%5Cfrac%7B%5Cpartial%20U%28%5Cmathbf%7Br%7D%29%7D%7B%5Cpartial%20%5Cmathbf%7Br%7D%28t%20%2B%20%5CDelta%20t%29%7D) is the negative gradiant of the potential energy ![pot_energy](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20U%28%5Cmathbf%7Br%7D%29) with respect to the position ![position](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Br%7D) at ![time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t%20%2B%20%5CDelta%20t). The gradient is calculated by taking the first-order partial derivative with respect to the position ![position_at_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Br%7D%28t%20%2B%20%5CDelta%20t%29)
   
      <br>

5. **Velocity Update (Full-step)**
   * Finally, the Veloctiy Verlet algorithm updates the velocities to the full-time step using the new accelerations obtained from the calculations in step 4:
  
      <br>
      <div align="center">
         <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%28t%20%2B%20%5CDelta%20t%29%20%3D%20%5Cmathbf%7Bv%7D%5Cleft%28t%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7B2%7D%5Cright%29%20%2B%20%5Cfrac%7B1%7D%7B2%7D%20%5Cmathbf%7Ba%7D%28t%20%2B%20%5CDelta%20t%29%20%5CDelta%20t", alt="velocity-verlet-algorithm-full-step-velocity-equation"/>
      </div>

   where:
   * ![full_step_velocity](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%28t%20%2B%20%5CDelta%20t%29) is the full step velocity at time ![full_step_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t%20%2B%20%5CDelta%20t)
   * ![half_step_velocity](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%5Cleft%28t%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7B2%7D%5Cright%29) is the velocity at the half-step ![half_step_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7B2%7D)
   * ![full_step_acceleration](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Ba%7D%28t%20%2B%20%5CDelta%20t%29%20%3D%20%5Cfrac%7B%5Cmathbf%7BF%7D%28t%20%2B%20%5CDelta%20t%29%7D%7Bm%7D) is the acceleration at time ![full_step_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t%20%2B%20%5CDelta%20t), where ![full_step_force](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%20%2B%20%5CDelta%20t%29) is the force acting on the particle/atom (calculated in step 4) at time ![full_step_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t%20%2B%20%5CDelta%20t), and *m* is the mass of the particle/atom

      <br>

In MD simulations, because the particles/atoms have a velocity, it is not sufficient to describe the total energy of the system in terms of just potential energy. Due to the motion of atoms/particles, to calculate the total system energy we must also take into account **kinetic energy**. Therefore, total system energy in MD, at time *t*, is given by:

   <br>
   <div align="center">
      <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20E%28t%29%20%3D%20P_%7BE%7D%28t%29%20%2B%20K_%7BE%7D%28t%29", alt="md-total-system-energy-equation"/>
   </div>

where:
* ![total)_system_energy](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20E%28t%29) is the total system energy at time *t*
* ![system_pot_energy](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20P_%7BE%7D%28t%29) is the system potential energy at time *t*
* ![system_kinetic_energy](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20K_%7BE%7D%28t%29) is the system kinetic energy at time *t*

The system kinetic energy at time *t* is calculated using the velocities obtained by the Velocity Verlet algorithm:

   <br>
   <div align="center">
      <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20K_%7BE%7D%28t%29%20%3D%20%5Csum_%7Bi%3D1%7D%5EN%20%5Cfrac%7B1%7D%7B2%7D%20m_i%20v_i%5E2%28t%29", alt="md-system-kinetic-energy-equation"/>
   </div>

   where:
   * ![system_kinetic_energy](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20KE_%7BE%7D%28t%29) is the total kinetic energy of all particles/atoms in the system at time *t*
   * ![atom_summation](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Csum_%7Bi%3D1%7D%5EN) is the summation of all atoms/particles in the system
   * ![atom_mass](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20m_i) is the mass of atom *i*
   * ![atom_velocity_squared](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20v_i%5E2%28t%29) is the squared velocity of of atom *i* at time *t*

## 5.2. Simulation Commands

### 5.2.1. Molecular Dynamics (MD) Simulation
1. The MD simulation parameter file (`md-charmm.mdp`) uses velocity-rescaling temperature coupling as the thermostat and stochastic cell rescaling as the barostat. It also specifies a 1 ns MD simulation
2. Navigate to `1fjs-protein/data/input/topology/protein` and run:
    * `gmx grompp -f ../../../../molecular-dynamics/data/input/md-charmm.mdp -c ../../../../npt-equilibration/data/processed/npt.gro -t ../../../../npt-equilibration/data/processed/npt.cpt -p topol.top -o md.tpr`
    * `mv md.tpr ../../../../molecular-dynamics/data/processed`
    * `mv mdout.mdp ../../../../molecular-dynamics/data/processed`
3. Navigate to `1fjs-protein/molecular-dynamics/data/processed`
    * `gmx mdrun -ntmpi 1 -v -deffnm md`

### 5.2.2. Molecular Dynamics (MD) Simulation Post-Processing and VMD Processing
1. For post-processing MD simulation analysis, the first GROMACS tool we will use is `gmx trjconv` which is used to manipulate trajectory files e.g. strip out coordinates, correct for periodicity (i.e. fix any parts of the protein that have extended beyond a simulation box with PBC on one side, and re-entered on the opposite side) or manually alter the trajectory (time units, frame frequency etc.)
2. The GROMACS documentation has a [suggested workflow](https://manual.gromacs.org/2021/user-guide/terminology.html?highlight=periodic%20boundary) for `gmx trjconv` to fix periodicity effects
3. We will pipe in the `1\n1\n"` string to select atom 'Group 1' ('Protein') to be centered and atom 'Group 1' ('Protein') to be output to the `center.xtc` file. If we omit this piping GROMACS will prompt us for centering and output groups
4. We use the `-center` flag to center the chosen group of atoms and the `-pbc mol` flag corrects for periodic boundary conditions (PBC). The `mol` option means that whole molecules are kept together when correcting for PBC
5. Navigate to `1fjs-protein/molecular-dynamics/data/processed` and run:
   * `printf "1\n1\n" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol`
6. Consolidate files into VMD directory:
	* `mkdir ../vmd && cp md_center.xtc ../vmd && cd../vmd`
	* `cp ../../../protein/data/topology/simulation-parameters/1fjs_newbox.gro .`
7. In directory `fjs-protein/molecular-dynamics/data/vmd` open MD simulation in VMD and modify graphical representation once open run:
	* `vmd 1fjs_newbox.gro md_center.xtc`
8. Render image via File > Renderer > Render the current scene using > Tachyon > Start Rendering
9. Render MP4 via Extensions > Visualizations > Movie Maker then Renderer > Tachyon and Movie Settings > Trajectory (Un-tick 'Delete image files')
10. Once `.ppm` files are generated, navigate to `md-intro-tutorial/scripts` and run:
    * `./vmd_movie_processor.sh <vmd_directory_absolute_path>`

### 5.2.3. Molecular Dynamics (MD) Simulation Data Analysis
**Protein-Periodic Boundary Minimum Distance**
1. Another GROMACS MD simulation post-processing tool is `gmx mindist`. This calculates the minimum distance between a given group of atoms and its periodic image (i.e. the distance to the periodic boundary). Here we pipe in `1` for the 'Protein' group:
	* `printf "1\n" | gmx mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg`
2. The distance between the protein and its periodic image should not be smaller than the cut-off distance used to describe non-bonded interactions (in this case, 1.2 nm as per `md-charmm.mdp`)

<div align="center">
  <img src="https://github.com/c-vandenberg/gromacs-tutorials/assets/60201356/87b89839-7157-440b-aeb9-520224db3a2a" alt="" width="protein-periodic-boundary-distance-vs-time">
</div>
   
3. From the data analysis in Python & JupyterLabs we can see that the distance doesn't drop below 1.4 nm
<br>

**Root Mean Square Deviation (RMSD)**
1. The post-processing tool `gmx rms` calculates the root mean square deviation (RMSD) of a group of atoms from a reference structure over time
2. RMSD is a measure of the average distance between the atoms of superimposed proteins (usually the backbone of atoms). It is commonly used to assess the similarity between protein structures or to monitor structural changes over time in MD simulations
3. Here the reference structure is the backbone ('Group 4') of the energy minimized topology (`em.tpr`):
    * `printf '4\n1\n' | gmx rms -s ../../../energy-minimization/data/processed/em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns -xvg none`
  
<div align="center">
  <img src="https://github.com/c-vandenberg/gromacs-tutorials/assets/60201356/7dba8748-baf7-4932-af10-edc40eaf8827" alt="rmsd-vs-time" width="">
</div>
<br>

**Radius of Gyration (Rg)**
1. The post-processing tool `gmx gyrate` calculates the radius of gyration (Rg) of a specified group of atoms over time
2. The Rg of a protein is a measure of its compactness. If a protein is stably folded, it will likely maintain a relatively steady value of Rg. If a protein unfolds, its Rg will change over time
3. Again we pipe in `1` for the 'Protein' group:
	* `echo "1" | gmx gyrate -f md_center.xtc -s md.tpr -o gyrate.xvg -xvg none`

<div align="center">
  <img src="https://github.com/c-vandenberg/gromacs-tutorials/assets/60201356/e2c0a6ec-380a-489b-a8d4-6205220977c9" alt="radius-gyration-vs-time" width="">
</div>
<br>

**Index File**
1. The post-processing tool `gmx make_ndx` is used to create and manipulate index groups. Index groups are sets of atoms that are used for various analyses and operations within GROMACS
2. The string `splitch 1` that is piped in splits chain 1 into separate index groups based on criteria such as residues:
	* `printf "splitch 1\nq\n" | gmx make_ndx -f md.tpr -o`
3. We can now calculate the hydrogen bonds between the two protein chains using `gmx hbond`, with the `-num` flag specifying GROMACS to output the number of hydrogen bonds as a function of time

<div align="center">
  <img src="https://github.com/c-vandenberg/gromacs-tutorials/assets/60201356/7be8a026-3f86-43fb-815c-d28ef28e4cea" alt="chain-1-chain-2-h-bonds-vs-time" width="">
</div>
<br>
	
**Report Methods**
1. Once we have run the simulation, it is good practice to report what type of simulation we have performed, as well as the basic system information
2. This can be achieved with the `gmx report-methods` GROMACS tool:
	* `gmx report-methods -s md.tpr`


 ## References
**[1]** Jensen, F. (2017) ‘15. 2 Time-Dependent Methods’, in *Introduction to Computational Chemistry*. 3rd edn. Newark: John Wiley & Sons, Incorporated, pp. 474–478. 
