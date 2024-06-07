# 3) NVT Ensemble Equilibration

## Introduction

Energy minimization has given us a reasonable starting structure in terms of protein geometry and solvent orientation. 

However, before the MD simulation we must equilibrate the solvent and ions around the protein. If we were to attempt an MD simulation without this, the system may collapse. This is because the solvent is mostly optimised within itself, not necessarily with the solute, and the ions have been randomly placed.

We therefore must first run an **NVT ensemble equilibration**. The NVT ensemble is a statistical ensemble in which the **N**umber of particles, **V**olume & **T**emperature of the system is kept constant. It is commonly used in MD simulations to **model systems at a fixed temperature and volume.**

In the context of equilibration, NVT equilibration equilibrates the system at constant volume & temperature, using a **thermostat** (in our case a **modified Berendsen thermostat**) to control system temperature.

NVT equilibration ensures that the kinetic energy distribution of the system corresponds to the desired temperature, allowing the system to reach a stable temperature.

## Simulation Commands

### 3.1) Equilibration Run - Temperature (NVT Ensemble Equilibration)
1. Within the `topol.top` file, there are `porse.itp` position restraint files that apply a position restraining force on the heavy atoms of the protein (i.e. anything that is not H). Movement is permitted, but only after overcoming a substantial energy penalty
2. These position restraints allow us to relax the solvent and ions around the protein, without the added variable of structural changes in the protein
3. To use position restraints, we have added the following line to the NVT equilibration simulation parameter file (`nvt-charmm.mdp`):
    * `define                  = -DPOSRES  ; position restrain the protein`
4. Equilibration of temperature is conducted under an NVT ensemble. This ensemble is also called the 'isothermal-isochoric' ensemble
5. The NVT equilibration simulation parameter file (`nvt-charmm.mdp`) uses a modified Berendsen thermostat to control temperature, and specifies a 100 ps NVT equilibration
6. Navigate to `1fjs-protein/protein/topology/protein` and run:
    * `gmx grompp -f ../../../../nvt-equilibration/data/input/nvt-charmm.mdp -c ../../../../energy-minimization/data/processed/em.gro -r ../../../../energy-minimization/data/processed/em.gro -p topol.top -o nvt.tpr`
    * `mv nvt.tpr ../../../../nvt-equilibration/data/processed`
    * `mv mdout.mdp ../../../../nvt-equilibration/data/processed`
7. Navigate to `1fjs-protein/nvt-equilibration/data/processed` and run:
    * `gmx mdrun -ntmpi 1 -v -deffnm nvt`
8. The `-ntmpi 1` flag specifies the number of Message Passing Interface (MPI) ranks to use. Here we are using only one, so we are not taking advantage of parallel execution. If we were to specify 4 MPI ranks for example (`-ntmpi 4`), if on a single computer (as is the case here), the simulation workload could be distributed between separate CPU cores

### 3.2) NVT Ensemble Equilibration Data Analysis
1. Navigate to `1fjs-protein/nvt-equilibration/data-analysis`
2. Analyse temperature & time data via Python & JupyterLab:
	* `echo "Temperature" | gmx energy -f ../data/processed/nvt.edr -o temperature.xvg -xvg none -b 20`

<br>
<div align="center">
  <img src="https://github.com/c-vandenberg/gromacs-tutorials/assets/60201356/b7bf02ad-1f85-4633-9538-912eb6b3f781" alt="temp-vs-time" width="">
</div>
