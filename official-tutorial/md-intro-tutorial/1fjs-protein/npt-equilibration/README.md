# 4) NPT Ensemble Equilibration

## Introduction

Using NVT equilibration we have stabilised the temperature of the system. However, prior to the MD simulation, we must also stabilise the pressure (and thus the density) of the system.

We therefore must first run an **NPT ensemble equilibration**. The NPT ensemble is a statistical ensemble in which the **N**umber of particles, **P**ressure & **T**emperature of the system is kept constant. It is commonly used in MD simulations to **model systems at a fixed temperature and pressure.**

In the context of equilibration, NPT equilibration equilibrates the system at constant volume & pressure, using a **barostat** (in our case the Berendsen barostat) to control system pressure.

NPT equilibration allows the system volume to adjust to achieve the target pressure, thus allowing the system to react a realistic, stable density before any MD production runs.

## Simulation Commands

### 4.1) Equilibration Run - Pressure (NPT Ensemble Equilibration)
1. Equilibration of pressure is conducted under an NPT ensemble. This ensemble is also called the 'isothermal-isobaric' ensemble, and most closely resembles experimental conditions
2. The NPT equilibration simulation parameter file (`npt-charmm.mdp`) uses the Berendsen barostat to control pressure, and specifies a 100 ps NPT equilibration
3. Navigate to `1fjs-protein/data/input/topology/protein`
    * `gmx grompp -f ../../../../npt-equilibration/data/input/npt-charmm.mdp -c ../../../../energy-minimization/data/processed/em.gro -r ../../../../energy-minimization/data/processed/em.gro -p topol.top -o npt.tpr`
    * `mv npt.tpr ../../../../npt-equilibration/data/processed`
    * `mv mdout.mdp ../../../../npt-equilibration/data/processed`
4. Navigate to `1fjs-protein/nvt-equilibration/data/processed`
    * `gmx mdrun -ntmpi 1 -v -deffnm npt`

### 4.2) NPT Ensemble Equilibration Data Analysis
1. Navigate to `1fjs-protein/nvt-equilibration/data-analysis`
2. Extract pressure & time data
	* `echo "Pressure" | gmx energy -f ../data/processed/npt.edr -o pressure.xvg -xvg none`

<br>
<div align="center">
  <img src="https://github.com/c-vandenberg/gromacs-tutorials/assets/60201356/e152716f-0775-484d-ad28-2087dbb1887e" alt="pressure-vs-time" width="">
</div>
<br>

3. Extract density & time data
	* `echo "Density" | gmx energy -f ../data/processed/npt.edr -o density.xvg -xvg none`

<br>
<div align="center">
  <img src="https://github.com/c-vandenberg/gromacs-tutorials/assets/60201356/296cd152-e94e-45eb-a63d-f003702d6fcc" alt="density-vs-time" width="">
</div>
<br>

4. The average system density is close to the experimental value of 100 kg/m<sup>3</sup> and the expected density of the TIP3P model of 1001 kg/m<sup>3</sup>. The density values are also very stable over time, indicating that the system is now well-equilibrated with respect to pressure & density
