# 5) Molecular Dynamics Simulation

## Introduction

After the NVT & NPT equilibration runs, the system is now well-equilibrated at the desired temperature (NVT equilibration) and pressure/density (NPT equilibration). We are now ready to release the position restraints and run the production molecular dynamics simulation for data collection.

A molecular dynamics (MD) simulation predicts the motion of atoms/molecules at a given temperature for a given simulation time length. 
* The basic idea is that by heating the system to a certain temperature, molecules in the system have sufficient energy to cross potential energy barriers and move across the PES to adopt new energy minimum structures/conformations. 
* This is similar to **Monte Carlo Simulation** energy minimization, with the main difference being we are generating structures/conformations sequentially in time (i.e. dynamics) in MD.

At the heart of MD simulations are **Newton's equations of motions**, which describe how the positions & velocities of particles evolve over time. MD simulations make use of **integration algorithms** such as the **Velocity Verlet algorithm** to integrate and solve Newton's equations of motion.

The general steps of the Velocity Verlet algorithm are:
1. **Initial Position & Velocity**
    * Each particle/atom in the molecule is given an initial position & velocity
    * Using the force field (e.g. CHARMM 27), the forces acting on each atom at this initial position is calculated (sum of the bonding & non-bonding interactions)
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
    * ![acceleration](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Ba%7D%28t%29%20%3D%20%5Cfrac%7B%5Cmathbf%7BF%7D%28t%29%7D%7Bm%7D) is the acceleration at time ![time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t), ![force_at_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%29) is the force acting on the particle/atom at time ![time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t) and ![mass](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20m) is the mass of the particle

<br>

3. **Velocity Update (Half-step)**
    * The Velocity Verlet algorithm calculates the velocities of each particle/atom in the molecule at the half-time step:

    <br>
    <div align="center">
      <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%5Cleft%28t%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7B2%7D%5Cright%29%20%3D%20%5Cmathbf%7Bv%7D%28t%29%20%2B%20%5Cfrac%7B1%7D%7B2%7D%20%5Cmathbf%7Ba%7D%28t%29%20%5CDelta%20t", alt="velocity-verlet-algorithm-half-step-velocity-update-equation"/>
    </div>


    where:
    * ![half_step_velocity](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%5Cleft%28t%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7B2%7D%5Cright%29) is the velocity at the half-step ![half_step_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t%20%2B%20%5Cfrac%7B%5CDelta%20t%7D%7B2%7D)
    * ![current_velocity](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bv%7D%28t%29) is the velocity at time *t*
    * ![half_step_vacceleration](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Ba%7D%28t%29%20%3D%20%5Cfrac%7B%5Cmathbf%7BF%7D%28t%29%7D%7Bm%7D) is the acceleration at time ![time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t), ![force_at_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%29) is the force acting on the particle/atom at time ![time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t) and ![mass](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20m) is the mass of the particle

    <br>

4. **Force Calculation**
   * The Velocity Verlet algorithm then calculates the forces **F** *(t+Δt)*: at the new position *r(t+Δt)*:

    <br>
    <div align="center">
      <img src="https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%20%2B%20%5CDelta%20t%29%20%3D%20-%5Cfrac%7B%5Cpartial%20U%28%5Cmathbf%7Br%7D%29%7D%7B%5Cpartial%20%5Cmathbf%7Br%7D%28t%20%2B%20%5CDelta%20t%29%7D", alt="velocity-verlet-algorithm-half-step-force-calculation-equation"/>
    </div>

    where:
    * ![force_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7BF%7D%28t%20%2B%20%5CDelta%20t%29) is the force at time ![current_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Bt%7D%28t%20%2B%20%5CDelta%20t%29)
    * ![negative_pot_energy_gradient](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20-%5Cfrac%7B%5Cpartial%20U%28%5Cmathbf%7Br%7D%29%7D%7B%5Cpartial%20%5Cmathbf%7Br%7D%28t%20%2B%20%5CDelta%20t%29%7D) is the negative gradiant of the potential energy ![pot_energy](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20U%28%5Cmathbf%7Br%7D%29) with respect to the position ![position](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Br%7D) at ![time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20t%20%2B%20%5CDelta%20t). The gradient is calculated by taking the first-order partial derivative with respect to the position ![position_at_time](https://latex.codecogs.com/svg.latex?%5Ccolor%7Bwhite%7D%20%5Cmathbf%7Br%7D%28t%20%2B%20%5CDelta%20t%29)

