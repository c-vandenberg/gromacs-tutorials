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
