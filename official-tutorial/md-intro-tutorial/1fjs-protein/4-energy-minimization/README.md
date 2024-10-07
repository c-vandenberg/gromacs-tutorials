# 2 Energy Minimization

## 2.1 Introduction

We have now defined the protein topology, solvent topology and simulation box of our system. We have also ensured our solvated system is electronically neutral. 

However, before we can begin molecular dynamics, we must ensure that the system has no steric clashes or inappropriate geometry. This is achieved via energy minimization. We will use 'steepest descent', which is one of the simplest energy minimization algorithms.

**Energy Minimization by steepest descent works** as follows **<sup>1</sup>**:
1. **Initialization** 
    * **Input System Configuration**: Start with initial system configuration/topology
    * **Potential Energy Calculation**: Calculate the initial potential energy, *U* of the system using a defined force field. In our case it is the CHARMM27 all-atom force field
2. **Calculate Forces**
    * **Gradient Calculation**: Calculate the **gradient of the potential energy surface (PES)** with respect to the positions of the atoms in the current structure. The gradient points towards the direction of **most increasing energy**
      * The gradient is calculated by taking the first-order partial derivative of the potential energy, *U* of each atom with respect to its coordinates
      * Therefore, the gradient of the PES is a vector of these first-order partial derivatives and represents the forces acting on each atom
	<br>

 	$$\nabla_{\mathbf{}} U(\mathbf{r_i}) = \left( \frac{\partial U}{\partial x_i}, \frac{\partial U}{\partial y_i}, \frac{\partial U}{\partial z_i} \right)$$

	where:
	* $`i`$ is the current iteration

3. **Atom Displacement**
    * **Displacement Size Determination**: A displacement size scalar value, $`\gamma`$ is chosen to determine how far each atom will move in the direction of **most decreasing energy (the negative gradient)** in a single iteration
    * **Update Positions**: Update geometric coordinates of each atom in the molecule to get new positions, $`r_{i + 1}`$
      
	<br>
 	$$r_{i+1} = r_i - \gamma_i \nabla U(r_i)$$

  	where:
   	* $`r_{i + 1}`$ are the new x, y and z geometric coordinates
   	* $`r_{i}`$ are the old x, y and z geometric coordinates
   	* $`- \gamma_i \nabla U(r_i)`$ is negative gradient, multiplied by our displacement size scalar value. This represents how far we are displacing the atoms in the direction of most decreasing energy

5. **Recalculate Potential Energy**
    * **Potential Energy Calculation**: Calculate the new potential energy, *U* of the system using a defined force field
6. **Convergence Check**
    * For each of the convergence checks, the predfined threshold will be 0 with a certain decimal place tolerance
    	* **Energy Difference Check**: Calculate the change in system potential energy, *U* between the current and previous iteration. If it is below the predefined threshold, this indicates that the system has reached a local energy mimimum
       
		<br>
		$$|U(r_i) - U(r_{i+1})| \approx 0$$
		<br>
  
    	* **Gradient Magnitude Check**: Evaluate the magnitude of the current PES gradient at position *r*. If it is below the predefined threshold, this indicates that the system has reached a local energy minimum
       
		<br>
		$$||\nabla U(r_i)|| \approx 0$$
		<br>
  
		* **Atom Displacement Magnitude Check**: Evaluate the magnitude of atom displacement between the current iteration and the previous iteration. If it is below the predefined threshold, this indicates that the system has reached a local energy minimum

		<br>
  		$|| r_i - r_{i+1} || \approx 0$
		<br>
  
    * If all three of these convergence criteria are met, continue to step 8 and output final molecular structure
7. **Iteration Limit Check**: If number of iterations has reached a predefined limit, continue to step 8 and output final molecular structure
8. **Repeat**
    * **Iterative Process**: Else, repeat steps 2 to 6 until the convergence criteria are met, or the number of iterations has reached the predefined limit
9. **Output**
    * **Final Molecular Structure & Minimized Energy**: Once convergence criteria is met or iteration limit is reached, the final positions of the atoms represent a configuration with minimized potential energy

Steepest descent is a very fast approach to energy minimization and can be used to treat very large systems of millions of atoms. 

However, the main drawback of steepest descent is that it can only find the nearest local energy minimum. This is because it cannot overcome energy barriers within local energy minima as it only moves across the PES in the direction of a negative gradient. 

Therefore, steepest descent cannot climb out of the nearest local energy minimum, further explore the PES and find even lower energy minima. To do this, we would need an algorithm that was capable of **uphill steps**, such as **Monte Carlo Simulation**.

## 2.2 Simulation Commands

### 2.2.1 Energy minimisation
1. Before we can begin molecular dynamics, we must ensure that the solvated, electronically neutral system has no steric clashes or inappropriate geometry via energy minimization
2. Navigate to `1fjs-protein/energy-minimization/data/processed`
3. We again use `gmx grompp` to generate the atomic-level binary file (`em.tpr`) that GROMACS can use to run energy minimization:
	* `gmx grompp -f  ../input/energy-minimization/emin-charmm.mdp -c ../input/topology/solvent/1fjs_solv_ions.gro -p ../input/topology/protein/topol.top -o em.tpr`
4. The next GROMACS tool we use is `gmx mdrun` which is used for running simulations or computations from an atomic-level input binary file (`.tpr`) file
5. We specify GROMACS to use verbose terminal output (`-v`), and we specify that the input and all output files have the base name `em` (`-deffnm em`)
	* `gmx mdrun -v -deffnm em`
6. To ensure the energy minimization run was successful, in the end output, the potential energy (`epot`) must be negative & on the order of 1e5, and the maximum force `Fmax` < `emtol` (defined in `emin-charmm.mdp`)

### 2.2.2 Energy minimisation Data Analysis
1. Navigate to `1fjs-protein/energy-minimization/data-analysis`
2. The GROMACS tool `gmx energy` is used for extracting and analysing energy terms from an energy (`.edr`) file 
3. Pipe in `Potential` string into `gmx energy` command. Again this isn't needed as GROMACS will prompt you for it if you run it without it
4. Analyse potential energy & time data via Python & JupyterLab:
	* `printf "Potential\n0\n" | gmx energy -f ../data/processed/em.edr -o potential.xvg -xvg none`

<div align="center">
  <img src="https://github.com/c-vandenberg/gromacs-tutorials/assets/60201356/446530fd-7aac-4d4a-9087-8578fa5b4c78" alt="pot-energy-vs-time" width="">
</div>

## References
**[1]** Jensen, F. (2017) ‘13. 2 Optimizing General Functions: Finding Minima’, in *Introduction to Computational Chemistry*. 3rd edn. Newark: John Wiley & Sons, Incorporated, pp. 407–408. 
    
