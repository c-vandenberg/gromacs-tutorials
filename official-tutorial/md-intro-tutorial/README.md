# GROMACS Molecular Dynamics Introduction Tutorial

This GROMACS introductory project is based on the [official GROMACS introduction to molecular dynamics tutorial](https://tutorials.gromacs.org/md-intro-tutorial.html), adapted slightly in directory structure, molecular simulation trajectory analysis and Python/JupyterLab data analysis techniques.

The project involves setting up the topology, energy minimization, NVT & NPT equilibration and finally molecular dynamics simulation of the small protein **Coagulation Factor Xa** in ionized water.

## Introduction

Coagulation Factor X is a serine protease that sits at a pivotal point in the coagulation cascade. It has a role in the three major pathways of the cascade; the intrinsic, extrinsic and common pathways.<sup>1</sup> 

Factor X is a proenzyme and is activated into Factor Xa by hydrolysis, which can occur via the extrinsic or intrinsic pathways<sup>1</sup>.

The 3D structure of Factor Xa is available from the protein data bank (PDB), with the PDB code **1FJS**. The [crystal structure](https://www.rcsb.org/3d-view/1FJS/1) can be downloaded as a PDB file.

This order in which the simulation system was build, minimised, equilibrated and the MD simulation run with regards to the directories is as follows:
1. Generating Topologies, Defining Simulation Box & Solvating System (`protein` directory)
2. Energy Minimization (`1fjs-protein/energy-minimization` directory)
3. NVT Ensemble Equilibration (`1fjs-protein/nvt-equilibration` directory)
4. NPT Ensemble Equilibration (`1fjs-protein/npt-equilibration` directory)
5. Molecular Dynamics Simulation (`1fjs-protein/molecular-dynamics` directory)

## References
[1] Camire, R.M. (2021) 'Blood coagulation factor X: Molecular biology, inherited disease, and engineered therapeutics', *Journal of Thrombosis and Thrombolysis*, 52(2), pp. 383–390.<br>
