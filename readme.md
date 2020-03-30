# Molecular Dynamics of Liquid Argon

ICTP - CMP Diploma 2019-2020

Numerical methods II - Final Project

Jesus Andres Espinoza-Valverde

## Folder Organization

Our package contains the following directories:

- **sources**: This Folder contains the source code.
- **input**: Here we store the input files.
- **results**: Within this folder we store the raw output of the simulation.
- **executable**: This folder is meant to contain our executable binary.
- **scripts**: Here we have helper bash scripts.
- **data_analysis**: This folder contains a python script that we used to analyze the simulation data.
- **documents**: Here we find references and important documents.

## Source Code Files

The source code is composed by the following files: 

- **md.f90**: The main program, it contains the main molecular dynamics loop.
- **physics.f90**: This file contains physics-related modules and routines.
- **io.f90**: This file contains input/output and file management modules and routines.
- **utils.f90**: Here we find miscellaneous modules and routines.  

## Simulation Modules

#### MODULE Physical Constants (physconst)

**Description**: This module holds the needed physical constants used in the simulation.

#### MODULE Molecular Dynamics System (mdsys)

**Description**: This module holds the complete system information.

#### MODULE High precision kinds

**Description**: This module contains redefinitions for precision of floating point number and length file size names.

#### MODULE Utilities (utils)

**Description**: This module contains helper and miscellaneous routines.

**Routines**:

- **pbc routine (pbc)**: applies minimum image convention.
- **Box-Muller Method (box_muller_method)**: Generates Gaussian random numbers.

#### MODULE Physics

**Description**: This module contains all physically-related routines.

**Routines:**

- **Get kinetic energy routine (getekin)**: Computes total kinetic energy of the system.
- **Get temperature routine (gettemp)**: Computes the temperature of the system.
- **Force routine (force)**: Computes the temperature of the system.
- **Velocity Verlet routine (velverlet)**: Updates velocities and p via Verlet Algorithm.
- **Poors-man Thermostat routine (thermostat)**: Rescales temperature of the system.
- **Maxwell Bolzmann velocity initializator (MaxBoltz_Dist_vel_init)**: Initialize velocities according to MB velocity distribution at a given temperature.
- **FCC lattice positions iniatializator (fcc_lattice_positions_init):** Initialize positions in a FCC lattice.
- **Force to Zero routine (force_to_zero)**: Resets Force values.
-  **Get distances routine (get_distances):** Computes distances between the particles in the system.

#### MODULE Input/Output

**Description**: This module is in charge of the reading and writing files with simulation data.

**Routines**:

- **ioopen routine (ioopen):** Opens all the needed files.
- **ioclose routine (ioclose)**: Closes all the needed files.
- **output routine (output):** Writes simulation data into files.
- **output temps routine (output_temps):** writes temperature data only.

### Animation

A cool video of this simulation using VMD software can be found here:

https://youtu.be/M-hA6Vg5Mno