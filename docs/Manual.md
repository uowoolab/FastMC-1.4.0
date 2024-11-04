# FastMC 1.4.0 Manual

## Table of Contents
1. [Introduction](#introduction)
2. [GCMC Moves](#gcmc-moves)
    1. [Main Moves](#main-moves)
    2. [Sub Moves](#sub-moves)
3. [Input Files](#input-files)
    1. [FIELD File](#field-file)
    2. [CONFIG File](#config-file)
    3. [CONTROL File](#control-file)
4. [Output Files](#output-files)
    1. [OUTPUT File](#output-file)
    2. [REVCON File](#revcon-file)
    3. [REVIVE File](#revive-file)
    4. [History File (his.xyz)](#history-file)
    5. [Convergence File (numguests.out)](#convergence-file)
5. [How to Restart Simulation](#how-to-restart-a-simultion)

---

## Introduction

**FastMC** is a Grand-Canonical Monte Carlo (GCMC) simulation code written in Fortran. Optimized for efficiency and speed, it provides users with accurate adsorption information for periodic crystalline structures such as MOFs, COFs, and zeolites.

The code is based on and utilizes functions from the molecular dynamics code **DL POLY Classics 1.10**. FastMC leverages DL POLY’s energy subroutines to calculate configurations. Due to this structural basis, FastMC's input and output files share similar characteristics with DL POLY's file formats.

## GCMC Moves

This repository contains the **1.4.0 stable release** of FastMC. In this version, FastMC utilizes a total of three main moves and two sub-moves:

### Main Moves

- **Insertion**: Inserts the guest molecule at a random location in the framework.
- **Deletion**: Deletes a guest molecule at random.
- **Displacement**: Uses two sub-moves to translate and then rotate a guest molecule at random.

### Sub Moves

- **Rotation**: Rotates the guest molecule by a random degree around its center of mass.
- **Translation**: Translates the guest molecule from its center of mass to a random location within the framework.

### Note

In the **OUTPUT** file, you may encounter the terms *jumping*, *flexing*, *swapping*, and *multi-jump* moves. Although these can be specified in the **CONTROL** file, their subroutines are not implemented in FastMC 1.4.0. Attempting to use these parameters will cause your simulation to fail.

## Input Files

### FIELD File

Similarly to the DL POLY FIELD file, the FIELD file specifies force field information for the guest and framework molecules. Below is an example FIELD file for the adsorption of N2 (using the N2-NIMF force field) in the VIXGAM framework (from the CoRE database).

~~~
Absorbtion of N2 in VIXGAM
UNITS   kcal

molecular types 2
&guest N2-NIMF: Nitrogen In MoFs (NIMF): DOI
NUMMOLS 0
ATOMS 3
Nx        14.006700    -0.482000   -0.550000     0.000000     0.000000
                                        ...
COM        0.000000     0.964000    0.000000     0.000000     0.000000
rigid 1
3 1 2 3
finish

Framework
NUMMOLS 36
ATOMS 28
Zr        91.224000     1.65522000000000      1      1
                        ...
O         15.999400    -0.69296200000000      1      1
finish

VDW 15
Zr     C      lj 0.085118 3.107050
                ...
Nx     Nx     lj 0.079420 2.454910
close
~~~

#### General Information

The first section of the FIELD file specifies basic simulation parameters.

- **Header 1**  
  FIELD file header: the title of the simulation. Any string is accepted as a description.

- **Header 2**  
  UNITS n
    
    Defines the units used in the simulation. Choose from the following options:
  
  - **eV**
  - **kcal** (recommended for GCMC)
  - **kJ**
  - **K**
  - **internal**

  **Note**: The recommended unit for GCMC simulations is **kcal** due to its compatibility with FastMC’s energy and force field calculations.

#### Molecular Information

This section specifies the molecular details for each type of molecule in the system. Each guest molecule and framework molecule has its own molecular section, with details depending on the specifics of your system.

##### Guest Molecule

- **molecular types n**  
  `n` represents the number of distinct molecule types, including guest molecules and the framework.

- **&guest**  
  Provides information about the guest molecule, such as its identity, name, and DOI reference.

- **NUMMOL n**  
  `n` is the number of times this guest molecule appears in the CONFIG file. For an empty framework, set `n = 0`. If using a REVCON file, specify the number of molecules in the simulated system.

- **ATOMS n**  
  Specifies the number of atoms in the guest molecule. After this keyword, include the guest-specific information as a table (each row being an atomic site):

    ```
    AS_label  AS_mass       AS_charge   AS_coords_x  AS_coords_y  AS_coords_z
    Ex:
    Cx        12.010700     0.651200    0.000000     0.000000     0.000000
    ```

    **Note**: Entries are order-sensitive. For neutral guest molecule sites, specify a charge value of zero.

- **rigid n**

  `n` represents the number of rigid units in the molecule. The entries following this keyword define the connectivity for each rigid unit, in the order specified in the **ATOMS** section.
  
  - **m int int int...**
  
    - `m` is the number of sites in the rigid unit.
    - The sequence of integers (`int int int...`) following `m` specifies the index of each site, corresponding to the order in the **ATOMS** section.

    **Example for CO2**:
    ```
    rigid 1
    3 1 2 3
    ```
    In this example, `rigid 1` defines one rigid unit with three sites (1, 2, and 3) forming the rigid connectivity.


- **finish**  
  This keyword terminates the guest molecule section. It is required for the section to complete successfully.

#### Framework Molecule

Similar to the guest molecule section, the framework section defines the structure and composition of the framework. This section begins with the keyword **Framework**.

- **NUMMOL n**  
  `n` represents the number of times the unit cell is replicated. For example, in a simulation cell of size 2x2x2, the unit cell is replicated 8 times.

- **ATOMS n**  
  `n` is the number of atoms in the unit cell. Multiplying `NUMMOL` and `ATOMS` gives the total number of atoms in the supercell, which should match the atom count in the CONFIG file header. The following information provides site details for atoms in the unit cell (atomic site = AS):

    ~~~
        AS_label    AS_mass      AS_charge    Repeat    AS_Frozen
    Ex:
        Zr          91.224000    1.65522      1         1
    ~~~

    **Note**: `Repeat` and `AS_Frozen` should be set to **1**, as GCMC simulations using FastMC assume the framework is periodic and rigid, even if the supercell is set to 1x1x1.

- **finish**  
  This keyword must be included to properly close the framework section. Omitting `finish` will result in an error, as it is required to complete the FIELD file's framework definition.


#### Lennard-Jones Parameter Section

In this section of the FIELD file, the pairwise mixed Lennard-Jones (LJ) parameters are defined for each combination of framework-framework, guest-framework, and guest-guest atoms. You may apply any mixing rules to specify the Lennard-Jones parameters according to your requirements.

- **VDW n**  
  `n` represents the number of pairwise LJ parameters in the system. After this keyword, provide a table of mixed parameters. Each entry should include the following fields (atomic_site = AS):

    ~~~
        AS_label_1    AS_label_2    Type_of_interaction    Epsilon      Sigma   
    Ex:
        Zn            Zn            lj                     0.124000     2.461600
    ~~~

- **close**  
This keyword terminates the FIELD file subroutine parser. It is required to complete the parameter section. If `close` is omitted, an error will be raised.

### CONTROL File

The CONTROL file contains all simulation variables for the GCMC simulation. The format remains similar to the DL POLY CONTROL file, but the keywords are tailored specifically for GCMC applications. Below is an example of keywords used for a GCMC run.

~~~
GCMC Run
temperature  298.0
steps    10000
equilibration    -10000
max guest atoms      3000
max framework atoms      15000
cutoff          4.000000 angstrom
delr            1.0 angstrom
ewald precision  1d-6
numguests 1000
history 100

&guest 1
  pressure (bar) 0.994531
  probability 2
  1  1
  2  2  3
&end

move insert 0.333000
move delete 0.333000
move displace 0.333000
move swap 0.000000
move switch 0.000000
move jump 0.000000

grid factors 2 2 2

grid spacing 0.150000

finish
~~~

The CONTROL file begins with an identifier line titled `(GCMC Run)`. Descriptions for all possible keywords are found below.

- **temperature n**  
  `n` represents the temperature of the simulation.

- **steps n**  
  `n` is the number of production steps in the GCMC simulation. During this phase, after reaching equilibrium, each step records essential data on properties like adsorption and energy for analysis.

- **equilibration n**  
  `n` is the number of steps the GCMC simulation takes to equilibrate the system. During equilibration, properties like particle density and energy gradually stabilize, moving away from initial conditions to reflect a steady-state behavior.

  **Note**: In some examples, you may encounter a negative sign before `n` for production and equilibration steps. A negative `-n` specifies cycles instead of individual steps. A cycle refers to a sequence of moves (such as insertion, deletion, and displacement) that together represent a complete sampling round of the system.

- **max guest atoms n**  
  `n` is the maximum number of guest atoms in the system.

- **max framework atoms n**  
  `n` represents the maximum number of framework atoms.

- **cutoff n angstrom**  
  `n` specifies the cutoff distance in angstroms, defining the maximum interaction distance for non-bonded forces.

- **delr n angstrom**  
  `n` sets the width of the Verlet neighbor list shell in angstroms.

- **ewald precision n**  
  `n` represents the precision for Ewald summation, determining the accuracy for long-range electrostatic interactions.

- **numguests n**  
  `n` is the maximum number of guest molecules in the simulation cell at any time.

- **history n**  
  `n` sets the frequency of configuration saves in the history file. Set to 0 to disable, as frequent saves can create large files.

- **&guest n**  
  `n` is the index for guest information from the FIELD file (more details provided later). Must be terminated with `&end`.

- **move insert n**  
  - **move delete n**
  - **move displace n**
  - **move rotate n**
  - **move translate n**  

    Each `n` represents the probability of using the respective move during GCMC.

- **accept displace n**  
  - **accept rotate n**
  - **accept translate n**  

    Each `n` represents the probability of accepting a move during GCMC.

- **restart**  
  Restarts the job from the endpoint of a previous run.

- **averaging window n**  
  `n` specifies the number of steps over which properties are averaged to smooth fluctuations in recorded data.

#### Probability Section

The probability section of the CONTROL input file, located under the **&guest** keyword, specifies the number of probability plots and the atoms included in each plot.

- **probability n**  
  `n` is the number of probability plots to be generated.

  Each probability plot entry follows this format:

  ~~~
  n  index
  n  index  index
  ~~~

- `n` specifies the number of sites included in each probability plot.
- `index` refers to the order of the sites as represented in the **FIELD** file.

**Example**:

For example, here is the FIELD file entry for a CO2 guest molecule:
~~~

    &guest CO2: Carbon Dioxide for Adsorption in Zeolites: http://dx.doi.org/10.1021/jp810871f
    NUMMOLS 0
    ATOMS 3
    Cx        12.010700     0.651200    0.000000     0.000000     0.000000
    Ox        15.999400    -0.325600    1.149000     0.000000     0.000000
    Ox        15.999400    -0.325600   -1.149000     0.000000     0.000000
    rigid 1
    3 1 2 3
    finish
~~~

Below are the specifications for the probability plot of the same CO2 guest molecule. The first line, `probability 2`, indicates that two probability plots will be generated.

- The **first probability plot** includes only 1 site, corresponding to the **Cx** atom specified in the FIELD file.
- The **second probability plot** includes 2 sites, corresponding to the two **Ox** atoms from the FIELD file.

~~~
    probability 2
    1  1
    2  2  3
~~~

#### Guest Specific Moves
It is worth noting that users can specify guest-specific moves and acceptance criteria, which becomes particularly useful in multi-component simulations. Here is an example of **&guest** inputs for a binary simulation. 

**Note**: The acceptance rate for the rotation move has not been implemented for guest-specific configurations.

**Example:**
~~~
&guest 1
  pressure (bar) 0.50

  move insert 0.25
  move delete 0.25
  move displace 0.50

  accep displace  0.25
  accep translate 0.75

  delr = 0.5 angstrom

  probability 2
  1  1
  2  2  3
&end
&guest 2
  pressure (bar) 0.50

  move insert 0.75000
  move delete 0.10000
  move displace 0.15

  accept displace  0.3
  accept translate 0.7

  delr = 1.2 angstrom

  probability 1
  2  1  2
&end
~~~

In this example:
- **&guest 1** and **&guest 2** sections specify individual guest properties for a binary simulation.
- Each **&guest** entry includes:
  - **pressure**: Operating pressure for the guest molecule.
  - **move insert/delete/displace**: Probabilities for inserting, deleting, and displacing the guest.
  - **accept displace/translate**: Acceptance probabilities for displacing and translating the guest molecule.
  - **delr**: Specifies the interaction range in angstroms.
  - **probability**: Defines the number of probability plots and the sites included in each plot.

### CONFIG File

The CONFIG file is similar to the DL POLY version, containing the lattice vectors, atom labels, and coordinates of the supercell framework structure. Below is an example of the start of a CONFIG file for the unit cell of calf-20.

~~~
calf-20
         0         3        44
      8.913800000000      0.000000000000      0.000000000000
      0.000000000000      9.693500000000      0.000000000000
     -4.141708950172      0.000000000000      8.531407617273
Zn             1
     -0.241297908345      0.559411885000      3.726433533149
Zn             2
      7.084243433260      5.406161885000      0.539270275488
~~~


#### CONFIG File Variables

- **header**  
  The title line, specifying the name of the periodic structure.

- **System information**

  - **levcfg n**  
    `n` specifies what data is included in the CONFIG file. For simplicity, only coordinates are included in these simulations.

    - **0**: Coordinates included
    - **1**: Coordinates and velocities included
    - **2**: Coordinates, velocities, and forces included

  - **imcon n**  
    `n` defines the boundary conditions. Since a custom supercell is generated, this option is set to **0** (no periodic boundaries).

    Possible options:
    - **0**: No periodic boundaries
    - **1**: Cubic boundary conditions
    - **2**: Orthorhombic boundary conditions
    - **3**: Parallelepiped boundary conditions
    - **4**: Truncated octahedral boundary conditions
    - **5**: Rhombic dodecahedral boundary conditions
    - **6**: X-Y parallelogram boundary conditions (no periodicity in the z direction)
    - **7**: Hexagonal prism boundary conditions

  - **natms n**  
    `n` represents the total number of atoms in the system.

  - **engcfg n (optional)**  
    `n` specifies the configuration energy of the system.

- **Lattice parameters**  
  Defined as a matrix in the form:  
  ~~~
  [ a 0 0 ] [ 0 b 0 ] [ 0 0 c ]
  ~~~

Following the lattice parameters, each atomic site is specified with:
- **Label and index**: The atom label and its index on the same line.
- **Cartesian coordinates**: The x, y, and z coordinates for each site in the supercell.

This structure is repeated for each site in the supercell under P1 symmetry.

## Output Files

The main output file is saved in the directory where the input files were created. Additional output files are located in the **branch1** directory.

### OUTPUT File

The OUTPUT file contains all final simulation data in a human-readable format. It is divided into two sections:
- **Simulation initiation**: Contains the initial parameters and setup details.
- **GCMC simulation results**: Reports the outcomes, including the average number of guest molecules.

FastMC will provide the average guest molecule count, which can be converted to mmol/g using the following equation:

$$
\text{Uptake (mmol/g)} = \frac{\text{N}_{\text{guest, sc}}}{\text{uc count}} \times \frac{1000}{\text{MMOF}}
$$

Where:
- **Nguest, sc**: Average number of guest molecules in the simulation supercell (from FastMC OUTPUT file).
- **uc count**: Number of unit cells in the supercell (e.g., a 2 x 2 x 2 fold factor gives 8 repeated unit cells).
- **MMOF**: Molecular mass of the MOF unit cell (calculated MOF's unit cell composition).
- The factor `1000` converts mol to mmol.

### REVCON File

The REVCON file serves as the restart configuration file, containing the configuration from the simulation's last step. The REVCON file format matches that of the CONFIG file. To restart a simulation, move the REVCON file to the main directory and rename it to **CONFIG**.

### REVIVE File

The REVIVE file is a binary file containing accumulated statistical data from the simulation. To continue a simulation, move the REVIVE file to the main directory and rename it to **REVOLD**.

### History File

The **his.xyz** file is generated only if the `history` keyword is specified in the CONTROL file. This file records the guest configuration at each specified interval in XYZ format, capturing the system’s changes over time. It can be visualized using any XYZ-compatible application. Guest specific hsitory files will be generated if simulation contain more than one guest molecule.


### Convergence File

The **numguests.out** file is helpful for determining whether the simulation has converged. Multiple `numguests.out` files may be generated, each corresponding to a different guest molecule in the system. The columns in this file are as follows:

~~~
Simulation_step   Number_of_guests_at_step   Configuration_energy   Variation_in_energy_with_move   is_move_ins   is_move_del   is_move_dis   is_move_jmp   is_move_flx   is_move_swap
~~~

Plotting **Simulation_step** against **Number_of_guests_at_step** provides an indication of convergence. In a converged simulation, you should see a plateau over the last 10,000 steps (or equal to your windowed average), indicating stability in the guest molecule count.

## How To Restart a Simultion

To restart a simulation, you need a completed FASTMC directory. Below is an example of a directory structure from a unary CO2 absorption simulation within the calf-20 MOF framework:

~~~
[omarchand@wooki co2]$ tree .
|-- branch01
| |-- numguests.out
| |-- REVCON
| |-- REVIVE
| `-- runningstats.out
|-- CONFIG
|-- CONTROL
|-- FIELD
|-- jobcontrol.in
|-- OUTPUT
|-- prob_guest01_prob_01.cube
`-- prob_guest01_prob_02.cube
~~~

The key files required to restart a simulation are the REVCON and REVIVE files found in the branch01 directory.

To restart a simulation, move the **REVCON** and **REVIVE** files from the `branch01` directory to the main FASTMC directory where the OUTPUT file was created:

1. Change to the `branch01` directory:
```bash
cd branch01
```

2. Move REVCON and REVIVE to the main directory:
```bash
mv REVCON REVIVE ../
cd ../
```

Next, rename the REVCON file to CONFIG and the REVIVE file to REVOLD. You can do this manually or use the copy executable:

```bash
mv REVCON CONFIG
mv REVIVE REVOLD
```

Or, if you have DL POLY classic 1.10 installed on your machine
```bash
module load dl_poly/classic_1.10
copy
```

The copy executable renames the files and saves the original CONFIG file. The directory should now look like this:

~~~
[omarchand@wooki co2]$ tree .
|-- branch01
|   |-- numguests.out
|   `-- runningstats.out
|-- CONFIG
|-- CONFIG.OLD
|-- CONTROL
|-- FIELD
|-- OUTPUT
|-- prob_guest01_prob_01.cube
|-- prob_guest01_prob_02.cube
`-- REVOLD
~~~

### Updating the FIELD File

To update the FIELD file, first determine the number of guest molecules in the new CONFIG file using the grep command:

```bash
grep -o "Cx" CONFIG | wc -l
```

In this example, each CO2 molecule contains the Cx site once, so the output from grep directly reflects the number of CO2 molecules in the system. However, if your guest molecule contains multiple instances of a site, divide the result by the stoichiometric factor of the counted sites.

Update the FIELD file (NUMMOL input) to reflect the correct number of guest molecules:

~~~
nano FIELD

calf-20
UNITS kcal
molecular types 2
&guest CO2: Carbon Dioxide for Adsorption in Zeolites:
http://dx.doi.org/10.1021/jp810871f
NUMMOLS 90
ATOMS 3
Cx        12.010700   0.651200   0.000000   0.000000   0.000000
Ox        15.999400  -0.325600   1.149000   0.000000   0.000000
Ox        15.999400  -0.325600  -1.149000   0.000000   0.000000
rigid 1
3 1 2 3
finish
~~~

### Updating the CONTROL File

Add the restart keyword in the CONTROL file:

~~~
GCMC Run
restart
temperature 303.150000
steps 100000000
equilibration -1500000
max guest atoms 3000
max framework atoms 15000
...
~~~

After correctly setting up the restart job, you can adjust the CONTROL file as needed, such as reducing the number of steps for the new simulation. Once your modifications are complete, clean up unnecessary files from the previous simulation:

```bash
find . -maxdepth 1 ! -name 'CONFIG*' ! -name 'FIELD' ! -name 'CONTROL' ! -name 'REVOLD' ! -name '.' -exec rm -rf {} +
```

Finally, submit the restart job.