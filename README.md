# GCMC Code for Simulating Adsorption in Crystalline Frameworks

<center><pre>
███████╗ █████╗ ███████╗████████╗███╗   ███╗ ██████╗
██╔════╝██╔══██╗██╔════╝╚══██╔══╝████╗ ████║██╔════╝
█████╗  ███████║███████╗   ██║   ██╔████╔██║██║     
██╔══╝  ██╔══██║╚════██║   ██║   ██║╚██╔╝██║██║     
██║     ██║  ██║███████║   ██║   ██║ ╚═╝ ██║╚██████╗
╚═╝     ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚═╝     ╚═╝ ╚═════╝
</pre></center>

## Overview
FastMC is a Grand Canonical Monte Carlo (GCMC) program designed for simulating gas adsorption, selectivity, and diffusion in periodic crystalline structures. This includes metal-organic frameworks (MOFs), covalent organic frameworks (COFs), zeolitic imidazolate frameworks (ZIFs), and zeolites. FastMC is built on the subroutines from DL POLY Classic 1.10 and is written in Fortran.

## Dependencies
Before building FastMC, ensure that the following dependencies are installed on your machine:
- **GNU Compiler** 9.3.0
- **Intel Compilers** 2022.2
- **OpenMPI** 4.0.4 (for parallel builds)

## Installation
FastMC supports several build configurations, including serial, gfortran, and parallel modes.

### Serial Build
```bash
export PATH=/path/to/gnu9/bin:$PATH
export PATH=/path/to/intel_compilers/bin:$PATH

make clean
make -f Makefile serial
```

### Gfortran Build
```bash
export PATH=/path/to/gnu9/bin:$PATH

make clean
make -f Makefile gfortran
```

### Parallel Build
```bash
export PATH=/path/to/gnu9/bin:$PATH
export PATH=/path/to/openmpi4/bin:$PATH
export PATH=/path/to/intel_compilers/bin:$PATH

make clean
make -f Makefile parallel
```

## Running the Code
```bash
export PATH=/path/to/fastmc-1.4.0/bin:$PATH
cd /path/to/simulation

gcmc.x
```