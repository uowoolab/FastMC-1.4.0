# FastMC-1.4.0

# Running

cd to/an/example/dir
./call/path/to/gcmc.x

# Build from source

Installation options;

### serial

module load gnu9/9.3.0
module load intel_compilers/2022.2

make clean
make -f Makefile serial

### gfortran

module load gnu9/9.3.0

make clean
make -f Makefile gfortran

### parallel

module load gnu9/9.3.0
module load openmpi4/4.0.4
module load intel_compilers/2022.2

make clean
make -f Makefile parallel



