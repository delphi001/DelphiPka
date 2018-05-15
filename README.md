# DelphiPka v2.3
Delphi-PKA is a DelPhi-based C++ program, allowing to predict pKa's for ionizable groups in proteins, RNA and DNA.
Unique approach stems from:
1. Use gaussian-based smooth function to mimic conformational changes associated with ionization changes.
2. Calculate the electrostatic energy without defining molecular surface.

Our web server: <http://compbio.clemson.edu/pka_webserver/>

For questions and help, visit <http://compbio.clemson.edu/forum/> or email to <delphi@g.clemson.edu>

## Running DelphiPka on a PC
### Compile DelphiPka
To compile the source code, you will need following library:
* [boost library](http://www.boost.org)
* [GSL library](http://www.gnu.org/software/gsl/) (GNU Scientific Library)
* [MPI library (Multi cores support)](http://www.open-mpi.org) (open-mpi v1.6.4 recommended)

C++ compiler needs to be 4.4 or higher, which includes C++11 feature.

With all the required above, run:
```bash
make
```

### For the users that do not want to compile MPI version:
Edit `src\delphiPKa\prime_environment.h`:

And delete two lines `#define MPI_PARALLEL` and `#include <mpi.h>`.

Then change the `CC` entry in `Makefile` to `g++` or `clang++`

If you want to build the static version, you can modify the Makefile,
and add `-static` in `CFLAGS` entry and `LDFLAGS` entry (only support without MPI).

## Running DelphiPka on Palmetto HPC
### With PBS script (Recommend)
```bash
# login to the palmetto HPC login node
git clone https://github.com/delphi001/DelphiPka.git
# upload your pdb file
cd DelphiPka
# modify the run.prm and sample.pbs
qsub sample.pbs
```
The pbs script will claim the computing node, setup the environment, build the program, and run the delphipka automatically.
This script will run with 8 core by default. To change the cpu number, and memory size, you can modify the sample.pbs file.

### Build the program on Palmetto HPC, and run it manually
You need to load the following modules on Palmetto first:
* gcc/4.8.1
* openmpi/1.6.4
* gsl/1.16
```
# claim a computing node
# Recommend flag: select=1:ncpus=8:mem=40gb:mpiprocs=8:interconnect=fdr,walltime=72:00:00
module purge
module add gcc/4.8.1 openmpi/1.6.4 gsl/1.16
git clone https://github.com/delphi001/DelphiPka.git
cd DelphiPka
make -j8

mpirun --mca btl openib,self --mca btl_openib_warn_nonexistent_if 0 -np 8 bin/delphiPKa run.prm
```

## How to change config file
Steps:
1. Open the run.prm and edit pdb name entry, charge and radius parameter entries with your desire. Currently supports amber, charmm22 and parse parameters.
2. For other entries, ref the [manual](http://compbio.clemson.edu/pka_webserver/assets/manual/DelPhiPKa_User_Manual.html) and edit with your own desire.
3. Change the sample.pbs to your desired job name. And qsub the pbs script.
