# DelphiPka

### Compile DelphiPka
To compile the source code, you will need following library:
1. [boost library](http://www.boost.org)
2. [GSL library](http://www.gnu.org/software/gsl/) (GNU Scientific Library)
3. [MPI library](http://www.open-mpi.org) (open-mpi recommended)

If you do not want to compile MPI version,
edit prime_environment.h in delphiPKa folder.
And delete two lines `#define MPI_PARALLEL` and `#include <mpi.h>`.
Then change the CC entry in Makefile to g++ or clang++

C++ compiler needs to be 4.4 or higher, which includes C++11 feature.

With all the required above, run:
```bash
make
```

Good luck!

### Run DelphiPka

This is a brief intro to run DelphiPKA on Palmetto

You need to load the following modules on Palmetto first:
gcc/4.8.1
openmpi/1.6.4
gsl/1.16

To run the program, you need:
1. PDB file
2. run.prm file (A sample can be found in current directory)
3. PBS script (A sample can be found in current directory)

Steps:
1. Copy the run.prm and PBS script with your PDB file in any working directory.

2. Open the run.prm and edit pdb name entry, charge and radius parameter entries with your desire. Currently supports amber, charmm22 and parse parameters. 

3. For other entries, ref the manual and edit with your own desire.

4. Change the sample.pbs to your desired job name. And qsub the pbs script.
