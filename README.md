# DelphiPka v2.3
Delphi-PKA is a DelPhi-based C++ program, allowing to predict pKa's for ionizable groups in proteins, RNA and DNA.
Unique approach stems from:
1. Use gaussian-based smooth function to mimic conformational changes associated with ionization changes.
2. Calculate the electrostatic energy without defining molecular surface.

Our web server: <http://compbio.clemson.edu/pka_webserver/>

For questions and help, visit <http://compbio.clemson.edu/forum/> or email to <delphi@g.clemson.edu>

### Compile DelphiPka
To compile the source code, you will need following library:
* [boost library](http://www.boost.org)
* [GSL library](http://www.gnu.org/software/gsl/) (GNU Scientific Library)
* [MPI library](http://www.open-mpi.org) (open-mpi v1.6.4 recommended)

C++ compiler needs to be 4.4 or higher, which includes C++11 feature.

With all the required above, run:
```bash
make
```

#### For the users that do not want to compile MPI version:
Edit `src\delphiPKa\prime_environment.h`:

And delete two lines `#define MPI_PARALLEL` and `#include <mpi.h>`.

Then change the `CC` entry in `Makefile` to `g++` or `clang++`

### Run DelphiPka on Palmetto HPC

You need to load the following modules on Palmetto first:
* gcc/4.8.1
* openmpi/1.6.4
* gsl/1.16

To run the program, you need:
1. PDB file
2. run.prm file (A sample can be found in current directory)
3. PBS script (A sample can be found in current directory)

Steps:
1. Open the run.prm and edit pdb name entry, charge and radius parameter entries with your desire. Currently supports amber, charmm22 and parse parameters.
2. For other entries, ref the [manual](http://compbio.clemson.edu/pka_webserver/assets/manual/DelPhiPKa_User_Manual.html) and edit with your own desire.
3. Change the sample.pbs to your desired job name. And qsub the pbs script.
