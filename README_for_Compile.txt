
To compile the source code, you will need following library:

1.boost library
http://www.boost.org

2. GSL library (GNU Scientific Library)
http://www.gnu.org/software/gsl/

3. MPI library (open-mpi recommended)
http://www.open-mpi.org

If you do not want to compile MPI version, edit prime_environment.h in delphiPKa folder.
And delete two lines "#define MPI_PARALLEL" and "#include <mpi.h>"
And change the CC entry in Makefile to g++ or c++

C++ compiler needs to be 4.4 or higher, which includes C++11 feature.

With all the required above, run:

make




Good luck!



