//
//  environment.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*
 *  global variables
 *  physical and chemic constants
 *  MPI define is in this file and could be commented out if user want to compile WITHOUT MPI library.
 *
 */

#ifndef PRIME_ENVIRONMENT_H_
#define PRIME_ENVIRONMENT_H_


//#define PRIME_DEBUG

//#define DELPHI_OUTPUT
//#define DEBUG_ENERGY_EXPORT
//#define PRINT_CHARGED_PQR


#define MPI_PARALLEL    // This two lines need to be comment out if compiling the code WITHOUT MPI
#include <mpi.h>        // This two lines need to be comment out if compiling the code WITHOUT MPI

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <ctime>
#include <map>
#include <unordered_map>



using namespace std;

const float PI = 3.14159;
const float AbsoluteZero = -273.15; /// temperature of absolute zero
const float AtomicUnitCrg = 1.602176487e-19; ///< e, Coulomb C
const float BoltzmannConstant = 1.3806504e-23; ///< k, Joule per Kelvin JK^(-1)
const float EulerNum  = 2.718281828;
const float KCAL2KT   = 1.688;
const float ROOMTEMPER = 298.15;
const float VacuumPermittivity = 8.8541878176e-12; ///< e0, farads per meter Fm^(-1)
const float coval_rad_H = 0.300;  /// covalent radii for Hydrogen atom
const float coval_rad_C = 0.772;  /// covalent radii for Carbon atom
const float coval_rad_N = 0.700;  /// covalent radii for Nitrogen atom
const float coval_rad_O = 0.660;  /// covalent radii for Oxygen atom



#endif // PRIME_ENVIRONMENT_H_
