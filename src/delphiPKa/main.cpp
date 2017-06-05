//
//  main.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/**
 *  @file main.cpp
 *  @author Lin Wang, lwang3@g.clemson.edu or wonglynn2004@gmail.com
 *
 *  DelPhiPKa is a pKa calculator written in C++, which takes a PDB as input and protonate the structure,
 *  calculate the energies for each ionizable residue in protonated and de-protonated state, then generate
 *  energy matrix. And with all energy terms, it will calculate the Boltzmann distribution and furhter to 
 *  calculate the pKa value by drawing a titration curve from pH 0 to 14. The energy calculator is using
 *  DelPhi C++ and paralleled with OpenMPI library.
 *
 *  ====================================================================================================================
 *
 *  Several libraries are required to compile:
 *  Boost library for DelPhi C++
 *  OpenMPI library for pKa code if you want it work paralleled. (It can be turned off by modify the prime_environment.h)
 *  GSL - GNU Scientific Library for fitting. 
 *  The code contains some C++11 feature and syntax, so make sure to add -std=c++11 flag in makefile
 *
 *  ====================================================================================================================
 *
 *  [Changelog]
 *
 *  - ver 1.0
 *  This is the first public release version in 7/2015 as version 1.0
 *  It contains two components (as two folders). One is delphi c++ and the other is delphiPKa
 *  delphi folder can be updated with newer version delphi c++ by replacing the sub-folders EXCEPT app folder,
 *  which contains the API that works only for delphiPKa and this file is not contained in delphi c++ code. 
 *  As now, DelPhiPKa v1.0 use delphicpp.r67 source code. For further development, mark the delphicpp version 
 *  in this Changelog if delphicpp component has been updated or modified.
 *
 */
 
 


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>


#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"
#include "readparam.h"
#include "topology_import.h"
#include "pdb_import.h"
#include "orbt_type.h"
#include "place_hydrogen.h"
#include "clustering.h"
#include "network.h"
#include "energy.h"
#include "titration.h"
#include "titration2.h"




using namespace std;


void banner() {
    
    cout << "   _________________________    DelPhiPKa   ver 1.0    _________________________ " << endl;
    cout << "  |                                                                             |" << endl;
    cout << "  |   DelPhiPKA is a DelPhi-based C++ program, allowing to predict pKa's for    |" << endl;
    cout << "  |   ionizable groups in proteins, RNA and DNA.  Unique approach stems from:   |" << endl;
    cout << "  |   a) Usage of gaussian-based smooth function to mimic conformational        |" << endl;
    cout << "  |      changes associated with ionization changes.                            |" << endl;
    cout << "  |   b) Calculate the electrostatic energy without defining molecular surface  |" << endl;
    cout << "  |                                                                             |" << endl;
    cout << "  |                         For questions and help, visit                       |" << endl;
    cout << "  |                       http://compbio.clemson.edu/forum/                     |" << endl;
    cout << "  |                        or email to delphi@g.clemson.edu                     |" << endl;
    cout << "  |                                                                             |" << endl;
    cout << "  |                          Developed by Lin Wang, 2015                        |" << endl;
    cout << "  |              ( Dr. Emil Alexov's group at Clemson University )              |" << endl;
    cout << "  |                                                                             |" << endl;
    cout << "  |_________________________    DelPhiPKa   ver 1.0    _________________________|" << endl;

    cout << endl;cout << endl;
}

/*  For invalid input */

void invalid_input() {
    
    cout << "[Error] Invalid input. Please follow the correct execution below." << endl;

    cout << "                   delphiPKa run.prm              " << endl;

    cout << "  Program now exit.." << endl;
}


/*  The main function of DelPhiPKa  */
/*  MPI is merged in single set of code. Thus, MPI is initialed at this point. Anytime later, use
 *
 *  #ifdef MPI_PARALLEL
 *  if(id==0) {
 *  #endif
 *
 *  to perform operations on master node.
 */
 

int main(int argc, const char * argv[]) {

    
#ifdef MPI_PARALLEL
    
    MPI::Init();

    int id = MPI::COMM_WORLD.Get_rank();
    
#endif

    
    
#ifdef MPI_PARALLEL
    if(id==0) {
#endif
        banner();
        
#ifdef MPI_PARALLEL
    }
#endif

            
    if(argc!=2) {
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
            invalid_input();
        
#ifdef MPI_PARALLEL
        }
        MPI::Finalize();
#endif
        exit(0);
    }
    
//  smart pointer, use "-std=c++11" flag to compile  //
    
    shared_ptr<DATA_STORE> pData ( new DATA_STORE() ); // the shared variables, refer to data_store.h
    pData->paramfile = argv[1];
    
    unique_ptr<CReadParam> pReadParam ( new CReadParam(pData)); // read the runtime paramters. default is run.prm.
    pReadParam->run();

    unique_ptr<CTologyImport> pTologyImport ( new CTologyImport(pData) ); // read the topology parameter file.
    pTologyImport->run();

    unique_ptr<CPdbImport> pPdbImport ( new CPdbImport(pData) ); // read PDB file
    pPdbImport->run();

    if(pData->bDoProton){
        unique_ptr<CPlaceHydrogen> pPlaceHydrogen ( new CPlaceHydrogen(pData) );
        pPlaceHydrogen->run();
        if(pData->bOutPQRtopo)     pPlaceHydrogen->output_newPDB();
    }

    if(!pData->bDoEnergy && !pData->bDoPka)
    {
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
        cout << endl;
        cout << "Energy and pKa's calculation are skipped upon request, program is about to exit.  " << endl;
        cout << endl;
#ifdef MPI_PARALLEL
        }
        MPI::Finalize();
#endif
                
        exit(0);
    }
    
    unique_ptr<CEnergy> pEnegy ( new CEnergy(pData) ); // energy class
    if(!pData->bDoEnergy){
        
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
        cout << endl;
        cout << "Energy calcultion with DelPhi is skipped by request ." << endl;
        cout << "Energy table is read from previous logfiles : energies.txt and pairwise.txt . "      << endl;
            
#ifdef MPI_PARALLEL
        }
#endif
        pEnegy->readFromFile();
    }
        

    if(pData->bDoEnergy){
    
        pEnegy->run();
        
    }
    
    
    if(pData->bDoPka) {

        //  The following is to evoke dynamic network algorithm, and now is primary approach
        
        unique_ptr<CNetwork> pNetwork ( new CNetwork(pData) );
        pNetwork->run();
        
        unique_ptr<CTitration2>  pTitration ( new CTitration2(pData) );
        pTitration->run();

        
        //  The following is to evoke kmean++ algorithm, which is optional, not primary approach
    /*

        unique_ptr<CClustering> pClustering ( new CClustering(pData) );
        pClustering->run();
        
        unique_ptr<CTitration>  pTitration ( new CTitration(pData) );
        pTitration->run();
    */
        
    }
        

    
    pData.reset();
    

#ifdef MPI_PARALLEL
    
    MPI::Finalize();
    
#endif
    
    return 0;
}
