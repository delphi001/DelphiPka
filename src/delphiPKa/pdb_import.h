//
//  pdb_import.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*
 *  This class is used to read PDB file.
 *  Special naming rules are needed to include in pdb_import.cpp file, 
 *  so that delphiPKa can recgonize the non-standard residue naming.
 */



#ifndef PDB_IMPORT_
#define PDB_IMPORT_

#include <iostream>
#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"
using namespace std;

class CPdbImport {
private:
    
    /* From Data_STORE */
    
    string& pdb_input;
    vector<PDBFORM>& PDB;
    mapPDBFORM& mapPDB;
    vector<string>& strHETATM;
    vector<string>& strHETATM2;
    vector<float>&  vecCrgHETATM;
    bool& bRemoveHETATM;
    bool& bRemoveWater;
    
    /* End */

    int id;
    string line;
    string key;
    string strCrgHETATM;
    PDBFORM pdb_temp;
    
public:
    CPdbImport(shared_ptr<DATA_STORE> pData) :
    pdb_input(pData->pdb_input),
    PDB(pData->PDB),
    mapPDB(pData->mapPDB),
    strHETATM(pData->strHETATM),
    strHETATM2(pData->strHETATM2),
    vecCrgHETATM(pData->vecCrgHETATM),
    bRemoveHETATM(pData->bRemoveHETATM),
    bRemoveWater(pData->bRemoveWater)
    
    {};

    
    void run();
    
    ~CPdbImport(){};
};


#endif // PDB_IMPORT_