//
//  place_hydrogen.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*
 *  With topology and atomic orbitals, this class is used to protonate the PDB by adding hydrogens.
 *  PQR file will be generated for the rest of subroutines.
 */

#ifndef PLACE_HYDROGEN_
#define PLACE_HYDROGEN_

#include <iostream>
#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"
using namespace std;

class CPlaceHydrogen {
private:
    
    /* From Data_STORE */
    
    
    const string& pdb_input;
    const bool& bRemoveHETATM;
    const bool& bOutPQRtopo;
    
    mapTopology& map;
    mapPDBFORM& mapPDB;
    
    vector<PDBFORM>& PDB;
    vector<PDBFORM>& newPDB;
    vector<string>&  strHETATM;
    
    unordered_map<string, float>& sizmap;
    unordered_map<string, float>& crgmap;
    
    /* End */

    int id;
    int found_bond_num;
    int num_of_H_to_add;
    PDBFORM atom_to_insert;
    TOPOLOGY value_found;
    
    nVECTOR Vector0, Vector1, Vector2, Vector3, Vector4, Vector5;
    
    float BondLength, BondAngle_sp3=PI*120/180, BondAngle_sp2=PI*109/180;
    float TorsionAngle_sp2 = PI, TorsionAngle_sp3 = PI/3.;
    
    string tpl_key;
    string pdb_search_key;

    vector<PDBFORM>::iterator it_pdb;
    
public:
    CPlaceHydrogen(shared_ptr<DATA_STORE> pData) :
    map(pData->map),
    mapPDB(pData->mapPDB),
    PDB(pData->PDB),
    newPDB(pData->newPDB),
    sizmap(pData->sizmap),
    crgmap(pData->crgmap),
    strHETATM(pData->strHETATM),
    pdb_input(pData->pdb_input),
    bOutPQRtopo(pData->bOutPQRtopo),
    bRemoveHETATM(pData->bRemoveHETATM)
    
    
    {};

    
    void run();
    
    void output_newPDB();
    
    ~CPlaceHydrogen() {};
    
};


#endif // PLACE_HYDROGEN_
