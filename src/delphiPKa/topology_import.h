//
//  topology_import.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*
 *  This class is used to read topology file.
 *  Topology file contains bond connectivity, hydrogen information, reference pKa values and etc.
 *  The format of topology file is essential for correctly reading. DO NOT change the format.
 *  If new information is needed, follow the current format.
 *
 */

#ifndef TOPOLOGY_IMPORT_
#define TOPOLOGY_IMPORT_

#include <iostream>
#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"


class CTologyImport {
private:
    
    /* From Data_STORE */
    
    string& parameter_input;
    string& siz_file;
    string& crg_file;
    string& HofGLU;
    string& HofASP;
    
    unordered_map<string, float>& sizmap;
    unordered_map<string, float>& crgmap;
    unordered_map<string, float>& pkamap;
    
    mapTopology& topology_map;
    
    /* End */
    
    int id;
    TOPOLOGY topology_temp;
    string line;
    string key;
    string linetmp;
    float val;
    
    void string_to_topology(const string& line, TOPOLOGY& topology_temp);
    void topology_import();
    void chargeinfo_import();
    void sizeinfo_import();
    
public:
    
    CTologyImport(shared_ptr<DATA_STORE> pData):
    parameter_input( pData->parameter_input ),
    siz_file(pData->siz_file),
    crg_file(pData->crg_file),
    HofGLU(pData->HofGLU),
    HofASP(pData->HofASP),
    sizmap(pData->sizmap),
    crgmap(pData->crgmap),
    pkamap(pData->pkamap),
    topology_map( pData->map )
    
    
    {};
    
    void run();
    
    
    
    ~CTologyImport(){};
};

#endif // TOPOLOGY_IMPORT_

