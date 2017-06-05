//
//  network.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*
 *  This class is used to partition the system with network algorithm
 */

#ifndef NETWORK_H_
#define NETWORK_H_

#include <iostream>
#include <fstream>
#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"
using namespace std;

class CNetwork {
private:
    
    /* From Data_STORE */
    
    float& clusterThreshold;
    
    const vector<PDBFORM>& newPDB;
    
    vector<vector<int> >& vecCluster;
    vector<vector<string> >& vecClusterChID;
    
    
    /* END */
    
    vector<PDBFORM> vecPDB;
    
    void importcenter();
    
    void grouping();
    
    float dist2(PDBFORM a, PDBFORM b);
    
    
public:
    CNetwork(shared_ptr<DATA_STORE> pData) :
    
    newPDB(pData->newPDB),
    vecCluster(pData->vecCluster),
    vecClusterChID(pData->vecClusterChID),
    clusterThreshold(pData->clusterThreshold)
    {};
    
    void run();
    
    
    ~CNetwork() {}
};


#endif // NETWORK_H_
