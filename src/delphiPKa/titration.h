//
//  titration.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//
//  This titration is used for kmean++ clustring algorithm which default is optional, primary clustering is dynamic network in network.cpp
//

#ifndef __TITRATION__
#define __TITRATION__

#include <iostream>
#include <cmath>
#include <deque>
#include <bitset>
#include <gsl/gsl_fit.h>
#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"

using namespace std;

class CTitration {
private:
    
    /* From Data_STORE */
    
    
    const int& n_cluster;
    const vector<vector<int> >&         vecCluster;
    const vector<vector<int> >&         vecIonRes;
    const vector<PDBFORM>&              newPDB;
    const unordered_map<string, float>& pkamap;
    
    const vector<vector<Energy> >&   EnergyPair;
    const vector<float>&             EnergyPolarCrg;
    const vector<float>&             EnergyPolarNeu;
    const vector<float>&             EnergyRxnCrg;
    const vector<float>&             EnergyRxnNeu;
    
    const float& pH_initial;
    const float& pH_end;
    const float& pH_step;
    
    /* End */

    vector<vector<int> >              vecState;
    vector<vector<int> >              vecCluster1;
    vector<vector<respair> >          vecPair;
    vector<vector<float> >            vec2dProb;
    
    void genState(int&);
    void genCluster1();
    void linearReg();
    float genEnergy0(int&, float&);
    

    
public:
    CTitration(shared_ptr<DATA_STORE> pData) :
    
    n_cluster(pData->n_cluster),
    newPDB(pData->newPDB),
    pkamap(pData->pkamap),
    vecCluster(pData->vecCluster),
    vecIonRes (pData->vecIonRes),
    
    EnergyPair(pData->EnergyPair),
    EnergyPolarCrg(pData->EnergyPolarCrg),
    EnergyPolarNeu(pData->EnergyPolarNeu),
    EnergyRxnCrg(pData->EnergyRxnCrg),
    EnergyRxnNeu(pData->EnergyRxnNeu),
    
    pH_initial(pData->pH_initial),
    pH_end(pData->pH_end),
    pH_step(pData->pH_step)
    
    {};
    
    void run();
    
    
    ~CTitration() {};
};


#endif // __TITRATION__
