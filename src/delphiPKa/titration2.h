//
//  titration2.h
//  DelPhiPKA
//
//  Created by Lin Wang on 12/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//
//  This titration2.h and titration2.cpp are providing for dynamic network algorithm in network.cpp, it is now the primary approach.
//

#ifndef __TITRATION2_H_
#define __TITRATION2_H_

#include <iostream>
#include <cmath>
#include <deque>
#include <bitset>
#include <gsl/gsl_fit.h>
#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"

using namespace std;

class CTitration2 {
private:

    /* From Data_STORE */

    const string                       &pdb_input;
    const int                          &n_cluster;
    const vector<vector<int> >         &vecCluster;
    const vector<vector<string> >      &vecClusterChID;
    const vector<vector<int> >         &vecIonRes;
    const vector<string>               &strHETATM;
    vector<PDBFORM>                    &newPDB;
    const unordered_map<string, float> &pkamap;

    const vector<vector<Energy> > &EnergyPair;
    const vector<float>           &EnergyPolarCrg;
    const vector<float>           &EnergyPolarNeu;
    const vector<float>           &EnergyRxnCrg;
    const vector<float>           &EnergyRxnNeu;


    const float &pH_initial;
    const float &pH_end;
    const float &pH_step;
    const float &fGivenPhVal;
    const bool  &bOutPQRpka;
    const bool  &bRemoveHETATM;
    const bool  &bHETATMinPQR;

    unordered_map<string, float> &sizmap;
    unordered_map<string, float> &crgmap;


    /* End */

    vector<vector<int> >     vecState;
    vector<vector<int> >     vecCluster1;
    vector<vector<respair> > vecPair;
    vector<vector<float> >   vec2dProb;

    void genState(int &);

    void genCluster1();

    void linearReg();

    void output_PQR();

    float genEnergy0(int &, float &);


public:
    CTitration2(shared_ptr<DATA_STORE> pData) :

            pdb_input(pData->pdb_input),
            n_cluster(pData->n_cluster),
            newPDB(pData->newPDB),
            pkamap(pData->pkamap),
            vecCluster(pData->vecCluster),
            vecClusterChID(pData->vecClusterChID),
            vecIonRes(pData->vecIonRes),
            strHETATM(pData->strHETATM),

            sizmap(pData->sizmap),
            crgmap(pData->crgmap),

            EnergyPair(pData->EnergyPair),
            EnergyPolarCrg(pData->EnergyPolarCrg),
            EnergyPolarNeu(pData->EnergyPolarNeu),
            EnergyRxnCrg(pData->EnergyRxnCrg),
            EnergyRxnNeu(pData->EnergyRxnNeu),

            bOutPQRpka(pData->bOutPQRpka),
            bRemoveHETATM(pData->bRemoveHETATM),
            bHETATMinPQR(pData->bHETATMinPQR),
            pH_initial(pData->pH_initial),
            pH_end(pData->pH_end),
            pH_step(pData->pH_step),
            fGivenPhVal(pData->fGivenPhVal) {};

    void run();


    ~CTitration2() {};
};

#endif // __TITRATION2_H_
