//
//  readparam.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*
 *  This class is used to read paramter file, force-field files (crg and siz files)
 *  The force-field files are the same as DelPhi c++ use
 */

#ifndef READPARAM_
#define READPARAM_

#include <iostream>
#include <fstream>
#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"

class CReadParam {
private:
    
    /* From Data_STORE */
    
    string& paramfile;
    string& pdbname;
    string& parameter_input;
    string& siz_file;
    string& crg_file;
    string& HofGLU;
    string& HofASP;
    
    bool& bHETATMinPQR;
    bool& bClusterAuto;
    bool& bRemoveHETATM;
    bool& bRemoveWater;
    bool& bOutPQRtopo;
    bool& bOutPQRpka;
    
    bool& bDoProton;
    bool& bDoEnergy;
    bool& bDoPka;
    int& n_cluster;
    int& iGaussian;
    float& fSigma;
    float& fSrfcut;
    float& indi;
    float& exdi;
    float& ionrad;
    float& prbrad;
    float& maxc;
    float& clusterThreshold;
    float& fGivenPhVal;
    
    float& pH_initial;
    float& pH_end;
    float& pH_step;
    
    /*  END */
    
    string line;
    
public:
    
    CReadParam(shared_ptr<DATA_STORE> pData) :
    paramfile(pData->paramfile),
    pdbname(pData->pdb_input),
    parameter_input(pData->parameter_input),
    siz_file(pData->siz_file),
    crg_file(pData->crg_file),
    HofGLU(pData->HofGLU),
    HofASP(pData->HofASP),
    iGaussian(pData->iGaussian),
    fSigma(pData->fSigma),
    fSrfcut(pData->fSrfcut),
    indi(pData->indi),
    exdi(pData->exdi),
    ionrad(pData->ionrad),
    prbrad(pData->prbrad),
    maxc(pData->maxc),
    clusterThreshold(pData->clusterThreshold),
    fGivenPhVal(pData->fGivenPhVal),
    bHETATMinPQR(pData->bHETATMinPQR),
    bRemoveHETATM(pData->bRemoveHETATM),
    bRemoveWater(pData->bRemoveWater),
    bOutPQRtopo(pData->bOutPQRtopo),
    bOutPQRpka(pData->bOutPQRpka),
    bClusterAuto(pData->bClusterAuto),
    bDoProton(pData->bDoProton),
    bDoEnergy(pData->bDoEnergy),
    bDoPka(pData->bDoPka),
    n_cluster(pData->n_cluster),
    pH_initial(pData->pH_initial),
    pH_end(pData->pH_end),
    pH_step(pData->pH_step)
    
    
    {};
    
    void run();
    
    ~CReadParam(){};
    
};

#endif // READPARAM__
