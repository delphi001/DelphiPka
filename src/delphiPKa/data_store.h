//
//  data_store.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//


/* This file is the data storage. The shared varaibles including parameters read from input are stored in this class
 *
 */

#ifndef DATA_STORE_H_
#define DATA_STORE_H_

#include "global_type.h"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
using namespace std;

class DATA_STORE {
    

public:
    DATA_STORE(){};
    
    string paramfile;
    string pdb_input;
    string parameter_input; // The file contains topology information (Residue structure, connector, etc)
    string siz_file;
    string crg_file;
    string HofGLU  = "OE2";
    string HofASP  = "OD2";
    bool bHETATMinPQR = false;
    bool bDoProton = true;
    bool bDoEnergy = true;
    bool bDoPka    = true;
    bool bClusterAuto = true;
    bool bRemoveHETATM= true;
    bool bRemoveWater = true;
    bool bOutPQRtopo = false;
    bool bOutPQRpka  = false;
    
    int n_cluster;
    int iGaussian = 1;
    float fSigma  = 0.70;
    float fSrfcut = 40.0;
    float indi    =  8.0;
    float exdi    = 80.0;
    float ionrad  =  2.0;
    float prbrad  =  1.4;
    float maxc    =  0.001;
    float clusterThreshold = 15.0;
    float fGivenPhVal=  7.0;
    
    float pH_initial =  0.0;
    float pH_end     = 14.0;
    float pH_step    =  1.0;
    
    unordered_map<string, float> sizmap;
    unordered_map<string, float> crgmap;
    unordered_map<string, float> pkamap;
    
    
           // The input pdb file.

    mapTopology map;
    
    mapPDBFORM mapPDB;
    
    vector<PDBFORM> PDB; // The reserve operation is done in pdbimport class.
    
    vector<PDBFORM> newPDB;
    
    vector<string>  strHETATM;
    
    vector<string>  strHETATM2;
    
    vector<float>   vecCrgHETATM;
    
    vector<vector<int> > vecCluster;
    
    vector<vector<string> > vecClusterChID;

    vector<vector<int> > vecIonRes;
    
    vector<vector<Energy> > EnergyPair;
    
    vector<float> EnergyPolarCrg;
    
    vector<float> EnergyPolarNeu;
    
    vector<float> EnergyRxnCrg;
    
    vector<float> EnergyRxnNeu;
    
    ~DATA_STORE(){};

    
};

#endif // DATA_STORE_H_
