//
//  energy.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*  This class is used to do energy calculations.
 *  DelPhi is implemented in this class, and use MPI to do parallelization
 */


#ifndef ENERGY_H_
#define ENERGY_H_

#include <iostream>
#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"

#include "../delphi/delphi/delphi_constants.h"
#include "../delphi/io/io.h"


using namespace std;

class CEnergy {
private:
    
    /* From Data_STORE */
    
    const string& pdb_input;
    const float& indi;
    const float& exdi;
    const float& ionrad;
    const float& prbrad;
    const float& maxc;
    
    const int&   iGaussian;
    const float& fSigma;
    const float& fSrfcut;
    
    vector<PDBFORM>&           newPDB;
    vector<vector<int> >&      vecIonRes;
    vector<vector<Energy> >&   EnergyPair;
    vector<float>&             EnergyPolarCrg;
    vector<float>&             EnergyPolarNeu;
    vector<float>&             EnergyRxnCrg;
    vector<float>&             EnergyRxnNeu;
    
    const bool&                bHETATMinPQR;
    const vector<string>&      strHETATM;
    const vector<string>&      strHETATM2;
    const vector<float>&       vecCrgHETATM;
    const unordered_map<string, float>& sizmap;
    const unordered_map<string, float>& crgmap;
    
    /* End */
    
    int id, numprocs, MPI_START, MPI_END;
    string frcnam;
    string strFRCOutFile;
    vector<PDBFORM>::iterator it_pdb;
    vector<string> PQR_to_DelPhi;
    vector<string> singleRes_to_DelPhi;
    vector<bool>   vecCrgStat;
    vector<bool>   resInBox;
    vector<bool>   itr;
    SGrid<double>  center0, center1;
    
    vector<double> vecGridPotential;
	vector<float> globalEnergyVec;

    /*
     *  make the ionizable residues in their charged or neutral states depends on the passing variable
     *  Different naming rules are applied here. For charged residue naming, "ASP", "GLU", "HIS", "ARG",
     *  "LYS" are used. Corresponding residues in neutral is labled as "AS0", "GL0", "HI0", "AR0", "LY0".
     *
     */
    void charge(bool&); // energy_misc.cpp
    
    
    
    /*
     *  generate ionizable residue array, the array is used to further distribute with MPI
     */
    void ionizableArr();
    
    
    /*
     *  use to force each residue to be uncharged
     */
    void ResetCrgStat();
    
    
    /*
     *  evoke the homogeneous dielectric model in delphicpp.
     *  since gaussian flag is set to be true by default. this option needs some modification before it works
     *  only for testing purpose, otherwise do not use this function.
     */
    void HomoDielec(shared_ptr<SPrime> param);
    
    
    
    /*
     *  evoke the gaussian dielectric model in delphicpp.
     */
    void GaussDielec(shared_ptr<SPrime> param);
    
    
    
    /*
     *  generate residue side-chain in free state for desolvation energy calculation
     */
    void resSolv(shared_ptr<SPrime> param);
    
    
    /*
     *  set up the first step in a three focusing run, with scale = 1
     */
    void setbase(shared_ptr<SPrime> param, SGrid<double>& center1);
    
    
    
    /*
     *  set up the second and third steps in a three focusing run, with scale = 2 and 4
     */
    void setfocus(shared_ptr<SPrime> param);
    
    
    
    /*
     *  update potential map after each step during a focusing run. only potentials within the newer box will be updated.
     */
    void updateGridPotential(shared_ptr<SPrime> param);
    
    
    
    
    /*
     *  followed with setbase() and setfocus. call delphi class to run the focusing.
     */
    bool runFocus(shared_ptr<SPrime> param, SGrid<double>& center1);


#ifdef MPI_PARALLEL
    /*
     *  synchronize the energy map if MPI is applied
     */
    void genEnergyMap(const int&, float&, const bool&, int&, float *);
#endif
    
    
    
#ifndef MPI_PARALLEL
    /*
     *  synchronize the energy map if MPI is not applied
     */
    void genEnergyMap(const int&, float&, const bool&);
#endif
    
    
public:
    CEnergy(shared_ptr<DATA_STORE> pData) :
    pdb_input(pData->pdb_input),
    indi(pData->indi),
    exdi(pData->exdi),
    ionrad(pData->ionrad),
    prbrad(pData->prbrad),
    maxc(pData->maxc),
    iGaussian(pData->iGaussian),
    fSigma(pData->fSigma),
    fSrfcut(pData->fSrfcut),
    newPDB(pData->newPDB),
    sizmap(pData->sizmap),
    crgmap(pData->crgmap),
    vecIonRes(pData->vecIonRes),
    bHETATMinPQR(pData->bHETATMinPQR),
    strHETATM(pData->strHETATM),
    strHETATM2(pData->strHETATM2),
    vecCrgHETATM(pData->vecCrgHETATM),
    EnergyPair(pData->EnergyPair),
    EnergyPolarCrg(pData->EnergyPolarCrg),
    EnergyPolarNeu(pData->EnergyPolarNeu),
    EnergyRxnCrg(pData->EnergyRxnCrg),
    EnergyRxnNeu(pData->EnergyRxnNeu)
    
    {};
    
    void run();
    
    
    /*
     *  This funciton is used to read energy matrix and terms from text files: pairwise.txt and energy.txt generated by energy class
     *  This is used in case that energy class is set to be skipped but the rest of program needs to run. 
     *  Before the clustering class, call this function and the required energy terms are read into the data store class.
     *  Do not use this function if the energy class is set to be run, otherwise it will slow down the runtime.
     */
    void readFromFile();
    
    ~CEnergy() {};
};

#endif // ENERGY_H_
