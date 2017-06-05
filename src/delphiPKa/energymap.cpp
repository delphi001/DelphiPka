//
//  energymap.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "energy.h"


/*
 *  synchronize the energy map.
 *  when MPI is applied, the sync stating point from each array is defined as MPI_START
 */

#ifdef MPI_PARALLEL
void CEnergy :: genEnergyMap(const int& k, float& ergpolar, const bool& bNeutral, int& MPI_START, float * globalEnergyVec_partial)
#endif

#ifndef MPI_PARALLEL
void CEnergy :: genEnergyMap(const int& k, float& ergpolar, const bool& bNeutral)
#endif

{
    int i, j, sizei,resnum1, resnum2;
    string key, resnam2, strtmp;
    float fCharge, enrgCharge, enrgNeutral;
    
    resnum1 = newPDB[vecIonRes[k][0]].res_num;

    for (i=0;i<vecIonRes.size();i++) {
		sizei = vecIonRes.size();
        if (i == k) {
            
            if(bNeutral) {
                
#ifndef MPI_PARALLEL
                EnergyPair[2*k+1][2*i  ].resnum1 = resnum1;
                EnergyPair[2*k+1][2*i  ].resnum2 = resnum1;
                EnergyPair[2*k+1][2*i+1].resnum1 = resnum1;
                EnergyPair[2*k+1][2*i+1].resnum2 = resnum1;
                
                EnergyPair[2*k+1][2*i  ].fEnergy = 0.0;
                EnergyPair[2*k+1][2*i+1].fEnergy = 0.0;
#endif

#ifdef MPI_PARALLEL
				globalEnergyVec_partial[4*(k-MPI_START)*sizei + 4*i + 0] = 0.0;
				globalEnergyVec_partial[4*(k-MPI_START)*sizei + 4*i + 1] = 0.0;
#endif
				
            }
            else if(!bNeutral) {
                
#ifndef MPI_PARALLEL
                EnergyPair[2*k][2*i  ].resnum1 = resnum1;
                EnergyPair[2*k][2*i  ].resnum2 = resnum1;
                EnergyPair[2*k][2*i+1].resnum1 = resnum1;
                EnergyPair[2*k][2*i+1].resnum2 = resnum1;
                
                EnergyPair[2*k][2*i  ].fEnergy = 0.0;
                EnergyPair[2*k][2*i+1].fEnergy = 0.0;
#endif

#ifdef MPI_PARALLEL
				globalEnergyVec_partial[4*(k-MPI_START)*sizei + 4*i + 2] = 0.0;
				globalEnergyVec_partial[4*(k-MPI_START)*sizei + 4*i + 3] = 0.0;
#endif
				
            }
            
        }
        
        else {
            
            resnum2 = newPDB[vecIonRes[i][0]].res_num;
            resnam2 = newPDB[vecIonRes[i][0]].res_name;
            
            enrgCharge  = 0.0;
            enrgNeutral = 0.0;
            
            for(j=0;j<vecIonRes[i].size();j++) {
             
                // Calculate the charged state
                
                key         = resnam2 + " " + newPDB[vecIonRes[i][j]].atom_name;
                
                if (crgmap.find(key) != crgmap.end()) {
                    fCharge = crgmap.find(key)->second ;
                }
                else {
                    fCharge = 0.0;
                }
                
                
                enrgCharge += fCharge * vecGridPotential[vecIonRes[i][j]];
                
                
                
                // Calculate the neutral state
                
                if (resnam2 == "ASP")   strtmp = "AS0";
                if (resnam2 == "GLU")   strtmp = "GL0";
                if (resnam2 == "ARG")   strtmp = "AR0";
                if (resnam2 == "HIS")   strtmp = "HI0";
                if (resnam2 == "LYS")   strtmp = "LY0";
                if (resnam2 == "A")     strtmp = "A0";
                if (resnam2 == "C")     strtmp = "C0";
				if (resnam2 == "DA")	strtmp = "DA0";
				if (resnam2 == "DC")	strtmp = "DC0";
                
                key          = strtmp + " " + newPDB[vecIonRes[i][j]].atom_name;
                
                if (crgmap.find(key) != crgmap.end()) {
                    fCharge = crgmap.find(key)->second ;
                }
                else {
                    fCharge = 0.0;
                }
                
                enrgNeutral += fCharge * vecGridPotential[vecIonRes[i][j]] ;
                
                
            }

            
            if(!bNeutral) {
                
#ifndef MPI_PARALLEL
                EnergyPair[2*k][2*i  ].resnum1 = resnum1;
                EnergyPair[2*k][2*i  ].resnum2 = resnum2;
                EnergyPair[2*k][2*i+1].resnum1 = resnum1;
                EnergyPair[2*k][2*i+1].resnum2 = resnum2;
                
                EnergyPair[2*k][2*i  ].fEnergy = enrgCharge;
                EnergyPair[2*k][2*i+1].fEnergy = enrgNeutral;
#endif

#ifdef MPI_PARALLEL
				globalEnergyVec_partial[4*(k-MPI_START)*sizei + 4*i + 0] = enrgCharge;
				globalEnergyVec_partial[4*(k-MPI_START)*sizei + 4*i + 1] = enrgNeutral;
#endif
            }
            
            else if (bNeutral) {
                
#ifndef MPI_PARALLEL
                EnergyPair[2*k+1][2*i  ].resnum1 = resnum1;
                EnergyPair[2*k+1][2*i  ].resnum2 = resnum2;
                EnergyPair[2*k+1][2*i+1].resnum1 = resnum1;
                EnergyPair[2*k+1][2*i+1].resnum2 = resnum2;

				EnergyPair[2*k+1][2*i  ].fEnergy = enrgCharge;
				EnergyPair[2*k+1][2*i+1].fEnergy = enrgNeutral;
#endif				

#ifdef MPI_PARALLEL
				globalEnergyVec_partial[4*(k-MPI_START)*sizei + 4*i + 2] = enrgCharge;
				globalEnergyVec_partial[4*(k-MPI_START)*sizei + 4*i + 3] = enrgNeutral;
#endif
            }

            
        }
        
        
    }

   ///////////// Calculate the Polar Energy Component /////////////
    
    ergpolar = 0;

    for (int itr=0; itr<newPDB.size();itr++) {

        if (newPDB[itr].bIonizable && newPDB[itr].conf == "BK") {
            if (newPDB[itr].res_name == "ASP")   strtmp = "AS0";
            if (newPDB[itr].res_name == "GLU")   strtmp = "GL0";
            if (newPDB[itr].res_name == "ARG")   strtmp = "AR0";
            if (newPDB[itr].res_name == "HIS")   strtmp = "HI0";
            if (newPDB[itr].res_name == "LYS")   strtmp = "LY0";
            if (newPDB[itr].res_name == "A")     strtmp = "A0";
            if (newPDB[itr].res_name == "C")     strtmp = "C0";
			if (newPDB[itr].res_name == "DA")    strtmp = "DA0";
			if (newPDB[itr].res_name == "DC")    strtmp = "DC0";
            
            key       = strtmp + " " + newPDB[itr].atom_name;
            
            if (crgmap.find(key) != crgmap.end()) {
                fCharge = crgmap.find(key)->second ;
            }
            else {
                fCharge = 0.0;
            }
            
            ergpolar += vecGridPotential[itr] * fCharge;

        }
        
        else if (!newPDB[itr].bIonizable) {
            
            key       = key_to_hashmapQR(newPDB[itr]);
            
            if (crgmap.find(key) != crgmap.end()) {
                fCharge = crgmap.find(key)->second ;
            }
            else {
                fCharge = 0.0;
            }
            
            ergpolar += vecGridPotential[itr] * fCharge;
            
        }
        
    }
    
    int newPDBSIZE = newPDB.size();
    
    for (int itr=0; itr<vecCrgHETATM.size(); itr++) {
        
        ergpolar += vecGridPotential[newPDBSIZE+itr] * vecCrgHETATM[itr];
        
//        cout << vecGridPotential[newPDBSIZE+itr] << endl;

    }
//    cout << endl;
    
}