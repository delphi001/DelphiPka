/*
 * energy_clb.cpp
 *
 * Author: Lin Wang, lwang3@g.clemson.edu
 *
 * This file is function of coulombic energy calculation in linear case. 
 * It's called by energy_run function if linear case evokes and media number is only 1.
 * 
 */

#include "energy.h"


void CDelphiEnergy::energy_clb(delphi_real& fEnergy_Coulombic)
{
    int i, j;
	delphi_real dx, dy, dz;
    delphi_real fEnergy_Temp = 0.0;
	delphi_real fEnergy_Coulombic_Temp = 0.0;
    delphi_real fDistance = 0.0;
 
// +++++++ Define these variables in this scope for OpenMP ++++++ //
        
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);

    for(j=0;j<iCrgGridNum;j++)
    {
        
        prgfgCrgPoseA_nX[j] = prgfgCrgPoseA[j].nX;
        prgfgCrgPoseA_nY[j] = prgfgCrgPoseA[j].nY;
        prgfgCrgPoseA_nZ[j] = prgfgCrgPoseA[j].nZ;
        prggvAtomicCrg_nValue[j] = prggvAtomicCrg[j].nValue;
    }
    
// ++++++++ Starting parallel computing and evoke OpenMP +++++++++++ //
// Refer OpenMP website for more details about tutorials and manuals //
// http://openmp.org/wp/openmp-specifications 						 //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

#ifdef PARALLEL_OMP

   #pragma omp parallel shared(prgfgCrgPoseA_nX,prgfgCrgPoseA_nY,prgfgCrgPoseA_nZ,prggvAtomicCrg_nValue) private(j,dx,dy,dz,fDistance,fEnergy_Temp)
{	
	#pragma omp for reduction( + : fEnergy_Coulombic_Temp )

#endif
	
	for(i=0;i<iCrgGridNum-1;i++){
	
        fEnergy_Temp = 0.0;

		for(j=i+1;j<iCrgGridNum;j++){
			
			// ++++++ Calculation dot product of the vector in x,y,z diction ++++ //
           
            dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
            
            dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
            
            dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
            
            fDistance = sqrt(dx*dx + dy*dy + dz*dz);
            
            fEnergy_Temp += prggvAtomicCrg_nValue[j] / fDistance;
        
        }


        fEnergy_Coulombic_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp;
        
	}

#ifdef PARALLEL_OMP
	
}	// end of #pragma omp parallel

#endif

	fEnergy_Coulombic = fEnergy_Coulombic_Temp;

// +++++++++++++++++++++ Clean memory occupied by these vectors ++++++++++++++++ //

    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);
}
