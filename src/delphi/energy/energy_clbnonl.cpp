/*
 * energy_clbnonl.cpp
 *
 * Author: Lin Wang, lwang3@g.clemson.edu
 *
 * This file is function of coulombic energy calculation in non-linear case. 
 * It's called by energy_run function if non-linear case evokes.
 * 
 */

#include "energy.h"

void CDelphiEnergy::energy_clbnonl(delphi_real& fEnergy_Coulombic, delphi_real& fEnergy_SolvToChgIn, int& iGridOutput)
{
	int i, j, n;
	delphi_real dx, dy, dz;
	delphi_real fEnergy_Temp1, fEnergy_Temp2, fDistance = 0.0;
    delphi_real fEnergy_Temp = 0.0, fEnergy_SolvToChgIn_Temp=0.0;
	delphi_real c = 0.0006023;

    
    // +++++++ Define these variables in this scope for OpenMP ++++++ //
    
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> sout_nX(iGridOutput);
    vector<delphi_real> sout_nY(iGridOutput);
    vector<delphi_real> sout_nZ(iGridOutput);
    vector<delphi_real> sout_nValue(iGridOutput);
    vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);
    vector<delphi_real> prgfAtomEps_Val(iCrgGridNum);
    
    
    for(j=0;j<iCrgGridNum;j++)
    {
        
        prgfgCrgPoseA_nX[j] = prgfgCrgPoseA[j].nX;
        prgfgCrgPoseA_nY[j] = prgfgCrgPoseA[j].nY;
        prgfgCrgPoseA_nZ[j] = prgfgCrgPoseA[j].nZ;
        prggvAtomicCrg_nValue[j] = prggvAtomicCrg[j].nValue;
        prgfAtomEps_Val[j] = prgfAtomEps[j];
    }
    for(n=0;n<iGridOutput;n++)
    {
        sout_nX[n] = sout[n].nGrid.nX;
        sout_nY[n] = sout[n].nGrid.nY;
        sout_nZ[n] = sout[n].nGrid.nZ;
        sout_nValue[n] = sout[n].nValue;
    }
    
  // ++++++++ Starting parallel computing and evoke OpenMP +++++++++++ //
  // Refer OpenMP website for more details about tutorials and manuals //
  // http://openmp.org/wp/openmp-specifications 					   //
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

#ifdef PARALLEL_OMP

#pragma omp parallel shared(prgfgCrgPoseA_nX,prgfgCrgPoseA_nY,prgfgCrgPoseA_nZ,sout_nX,sout_nY,sout_nZ,sout_nValue,prggvAtomicCrg_nValue,prgfAtomEps_Val) private(j,n,dx,dy,dz,fDistance,fEnergy_Temp1,fEnergy_Temp2)
{
	#pragma omp for reduction( + : fEnergy_Temp, fEnergy_SolvToChgIn_Temp)

#endif
	
	for(i=0;i<iCrgGridNum;i++){
		fEnergy_Temp1 = 0.0;
		fEnergy_Temp2 = 0.0;
		for(j=0;j<iCrgGridNum;j++){
			if(i!=j){
                
                dx = prgfgCrgPoseA_nX[i] - prgfgCrgPoseA_nX[j];
                
                dy = prgfgCrgPoseA_nY[i] - prgfgCrgPoseA_nY[j];
                
                dz = prgfgCrgPoseA_nZ[i] - prgfgCrgPoseA_nZ[j];
                
                fDistance = sqrt(dx*dx + dy*dy + dz*dz);
                
                fEnergy_Temp1 += prggvAtomicCrg_nValue[j] / fDistance;
			}
		}

		fEnergy_Temp += prggvAtomicCrg_nValue[i] * fEnergy_Temp1 / prgfAtomEps_Val[i];

		// Calculation for the solvent contribution
	    // Array sout is declared in pointers module and allocated in nlener subroutine

		if(fIonStrength > 1.0e-6){
			for(n=0;n<iGridOutput;n++){
             
                dx = prgfgCrgPoseA_nX[i] - sout_nX[n];

                dy = prgfgCrgPoseA_nY[i] - sout_nY[n];
               
                dz = prgfgCrgPoseA_nZ[i] - sout_nZ[n];
                
                fDistance = sqrt(dx*dx + dy*dy + dz*dz);
                
                fEnergy_Temp2 += sout_nValue[n]/ fDistance;
			}

			fEnergy_SolvToChgIn_Temp += fEnergy_Temp2 * prggvAtomicCrg_nValue[i];
		}
	}
	
#ifdef PARALLEL_OMP

}	// end of #pragma omp parallel

#endif

	fEnergy_Coulombic = fEnergy_Temp / 2.0;
	fEnergy_SolvToChgIn = fEnergy_SolvToChgIn_Temp * c / (2.0 * fEpsOut);


// +++++++++++++++++++++ Clean memory occupied by these vectors ++++++++++++++++ //

    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);
    vector<delphi_real>().swap(sout_nX);
    vector<delphi_real>().swap(sout_nY);
    vector<delphi_real>().swap(sout_nZ);
    vector<delphi_real>().swap(sout_nValue);
    vector<delphi_real>().swap(prgfAtomEps_Val);
    

}

