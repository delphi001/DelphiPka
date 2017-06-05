/*
 * energy_nl.cpp
 *
 *  Created on: Feb 23, 2014
 *      Author: Lin Wang
 */

#include "energy.h"

void CDelphiEnergy::energy_nonl(delphi_real& fEnergy_Nonlinear, int& iGridOutput)
{
	int i, j, k, iw, n=0;
    delphi_real fEnergy_Solvation, fEnergy_Osmetic;
	delphi_real gX, gY, gZ;
	delphi_real fCarica, fCutedges_i, fCutedges_j, fCutedges_k, goff;
	delphi_real dhi1, dhi2, dhi3, dhi4, dhi5, fCarica_Temp, fPhi, c;
	SGridValue<delphi_real> fGridValue_Temp;

	fEnergy_Solvation=0.0; 
	fEnergy_Osmetic=0.0;

	goff = (iGrid + 1.0) / 2.0;
	gX = (-goff / fScale) + fgBoxCenter.nX;
	gY = (-goff / fScale) + fgBoxCenter.nY;
	gZ = (-goff / fScale) + fgBoxCenter.nZ;
	c = fScale*fScale*fScale;

	dhi5=-fTaylorCoeff5/6.0;
	dhi4=-fTaylorCoeff4/5.0;
	dhi3=-fTaylorCoeff3/4.0;
	dhi2=-fTaylorCoeff2/3.0;
	dhi1=-fTaylorCoeff1/2.0;


		for(k=0;k<iGrid;k++){
			fCutedges_k = 1.0;
			if(k==0 || k==iGrid-1){
				fCutedges_k = 0.5;
			}
			for(j=0;j<iGrid;j++){
				fCutedges_j = fCutedges_k;
				if(j==0 || j==iGrid-1){
					fCutedges_j = fCutedges_k * 0.5;
				}
				for(i=0;i<iGrid;i++){
					fCutedges_i = fCutedges_j;
					if(i==0 || i==iGrid-1){
						fCutedges_i = fCutedges_j * 0.5;
					}

					iw = k*iGrid*iGrid+j*iGrid+i;
					if(prgbDielecMap[iw]){
						fPhi = prgfPhimap[iw];


						fCarica_Temp = fTaylorCoeff5 * fPhi + fTaylorCoeff4;
						fCarica_Temp = fCarica_Temp * fPhi + fTaylorCoeff3;
						fCarica_Temp = fCarica_Temp * fPhi + fTaylorCoeff2;
						fCarica_Temp = fCarica_Temp * fPhi + fTaylorCoeff1;
						fCarica = fCutedges_i * fCarica_Temp * fPhi / c;

						fCarica_Temp = dhi5 * fPhi + dhi4;
						fCarica_Temp = fCarica_Temp * fPhi + dhi3;
						fCarica_Temp = fCarica_Temp * fPhi + dhi2;
						fCarica_Temp = fCarica_Temp * fPhi + dhi1;
						fEnergy_Osmetic = fEnergy_Osmetic + fCutedges_i * fCarica_Temp * fPhi * fPhi;


						if(bIonsEng){

							fGridValue_Temp.nGrid.nX = i + gX + 1;
							fGridValue_Temp.nGrid.nY = j + gY + 1;
							fGridValue_Temp.nGrid.nZ = k + gZ + 1;
							fGridValue_Temp.nValue = fCarica;
							sout.push_back(fGridValue_Temp);

						}
						fEnergy_Solvation -= fCarica*fPhi;
					}
				}
			}
		}


	fEnergy_Solvation = fEnergy_Solvation * 0.5 * 0.0006023;
	fEnergy_Osmetic = -fEnergy_Osmetic * 0.0006023 / c;

	iGridOutput = sout.size();


	cout << " rho*phi/2 term in solution       :               " << setw(8) << right << -fEnergy_Solvation << "  kt" << endl;
	cout << " osmotic pressure term            :               " << setw(8) << right << -fEnergy_Osmetic << "  kt" << endl;


	fEnergy_Nonlinear = fEnergy_Nonlinear + fEnergy_Osmetic +fEnergy_Solvation;

}

