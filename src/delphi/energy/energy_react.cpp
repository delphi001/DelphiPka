/*
 * energy_react.cpp
 *
 * Author: Lin Wang, lwang3@g.clemson.edu
 *
 * This file is function of reaction field energy (solvation) calculation.
 * It's called by energy_run function if energy function solvation(s or sol) flag set to be TRUE.
 * 
 */

#include "energy.h"

void CDelphiEnergy::energy_react(delphi_real& fEnergy_Solvation, delphi_real& fEnergy_AnalySurf, int& iisitpot){


	bool bRadiusWarn = false;
	int ix, iy, iz;
	int i, j, ii, jj, qq;
	delphi_real fRadius, fEnergy_SelfReaction, fCost, fVal;
	delphi_real fEnergy_Temp=0, fEnergy_TotalCharge=0;
	delphi_real dist, fEnergy_Temp1, fEnergy_Temp2, fEnergy_Temp3, fFact, ptemp;
	delphi_real temp, temp1, temp2, temp3, spt1, fConstSixth;
	delphi_real dx, dy, dz;

	vector<delphi_real> spdiv(iTotalBdyGridNum), spot(iTotalBdyGridNum), schrg_omp(iTotalBdyGridNum);
	vector < SGridValue<delphi_real> > cgrid(iTotalBdyGridNum);
    
    SGrid<int> ixyz;
	
	fFact = 0.9549296586/(2.0*fScale*fEPKT);
	fConstSixth = 1.0/6.0;
    fEnergy_Temp1=0.0;
    fEnergy_SelfReaction=0.0;

// +++++++ Define these variables in this scope for OpenMP ++++++ //

    vector<delphi_real> prgfgSurfCrgA_nX(iTotalBdyGridNum);
    vector<delphi_real> prgfgSurfCrgA_nY(iTotalBdyGridNum);
    vector<delphi_real> prgfgSurfCrgA_nZ(iTotalBdyGridNum);
    vector<delphi_real> prgfgCrgPoseA_nX(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nY(iCrgGridNum);
    vector<delphi_real> prgfgCrgPoseA_nZ(iCrgGridNum);
    vector<delphi_real> prggvAtomicCrg_nValue(iCrgGridNum);
    
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    
	for(i=0;i<iTotalBdyGridNum;i++){

		ix = prgigBndyGrid[i].nX;
		iy = prgigBndyGrid[i].nY;
		iz = prgigBndyGrid[i].nZ;
        

		temp1 = prgfPhimap[(iz-1)*iGrid*iGrid+(iy-1)*iGrid+ix]+prgfPhimap[(iz-1)*iGrid*iGrid+(iy-1)*iGrid+(ix-2)];
		temp2 = prgfPhimap[(iz-1)*iGrid*iGrid+(iy)*iGrid+(ix-1)]+prgfPhimap[(iz-1)*iGrid*iGrid+(iy-2)*iGrid+(ix-1)];
		temp3 = prgfPhimap[(iz)*iGrid*iGrid+(iy-1)*iGrid+(ix-1)]+prgfPhimap[(iz-2)*iGrid*iGrid+(iy-1)*iGrid+(ix-1)];
		temp = prgfPhimap[(iz-1)*iGrid*iGrid+(iy-1)*iGrid+(ix-1)]-(temp1+temp2+temp3)*fConstSixth;

		spdiv[i] = temp;
        cgrid[i].nGrid.nX = ix;
        cgrid[i].nGrid.nY = iy;
        cgrid[i].nGrid.nZ = iz;
        cgrid[i].nValue = temp;
        

	}

	if(iCrgBdyGrid!=0){
		
 
/* The folloing Loop for variable cgrid moved to above and vector push_back has discarded for higher performance.

        SGridValue<delphi_real> fCGridValue_Temp;
		for(i=0;i<iTotalBdyGridNum;i++){
			fCGridValue_Temp.nGrid.nX = prgigBndyGrid[i].nX;
			fCGridValue_Temp.nGrid.nY = prgigBndyGrid[i].nY;
			fCGridValue_Temp.nGrid.nZ = prgigBndyGrid[i].nZ;
			fCGridValue_Temp.nValue = spdiv[i];
			cgrid.push_back(fCGridValue_Temp);
		}*/
		
		for(i=0;i<iCrgBdyGrid;i++){
			ix = prgdgvCrgBndyGrid[i].fgCoord.nX;
			iy = prgdgvCrgBndyGrid[i].fgCoord.nY;
			iz = prgdgvCrgBndyGrid[i].fgCoord.nZ;
			fVal = prgdgvCrgBndyGrid[i].fVal1;

			for(j=0;j<iTotalBdyGridNum;j++){
				if(cgrid[j].nGrid.nX == ix && cgrid[j].nGrid.nY == iy && cgrid[j].nGrid.nZ == iz)
                {
					cgrid[j].nValue = cgrid[j].nValue - fVal;
				}
			}
		}

		for(i=0;i<iTotalBdyGridNum;i++){
			spdiv[i] = cgrid[i].nValue;
		}

        vector < SGridValue<delphi_real> >().swap(cgrid);
	}

	fEnergy_Temp1=0.0;
    schrg.reserve(iTotalBdyGridNum);
    
	for(i=0;i<iTotalBdyGridNum;i++){
		temp=spdiv[i]*fFact;
		schrg[i] = temp;
		schrg_omp[i] = temp;
		fEnergy_Temp1=fEnergy_Temp1+temp;

		} // omp

	if(bReactFieldlnFRC || bSolvEng || bNonlinearEng || bSurfEngOut || bSurfCrgOut){


		if(bSolvEng || bNonlinearEng){

			fEnergy_Temp=0.0; fEnergy_Temp2=0.0; fEnergy_Temp3=0.0; fEnergy_TotalCharge=0.0; fCost=0.0;
			if(iMediaNum>0){
				for(i=0;i<iCrgGridNum;i++){
					fRadius = 1.0; // fortran code, radius = radpolext, which is removed
					ii = prgiCrgAt[i];
					if(ii>0 && ii<=iAtomNum){
						fRadius = prgapAtomPdb[ii-1].getRadius();
						fCost = 1.0/(fEPKT*prgfAtomEps[i])-1.0;

						if(prgfAtomEps[i]<=0){
							cout << " atmeps error  " << i << "  " << ii << endl;
						}

						if(fRadius<=0){
							cout << " charged atom number " << setw(4) << fixed << right << ii << "  radius changed from zero to " << 1.0 << endl;
							fRadius = 1.0;
							bRadiusWarn=true;
						}

					}

					fEnergy_SelfReaction = fEnergy_SelfReaction + 0.5*prggvAtomicCrg[i].nValue*prggvAtomicCrg[i].nValue*fCost/fRadius; // not very clear

				}

				if(bRadiusWarn){
					CReacFieldEnergyOverEST waring;
				}

				fEnergy_SelfReaction = fEnergy_SelfReaction * fEPKT;

				cout << " self-reaction field energy       :               " << setw(8) << right << fEnergy_SelfReaction << " kt" << endl;

				if(bEngOut){
					ofstream ofEnergyFile;
					ofEnergyFile.open(strEnergyFile,std::fstream::app);
					ofEnergyFile << "self-reaction field energy     :  " << fEnergy_SelfReaction << " kt \n";
				    ofEnergyFile.close();
				}
			}
                
          	// +++++++++++++ Next TWO LOOPS for OpenMP ++++++++++ //
                
                for(i=0;i<iTotalBdyGridNum;i++)
                {
                    prgfgSurfCrgA_nX[i] = prgfgSurfCrgA[i].nX;
                    prgfgSurfCrgA_nY[i] = prgfgSurfCrgA[i].nY;
                    prgfgSurfCrgA_nZ[i] = prgfgSurfCrgA[i].nZ;
                    
                }
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

			#pragma omp parallel shared(prgfgSurfCrgA_nX,prgfgSurfCrgA_nY,prgfgSurfCrgA_nZ,prgfgCrgPoseA_nX,prgfgCrgPoseA_nY,prgfgCrgPoseA_nZ,prggvAtomicCrg_nValue,spot,schrg_omp) private(j,dx,dy,dz,dist,ptemp)
			{
			
			#pragma omp for reduction( + : fEnergy_TotalCharge, fEnergy_Temp)
			
#endif

				for(i=0;i<iTotalBdyGridNum;i++){
					ptemp = 0.0;
	
					for(j=0;j<iCrgGridNum;j++){
                        dx = prgfgSurfCrgA_nX[i] - prgfgCrgPoseA_nX[j];
                        dy = prgfgSurfCrgA_nY[i] - prgfgCrgPoseA_nY[j];
                        dz = prgfgSurfCrgA_nZ[i] - prgfgCrgPoseA_nZ[j];
                        
                        dist = sqrt(dx*dx + dy*dy + dz*dz);
                        ptemp = ptemp + prggvAtomicCrg_nValue[j]/dist;
					} // vectorization
					spot[i] = ptemp;
					fEnergy_Temp += spot[i]*schrg_omp[i];
					fEnergy_TotalCharge += schrg_omp[i];
				}

#ifdef PARALLEL_OMP

			}

#endif

				fEnergy_Solvation = fEnergy_Temp*fEPKT/2.0;

				cout << " total s.charge,no epsin carrying :               " << setw(8) << right << fEnergy_TotalCharge << endl;
				if(iTotalBdyGridNum==0 || iGrid<5){
					cout << " Midpoints are out side the cube and delphi cannot determine the molecular surface." << endl;
					cout << " Please enlarge the gsize or decrease the perfil value." << endl;
				}
				cout << " corrected reaction field energy  :               " << setw(8) << right << fEnergy_Solvation << " kt" << endl;
				cout << " total reaction field energy      :               " << setw(8) << right << (fEnergy_SelfReaction+fEnergy_Solvation) <<" kt" << endl;

				ergr = fEnergy_SelfReaction+fEnergy_Solvation;

				if(bEngOut){
					ofstream ofEnergyFile;
					ofEnergyFile.open(strEnergyFile,std::fstream::app);
					ofEnergyFile << "corrected reaction field energy : " << fEnergy_Solvation << " kt\n";
					ofEnergyFile << "total reaction field energy     : " << (fEnergy_SelfReaction+fEnergy_Solvation) <<" kt\n";
					ofEnergyFile.close();
				}
		}

		if(bSurfCrgOut){
			cout << " writing surface charge file      :               " << strScrgFile << endl;
			ofstream ofScrgFile;
			ofScrgFile.open(strScrgFile);
			// iSurfCrgFormatOut -> scrgfrm; integer Format of file scrgnam, = 0 if unknown format, = 1 if �PDB�.

			if(iSurfCrgFormatOut==1 || iSurfCrgFormatOut==2 ){
				ofScrgFile << "DELPHI FORMAT PDB" << endl;
				ofScrgFile << "FORMAT NUMBER = " << iSurfCrgFormatOut << endl;
				ofScrgFile << "       bgp#  atom SC   res#      pos                               scrg           surf energy" << endl;

			}
/*
    eBuffz function has removed from C++ DelPhi
			if(bBuffz){
				lim_min = 2+ieBuffz.nMin;
				lim_max = iGrid-1-ieBuffz.nMax;
			}
 */
			for(i=0;i<iTotalBdyGridNum;i++){
					ixyz = prgigBndyGrid[i];

/*
    eBuffz function has removed from C++ DelPhi
                if(bBuffz){
						ido = 1;
						if(optORLT<int>(ixyz,lim_min) || optORGT<int>(ixyz,lim_max)){
							ido = 0;
						}
						if(ido==0) continue;
					}
*/
					fEnergy_Temp1 = schrg[i];
					spt1 = spot[i]*fEnergy_Temp1*fEPKT/2.0;

					if(iSurfCrgFormatOut==0){
						ofScrgFile << setw(8) << fixed << right << i+1 << "  " << prgfgSurfCrgA[i].nX << "  " << prgfgSurfCrgA[i].nY << "  " << prgfgSurfCrgA[i].nZ << "  " << fEnergy_Temp1 << endl;
					}

                    if(iSurfCrgFormatOut==1 || iSurfCrgFormatOut==2){
                  		jj = atsurf[i];	   //from surface construction class
						qq = atoi(prgapAtomPdb[jj-1].getAtInf().substr(11,4).c_str());
						ofScrgFile << "ATOM  " << setw(5) << fixed << right << i+1 << " " <<  setw(5) << fixed << right << jj << " SC " << setw(6) << fixed << right << qq << "      " << setw(5) << fixed << right << prgfgSurfCrgA[i].nX << "  " << setw(5) << fixed << right << prgfgSurfCrgA[i].nY << "  " << setw(5) << fixed << right << prgfgSurfCrgA[i].nZ << "  " << setw(5) << fixed << right << scientific << fEnergy_Temp1 << "  " << setw(5) << fixed << right << scientific << spt1 << endl;
					}
			}
		}

		if(bSurfEngOut){
			cout << " writing surface energy file      :               surfen.dat" << endl;
			ofstream surfen;
			surfen.open("surfen.dat");
			fEnergy_Temp2=0.0;

			for(i=0;i<iTotalBdyGridNum;i++)
			{
				
                fEnergy_Temp1 = 0.000;
            
				surfen << setw(8) << fixed << right << i+1 << "  " << prgfgSurfCrgA[i].nX << "  " << prgfgSurfCrgA[i].nY << "  " << prgfgSurfCrgA[i].nZ << "  " << fEnergy_Temp1 << endl;

			}

			surfen.close();
		}

	}
 
// +++++++++++++++++++++ Clean memory occupied by these vectors ++++++++++++++++ //
   
    vector<delphi_real>().swap(prgfgSurfCrgA_nX);
    vector<delphi_real>().swap(prgfgSurfCrgA_nY);
    vector<delphi_real>().swap(prgfgSurfCrgA_nZ);
    vector<delphi_real>().swap(prgfgCrgPoseA_nX);
    vector<delphi_real>().swap(prgfgCrgPoseA_nY);
    vector<delphi_real>().swap(prgfgCrgPoseA_nZ);
    vector<delphi_real>().swap(prggvAtomicCrg_nValue);
    
    vector<delphi_real>().swap(spdiv);
    vector<delphi_real>().swap(spot);
	vector<delphi_real>().swap(schrg_omp);

}

