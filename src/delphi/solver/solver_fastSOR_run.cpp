/*
 * solver_fastSOR_run.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

void CDelphiFastSOR::run()
{
    debug_solver=false;
    const char * infoString = " Info> ";
    const char * timeString = " Time> ";
    const char * bndcString = " Bndc> ";
    size_t MAXWIDTH = 45;

    if(debug_solver) cout << "###### start run of solver ### " << endl;
	if (debug_solver) cout << "###: fsrfdens: " << fsrfdens << endl;
	if (debug_solver) cout << "###: inhomo: " << inhomo << endl;

    //ARGO modification made to accommodate iConvlute module

    if (inhomo == 1) fIonStrength=0;

    if( (iGaussian==1 && inhomo==0 ) ||  (iConvolute != 0 && inhomo==0) ) //reset some vectors for 2nd run of Gaussian or convolution
    {
        if (prgfGridCrg.size()>0 ) {
            if(debug_solver)cout << "prgfGridCrg.size(): " << prgfGridCrg.size() << endl;
            vector<delphi_real>().swap(prgfGridCrg) ;//gchrg
        }
        if (prgigGridCrgPose.size()>0 ){
            if(debug_solver)cout << "prgigGridCrgPose.size(): " << prgigGridCrgPose.size() << endl;
            vector< SGrid<delphi_integer> >().swap(prgigGridCrgPose) ; //gchrgp
        }

      //vector<delphi_integer> ibndx,ibndy,ibndz;
        if (ibndx.size()>0 ){
            if(debug_solver)cout << "ibndx.size(): " << ibndx.size() << endl;
            vector<delphi_integer>().swap(ibndx) ; //ibndx
        }
        if (ibndy.size()>0 ){
            if(debug_solver)cout << "ibndy.size(): " << ibndy.size() << endl;
            vector<delphi_integer>().swap(ibndy) ; //ibndy
        }
        if (ibndz.size()>0 ){
            if(debug_solver)cout << "ibndz.size(): " << ibndz.size() << endl;
            vector<delphi_integer>().swap(ibndz) ; //ibndz
        }
        if (phimap1.size()>0 ){
            if(debug_solver)cout << "phimap1.size(): " << phimap1.size() << endl;
            vector<delphi_real>().swap(phimap1) ; //phimap1
        }
        if (phimap2.size()>0 ){
            if(debug_solver)cout << "phimap2.size(): " << phimap2.size() << endl;
            vector<delphi_real>().swap(phimap2) ; //phimap2
        }

        if (prgiBndyDielecIndex.size()>0 ){
            if(debug_solver)cout << "prgiBndyDielecIndex.size(): " << prgiBndyDielecIndex.size() << endl;
            vector<delphi_integer>().swap(prgiBndyDielecIndex) ; //prgiBndyDielecIndex
        }
        if (prgfBndyDielec.size()>0 ){
            if(debug_solver)cout << "prgfBndyDielec.size(): " << prgfBndyDielec.size() << endl;
            vector< vector<delphi_real> >().swap(prgfBndyDielec) ; //prgfBndyDielec
        }
              //vector< SDoubleGridValue >& prgdgvCrgBndyGrid;   // cgbp(ibc)
        if (prgdgvCrgBndyGrid.size()>0 ){
            if(debug_solver)cout << "prgdgvCrgBndyGrid.size(): " << prgdgvCrgBndyGrid.size() << endl;
            vector< SDoubleGridValue >().swap(prgdgvCrgBndyGrid) ; //prgdgvCrgBndyGrid
        }
        if(debug_solver)cout << "prgfPhiMap.size(): " << prgfPhiMap.size() << endl;

		if (gaussianBoundaryDensity.size()>0) {
			if (debug_solver)cout << "gaussianBoundaryDensity.size(): " << gaussianBoundaryDensity.size() << endl;
			vector< delphi_real >().swap(gaussianBoundaryDensity); //prgdgvCrgBndyGrid
		}

        iDielecBndyOdd=0;

    }
    if(debug_solver)cout << "iDielecBndyOdd: " << iDielecBndyOdd << endl;

    validateInput();

    setDielecBndySaltMap();

    setCrg();
#ifdef VERBOSE
    cout << timeString << "iepsmp to db, and charging done at ";
    pTimer->showElapse();
    cout << endl;
    cout << infoString << left << setw(MAXWIDTH) << "Number of grid points assigned charge" << " : " << iCrgedGridSum << endl;
#endif

    if (bEpsOut && iGaussian==0 && iConvolute==0) //----- write dielectric map
    {
        unique_ptr<CIO> pio(new CIO()); // smart unique_ptr
        // pio->writeEpsMap(iAtomNum,iObjectNum,iGrid,fScale,fgBoxCenter,prgigEpsMap,prgbDielecMap,strEpsFile);
        pio->writeHomoEpsMap(iGrid,repsout, repsin,fScale,fgBoxCenter,prgbDielecMap,strEpsFile);
        pio.reset();
    }

    if (bEpsOut && iGaussian==1 && iConvolute==0) //----- write dielectric map
    {
        unique_ptr<CIO> pio(new CIO()); // smart unique_ptr
        // pio->writeEpsMap(iAtomNum,iObjectNum,iGrid,fScale,fgBoxCenter,prgigEpsMap,prgbDielecMap,strEpsFile);
        pio->writeGaussEpsMap(iGrid,repsout, fScale,fgBoxCenter,gepsmp2,strEpsFile);
        pio.reset();
    }

    //----- the order of calculateRelaxFactor and setBndy cannnot be reverted! prgfPhiMap is used
    //      as temporary array in calculateRelaxFactor...
    if (bSpectralRadius)
    {
        fSpec = fSpectralRadius;
        cout << infoString << left << setw(MAXWIDTH) << "Using entered value for relaxation of" << " : " << fSpec << endl;
    }
    else
    {
        fSpec = calculateRelaxFactor();
    }

    int noit = (int)(7.8/log(1.0+sqrt(1-fSpec)));

    cout << left << setw(MAXWIDTH) << " Estimated iterations to convergence" << " : " << noit << endl;

    if (bAutoConverge) iLinIterateNum = noit;

    setBndy();
#ifdef VERBOSE
    cout << endl;
    cout << timeString << "Setup time was ";
    pTimer->showElapse();
    cout << endl;
#endif
    if (0.3 < (delphi_real)iCrgBndyGridNum/iBndyGridNum && bAutoConverge)
    {
        iLinIterateNum = iLinIterateNum*iCrgBndyGridNum/(0.3*iBndyGridNum);
        cout << left << setw(MAXWIDTH) << " Re-estimated iterations now " << " : " << iLinIterateNum << endl;
    }

#ifdef VERBOSE
    cout << timeString << "Now iterating on ";
    pTimer->showTime();
    cout << endl;
#endif

    if (0 == iNonIterateNum || fZero > fIonStrength || inhomo == 1)
    {
        if (0 >= iLinIterateNum) throw CZeorLinIterateNum(bAutoConverge,iLinIterateNum);
        itit();
    }
    else
    {
        if (50 < noit || 0 >= iLinIterateNum)
        {
            iLinIterateNum = noit/2;
            cout << left << setw(MAXWIDTH) << " Re-estimated iterations now " << " : " << iLinIterateNum << endl;
        }

        delphi_real qfact = abs(fNetCrg)*(delphi_real)iCrgBndyGridNum/iBndyGridNum;

        nitit(qfact);
    }
}
