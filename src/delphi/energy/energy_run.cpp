/*
 * energy_run.cpp
 *
 *  Created on: Feb 23, 2014
 *      Author: Lin Wang
 */

#include "energy.h"

#ifdef PARALLEL_OMP
#include <omp.h>
#endif


void CDelphiEnergy::run()
{
    bool ido;
    int i, iw, iGridOutput, iisitpot;
    debug_energy=false;

    delphi_real fEnergy_Grid=0.0;
    delphi_real fEnergy_Solvation=0.0;
    delphi_real fEnergy_AnalySurf=0.0;
    delphi_real fEnergy_Nonlinear=0.0;
    delphi_real fEnergy_Coulombic=0.0;
    delphi_real fEnergy_AnalyGrid=0.0;
    delphi_real fEnergy_SolvToChgIn=0.0;
    delphi_real fEnergy_SolvToChgOut = 0.0;  // Solvent contribution to fixed charges  outside the cube is zero.
    delphi_real fEnergy_Total=0.0;

    SGrid<int> ixyz;

    /*
    #ifdef PARALLEL_OMP

    	// +++++++++++++++++++++ Setup OpenMP  +++++++++++++++++++++ //

        // User should specify how many threads are going to set for OpenMP by "export OMP_NUM_THREADS=Num"
        // Otherwise, the program will automatically set the max of threads available
        // To get the best performance, set the threads equal to the number of cores on the node
    	// Get the number of max threads in this system

    	int iProcs, iThreads;
    	iProcs = omp_get_num_procs();
    	iThreads = omp_get_max_threads();
    	if(iThreads<iProcs)	CThreadsLessThanProcs waring(iProcs, iThreads);

    	cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
    	cout << " Energy Calculation with OpenMP running in " << iThreads << " Threads" << endl;
    	cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
    	cout << endl;

    	// Now set the number of threads
    	omp_set_num_threads(iThreads);

    	// Now set the OpenMP Timer
    	double omp_timer;
    	omp_timer = omp_get_wtime();

    #endif
    */
    // +++++++++++++++++++++ Analytic grid energy +++++++++++++++++++++ //

    if(bEngOut) 		// output energy calculation results "energy.dat"
    {
        ofstream ofEnergyFile;
        ofEnergyFile.open(strEnergyFile);
        ofEnergyFile.close();
    }

    if(bAnalyEng)
    {
        fEnergy_AnalyGrid=0.0;
        cout << " analytic grid energy is no longer available." << endl;


        if(bEngOut)
        {
            ofstream ofEnergyFile;
            ofEnergyFile.open(strEnergyFile,std::fstream::app);
            ofEnergyFile << " analytic grid energy is " << fEnergy_AnalyGrid << " kt.\n";
            ofEnergyFile.close();
        }
        exit (EXIT_FAILURE);
    }

    //  +++++++++++++++ Total grid energy ++++++++++++++++++++++++++++ //

    if(bGridEng)
    {
        fEnergy_Grid=0.0;
        lim_min = 2+ieBuffz.nMin;  //ieBuffz(pdc->getKey_constRef< SExtrema<int> >("bufz"))
        lim_max = iGrid-1-ieBuffz.nMax;	 //lim%min=2+bufz%min ; lim%max=igrid-1-bufz%max


        for(i=0; i<iCrgedGridB; i++)
        {
            ixyz = prgigGridCrgPose[i];


            ido = 1;
            if(optORLT<int>(ixyz,lim_min) || optORGT<int>(ixyz,lim_max))
            {
                ido = 0;
            }
            if(ido)
            {
                iw = ((ixyz.nZ-1)*iGrid*iGrid)+(ixyz.nY-1)*iGrid+(ixyz.nX-1);  // i, j, k -> k, j, i (z, y, x)

                fEnergy_Grid = fEnergy_Grid + prgfPhimap[iw]*prgfGridCrg[i];

                //cout << "LinLi: in grid energy: " << setw(10)<< i << setw(20) << fixed << setprecision(7) << prgfPhimap[iw] << setw(20) << prgfGridCrg[i] << endl;

            }


        }

        fEnergy_Grid = fEnergy_Grid / 2.0;

        if(debug_energy)cout << "iGaussian,inhomo,bSolvEng: " << iGaussian << " " << inhomo << " " << bSolvEng << endl;
        if(iGaussian==1&&inhomo==1&&bSolvEng)
        {
            ergsgaussian=fEnergy_Grid;
            if(debug_energy)cout << "ergsgaussian: " << ergsgaussian << endl;
        }
        else
        {

            cout << endl << " total grid energy                :               " << setw(8) << right << fEnergy_Grid << " kt" << endl;
        }
        if(bEngOut)
        {
            ofstream ofEnergyFile;
            ofEnergyFile.open(strEnergyFile,std::fstream::app);
            ofEnergyFile << "total grid energy		:   " << setw(8) << fEnergy_Grid << " kt \n";
            ofEnergyFile.close();
        }

    }

    //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    if(iGaussian==1&&inhomo==0&&bSolvEng)
    {
        //write(6,"(a,f20.4,a)"),        ' corrected reaction field energy :',ergg-ergsgaussian,' kt'
        cout << " corrected reaction field energy  :               " << setw(8) << right << fEnergy_Grid-ergsgaussian << " kt" << endl;
        //ergs=fEnergy_Grid-ergsgaussian;
        ergs=fEnergy_Grid-ergsgaussian;
    }


    //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    if(iGaussian==0)   //LinLi:if iGaussian==1, skip other energy terms
    {
        //goto 1212
        if(bGridEng && bAnalyEng)
        {
            cout << " difference energy, in kt, is  " << fEnergy_Grid-fEnergy_AnalyGrid << endl;
            cout << " difference energy, in kcals, is  " << (fEnergy_Grid-fEnergy_AnalyGrid)*0.6 << endl;
        }

        //  +++++++++++++++ For polarized membrane (not supported and not tested) +++++ //


        if(iBndyType == 5)
        {
            cout << " WARNING!!!Not completely tested routine for polarized membrane!!" << endl;
            exit (EXIT_FAILURE);

            if(bNonlinearEng || bAnalySurfEng)
            {
                cout << "This option is not yet working with fixed potential difference!" << endl;

                fEnergy_AnalySurf=0.0;
                fEnergy_Nonlinear=0.0;

            }

            /*
                Potential drop polarized membrane has removed from C++ DelPhi

            		deltaphi=0.0;
            		zfieldown=0.0;
            		zfieldup=0.0;

            		ofstream fields;
            		fields.open("fields.txt");
            		for(ix=1; ix<iGrid-1; ix++){
            			for(iy=1; iy<iGrid-1;iy++){
            				phiup = vector<delphi_real> prgfPhimap(ix,iy,iGrid);
            				phidown = vector<delphi_real> prgfPhimap(ix,iy,0);
            				deltaphi = deltaphi + (phiup - vector<delphi_real> prgfPhimap(ix,iy,iGrid-2));
            				zfieldup = zfieldup - (phiup - vector<delphi_real> prgfPhimap(ix,iy,iGrid-2))*fScale;
            //			averaging over upper and lower cube faces in order to cancel numerical errors.	//
            				deltaphi = deltaphi - (phidown - vector<delphi_real> prgfPhimap(ix,iy,1));
            				zfieldown = zfieldown + (phidown - vector<delphi_real> prgfPhimap(ix,iy,1))*fScale;
            			}
            		}

            		enface = deltaphi * fgPotentialDrop.nZ * fEpsOut; * ((iGrid-1.0)/((iGrid-2.0)*(iGrid-2.0))/(4.0*fPi*fScale));// unclear operation here!!
            		cout << "Energy contribution from voltage drop = " << enface << " kt" << endl;
            		fields.close();

            		ofstream potcen;
            		potcen.open("potcen.txt");
            		cout << "fieldup medio: " << zfieldup/((iGrid - 2)*(iGrid - 2))<< endl;
            		cout << "fieldown medio: " << zfieldown/((iGrid - 2)*(iGrid - 2)) << endl;

            		for(iz=0;iz<=iGrid;iz++){

            			iw = (iGrid+1)/2*iGrid*iGrid + (iGrid+1)/2*iGrid + iz;
            			potcen << iz << "   " << prgfPhimap[iw] << endl;
            		}
            		potcen.close();
            */
        }


        //  ++++++++++++++++++ Reaction field energy calculation ++++++++++++++++++++++++++ //


        if(bReactFieldlnFRC || bSolvEng || bNonlinearEng || bAnalySurfEng || bSurfEngOut || bSurfCrgOut)
        {
            // fEnergy_SolvToChgIn = interaction energy of the solvent and the fixed charges. //

            fEnergy_Solvation = 0.0;
            fEnergy_AnalySurf = 0.0;
            fEnergy_Nonlinear=0.0;
            fEnergy_SolvToChgIn = 0.0;

            if(bPotentiallnSite)
            {
                iisitpot = 1;
            }
            else
            {
                iisitpot = 0;
            }
            energy_react(fEnergy_Solvation, fEnergy_AnalySurf, iisitpot);  // call reaction field function
        }

        //  ++++++++++++++++++ Coulombic energy calculation +++++++++++++++++++++++++++++++ //
    }//1212


    //if( bCoulombEng && ( !bIonsEng || !bNonlinearEng ) )
    if( bCoulombEng && ( !bIonsEng || !bNonlinearEng ) && !(iGaussian==1 && inhomo==1 ))
    {

        fEnergy_Coulombic = 0.0;

        if(bIonsEng)
        {

            fEnergy_SolvToChgIn = 0.0;

            energy_clbtotal(fEnergy_SolvToChgIn, fEnergy_Coulombic);	// call clbtotal function.

            {
                CIonicCalcUnderTest waring;
            }

            cout << " solvent contribution to fixed charges" << endl;

            cout << " respectively inside and outside the cube :               ";

            cout << setw(8) << right << fEnergy_SolvToChgIn << "  kt   " << fEnergy_SolvToChgOut << "  kt" << endl;  // where is ergestout??

            cout << " total ionic direct contribution  :               " << setw(8) << right << (fEnergy_SolvToChgIn+fEnergy_SolvToChgOut) << "  kt" << endl;

        }

        else
        {

            if(iMediaNum == 1)
            {

                energy_clb(fEnergy_Coulombic);	// call clb function.

                fEnergy_Coulombic = fEnergy_Coulombic / fEpsin;
            }

            else
            {

                energy_clbmedia(fEnergy_Coulombic);	// call clbmedia function.

            }
        }

        cout << " coulombic energy                 :               " << setw(8) << right << fEnergy_Coulombic << " kt" << endl;

        if(bEngOut)
        {

            ofstream ofEnergyFile;

            ofEnergyFile.open(strEnergyFile,std::fstream::app);

            ofEnergyFile << "total coulombic energy  :  " << fEnergy_Coulombic << " kt \n";

            ofEnergyFile.close();
        }
    }

    //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    if(bNonlinearEng)
    {
        energy_nonl(fEnergy_Nonlinear, iGridOutput);	// call nonlinear function.

        fEnergy_Coulombic = 0.0;
        fEnergy_SolvToChgIn = 0.0;

        if(bIonsEng)
        {

            energy_clbnonl(fEnergy_Coulombic, fEnergy_SolvToChgIn, iGridOutput);

            cout << " direct ionic contribution inside the box :       " << setw(8) << right << fEnergy_SolvToChgIn << " kt" << endl;

            cout << " coulombic energy                 :               " << setw(8) << right << fEnergy_Coulombic << " kt" << endl;

            if(bEngOut)
            {
                ofstream ofEnergyFile;

                ofEnergyFile.open(strEnergyFile,std::fstream::app);

                ofEnergyFile << "total coulombic energy :  " << setw(8) << fEnergy_Coulombic << " kt \n";

                ofEnergyFile.close();
            }

        }

        vector < SGridValue<delphi_real> >().swap(sout);
    }

    //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    if(iGaussian==0) //if iGaussian==1,skip these terms
    {
        //goto 1213

        if(bSolvEng && bIonsEng)
        {

            cout << " Energy arising from solvent and boundary pol.  " << setw(8) << right << (fEnergy_Nonlinear+fEnergy_Solvation+fEnergy_SolvToChgIn+fEnergy_SolvToChgOut) << " kt" << "\n";
        }

        //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //


        if(bNonlinearEng && bGridEng)
        {

            cout << " total non linear grid energy	  :               " << setw(8) << right << (fEnergy_Grid+fEnergy_Nonlinear) << " kt \n";
        }

        //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        fEnergy_Total = fEnergy_Nonlinear + fEnergy_Coulombic + fEnergy_Solvation + fEnergy_SolvToChgIn + fEnergy_SolvToChgOut;

        if(bSolvEng || bCoulombEng)
        {
            cout << " all required energy terms but grid and self_react.:  " << setw(8) << right << fEnergy_Total << " kt" << endl;

            if(bEngOut)
            {
                ofstream ofEnergyFile;

                ofEnergyFile.open(strEnergyFile,std::fstream::app);

                ofEnergyFile << "total required energy (everything calculated but grid and self_reaction energies: " << setw(8) << fEnergy_Total << " kt \n";

                ofEnergyFile.close();
            }
        }

        //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        if(bAnalySurfEng && bAnalyEng && bGridEng)
        {

            cout << " excess grid energy =               " << setw(8) << right  << (fEnergy_Grid - fEnergy_AnalySurf - fEnergy_AnalyGrid) << endl;
        }

        //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        cout << endl;

        //  ++++++++++++++++++ Write into data container ++++++++++++++++++++ //

        ergg    = fEnergy_Grid;
        ergc    = fEnergy_Coulombic;
        ergs    = fEnergy_Solvation;
        ergions = fEnergy_SolvToChgIn + fEnergy_SolvToChgOut;

    } //1213
    else
    {
        ergg    = fEnergy_Grid;
        ergc    = fEnergy_Coulombic;
        ergs    = fEnergy_Solvation;

        //cout << "ergg: " << ergg << " ergc: " << ergc << " ergs: " << ergc << endl;
    }


    /*
    #ifdef PARALLEL_OMP

    	cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
    	cout << " Energy calculation with OpenMP finishes in  " << (omp_get_wtime() - omp_timer) << " sec" << endl;
      	cout << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;

    #endif
    */

}
