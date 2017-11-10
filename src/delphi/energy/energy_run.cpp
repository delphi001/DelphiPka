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
    debug_energy=true;

    delphi_real fEnergy_Grid=0.0;
    delphi_real fEnergy_Solvation=0.0;
    delphi_real fEnergy_AnalySurf=0.0;
    delphi_real fEnergy_Nonlinear=0.0;
    delphi_real fEnergy_Coulombic=0.0;
    delphi_real fEnergy_AnalyGrid=0.0;
    delphi_real fEnergy_SolvToChgIn=0.0;
    delphi_real fEnergy_SolvToChgOut = 0.0;  // Solvent contribution to fixed charges  outside the cube is zero.
    delphi_real fEnergy_Total=0.0;

	if (inhomo == 1) fIonStrength = 0;

    SGrid<int> ixyz;

      infoString = " Info> ";
      timeString = " Time> ";
      enerString = " Energy> ";
      MAXWIDTH = 45;
      NUMWIDTH = 12;
      //ARGO: Remove (if) TRUE in the final version. Only for debugging here
      //      bool debug_energy=false;
      //

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
        cout << " WARNING !!! Analytic grid energy is no longer available." << endl;

        if(bEngOut)
        {
            ofstream ofEnergyFile;
            ofEnergyFile.open(strEnergyFile,std::fstream::app);
            ofEnergyFile << " Analytic grid energy is " << fEnergy_AnalyGrid << " kT.\n";
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

		if (debug_energy) cout << " iCrgedGridB: " << iCrgedGridB << endl;
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

              //  cout << "LinLi: in grid energy: " << setw(10)<< i << setw(10) << iw << setw(20) << fixed << setprecision(7) << prgfPhimap[iw] << setw(20) << prgfGridCrg[i] << endl;

            }


        }
		
        fEnergy_Grid = fEnergy_Grid / 2.0;

		if (debug_energy) cout << "fEnergy_Grid: " << fEnergy_Grid << endl;

        if(debug_energy) cout << " Gaus> iGaussian,inhomo,bSolvEng: " << iGaussian << " " << inhomo << " " << bSolvEng << endl;
        if(debug_energy) cout << " Conv> iConvolute,inhomo,bSolvEng: " << iConvolute << " " << inhomo << " " << bSolvEng << endl;

        //ARGO modification of IF condition
        if((iGaussian==1||iConvolute!=0)&&inhomo==1&&bSolvEng)
        {
            ergsgaussian=fEnergy_Grid;
            if(debug_energy) cout << "ergsgaussian: " << ergsgaussian << endl;
        }
        else
        {

            cout << enerString << left << setw(MAXWIDTH) << "Total grid energy" << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_Grid << " kT" << endl;
        }
        if(bEngOut)
        {
            ofstream ofEnergyFile;
            ofEnergyFile.open(strEnergyFile,std::fstream::app);
            ofEnergyFile << "Total grid energy		:   " << setw(NUMWIDTH) << setprecision(NUMPRECISION) << fixed << fEnergy_Grid << " kT \n";
            ofEnergyFile.close();
        }

    }

    //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    //ARGO modification of IF condition
    if((iGaussian==1||iConvolute!=0)&&inhomo==0&&bSolvEng)
    {
        //write(6,"(a,f20.4,a)"),        ' corrected reaction field energy :',ergg-ergsgaussian,' kT'
        cout << enerString << left << setw(MAXWIDTH) << "Corrected reaction field energy" << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_Grid-ergsgaussian << " kT" << endl;
        //ergs=fEnergy_Grid-ergsgaussian;
        ergs=fEnergy_Grid-ergsgaussian;
        if (debug_energy) cout << "Gaus> iGaussian,ergs = " << iGaussian << " " << ergs << endl;
        if (debug_energy) cout << "Conv> iConvolute,ergs = " << iConvolute << " " << ergs << endl;
    }


    //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    //ARGO modification of IF-condition
    if(iGaussian==0 && iConvolute==0)   //LinLi,Argo :if iGaussian==1 || iConvolute!=0, skip other energy terms
    {
        //goto 1212
        if(bGridEng && bAnalyEng)
        {
			#ifdef VERBOSE
            cout << " Difference energy, in kT, is  " << fEnergy_Grid-fEnergy_AnalyGrid << endl;
            cout << " Difference energy, in kcals, is  " << (fEnergy_Grid-fEnergy_AnalyGrid)*0.6 << endl;
			#endif
        }

        //  +++++++++++++++ For polarized membrane (not supported and not tested) +++++ //


        if(iBndyType == 5)
        {
            cout << " WARNING!!! Not completely tested routine for polarized membrane!!" << endl;
            exit (EXIT_FAILURE);

            if(bNonlinearEng || bAnalySurfEng)
            {
				#ifdef VERBOSE
                cout << " WARNING !!! This option is not yet working with fixed potential difference!" << endl;
				#endif

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
            		cout << "Energy contribution from voltage drop = " << enface << " kT" << endl;
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
            fEnergy_Nonlinear = 0.0;
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
    //ARGO modification of the IF-condition
    if( bCoulombEng && ( !bIonsEng || !bNonlinearEng ) && !((iGaussian==1 || iConvolute!=0) && inhomo==1 ))
    {

        fEnergy_Coulombic = 0.0;

        if(bIonsEng)
        {

            fEnergy_SolvToChgIn = 0.0;

            energy_clbtotal(fEnergy_SolvToChgIn, fEnergy_Coulombic);	// call clbtotal function.

            {
                CIonicCalcUnderTest waring;
            }
#ifdef VERBOSE
            cout << enerString << "Solvent contribution to fixed charges" << endl;

            cout << enerString << left << setw(MAXWIDTH) << "Respectively inside and outside the cube" << " : "
            << setw(NUMWIDTH) << right << fEnergy_SolvToChgIn << "  kT   " << setprecision(NUMPRECISION) << fixed << fEnergy_SolvToChgOut << "  kT" << endl;  // where is ergestout??

            cout << enerString << left << setw(MAXWIDTH) << "Total ionic direct contribution" << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << (fEnergy_SolvToChgIn+fEnergy_SolvToChgOut) << "  kT" << endl;
#endif
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

        cout << enerString << left << setw(MAXWIDTH) << "Coulombic energy" << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT" << endl;

        if(bEngOut)
        {

            ofstream ofEnergyFile;

            ofEnergyFile.open(strEnergyFile,std::fstream::app);

            ofEnergyFile << "Total coulombic energy  :  " << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT \n";

            ofEnergyFile.close();
        }
    }

    //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    if(bNonlinearEng)
    {
        energy_nonl(fEnergy_Nonlinear, iGridOutput);	// call nonlinear function.



        if(bIonsEng)
        {
			fEnergy_Coulombic = 0.0;
			fEnergy_SolvToChgIn = 0.0;

            energy_clbnonl(fEnergy_Coulombic, fEnergy_SolvToChgIn, iGridOutput);

            cout << enerString << left << setw(MAXWIDTH) <<"Direct ionic contribution inside the box" << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_SolvToChgIn << " kT" << endl;

            cout << enerString << left << setw(MAXWIDTH) <<"Coulombic energy" << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT" << endl;

            if(bEngOut)
            {
                ofstream ofEnergyFile;

                ofEnergyFile.open(strEnergyFile,std::fstream::app);

                ofEnergyFile << "Total coulombic energy :  " << setw(NUMWIDTH) << setprecision(NUMPRECISION) << fixed << fEnergy_Coulombic << " kT \n";

                ofEnergyFile.close();
            }

        }

        vector < SGridValue<delphi_real> >().swap(sout);
    }

    //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    //Argo modification of the IF-condition
    if(iGaussian==0 && iConvolute==0) //if iGaussian==1 || iConvolute != 0,skip these terms
    {
        //goto 1213

        if(bSolvEng && bIonsEng)
        {

            cout << " Energy arising from solvent and boundary pol.  " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << (fEnergy_Nonlinear+fEnergy_Solvation+fEnergy_SolvToChgIn+fEnergy_SolvToChgOut) << " kT" << "\n";
        }

        //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //


        if(bNonlinearEng && bGridEng)
        {

            cout << enerString << left << setw(MAXWIDTH) << "Total non linear grid energy" << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << (fEnergy_Grid+fEnergy_Nonlinear) << " kT \n";
        }

        //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        fEnergy_Total = fEnergy_Nonlinear + fEnergy_Coulombic + fEnergy_Solvation + fEnergy_SolvToChgIn + fEnergy_SolvToChgOut;

        if(bSolvEng || bCoulombEng)
        {
            cout << enerString << left << setw(MAXWIDTH) << "All required energy terms but grid energy    " << " : " << setw(NUMWIDTH) << right << setprecision(NUMPRECISION) << fixed << fEnergy_Total << " kT" << endl;

            if(bEngOut)
            {
                ofstream ofEnergyFile;

                ofEnergyFile.open(strEnergyFile,std::fstream::app);

                ofEnergyFile << "Total required energy (everything calculated but grid and self_reaction energies: " << setw(NUMWIDTH) << setprecision(NUMPRECISION) << fixed << fEnergy_Total << " kT \n";

                ofEnergyFile.close();
            }
        }

        //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        if(bAnalySurfEng && bAnalyEng && bGridEng)
        {

            cout << enerString << left << setw(MAXWIDTH) << "Excess grid energy" << " : " << setw(NUMWIDTH) << right  << setprecision(NUMPRECISION) << fixed << (fEnergy_Grid - fEnergy_AnalySurf - fEnergy_AnalyGrid) << endl;
        }

        //  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

        cout << endl;

        //  ++++++++++++++++++ Write into data container ++++++++++++++++++++ //

        ergg    = fEnergy_Grid;
        ergc    = fEnergy_Coulombic;
        ergs    = fEnergy_Solvation;
        ergions = fEnergy_SolvToChgIn + fEnergy_SolvToChgOut;

    } //1213
    else if (iGaussian==1 || iConvolute!=0)
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
