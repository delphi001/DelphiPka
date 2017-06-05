//
//  energy.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unordered_set>
#include <cstdio>
#include "energy.h"



using namespace std;

extern void runDelphi(shared_ptr<SPrime> param);


void CEnergy :: run()
{
    int i,j;
    bool bNeutral;
    int    currentResNum;
    string currentResNam;
    strFRCOutFile = "frclog";
    if(bHETATMinPQR) {
        vecGridPotential.resize(newPDB.size()+vecCrgHETATM.size());
    }
    else {
        vecGridPotential.resize(newPDB.size());
    }
    
    vecCrgStat.resize(newPDB.size());
    
    if(!bHETATMinPQR) {
        PQR_to_DelPhi.resize(newPDB.size());
    }
    else {
        PQR_to_DelPhi.resize(newPDB.size() + strHETATM.size());
    }
    
    float ergpolar_crg, ergpolar_neu;
    float ergg_prot, ergg_wat;
    
    
    
    
    /////// A smart pointer to create the parameters for DelPhi fun ////
    
    shared_ptr<SPrime> param ( new SPrime );
    
    param->pdbfile     = pdb_input;
    param->pdbformat   = CommPDB;
    param->indi        = indi;
    param->exdi        = exdi;
    param->prbrad      = prbrad;
    param->ionrad      = ionrad;
    param->fMaxc       = maxc;
    param->salt        = 0.15;
    param->bCommFrc    = true;
    param->bSolvEng    = false;
    
    
    
    ///////////////////////////////////////////////////////////////////
    
    
    ionizableArr(); // This function is to genereate a 2d array to store charged residue and their atoms.
    
    
    vecIonRes.erase(vecIonRes.begin()); // vecIonRes[0][n] array is EMPTY. Therefore, remove from the vector
    
    
#ifdef PRIME_DEBUG   // print out the 2d array of ionizable residue and their atom number (index) in the newpdb vector
    
    for(int j=0;j<vecIonRes.size();j++)
    {
        for(int k=0;k<vecIonRes[j].size();k++)
        {
            cout << vecIonRes[j][k] << " ";
        }
        cout << endl;
    }
    
#endif
    
    EnergyPolarCrg.resize(vecIonRes.size());
    EnergyPolarNeu.resize(vecIonRes.size());
    EnergyRxnCrg.resize(vecIonRes.size());
    EnergyRxnNeu.resize(vecIonRes.size());
    
    EnergyPair.resize(vecIonRes.size()*2);
    for(i=0;i<EnergyPair.size();i++) { EnergyPair[i].resize(vecIonRes.size()*2); }
    
    
#ifdef MPI_PARALLEL
    
    int globalSIZE = 4*vecIonRes.size()*vecIonRes.size();
    globalEnergyVec.resize(globalSIZE);
    
    for(i=0; i<vecIonRes.size(); i++) {
        for(j=0; j<vecIonRes.size(); j++) {
            EnergyPair[2*i][2*j].resnum1 = newPDB[vecIonRes[i][0]].res_num;
            EnergyPair[2*i][2*j].resnum2 = newPDB[vecIonRes[j][0]].res_num;
            
            EnergyPair[2*i][2*j + 1].resnum1 = newPDB[vecIonRes[i][0]].res_num;
            EnergyPair[2*i][2*j + 1].resnum2 = newPDB[vecIonRes[j][0]].res_num;
            
            EnergyPair[2*i + 1][2*j].resnum1 = newPDB[vecIonRes[i][0]].res_num;
            EnergyPair[2*i + 1][2*j].resnum2 = newPDB[vecIonRes[j][0]].res_num;
            
            EnergyPair[2*i + 1][2*j + 1].resnum1 = newPDB[vecIonRes[i][0]].res_num;
            EnergyPair[2*i + 1][2*j + 1].resnum2 = newPDB[vecIonRes[j][0]].res_num;
        }
    }
    
    
    int N = vecIonRes.size();
    
    id = MPI::COMM_WORLD.Get_rank();
    numprocs = MPI::COMM_WORLD.Get_size();
    
    int rank_size = N / numprocs;
    MPI_START = ( N / numprocs) * id;
    
    MPI::COMM_WORLD.Barrier();
    
    double MPI_Timer = MPI_Wtime();
    
    float EnergyRxnCrg_partial[rank_size+1];
    float EnergyRxnNeu_partial[rank_size+1];
    float EnergyPolarCrg_partial[rank_size+1];
    float EnergyPolarNeu_partial[rank_size+1];
    float globalEnergyVec_partial[(rank_size+1) * 4 * N];
    
    
    
    if (N % numprocs > id){
        MPI_START += id;
        MPI_END = MPI_START + ( N / numprocs) + 1;
    }
    else{
        MPI_START += N % numprocs;
        MPI_END = MPI_START + ( N / numprocs);
    }
    
#endif
    

#ifdef MPI_PARALLEL
    if(id==0) {
#endif
        cout << endl;
        cout << "Energy calculation with DelPhi starting... " << endl;
        cout << " progress details can be found in delphi_progress.txt file during calculation. " << endl;
        cout << " progress log files will be deleted afterwards. " << endl;


#ifdef MPI_PARALLEL
    }
#endif
    
    ofstream delphi_progress;
    delphi_progress.open("delphi_progress.txt");
    
    
#ifdef MPI_PARALLEL
    for(int k=MPI_START; k<MPI_END;k++) {
#endif
        
#ifndef MPI_PARALLEL
        
        clock_t begin_time = clock();

        for(int k=0; k<vecIonRes.size(); k++) {
#endif
            ResetCrgStat();
            
            for(j=0;j<vecIonRes[k].size();j++) {     ///  Charging residue Loop Start here
                
                vecCrgStat[vecIonRes[k][j]] = true;
                
                currentResNum = newPDB[vecIonRes[k][j]].res_num;
                currentResNam = newPDB[vecIonRes[k][j]].res_name;
                
                if (newPDB[vecIonRes[k][j]].atom_name == "CG" || newPDB[vecIonRes[k][j]].atom_name == "C1'") {
                    center1.nX = newPDB[vecIonRes[k][j]].coord.X;
                    center1.nY = newPDB[vecIonRes[k][j]].coord.Y;
                    center1.nZ = newPDB[vecIonRes[k][j]].coord.Z;
                }
            }
            
            bNeutral = false;
            
            charge(bNeutral);  //  This function is to charge the corresponding residue (all atoms belong to that residue)
            
            
            /////  Create a singleRes_to_DelPhi intermediate array to DelPhi for solvation energy in water //////
            
            for(i=0;i<newPDB.size();i++) {
                if (   newPDB[i].res_name == currentResNam
                    && newPDB[i].res_num  == currentResNum
                    && newPDB[i].conf != "BK") {
                    singleRes_to_DelPhi.push_back(PQR_to_DelPhi[i]);
                }
            }
            
            
            //////////////////////////////
            /*
             Starting the following, DelPhi starts and will have a base run, then focusing with scale=2 and scale=4.
             If molecule is small enough that when scale=2 the grid size is larger than base run, scale=2 will be skipped.
             If too small, scale=4 run will be skipped.
             */
            /////////////////////////////
            
            delphi_progress << "    Now Starting the DelPhi focusing with charged residue ";
            delphi_progress << currentResNam  << setw(4) << currentResNum << endl;
            
            param->ntimes = k;
            
            setbase(param, center1);
            
            runDelphi(param);
            
            vecGridPotential = param->vecGridPot;

            
            
            if ( runFocus(param, center1) )
            {
                
                param->bAcent = true;
                param->center[0] = center1.nX; param->center[1] = center1.nY; param->center[2] = center1.nZ;
                
                frcnam = "log.frc";
                setfocus(param);
                
                runDelphi(param);
                
                updateGridPotential(param);
                
                if ( runFocus(param, center1) ) {
                    
                    frcnam = "log.frc";
                    setfocus(param);
                    
                    runDelphi(param);
                    
                    updateGridPotential(param);
                    
                    ergg_prot = param->ergg;
                    
                }
            }
            
#ifdef MPI_PARALLEL
            genEnergyMap(k, ergpolar_crg, bNeutral, MPI_START, globalEnergyVec_partial);
#endif
            
#ifndef MPI_PARALLEL
            genEnergyMap(k, ergpolar_crg, bNeutral);
#endif
            
            resSolv(param);
            runDelphi(param);
            
            
#ifdef MPI_PARALLEL
            EnergyPolarCrg_partial[k-MPI_START] = ergpolar_crg;
#endif
            
#ifndef MPI_PARALLEL
            EnergyPolarCrg[k] = ergpolar_crg;
#endif
            
            ergg_wat = param->ergg; // cout << ergs1res << endl;
            
            
#ifdef MPI_PARALLEL
            EnergyRxnCrg_partial[k-MPI_START] = ergg_prot - ergg_wat;
#endif
            
#ifndef MPI_PARALLEL
            EnergyRxnCrg[k] = ergg_prot - ergg_wat;
#endif
            
            singleRes_to_DelPhi.clear();
            
            
            
            ///////////////////////////////  Second with uncharged residue    /////////////////////////////
            
            
            
            ResetCrgStat();
            
            
            for(j=0;j<vecIonRes[k].size();j++) {     ///  Charging residue Loop Start here
                
                vecCrgStat[vecIonRes[k][j]] = true;
                
                currentResNum = newPDB[vecIonRes[k][j]].res_num;
                currentResNam = newPDB[vecIonRes[k][j]].res_name;
                
                if (newPDB[vecIonRes[k][j]].atom_name == "CG" || newPDB[vecIonRes[k][j]].atom_name == "C1'") {
                    center1.nX = newPDB[vecIonRes[k][j]].coord.X;
                    center1.nY = newPDB[vecIonRes[k][j]].coord.Y;
                    center1.nZ = newPDB[vecIonRes[k][j]].coord.Z;
                }
            }
            
            bNeutral = true;
            
            charge(bNeutral);  //  This function is to charge the corresponding residue (all atoms belong to that residue)
            
            
            /////  Create a singleRes_to_DelPhi intermediate array to DelPhi for solvation energy in water //////
            
            for(i=0;i<newPDB.size();i++) {
                if (   newPDB[i].res_name == currentResNam
                    && newPDB[i].res_num  == currentResNum
                    && newPDB[i].conf != "BK") {
                    singleRes_to_DelPhi.push_back(PQR_to_DelPhi[i]);
                }
            }
            
            
            delphi_progress << "    Now Starting the DelPhi focusing with neutral residue ";
            delphi_progress << currentResNam  << setw(4) << currentResNum << endl;
            
            setbase(param, center1);
            
            runDelphi(param);
            
            vecGridPotential = param->vecGridPot;
            
            if ( runFocus(param, center1) )
            {
                
                param->bAcent = true;
                param->center[0] = center1.nX; param->center[1] = center1.nY; param->center[2] = center1.nZ;
                
                frcnam = "log.frc";
                setfocus(param);
                runDelphi(param);
                updateGridPotential(param);
                
                if ( runFocus(param, center1) ) {
                    
                    frcnam = "log.frc";
                    setfocus(param);
                    runDelphi(param);
                    updateGridPotential(param);
                    
                    ergg_prot = param->ergg;
                }
            }
            
#ifdef MPI_PARALLEL
            genEnergyMap(k, ergpolar_neu, bNeutral, MPI_START, globalEnergyVec_partial);
#endif
            
#ifndef MPI_PARALLEL
            genEnergyMap(k, ergpolar_neu, bNeutral);
#endif
            
            
            resSolv(param);
            runDelphi(param);
            
#ifdef MPI_PARALLEL
            EnergyPolarNeu_partial[k-MPI_START] = ergpolar_neu;
#endif
            
#ifndef MPI_PARALLEL
            EnergyPolarNeu[k] = ergpolar_neu;
#endif
            
            
            ergg_wat = param->ergg; // cout << ergs0res << endl;
            
            
            
#ifdef MPI_PARALLEL
            EnergyRxnNeu_partial[k-MPI_START] = ergg_prot - ergg_wat;
#endif
            
#ifndef MPI_PARALLEL
            EnergyRxnNeu[k] = ergg_prot - ergg_wat;
#endif
            singleRes_to_DelPhi.clear();
            
            
        }
        
        param.reset();
        
#ifndef MPI_PARALLEL
        clock_t end_time = clock();
        float time_elapsed = float(end_time - begin_time) / CLOCKS_PER_SEC ;
        
        cout << "Done." << endl;
        
        cout << endl;
        cout << " Delphi energy calculation finishes in  " << time_elapsed << " sec" << endl;
        cout << endl;
#endif
        
#ifdef MPI_PARALLEL
        
        delphi_progress.close();
        
        MPI::COMM_WORLD.Barrier();
        
        if(id==0)	cout << endl;
        
        int count = MPI_END - MPI_START;
        int rcount[numprocs];
        int displs[numprocs];
        
        MPI::COMM_WORLD.Gather ( &count, 1, MPI::INT, &rcount[0], 1, MPI::INT, 0 );
        
        if(id==0)
        {
            displs[0] = 0;
            for(i=1;i<numprocs;i++)
            {
                displs[i] = displs[i-1] + rcount[i-1];
            }
        }
        
        MPI::COMM_WORLD.Gatherv( &EnergyRxnCrg_partial[0], count, MPI::FLOAT, &EnergyRxnCrg[0], rcount, displs, MPI::FLOAT, 0 );
        MPI::COMM_WORLD.Gatherv( &EnergyRxnNeu_partial[0], count, MPI::FLOAT, &EnergyRxnNeu[0], rcount, displs, MPI::FLOAT, 0 );
        MPI::COMM_WORLD.Gatherv( &EnergyPolarCrg_partial[0], count, MPI::FLOAT, &EnergyPolarCrg[0], rcount, displs, MPI::FLOAT, 0 );
        MPI::COMM_WORLD.Gatherv( &EnergyPolarNeu_partial[0], count, MPI::FLOAT, &EnergyPolarNeu[0], rcount, displs, MPI::FLOAT, 0 );
        
        MPI::COMM_WORLD.Bcast( &EnergyPolarCrg[0], N, MPI::FLOAT, 0 );
        MPI::COMM_WORLD.Bcast( &EnergyPolarNeu[0], N, MPI::FLOAT, 0 );
        MPI::COMM_WORLD.Bcast( &EnergyRxnCrg[0], N, MPI::FLOAT, 0 );
        MPI::COMM_WORLD.Bcast( &EnergyRxnNeu[0], N, MPI::FLOAT, 0 );
                              
                              
        int count1 = ( MPI_END - MPI_START ) * N * 4;
        int rcount1[numprocs];
        int displs1[numprocs];
        
        MPI::COMM_WORLD.Gather ( &count1, 1, MPI::INT, &rcount1[0], 1, MPI::INT, 0 );
        
        if(id==0)
        {
            displs1[0] = 0;
            for(i=1; i<numprocs; i++)
            {
                displs1[i] = displs1[i-1] + rcount1[i-1];
            }
        }
        
        MPI::COMM_WORLD.Gatherv( &globalEnergyVec_partial[0], count1, MPI::FLOAT, &globalEnergyVec[0], rcount1, displs1, MPI::FLOAT, 0 );
        
        MPI::COMM_WORLD.Bcast( &globalEnergyVec[0], globalSIZE, MPI::FLOAT, 0 );
        
        

        for(i=0; i<vecIonRes.size(); i++) {
            for(j=0; j<vecIonRes.size(); j++) {
                EnergyPair[2*i  ][2*j  ].fEnergy = globalEnergyVec[ 4*N*i + j*4 + 0 ];
                EnergyPair[2*i][2*j + 1].fEnergy = globalEnergyVec[ 4*N*i + j*4 + 1 ];
                EnergyPair[2*i + 1][2*j].fEnergy = globalEnergyVec[ 4*N*i + j*4 + 2 ];
                EnergyPair[2*i + 1][2*j + 1].fEnergy = globalEnergyVec[ 4*N*i + j*4 + 3 ];
                
            }
        }

        
        
        
        if(id==0){
            
            cout << "Done. " << endl;
            
            cout << endl;
            cout << " Delphi energy calculation with MPI finishes in  " << (MPI_Wtime() - MPI_Timer) << " sec" << endl;
            cout << endl;
            
        }
        
#endif
        

            /// Average the pairwise energies ///
            
            for (i=0;i<EnergyPair.size();i++) {
                for (j=0;j<EnergyPair[i].size();j++) {
                    EnergyPair[i][j].fEnergy = ( EnergyPair[i][j].fEnergy + EnergyPair[j][i].fEnergy ) / 2.0;
                    EnergyPair[j][i].fEnergy = EnergyPair[i][j].fEnergy;
                }
            }
            ////////////////////////////////////
            
            
#ifdef DEBUG_ENERGY_EXPORT
            
            ofstream fullenergy;
            fullenergy.open("energyfull.txt");
            
            
            for (auto itr : EnergyPolarCrg)
                fullenergy << "EnergyPolarCrg :" << itr << endl;
            
            for (auto itr : EnergyPolarNeu)
                fullenergy << "EnergyPolarNeu :" << itr << endl;
            
            for (auto itr : EnergyRxnCrg)
                fullenergy << "EnergyRxnCrg   :" << itr << endl;
            
            for (auto itr : EnergyRxnNeu)
                fullenergy << "EnergyRxnNeu   :" << itr << endl;
            
            
            for (auto itr : EnergyPair) {
                fullenergy << "EnergyPair     :";
                for (auto itr2 : itr) {
                    fullenergy << itr2.resnum1 << " " << itr2.resnum2 << " " << itr2.fEnergy << " ";
                }
                fullenergy << endl;
            }
            
            fullenergy.close();
            
#endif

        
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
            
            ///  Creat the table contains E_pol and E_desol eneries. ///
            
            ofstream ergsOutput;
            ergsOutput.open("energies.txt");
            
            ergsOutput <<  " RESNAME             epol        desol            ** Energies are shown in kcal/mol " << endl;
            for (i=0;i<vecIonRes.size();i++) {
                
                string sign;
                string resnam = newPDB[vecIonRes[i][0]].res_name;
                
                if (resnam == "GLU" || resnam == "ASP") sign = "-";
                if (resnam == "ARG" || resnam == "HIS" || resnam == "LYS" || resnam == "A" || resnam == "C")  sign = "+";
                
                ergsOutput << newPDB[vecIonRes[i][0]].res_name + "0" + newPDB[vecIonRes[i][0]].chain_id;
                ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[i][0]].res_num; ergsOutput << setfill(' '); ergsOutput << "         ";
                
                ergsOutput << fixed << setw(8) << right << setprecision(4) << EnergyPolarNeu[i] * 0.59 << "    ";
                ergsOutput << fixed << setw(8) << right << setprecision(4) << EnergyRxnNeu[i]   * 0.59   << endl;
                // multiply 1.36/2.3, unit convert(kt->k/mol)
                
                ergsOutput << newPDB[vecIonRes[i][0]].res_name + sign + newPDB[vecIonRes[i][0]].chain_id;
                ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[i][0]].res_num; ergsOutput << setfill(' '); ergsOutput << "         ";
                
                ergsOutput << fixed << setw(8) << right << setprecision(4) << EnergyPolarCrg[i] * 0.59 << "    ";// multiply 1.36/2.3, unit convert(kt->k/mol)
                ergsOutput << fixed << setw(8) << right << setprecision(4) << EnergyRxnCrg[i]   * 0.59 << endl;
                
            }
            
            ergsOutput.close();
            
            
            //// Create the table contains Pairwise energy matrix. ////
            
            ergsOutput.open("pairwise.txt");
            
            ergsOutput << " RESNAME          RESNAME         E_PAIR(in kt) " << endl;
            
            for (i=0; i<vecIonRes.size(); i++) {
                for (j=0; j<vecIonRes.size(); j++) {
                    
                    string sign;
                    string resnam = newPDB[vecIonRes[i][0]].res_name;
                    
                    if (resnam == "GLU" || resnam == "ASP") sign = "-";
                    if (resnam == "ARG" || resnam == "HIS" || resnam == "LYS" || resnam == "A" || resnam == "C")  sign = "+";
                    
                    ergsOutput << newPDB[vecIonRes[i][0]].res_name + sign + newPDB[vecIonRes[i][0]].chain_id;
                    ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[i][0]].res_num; ergsOutput << setfill(' ');ergsOutput << "        ";
                    ergsOutput << newPDB[vecIonRes[j][0]].res_name + sign + newPDB[vecIonRes[i][0]].chain_id;
                    ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[j][0]].res_num; ergsOutput << setfill(' ');ergsOutput << "       ";
                    ergsOutput << fixed << setw(10) << right << setprecision(6) << EnergyPair[i*2][j*2].fEnergy << endl;
                    
                    
                    ergsOutput << newPDB[vecIonRes[i][0]].res_name + sign + newPDB[vecIonRes[i][0]].chain_id;
                    ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[i][0]].res_num; ergsOutput << setfill(' ');ergsOutput << "        ";
                    ergsOutput << newPDB[vecIonRes[j][0]].res_name + "0" + newPDB[vecIonRes[i][0]].chain_id;;
                    ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[j][0]].res_num; ergsOutput << setfill(' ');ergsOutput << "       ";
                    ergsOutput << fixed << setw(10) << right << setprecision(6) << EnergyPair[i*2][j*2+1].fEnergy << endl;
                    
                    ergsOutput << newPDB[vecIonRes[i][0]].res_name + "0" + newPDB[vecIonRes[i][0]].chain_id;
                    ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[i][0]].res_num; ergsOutput << setfill(' ');ergsOutput << "        ";
                    ergsOutput << newPDB[vecIonRes[j][0]].res_name + sign + newPDB[vecIonRes[i][0]].chain_id;
                    ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[j][0]].res_num; ergsOutput << setfill(' ');ergsOutput << "       ";
                    ergsOutput << fixed << setw(10) << right << setprecision(6) << EnergyPair[i*2+1][j*2].fEnergy << endl;
                    
                    ergsOutput << newPDB[vecIonRes[i][0]].res_name + "0" + newPDB[vecIonRes[i][0]].chain_id;
                    ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[i][0]].res_num; ergsOutput << setfill(' ');ergsOutput << "        ";
                    ergsOutput << newPDB[vecIonRes[j][0]].res_name + "0" + newPDB[vecIonRes[i][0]].chain_id;
                    ergsOutput << setfill('0') << setw(4) << newPDB[vecIonRes[j][0]].res_num; ergsOutput << setfill(' ');ergsOutput << "       ";
                    ergsOutput << fixed << setw(10) << right << setprecision(6) << EnergyPair[i*2+1][j*2+1].fEnergy << endl;
                    
                }
            }
            
            ergsOutput.close();
            
            remove("delphi_progress.txt");
            
#ifdef MPI_PARALLEL
        }
#endif
        //////////////////////////////////////////////////////
    }
