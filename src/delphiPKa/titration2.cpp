//
//  titration2.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 12/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "titration2.h"
#include <iomanip>
#include <fstream>
#include <unordered_set>

void CTitration2::run() {

    int   i, j, k, n, m1, m2, pos, nPH, iCounter;
    float Temper  = 298.15; // input temperature
    float fFactor = -KCAL2KT / (Temper / ROOMTEMPER);
    bool  bAcid;
    float fProb;
    float erg0, ergpolar, ergrxn, ergpair;

    // double ergtotal, expn, fTotalExpN;
    long double ergtotal, expn, fTotalExpN;

    respair         pairtmp;
    vector<respair> vecPairtmp;

    unordered_set<int> current_set;

    //  vector<double> vecExpN;
    vector<long double>   vecExpN;
    vector<float>         vecProb;
    vector<vector<bool> > vecBool2d;

#ifdef MPI_PARALLEL
    int id = MPI::COMM_WORLD.Get_rank();  // MPI rank

    if (id == 0) {

        cout << " Start titration prediction ...  " << endl;
        cout << endl;
    }
#endif


#ifndef MPI_PARALLEL

    cout << " Start titration prediction ...  " << endl; cout << endl;

#endif

    // Resive the prob table to x = ionres, y = nph .

    nPH = (pH_end - pH_initial) / pH_step + 1;
    vec2dProb.resize(vecIonRes.size());
    vecBool2d.resize(vecIonRes.size());

    for (i = 0; i < vec2dProb.size(); i++) {

        vec2dProb[i].resize(nPH);
        vecBool2d[i].resize(nPH);
        for (j = 0; j < vecBool2d[i].size(); j++) {
            vecBool2d[i][j] = false;
        }
    }

    //////////

    genCluster1();  // Generate 2d array vecCluster1 with new index for lookup in energy table

    /////////



    // genterate the pairwise interaction (2^n) within each cluster.
    for (int i = 0; i < vecCluster1.size(); i++) {

        vector<int> res_out_cluster;

        for (int j = 0; j < vecCluster1[i].size(); j++) {
            current_set.insert(vecCluster1[i][j]);
        }    // Create a hash_set to store the residue index in current cluster

        for (j = 0; j < vecIonRes.size(); j++) {
            if (current_set.find(j) == current_set.end()) res_out_cluster.push_back(j);
        }    // In the range of entire ionizable set, check each with the previous created hash_set, if not existed, push back to the vector.
        // The new vector that created contained all ionizable residue index not included in current cluster.

        for (int j = 0; j < vecCluster1[i].size(); j++) {

            for (int k = j + 1; k < vecCluster1[i].size(); k++) {
                pairtmp.res1 = vecCluster1[i][j];
                pairtmp.res2 = vecCluster1[i][k];

                vecPairtmp.push_back(pairtmp);
            } // For the residue index in the same cluster

            for (int k = 0; k < res_out_cluster.size(); k++) {
                pairtmp.res1 = vecCluster1[i][j];
                pairtmp.res2 = res_out_cluster[k];
                vecPairtmp.push_back(pairtmp);
            } // For the residue index out of the same cluster, search in vector res_out_cluster.

        }

        vecPair.push_back(vecPairtmp);
        vecPairtmp.clear();
        current_set.clear();
        res_out_cluster.clear();


    }
    // end of pairwise creation.

#ifdef MPI_PARALLEL

    /* MPI configuration */

    int N = vecCluster1.size();

    int numprocs  = MPI::COMM_WORLD.Get_size();
    int MPI_START, MPI_END;
    int rank_size = N / numprocs;
    MPI_START = (N / numprocs) * id;

    MPI::COMM_WORLD.Barrier();


    float vecProbMPI[N];
    bool  vecBoolMPI[N];

    float vecProbMPI_partial[rank_size + 1];
    float vecBoolMPI_partial[rank_size + 1];

    if (N % numprocs > id) {
        MPI_START += id;
        MPI_END = MPI_START + (N / numprocs) + 1;
    } else {
        MPI_START += N % numprocs;
        MPI_END = MPI_START + (N / numprocs);
    }

    double MPI_Timer = MPI_Wtime();

    /* MPI configure finish */
#endif

#ifndef MPI_PARALLEL
    clock_t begin_time = clock(); // Sequential version Timer
#endif


    float fPH = 0.0;
    for (fPH = pH_initial; fPH <= pH_end; fPH += pH_step) {

        iCounter = (fPH - pH_initial) / pH_step;


#ifndef MPI_PARALLEL
        for (i=0; i<vecCluster1.size(); i++) {  // Sequential verison loop
#endif


#ifdef MPI_PARALLEL
        for (i = MPI_START; i < MPI_END; i++) {  // MPI version loop
#endif
            n = vecCluster1[i].size();

            vecProb.resize(n);
            genState(n);


            vecExpN.resize(vecState.size());


            for (j = 0; j < vecState.size(); j++) {

                ergrxn = 0.0;
                ergpolar = 0.0;
                erg0     = 0.0;

                for (k = 0; k < vecState[j].size(); k++) {

                    if (1 == vecState[j][k]) {
                        m1 = vecCluster1[i][k];
                        ergpolar += EnergyPolarCrg[m1];
                        ergrxn += EnergyRxnCrg[m1];
                        erg0 += genEnergy0(vecCluster1[i][k], fPH);
                    }

                    if (0 == vecState[j][k]) {
                        m1 = vecCluster1[i][k];
                        ergpolar += EnergyPolarNeu[m1];
                        ergrxn += EnergyRxnNeu[m1];
                        erg0 += 0;
                    }

                }

                ergpair = 0.0;

                for (k = 0; k < vecPair[i].size(); k++) {

                    pos = find(vecCluster1[i].begin(), vecCluster1[i].end(), vecPair[i][k].res1) -
                          vecCluster1[i].begin();
                    if (pos != vecCluster1[i].size()) { // Found.
                        if (1 == vecState[j][pos]) { m1 = vecPair[i][k].res1 * 2; }
                        if (0 == vecState[j][pos]) { m1 = vecPair[i][k].res1 * 2 + 1; }
                    } else { m1 = vecPair[i][k].res1 * 2; }


                    pos = find(vecCluster1[i].begin(), vecCluster1[i].end(), vecPair[i][k].res2) -
                          vecCluster1[i].begin();

                    if (pos != vecCluster1[i].size()) { // Found.

                        if (1 == vecState[j][pos]) { m2 = vecPair[i][k].res2 * 2; }
                        if (0 == vecState[j][pos]) { m2 = vecPair[i][k].res2 * 2 + 1; }
                    } else {
                        if (vecBool2d[vecPair[i][k].res2][iCounter]
                            && vec2dProb[vecPair[i][k].res2][iCounter] > 0.5) { m2 = vecPair[i][k].res2 * 2; }

                        else {

                            string resnam = newPDB[vecIonRes[vecPair[i][k].res2][0]].res_name;

                            float pkaref = pkamap.find(resnam)->second;

                            if (resnam == "GLU" || resnam == "ASP" || resnam == "TYR" || resnam == "THR" || resnam == "SER" || resnam == "CYS") bAcid = true; else bAcid = false;

                            if (bAcid) {
                                if (fPH < pkaref) { m2 = vecPair[i][k].res2 * 2 + 1; }
                                if (fPH >= pkaref) { m2 = vecPair[i][k].res2 * 2; }
                            }

                            if (!bAcid) {
                                if (fPH < pkaref) { m2 = vecPair[i][k].res2 * 2; }
                                if (fPH >= pkaref) { m2 = vecPair[i][k].res2 * 2 + 1; }
                            }
                        }
                    }

                    //                cout << right << setw(8) << fixed << setprecision(4) << EnergyPair[m1][m2].fEnergy << "   ";

                    ergpair += EnergyPair[m1][m2].fEnergy;

                }

                ergtotal = (erg0 + ergrxn + ergpolar + ergpair);

                expn = exp(-ergtotal * KCAL2KT);

                vecExpN[j] = expn;

            }

            fTotalExpN = 0.0;

            for (j = 0; j < vecExpN.size(); j++) { fTotalExpN += vecExpN[j]; }

            for (j = 0; j < vecExpN.size(); j++) { vecExpN[j] /= fTotalExpN; }

            fProb = 0.0;

            for (j = 0; j < vecState.size(); j++) {
                if (1 == vecState[j][0]) {
                    fProb += vecExpN[j];
                }
            }
#ifdef MPI_PARALLEL
            vecProbMPI_partial[i - MPI_START] = fProb;
            vecBoolMPI_partial[i - MPI_START] = true;
            vecState.clear();
#endif


#ifndef MPI_PARALLEL

            vec2dProb[i][iCounter] = fProb;
            vecBool2d[i][iCounter] = true;
            vecState.clear();
#endif

        }

#ifdef MPI_PARALLEL

        MPI::COMM_WORLD.Barrier();

        int count = MPI_END - MPI_START;
        int rcount[numprocs];
        int displs[numprocs];

        MPI::COMM_WORLD.Gather(&count, 1, MPI::INT, &rcount[0], 1, MPI::INT, 0);

        if (id == 0) {
            displs[0] = 0;
            for (int m = 1; m < numprocs; m++) {
                displs[m] = displs[m - 1] + rcount[m - 1];
            }
        }

        MPI::COMM_WORLD.Gatherv(&vecProbMPI_partial[0], count, MPI::FLOAT, &vecProbMPI[0], rcount, displs, MPI::FLOAT,
                                0);

        MPI::COMM_WORLD.Gatherv(&vecBoolMPI_partial[0], count, MPI::BOOL, &vecBoolMPI[0], rcount, displs, MPI::BOOL, 0);


        if (id == 0) {
            for (i = 0; i < vecCluster1.size(); i++) {
                vec2dProb[i][iCounter] = vecProbMPI[i];
                vecBool2d[i][iCounter] = vecBoolMPI[i];
            }
        }
#endif

    }
    //   Titration from pH=0 to pH=14 is finished   //


    //   Save the titration information in titra.csv and [OPTIONAL] print onscreen if needed  //

#ifndef MPI_PARALLEL  // Sequential version Runtime report

    clock_t end_time = clock();
    float time_elapsed = float(end_time - begin_time) / CLOCKS_PER_SEC ;

    cout << endl;
    cout << " Titration finishes in  " << time_elapsed << " sec" << endl;
    cout << endl;
#endif

#ifdef MPI_PARALLEL   // MPI version Runtime report
    if (0 == id) {
        cout << endl;
        cout << " Titration finishes in  " << (MPI_Wtime() - MPI_Timer) << " sec" << endl;
        cout << endl;
    }
#endif


#ifdef MPI_PARALLEL
    if (id == 0) {
#endif

        ofstream probOutput;
        probOutput.open("titra.csv");
        //probOutput << "RES/pH     0.0      1.0      2.0      3.0      4.0      5.0      6.0      7.0      8.0      9.0      10.0     11.0     12.0     13.0     14.0" << endl;
        probOutput << "RES/pH  " << "  ";;
        for (fPH = pH_initial; fPH <= pH_end; fPH += pH_step) {
            probOutput << fixed << setw(5) << setprecision(4) << fPH << "   ";
        }
        probOutput << endl;

        for (i = 0; i < vec2dProb.size(); i++) {

            probOutput << newPDB[vecIonRes[i][0]].res_name;
            probOutput << setfill('0') << setw(4) << right << newPDB[vecIonRes[i][0]].res_num;
            probOutput << setfill(' ');
            probOutput << newPDB[vecIonRes[i][0]].chain_id << "  ";
            for (j = 0; j < vec2dProb[i].size(); j++) {
                probOutput << fixed << setw(5) << setprecision(4) << left << vec2dProb[i][j] << "   ";
            }
            probOutput << endl;
        }

        // Calculate the net charge as function of pH; added *2016-03-26* //
        vector<double> vecNetCrg;
        for (i = 0; i < vec2dProb[0].size(); i++) {
            double total_netCrg = 0.0;
            for (j = 0; j < vec2dProb.size(); j++) {
                if (newPDB[vecIonRes[j][0]].res_name == "A" ||
                    newPDB[vecIonRes[j][0]].res_name == "C" ||
                    newPDB[vecIonRes[j][0]].res_name == "DA" ||
                    newPDB[vecIonRes[j][0]].res_name == "DC") { total_netCrg += vec2dProb[j][i] * 1.0 + (-1.0); }

                else {
                    int sign = 1;
                    if (newPDB[vecIonRes[j][0]].res_name == "LYS" ||
                        newPDB[vecIonRes[j][0]].res_name == "ARG" ||
                        newPDB[vecIonRes[j][0]].res_name == "HIS")
                        sign = 1;
                    else if (newPDB[vecIonRes[j][0]].res_name == "ASP" ||
                             newPDB[vecIonRes[j][0]].res_name == "GLU" ||
                             newPDB[vecIonRes[j][0]].res_name == "TYR" ||
                             newPDB[vecIonRes[j][0]].res_name == "THR" ||
                             newPDB[vecIonRes[j][0]].res_name == "SER" ||
                             newPDB[vecIonRes[j][0]].res_name == "CYS")
                        sign = -1;
                    total_netCrg += vec2dProb[j][i] * sign;
                }
            }
            vecNetCrg.push_back(total_netCrg);
        }

        probOutput << "----------" << endl;
        probOutput << "NetCharge ";
        for (i = 0; i < vecNetCrg.size(); i++) {
            probOutput << fixed << setw(9) << setprecision(3) << left << vecNetCrg[i];
        }
        probOutput << endl;
        cout << " Titration table is saved in titra.csv file.   " << endl;
        cout << endl;

#ifdef PRIME_DEBUG

        cout << "RES/pH     0.0      1.0      2.0      3.0      4.0      5.0      6.0      7.0      8.0      9.0      10.0     11.0     12.0     13.0     14.0" << endl;
        for(i=0;i<vec2dProb.size();i++) {
            
            cout << newPDB[vecIonRes[i][0]].res_name;
            cout << setfill('0') << setw(4) << right << newPDB[vecIonRes[i][0]].res_num; cout << setfill(' ');
            cout << newPDB[vecIonRes[i][0]].chain_id << "  ";
            for (j=0;j<vec2dProb[i].size();j++) {
                cout << fixed << setw(5) << setprecision(4) << left << vec2dProb[i][j] << "   ";
            } cout <<endl;
        }
#endif

        if (bOutPQRpka) {

            output_PQR();

        }


        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        linearReg();      // Linear regression to fit the titration curve and calculate pKa value for each residue

        ////////////////////////////////////////////////////////////////////////////////////////////////////////


        cout << " pKa predication values are saved in pKa.csv file.  " << endl;
        cout << endl;

#ifdef MPI_PARALLEL
    }
#endif
}
