//
//  titration.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "titration.h"
#include <iomanip>
#include <fstream>

void CTitration::run() {

    int   i, j, k, n, m1, m2, pos, nPH, iCounter;
    float Temper  = 298.15; // input temperature
    float fFactor = -KCAL2KT / (Temper / ROOMTEMPER);

    bool                  bAcid;
    float                 erg0, ergpolar, ergrxn, ergpair;
    double                ergtotal, expn, fTotalExpN;
//    long double ergtotal, expn, fTotalExpN;
    respair               pairtmp;
    vector<respair>       vecPairtmp;
    vector<double>        vecExpN;
//    vector<long double> vecExpN;
    vector<float>         vecProb;
    vector<vector<bool> > vecBool2d;


    cout << " Start titration prediction ...  " << endl;
    cout << endl;
    clock_t begin_time = clock();

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

        for (int j = 0; j < vecCluster1[i].size(); j++) {
            for (int m = 0; m < vecCluster1.size(); m++) {
                if (i == m) {
                    for (int n = j; n < vecCluster1[m].size(); n++) {
                        if (vecCluster1[i][j] != vecCluster1[m][n]) {
                            pairtmp.res1 = vecCluster1[i][j];
                            pairtmp.res2 = vecCluster1[m][n];
                            vecPairtmp.push_back(pairtmp);
                        }
                    }
                } else {
                    for (int n = 0; n < vecCluster1[m].size(); n++) {
                        if (vecCluster1[i][j] != vecCluster1[m][n]) {
                            pairtmp.res1 = vecCluster1[i][j];
                            pairtmp.res2 = vecCluster1[m][n];
                            vecPairtmp.push_back(pairtmp);
                        }
                    }
                }
            }

        }

        vecPair.push_back(vecPairtmp);
        vecPairtmp.clear();

    }
    // end of pairwise creation.

    float fPH = 0.0;

    for (fPH = pH_initial; fPH <= pH_end; fPH += pH_step) {

        iCounter = (fPH - pH_initial) / pH_step;

        for (i = 0; i < vecCluster1.size(); i++) {
            n = vecCluster1[i].size();
            vecProb.resize(n);
            genState(n);


#ifdef PRIME_DEBUG

            for (j=0; j<vecPair[i].size(); j++) {cout << vecPair[i][j].res1 << " " << vecPair[i][j].res2 << "   ";}
            cout << endl;
#endif

            vecExpN.resize(vecState.size());

            for (j = 0; j < vecState.size(); j++) {

                ergrxn = 0.0;
                ergpolar = 0.0;
                erg0     = 0.0;

                for (k = 0; k < vecState[j].size(); k++) {

#ifdef PRIME_DEBUG
                    cout << vecCluster1[i][k] << ":" << vecState[j][k] << " ";
#endif

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

#ifdef PRIME_DEBUG
                cout << right << ergtotal << "  " << expn;
                cout << endl;
#endif

            }

            fTotalExpN = 0.0;
            for (j = 0; j < vecExpN.size(); j++) { fTotalExpN += vecExpN[j]; }

            for (j = 0; j < vecExpN.size(); j++) { vecExpN[j] /= fTotalExpN; }

#ifdef PRIME_DEBUG
            cout << " fTotalExpN = " << fTotalExpN << endl;
#endif

            fill(vecProb.begin(), vecProb.end(), 0);

            for (j = 0; j < vecState.size(); j++) {
                for (k = 0; k < vecState[j].size(); k++) {
                    if (1 == vecState[j][k]) {
                        vecProb[k] += vecExpN[j];
                    }
                }
            }

#ifdef PRIME_DEBUG
            for(auto itr : vecProb)
                cout << itr << "  ";
#endif


            for (j = 0; j < vecProb.size(); j++) {

                vec2dProb[vecCluster1[i][j]][iCounter] = vecProb[j];
                vecBool2d[vecCluster1[i][j]][iCounter] = true;
            }

            vecState.clear();
        }

    }
    cout << endl;


    clock_t end_time     = clock();
    float   time_elapsed = float(end_time - begin_time) / CLOCKS_PER_SEC;

    cout << endl;
    cout << " Titration finishes in  " << time_elapsed << " sec" << endl;
    cout << endl;


    ofstream probOutput;
    probOutput.open("titra.dat");
    //probOutput << "RES/pH     0.0      1.0      2.0      3.0      4.0      5.0      6.0      7.0      8.0      9.0      10.0     11.0     12.0     13.0     14.0" << endl;

    probOutput << "RES/pH  " << "  ";
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

    cout << endl;
    cout << " Titration table is saved in titra.txt file.   " << endl;
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


    linearReg();

    cout << " pKa predication values are saved in pKa.csv file.  " << endl;
    cout << endl;
}
