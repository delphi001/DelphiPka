//
//  titration_misc.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "titration2.h"
#include <iomanip>
#include <fstream>


void CTitration2::linearReg() {
    size_t i, j, n;
    float  fPKA;
    float  fVal1, fMin, fMax;
    double c0, c1, cov00, cov01, cov11, sumsq; // c1 is slope, c0 is intercept

    deque<float>   vecX, vecY;
    vector<size_t> pH_index;

    ofstream pkaOuptut;
    pkaOuptut.open("pKa.csv");
    pkaOuptut << " ResName           pKa    polar(charged) polar(neutral) de-solvation(charged) de-solvation(neutral)"
              << endl;

    for (i = 0; i < vec2dProb.size(); i++) {


        fMin = vec2dProb[i][0];
        fMax = vec2dProb[i][0];

        for (j = 0; j < vec2dProb[i].size(); j++) {
            fVal1 = vec2dProb[i][j];
            if (fVal1 > fMax) fMax = fVal1;
            if (fVal1 < fMin) fMin = fVal1;

            if (fVal1 >= 0.1 && fVal1 <= 0.9) {
                pH_index.push_back(j);
                vecX.push_back(pH_initial + j * pH_step);
                vecY.push_back(fVal1);
            }
        }

        if (vecX.empty()) {
            if (fMax < 0.1) fPKA = -1;
            if (fMin > 0.9) fPKA = 15;
        } else {

            if (vecX.front() != pH_initial) {
                vecY.push_front(vec2dProb[i][pH_index.front() - 1]);
                vecX.push_front(vecX.front() - pH_step);
            }

            if (vecX.back() != pH_end) {
                vecY.push_back(vec2dProb[i][pH_index.back() + 1]);
                vecX.push_back(vecX.back() + pH_step);
            }

            n = vecX.size();
            double x[n], y[n];

            for (j = 0; j < n; j++) {
                x[j] = vecX[j];
                y[j] = vecY[j];
            }

            gsl_fit_linear(x, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

            fPKA = (0.50 - c0) / c1;
        }

#ifdef PRIME_DEBUG
        for(int m=0;m<vecX.size(); m++) {
            cout << vecX[m] << ":" << vecY[m] << "   ";
        }
        cout << endl;
#endif

        pkaOuptut << newPDB[vecIonRes[i][0]].res_name;
        pkaOuptut << setfill('0') << setw(4) << right << newPDB[vecIonRes[i][0]].res_num;
        pkaOuptut << setfill(' ');
        pkaOuptut << newPDB[vecIonRes[i][0]].chain_id << "        ";

        string resname = newPDB[vecIonRes[i][0]].res_name;

        if (fPKA > 14.0) {
            if (resname == "ASP" || resname == "GLU" || resname == "TYR" || resname == "THR" || resname == "SER" || resname == "CYS")
                pkaOuptut << "  undetermined";
            else pkaOuptut << "  undetermined";
        } else if (fPKA < 0.0) {
            if (resname == "ARG" || resname == "HIS" || resname == "LYS")
                pkaOuptut << "  undetermined";
            else pkaOuptut << "  undetermined";
        } else pkaOuptut << fixed << setw(7) << setprecision(2) << right << fPKA;

        pkaOuptut << "    ";

        pkaOuptut << fixed << setw(11) << right << setprecision(4) << EnergyPolarCrg[i] * 0.59;
        pkaOuptut << fixed << setw(15) << right << setprecision(4) << EnergyPolarNeu[i] * 0.59;

        pkaOuptut << fixed << setw(20) << right << setprecision(4) << EnergyRxnCrg[i] * 0.59;
        pkaOuptut << fixed << setw(20) << right << setprecision(4) << EnergyRxnNeu[i] * 0.59;

        pkaOuptut << endl;


        vecX.clear();
        vecY.clear();
        pH_index.clear();
    }

}


void CTitration2::genState(int &n) {
    int         tmp;
    vector<int> vectemp;
    int         nState = pow(2, n);
    for (int    k      = 0; k < nState; k++) {

        string statebinary = bitset<32>(k).to_string();
        statebinary = statebinary.substr(32 - n, n);
        for (int j = 0; j < n; j++) {
            string s = statebinary.substr(j, 1);
            tmp      = stoi(s);
            vectemp.push_back(tmp);
        }
        vecState.push_back(vectemp);
        vectemp.clear();
    }
}


void CTitration2::genCluster1() {
    int            i, j, k;
    vector<int>    vecIonRes1d;
    vector<string> vecIonCh1d;
    vector<int>    vectmp;

    for (i = 0; i < vecIonRes.size(); i++) {
        vecIonRes1d.push_back(newPDB[vecIonRes[i][0]].res_num);
        vecIonCh1d.push_back(newPDB[vecIonRes[i][0]].chain_id);
    }


    for (i = 0; i < vecCluster.size(); i++) {
        for (j = 0; j < vecCluster[i].size(); j++) {
            for (k = 0; k < vecIonRes1d.size(); k++) {
                if (vecCluster[i][j] == vecIonRes1d[k] && vecClusterChID[i][j] == vecIonCh1d[k])
                    vectmp.push_back(k);
            }
        }
        vecCluster1.push_back(vectmp);
        vectmp.clear();
    }
}


float CTitration2::genEnergy0(int &n, float &fPH) {
    float  erg0   = 0, pkaref = 0;
    string resnam;
    float  iGamma = 0;

    resnam = newPDB[vecIonRes[n][0]].res_name;

    if ("ARG" == resnam || "HIS" == resnam || "LYS" == resnam || "A" == resnam || "C" == resnam || "DA" == resnam ||
        "DC" == resnam) {
        iGamma = 1.0;
        pkaref = pkamap.find(resnam)->second;
    }

    if ("ASP" == resnam || "GLU" == resnam || "TYR" == resnam || "THR" == resnam || "SER" == resnam || "CYS" == resnam) {
        iGamma = -1.0;
        pkaref = pkamap.find(resnam)->second;
    }

    erg0 = iGamma * 2.3 * (fPH - pkaref);

    return erg0;
}


void CTitration2::output_PQR() {

    ofstream outputPQR;

    string delimiter = ".";
    // string token = pdb_input.substr(0,pdb_input.find(delimiter));
    string token     = pdb_input.substr(0, pdb_input.find(".pdb"));
    string newPQR    = token + "_2.pqr";
    outputPQR.open(newPQR);

    cout << "Writing output PQR file due to pKa result at given pH value " << fGivenPhVal << endl;

    int atom_serial_num     = 1;

    float fSiz, fCrg, fProb = 0.0;

    string key;

    for (auto itr = newPDB.begin(); itr != newPDB.end(); itr++) {


        ///////////  retrieve the charge and size information from hashmap crgmap and sizmap and output pqr format  ///////

        int iCounter;
        if (itr->bIonizable) {
            for (int i = 0; i < vec2dProb.size(); i++) {
                if (newPDB[vecIonRes[i][0]].res_name == itr->res_name &&
                    newPDB[vecIonRes[i][0]].res_num == itr->res_num &&
                    newPDB[vecIonRes[i][0]].chain_id == itr->chain_id) {
                    // fProb = vec2dProb[i][fGivenPhVal];
                    iCounter = (fGivenPhVal - pH_initial ) / pH_step;
                    fProb = vec2dProb[i][iCounter];

                }
            }

            if (itr->res_name == "ASP") {
                if (fProb <= 0.50) key = "AS0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "GLU") {
                if (fProb <= 0.50) key = "GL0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "HIS") {
                if (fProb <= 0.50) key = "HI0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "ARG") {
                if (fProb <= 0.50) key = "AR0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "LYS") {
                if (fProb <= 0.50) key = "LY0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "TYR") {
                if (fProb <= 0.50) key = "TY0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "THR") {
                if (fProb <= 0.50) key = "TH0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "SER") {
                if (fProb <= 0.50) key = "SE0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "CYS") {
                if (fProb <= 0.50) key = "CY0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "DA") {
                if (fProb <= 0.50) key = "DA0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }

            if (itr->res_name == "DC") {
                if (fProb <= 0.50) key = "DC0 " + itr->atom_name;
                if (fProb > 0.50) key  = key_to_hashmapQR(*itr);
            }



        } else {
            key = key_to_hashmapQR(*itr);
        }

        if (crgmap.find(key) != crgmap.end()) {
            fCrg = crgmap.find(key)->second;
        } else fCrg = 0.0000;

        if (sizmap.find(key) != sizmap.end()) {
            fSiz = sizmap.find(key)->second;
        } else fSiz = 0.0000;

        ///////////////////   Crg and Siz retrieve done    /////////////////

        itr->fCharge = fCrg;

        itr->fRidus = fSiz;

        //////////////////    Output to the file in the pqr format  /////////////

        outputPQR << left << setw(6) << "ATOM"
                  << right << setw(5) << atom_serial_num << " "
                  << left << setw(4) << itr->atom_name << " "
                  << left << setw(3) << itr->res_name << " "
                  << right << setw(1) << itr->chain_id
                  << right << setw(4) << itr->res_num << "    "
                  << right << setw(8) << fixed << setprecision(3) << itr->coord.X
                  << right << setw(8) << fixed << setprecision(3) << itr->coord.Y
                  << right << setw(8) << fixed << setprecision(3) << itr->coord.Z
                  << right << setw(8) << fixed << setprecision(4) << fCrg
                  << right << setw(7) << fixed << setprecision(4) << fSiz
                  << endl;


        atom_serial_num++;
    }

    if (!bRemoveHETATM && bHETATMinPQR) {
        if (strHETATM.size() != 0) {

            outputPQR << "TER" << endl;

            for (int i = 0; i < strHETATM.size(); i++) {
                outputPQR << strHETATM[i] << endl;
            }
        }
    }

    outputPQR << "END" << endl;

    cout << "Done." << endl;

    cout << endl;


}
