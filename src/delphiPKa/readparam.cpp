//
//  readparam.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "readparam.h"


void CReadParam::run() {

#ifdef MPI_PARALLEL
    int id = MPI::COMM_WORLD.Get_rank();
#endif


    size_t pos;
    string strtmp, strtmp1, strtmp2;

    ifstream param(paramfile);

    if (!param) {
#ifdef MPI_PARALLEL
        if (id == 0) {
#endif

            cout << "[ERROR] The parameter file can not be found. " << endl;

#ifdef MPI_PARALLEL
        }
        MPI::Finalize();
#endif

        exit(0);
    }

#ifdef MPI_PARALLEL
    if (id == 0) {
#endif

        cout << " Reading paramter file from " << paramfile << " . " << endl;

#ifdef MPI_PARALLEL
    }
#endif

    while (getline(param, line)) {
        if (line.substr(0, 1) == "#")
            continue;

        remove_leadtail_whitespace_string(line);
        if (!line.empty()) {
            pos = line.find(":");
            if (pos == string::npos) continue;
            if (pos != string::npos) {

                strtmp1 = line.substr(0, pos - 1);
                strtmp2 = line.substr(pos + 1, line.length() - pos - 1);
                remove_leadtail_whitespace_string(strtmp1);
                remove_leadtail_whitespace_string(strtmp2);

                if (strtmp1 == "PDB file name") pdbname              = strtmp2;
                if (strtmp1 == "Charge parameter") crg_file          = strtmp2;
                if (strtmp1 == "Radius parameter") siz_file          = strtmp2;
                if (strtmp1 == "Topology parameter") parameter_input = strtmp2;
                if (strtmp1 == "Gaussian surface") iGaussian         = stoi(strtmp2);

                if (strtmp1 == "Kmean++ cluster number") {
                    if (strtmp2 == "AUTO") {
                        bClusterAuto = true;
                        n_cluster    = 0;
                    }
                    else {
                        bClusterAuto = false;
                        n_cluster = stoi(strtmp2);
                    }
                }

                if (strtmp1 == "HETATM in PQR format") {
                    strtmp                                       = strtmp2;
                    if (toupper(strtmp2[0]) == 'T') bHETATMinPQR = true;
                    else if (toupper(strtmp2[0]) == 'F') bHETATMinPQR = false;
                }

                if (strtmp1 == "Remove HETATM") {
                    strtmp                                        = strtmp2;
                    if (toupper(strtmp2[0]) == 'T') bRemoveHETATM = true;
                    else if (toupper(strtmp2[0]) == 'F') bRemoveHETATM = false;
                }

                if (strtmp1 == "Remove water molecule") {
                    strtmp                                       = strtmp2;
                    if (toupper(strtmp2[0]) == 'T') bRemoveWater = true;
                    else if (toupper(strtmp2[0]) == 'F') bRemoveWater = false;
                }

                if (strtmp1 == "Do Protonation") {
                    strtmp                                    = strtmp2;
                    if (toupper(strtmp2[0]) == 'T') bDoProton = true;
                    else if (toupper(strtmp2[0]) == 'F') bDoProton = false;
                }

                if (strtmp1 == "Do Energy Calculation") {
                    strtmp                                    = strtmp2;
                    if (toupper(strtmp2[0]) == 'T') bDoEnergy = true;
                    else if (toupper(strtmp2[0]) == 'F') bDoEnergy = false;
                }

                if (strtmp1 == "Do pKa's  Calculation") {
                    strtmp                                 = strtmp2;
                    if (toupper(strtmp2[0]) == 'T') bDoPka = true;
                    else if (toupper(strtmp2[0]) == 'F') bDoPka = false;
                }

                if (strtmp1 == "Output PQR file (With Topology)") {
                    strtmp                                      = strtmp2;
                    if (toupper(strtmp2[0]) == 'T') bOutPQRtopo = true;
                    else if (toupper(strtmp2[0]) == 'F') bOutPQRtopo = false;
                }

                if (strtmp1 == "Output PQR file (With pKa result)") {
                    if (toupper(strtmp2[0]) == 'T') bOutPQRpka = true;
                    else if (toupper(strtmp2[0]) == 'F') bOutPQRpka = false;
                }

                if (strtmp1 == "At given pH Value") fGivenPhVal = stof(strtmp2);

                if (strtmp1 == "Hydrogen of GLU Attached to Atom") HofGLU = strtmp2;
                if (strtmp1 == "Hydrogen of ASP Attached to Atom") HofASP = strtmp2;

                if (strtmp1 == "Internal Dielectric") indi = stof(strtmp2);
                if (strtmp1 == "External Dielectric") exdi = stof(strtmp2);
                if (strtmp1 == "Ion radius") ionrad        = stof(strtmp2);
                if (strtmp1 == "Probe radius") prbrad      = stof(strtmp2);
                if (strtmp1 ==
                    "Convergence threshold (MAXC)")
                    maxc             = stof(strtmp2);
                if (strtmp1 ==
                    "Cluster Delimitation Threshold (A)")
                    clusterThreshold = stof(strtmp2);
                if (strtmp1 ==
                    "Variance of Gaussian distribution")
                    fSigma           = stof(strtmp2);
                if (strtmp1 == "Surface cutoff") fSrfcut      = stof(strtmp2);
                if (strtmp1 == "pH End Value") pH_end         = stof(strtmp2);
                if (strtmp1 == "pH Initial Value") pH_initial = stof(strtmp2);
                if (strtmp1 == "pH Interval") pH_step         = stof(strtmp2);

                if (strtmp1 == "Salt Concentration") salt = stof(strtmp2);
            }

        }

    }

#ifdef MPI_PARALLEL
    if (id == 0) {
#endif

        cout << endl;

#ifdef MPI_PARALLEL
    }
#endif

}
