//
//  pdb_import.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "pdb_import.h"
#include <string>
#include <iostream>
#include <fstream>

void CPdbImport::run() {

#ifdef MPI_PARALLEL
    id                    = MPI::COMM_WORLD.Get_rank();
#endif

    ifstream pdb_sourcefile(pdb_input);

    if (!pdb_sourcefile) {
#ifdef MPI_PARALLEL
        if (id == 0) {
#endif

            cout << "[ERROR] The PDB file can not be found. " << endl;

#ifdef MPI_PARALLEL
        }
        MPI::Finalize();
#endif

        exit(0);
    }

    PDB.reserve(9999);


#ifdef MPI_PARALLEL
    if (id == 0) {
#endif

        cout << "Starting import PDB file from " << pdb_input << "  ..." << endl;

#ifdef MPI_PARALLEL
    }
#endif

    string secMod1, secMod2;
    bool   bHydrogenExist = false;

    while (getline(pdb_sourcefile, line)) {

        if (line.substr(0, 5) == "MODEL") {

#ifdef MPI_PARALLEL
            if (id == 0) {
#endif
                cout << endl;
                cout << "Warning!! PDB contains multiple MODELS, the first one will be used, the rest will be skipped. "
                     << endl;
                cout << endl;
#ifdef MPI_PARALLEL
            }
#endif
        }

        if (line.substr(0, 6) == "ENDMDL") {
            break;
        }

        if (line.substr(0, 4) == "ATOM") {
            pdb_temp.atom_name = line.substr(12, 4);
            pdb_temp.res_name  = line.substr(17, 3);
            pdb_temp.chain_id  = line.substr(21, 1);
            pdb_temp.res_num   = stoi(line.substr(22, 4)); // string type cast to integer
            pdb_temp.coord.X   = stof(line.substr(30, 8));
            pdb_temp.coord.Y   = stof(line.substr(38, 8));
            pdb_temp.coord.Z   = stof(line.substr(46, 8));

            secMod1 = line.substr(16, 1);
            secMod2 = line.substr(20, 1);

            remove_leadtail_whitespace_string(pdb_temp.atom_name);
            remove_leadtail_whitespace_string(pdb_temp.chain_id);
            remove_leadtail_whitespace_string(pdb_temp.res_name);

            if (pdb_temp.res_name == "ASD") pdb_temp.res_name = "ASP";
            if (pdb_temp.res_name == "ASH") pdb_temp.res_name = "ASP";
            if (pdb_temp.res_name == "AS0") pdb_temp.res_name = "ASP";
            if (pdb_temp.res_name == "GLH") pdb_temp.res_name = "GLU";
            if (pdb_temp.res_name == "GLD") pdb_temp.res_name = "GLU";
            if (pdb_temp.res_name == "GL0") pdb_temp.res_name = "GLU";
            if (pdb_temp.res_name == "HSE") pdb_temp.res_name = "HIS";
            if (pdb_temp.res_name == "HSD") pdb_temp.res_name = "HIS";
            if (pdb_temp.res_name == "HIE") pdb_temp.res_name = "HIS";
            if (pdb_temp.res_name == "HI0") pdb_temp.res_name = "HIS";
            if (pdb_temp.res_name == "LYD") pdb_temp.res_name = "LYS";
            if (pdb_temp.res_name == "LYE") pdb_temp.res_name = "LYS";
            if (pdb_temp.res_name == "LY0") pdb_temp.res_name = "LYS";
            if (pdb_temp.res_name == "ARD") pdb_temp.res_name = "ARG";
            if (pdb_temp.res_name == "AR0") pdb_temp.res_name = "ARG";


            if (secMod1 != "B" && secMod2 != "B") {

                if (pdb_temp.atom_name != "OXT") {

                    if (line.substr(12, 1) != "H" && line.substr(13, 1) != "H") {

                        PDB.push_back(pdb_temp);

                        key = pdb_temp.chain_id + " " + to_string(pdb_temp.res_num) + " " + pdb_temp.atom_name;

                        mapPDB.insert(pair<string, PDBFORM>(key, pdb_temp));
                    } else bHydrogenExist = true;

                }

            } else {

#ifdef MPI_PARALLEL
                if (id == 0) {
#endif
                    cout << "Warning!! PDB contains duplicate ATOM " << pdb_temp.atom_name << " "
                         << pdb_temp.res_name << " " << pdb_temp.res_num;

                    cout << " . Second ATOM will be deleted. " << endl;
#ifdef MPI_PARALLEL
                }
#endif

            }
        }

        if (!bRemoveHETATM) {
            if (line.substr(0, 6) == "HETATM") {
                if (bRemoveWater) {
                    if (line.substr(17, 3) != "HOH") {
                        strHETATM.push_back(line);
                        strCrgHETATM = line.substr(55, 7);
                        vecCrgHETATM.push_back(stof(strCrgHETATM));

                        line.replace(0, 6, "ATOM  ");
                        line.replace(55, 7, " 0.0000");
                        strHETATM2.push_back(line);

                    }
                } else {
                    strHETATM.push_back(line);
                    strCrgHETATM = line.substr(55, 7);
                    vecCrgHETATM.push_back(stof(strCrgHETATM));

                    line.replace(0, 6, "ATOM  ");
                    line.replace(55, 7, " 0.0000");
                    strHETATM2.push_back(line);
                }
            }
        }

    }

    if (bHydrogenExist) {
#ifdef MPI_PARALLEL
        if (id == 0) {
#endif
            cout << endl;
            cout << "Warning!! PDB contains hydrogen atoms, the existing hydrogen will be deleted. ";
            cout << endl;
#ifdef MPI_PARALLEL
        }
#endif
    }

    vector<PDBFORM>(PDB).swap(PDB); // Trim off the excess capacity of vector PDB.

#ifdef MPI_PARALLEL
    if (id == 0) {
#endif
        cout << "Done." << endl;
        cout << endl;

#ifdef MPI_PARALLEL
    }
#endif

}
