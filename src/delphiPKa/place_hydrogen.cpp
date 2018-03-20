//
//  place_hydrogen.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include "place_hydrogen.h"
#include "orbt_type.h"

using namespace std;

void CPlaceHydrogen::run() {

#ifdef MPI_PARALLEL
    id = MPI::COMM_WORLD.Get_rank();
#endif

    newPDB.reserve(9999);

    for (it_pdb = PDB.begin(); it_pdb != PDB.end(); ++it_pdb) {
        newPDB.push_back(*it_pdb);


        if (it_pdb->atom_name == "O") continue;

        if (it_pdb->atom_name == "N" &&
            it_pdb->res_num != 1) {  // when atom is N, it requires to search the connected atom in previous residue;

            if (it_pdb->res_name == "PRO") continue;   // N in Residue PROLINE doesn't have H connected. This skiped.

            Vector0 = it_pdb->coord;

            found_bond_num = 0;

            tpl_key = it_pdb->res_name + " " + it_pdb->atom_name;

            pdb_search_key = it_pdb->chain_id + " " + to_string(it_pdb->res_num - 1) + " " + "C";

            if (map.find(tpl_key) == map.end()) {
                cout << "[WARING] ATOM " << it_pdb->res_name << " " << it_pdb->res_num << " " << it_pdb->atom_name;
                cout << " is not found in topology. It may use non-standard notation and will be skipped. " << endl;
                continue;
            }

            if (mapPDB.find(pdb_search_key) != mapPDB.end()) {

                ++found_bond_num;

                Vector1 = mapPDB.find(pdb_search_key)->second.coord;
            }


            pdb_search_key = it_pdb->chain_id + " " + to_string(it_pdb->res_num) + " " + "CA";

            if (mapPDB.find(pdb_search_key) != mapPDB.end()) {
                ++found_bond_num;
                Vector2 = mapPDB.find(pdb_search_key)->second.coord;
            } else
                cout << "The search key " << pdb_search_key << " is not found in the mapPDB. " << endl;

            if (found_bond_num == 2) {
                BondLength = 1.5; /// ?????

                sp2_paint1(Vector0, Vector1, Vector2, Vector3, BondLength);

                atom_to_insert.res_name  = it_pdb->res_name;
                atom_to_insert.atom_name = ("H");
                atom_to_insert.chain_id  = it_pdb->chain_id;
                atom_to_insert.res_num   = it_pdb->res_num;
                atom_to_insert.coord     = Vector3;

                newPDB.push_back(atom_to_insert);

            }

        } else {  // when the case is not N(nitrogen), the search scope is in the same residue;


            Vector0         = it_pdb->coord;
            found_bond_num  = 0;
            num_of_H_to_add = 0;
            string tpl_keyword = it_pdb->res_name + " " + it_pdb->atom_name;

            if (map.find(tpl_keyword) == map.end()) {
                cout << "[WARING] ATOM " << it_pdb->res_name << " " << it_pdb->res_num << " " << it_pdb->atom_name;
                cout << " is not found in topology. It may use non-standard notation and will be skipped. " << endl;
                continue;
            }

            value_found = map.find(tpl_keyword)->second;

            if (value_found.atom.substr(0, 1) == "C") BondLength = coval_rad_C + coval_rad_H;
            if (value_found.atom.substr(0, 1) == "O") BondLength = coval_rad_O + coval_rad_H;
            if (value_found.atom.substr(0, 1) == "N") BondLength = coval_rad_N + coval_rad_H;


            if (!value_found.bondatom1.empty()) {
                if (value_found.bondatom1.substr(0, 1) == "H")
                    cout << "Fatal : The first connected atom should NOT be H." << endl;
                else {
                    pdb_search_key = it_pdb->chain_id + " " + to_string(it_pdb->res_num) + " " + value_found.bondatom1;
                    Vector1        = mapPDB.find(pdb_search_key)->second.coord;
                    found_bond_num++;
                }
            } else
                cout << "Fatal : The first connected atom should NOT be NULL." << endl;


            if (!value_found.bondatom2.empty()) {
                if (value_found.bondatom2.substr(0, 1) == "H") {
                    num_of_H_to_add++;
                } else {
                    pdb_search_key = it_pdb->chain_id + " " + to_string(it_pdb->res_num) + " " + value_found.bondatom2;
                    Vector2        = mapPDB.find(pdb_search_key)->second.coord;
                    found_bond_num++;
                }
            } else
                continue;


            if (!value_found.bondatom3.empty()) {
                if (value_found.bondatom3.substr(0, 1) == "H")
                    num_of_H_to_add++;
                else {
                    pdb_search_key = it_pdb->chain_id + " " + to_string(it_pdb->res_num) + " " + value_found.bondatom3;
                    Vector3        = mapPDB.find(pdb_search_key)->second.coord;
                    found_bond_num++;
                }

            }


            if (!value_found.bondatom4.empty()) {
                if (value_found.bondatom4.substr(0, 1) == "H")
                    num_of_H_to_add++;
                else {
                    pdb_search_key = it_pdb->chain_id + " " + to_string(it_pdb->res_num) + " " + value_found.bondatom4;
                    Vector4        = mapPDB.find(pdb_search_key)->second.coord;
                    found_bond_num++;
                }
            }


            if (value_found.orbital == "sp3" && found_bond_num == 1) {

                pdb_search_key = it_pdb->chain_id + " " + to_string(it_pdb->res_num) + " " + value_found.bondatom1;
                Vector1        = mapPDB.find(pdb_search_key)->second.coord;

                string next_key;
                next_key = it_pdb->res_name + " " + mapPDB.find(pdb_search_key)->second.atom_name;
                if (map.find(next_key)->second.bondatom3 == it_pdb->atom_name) {
                    pdb_search_key = it_pdb->chain_id + " " + to_string(it_pdb->res_num) + " " +
                                     map.find(next_key)->second.bondatom2;
                    Vector2        = mapPDB.find(pdb_search_key)->second.coord;
                } else if (map.find(next_key)->second.bondatom2 == it_pdb->atom_name) {
                    pdb_search_key = it_pdb->chain_id + " " + to_string(it_pdb->res_num) + " " +
                                     map.find(next_key)->second.bondatom1;
                    Vector2        = mapPDB.find(pdb_search_key)->second.coord;
                }

                sp3_paint3(Vector0, Vector1, Vector2, Vector3, Vector4, Vector5,
                           BondLength, BondAngle_sp3, TorsionAngle_sp3);

                if (num_of_H_to_add == 3) {
                    atom_to_insert.res_name = it_pdb->res_name;
                    atom_to_insert.chain_id = it_pdb->chain_id;
                    atom_to_insert.res_num  = it_pdb->res_num;

                    atom_to_insert.atom_name = value_found.bondatom2;
                    atom_to_insert.coord     = Vector3;
                    newPDB.push_back(atom_to_insert);

                    atom_to_insert.atom_name = value_found.bondatom3;
                    atom_to_insert.coord     = Vector4;
                    newPDB.push_back(atom_to_insert);

                    atom_to_insert.atom_name = value_found.bondatom4;
                    atom_to_insert.coord     = Vector5;
                    newPDB.push_back(atom_to_insert);
                } else if (num_of_H_to_add ==
                           1) {   // It's for the case SER OG/THR OG1/TYR OH which is sp3 orbital and one H to predict
                    atom_to_insert.res_name = it_pdb->res_name;
                    atom_to_insert.chain_id = it_pdb->chain_id;
                    atom_to_insert.res_num  = it_pdb->res_num;

                    atom_to_insert.atom_name = value_found.bondatom2;
                    atom_to_insert.coord     = Vector3;
                    newPDB.push_back(atom_to_insert);
                }
            }

            if (value_found.orbital == "sp3" && found_bond_num == 3 && num_of_H_to_add == 1) {
                sp3_paint1(Vector0, Vector1, Vector2, Vector3, Vector4, BondLength);

                atom_to_insert.res_name  = it_pdb->res_name;
                atom_to_insert.chain_id  = it_pdb->chain_id;
                atom_to_insert.res_num   = it_pdb->res_num;
                atom_to_insert.atom_name = value_found.bondatom4;
                atom_to_insert.coord     = Vector4;
                newPDB.push_back(atom_to_insert);
            }

            if (value_found.orbital == "sp3" && found_bond_num == 2 && num_of_H_to_add == 2) {
                sp3_paint2(Vector0, Vector1, Vector2, Vector3, Vector4, BondLength, BondAngle_sp3);

                atom_to_insert.res_name = it_pdb->res_name;
                atom_to_insert.chain_id = it_pdb->chain_id;
                atom_to_insert.res_num  = it_pdb->res_num;

                atom_to_insert.atom_name = value_found.bondatom3;
                atom_to_insert.coord     = Vector3;
                newPDB.push_back(atom_to_insert);

                atom_to_insert.atom_name = value_found.bondatom4;
                atom_to_insert.coord     = Vector4;
                newPDB.push_back(atom_to_insert);


            }


            if (value_found.orbital == "sp2" && found_bond_num == 1) { //
                pdb_search_key = it_pdb->chain_id + " " +
                                 to_string(it_pdb->res_num) + " " +
                                 value_found.bondatom1;

                Vector1 = mapPDB.find(pdb_search_key)->second.coord;

                string next_key;
                next_key = it_pdb->res_name + " " +
                           mapPDB.find(pdb_search_key)->second.atom_name;

                if (map.find(next_key)->second.bondatom3 == it_pdb->atom_name) {
                    pdb_search_key = it_pdb->chain_id + " " +
                                     to_string(it_pdb->res_num) + " " +
                                     map.find(next_key)->second.bondatom2;

                    Vector2 = mapPDB.find(pdb_search_key)->second.coord;

                } else if (map.find(next_key)->second.bondatom2 == it_pdb->atom_name) {
                    pdb_search_key = it_pdb->chain_id + " " +
                                     to_string(it_pdb->res_num) + " " +
                                     map.find(next_key)->second.bondatom1;

                    Vector2 = mapPDB.find(pdb_search_key)->second.coord;

                }

                sp2_paint2(Vector0, Vector1, Vector2, Vector3, Vector4, BondLength, BondAngle_sp2, TorsionAngle_sp2);


                if (num_of_H_to_add == 2) {


                    atom_to_insert.res_name = it_pdb->res_name;
                    atom_to_insert.chain_id = it_pdb->chain_id;
                    atom_to_insert.res_num  = it_pdb->res_num;

                    atom_to_insert.atom_name = value_found.bondatom2;
                    atom_to_insert.coord     = Vector3;
                    newPDB.push_back(atom_to_insert);

                    atom_to_insert.atom_name = value_found.bondatom3;
                    atom_to_insert.coord     = Vector4;
                    newPDB.push_back(atom_to_insert);
                } else if (num_of_H_to_add == 1) {  // It's for case CYS SG, the orbital is sp2 and one H to predict
                    atom_to_insert.res_name = it_pdb->res_name;
                    atom_to_insert.chain_id = it_pdb->chain_id;
                    atom_to_insert.res_num  = it_pdb->res_num;

                    atom_to_insert.atom_name = value_found.bondatom2;
                    atom_to_insert.coord     = Vector3;
                    newPDB.push_back(atom_to_insert);
                }


            }

            if (value_found.orbital == "sp2" && found_bond_num == 2 && num_of_H_to_add == 1) {
                sp2_paint1(Vector0, Vector1, Vector2, Vector3, BondLength);

                atom_to_insert.res_name  = it_pdb->res_name;
                atom_to_insert.chain_id  = it_pdb->chain_id;
                atom_to_insert.res_num   = it_pdb->res_num;
                atom_to_insert.atom_name = value_found.bondatom3;
                atom_to_insert.coord     = Vector3;
                newPDB.push_back(atom_to_insert);
            }


        }


    }

    vector<PDBFORM>(newPDB).swap(newPDB);  // Trim off the excess capacity of vector newPDB.

}


void CPlaceHydrogen::output_newPDB() {

#ifdef MPI_PARALLEL
    id = MPI::COMM_WORLD.Get_rank();
#endif

    ofstream output_PDB;

#ifdef MPI_PARALLEL
    if (id == 0) {
#endif
        string delimiter = ".";
        // string token = pdb_input.substr(0,pdb_input.find(delimiter));
        string token     = pdb_input.substr(0, pdb_input.find(".pdb"));
        string newPQR    = token + "_1.pqr";
        output_PDB.open(newPQR);

        cout << "Writing output PQR file....." << endl;

#ifdef MPI_PARALLEL
    }
#endif

    int atom_serial_num = 1;


    ///////////  retrieve the charge and size information from hashmap crgmap and sizmap and output pqr format  ///////

    float fSiz, fCrg;
    for (it_pdb = newPDB.begin(); it_pdb != newPDB.end(); ++it_pdb) {

        string key = key_to_hashmapQR(*it_pdb);

        bool bOutput_neutral = false;   //  If this boolean is TRUE, the output PQR file is assigned to be all neutral state.
        if (bOutput_neutral == true) {
            if (it_pdb->res_name == "ASP") key = "AS0 " + it_pdb->atom_name;
            if (it_pdb->res_name == "GLU") key = "GL0 " + it_pdb->atom_name;
            if (it_pdb->res_name == "HIS") key = "HI0 " + it_pdb->atom_name;
            if (it_pdb->res_name == "ARG") key = "AR0 " + it_pdb->atom_name;
            if (it_pdb->res_name == "LYS") key = "LY0 " + it_pdb->atom_name;
            if (bCalMoreRes){
                if (it_pdb->res_name == "TYR") key = "TY0 " + it_pdb->atom_name;
                if (it_pdb->res_name == "THR") key = "TH0 " + it_pdb->atom_name;
                if (it_pdb->res_name == "SER") key = "SE0 " + it_pdb->atom_name;
                if (it_pdb->res_name == "CYS") key = "CY0 " + it_pdb->atom_name;
            }
        }

        if (crgmap.find(key) != crgmap.end()) {
            fCrg = crgmap.find(key)->second;
        } else {
#ifdef MPI_PARALLEL
            if (id == 0) {
#endif
                cout << "The charge for atom " << it_pdb->atom_name << " in residue " << it_pdb->res_name << " "
                     << it_pdb->res_num << " is not found. The charge will be set to 0.0 " << endl;

#ifdef MPI_PARALLEL
            }
#endif
            fCrg = 0.0000;
        }

        if (sizmap.find(key) != sizmap.end()) {
            fSiz = sizmap.find(key)->second;
        } else {
#ifdef MPI_PARALLEL
            if (id == 0) {
#endif
                cout << "The radii  for atom " << it_pdb->atom_name << " in residue " << it_pdb->res_name << " "
                     << it_pdb->res_num << " is not found. The radii  will be set to 0.0 " << endl;
#ifdef MPI_PARALLEL
            }
#endif
            fSiz = 0.0000;
        }

        ///////////////////   Crg and Siz retrieve done    /////////////////

        it_pdb->fCharge = fCrg;

        it_pdb->fRidus = fSiz;


        ///////////////////   Label ionizable residue    /////////////////

        it_pdb->bIonizable =
                it_pdb->res_name == "HIS" ||
                it_pdb->res_name == "LYS" ||
                it_pdb->res_name == "ARG" ||
                it_pdb->res_name == "ASP" ||
                it_pdb->res_name == "GLU" ||
                (bCalMoreRes &&
                 (it_pdb->res_name == "TYR" ||
                  it_pdb->res_name == "THR" ||
                  it_pdb->res_name == "SER" ||
                  it_pdb->res_name == "CYS")
                ) ||
                it_pdb->res_name == "A" ||
                it_pdb->res_name == "C" ||
                it_pdb->res_name == "DA" ||
                it_pdb->res_name == "DC";

        ///////////////////   Label backbone and sidechain   /////////////////

        key = key_to_hashmapQR(*it_pdb);

        it_pdb->conf = map.find(key)->second.conf;


        //////////////////    Output to the file in the pqr format  /////////////

#ifdef MPI_PARALLEL
        if (id == 0) {
#endif

            output_PDB << left << setw(6) << "ATOM"
                       << right << setw(5) << atom_serial_num << " "
                       << left << setw(4) << it_pdb->atom_name << " "
                       << left << setw(3) << it_pdb->res_name << " "
                       << right << setw(1) << it_pdb->chain_id
                       << right << setw(4) << it_pdb->res_num << "    "
                       << right << setw(8) << fixed << setprecision(3) << it_pdb->coord.X
                       << right << setw(8) << fixed << setprecision(3) << it_pdb->coord.Y
                       << right << setw(8) << fixed << setprecision(3) << it_pdb->coord.Z
                       << right << setw(8) << fixed << setprecision(4) << fCrg
                       << right << setw(7) << fixed << setprecision(4) << fSiz
                       << endl;

#ifdef MPI_PARALLEL
        }
#endif
        atom_serial_num++;
    }


#ifdef MPI_PARALLEL
    if (id == 0) {
#endif
        if (!bRemoveHETATM) {
            if (strHETATM.size() != 0) {

                output_PDB << "TER" << endl;

                for (int i = 0; i < strHETATM.size(); i++) {
                    output_PDB << strHETATM[i] << endl;
                }
            }
        }

        output_PDB << "END" << endl;

        cout << "Done. " << endl;

#ifdef MPI_PARALLEL
    }
#endif

}
