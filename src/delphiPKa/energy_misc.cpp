//
//  energy_misc.cpp
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
#include "energy.h"


using namespace std;


/*
 *  generate ionizable residue array, the array is used to further distribute with MPI
 */

void CEnergy:: ionizableArr()
{
    int i=0,j=0,k;
    int tmpResnum = 0;
    vector<int> vectmp;

    for(k=0;k<newPDB.size();++k)
    {
        if (newPDB[k].bIonizable && newPDB[k].conf != "BK" ) {
            if (tmpResnum != newPDB[k].res_num) {
                vecIonRes.push_back(vectmp);
                vectmp.clear();
            }
            tmpResnum = newPDB[k].res_num;
            vectmp.push_back(k);
        }
        else {
            continue;
        }
    }
    vecIonRes.push_back(vectmp);

}




/*
 *  use to force each residue to be uncharged
 */

void CEnergy::ResetCrgStat()
{
    int i;
    
    for(i=0;i<newPDB.size();++i) {
        vecCrgStat[i] = false;
    }
}




/*
 *  make the ionizable residues in their charged or neutral states depends on the passing variable
 *  Different naming rules are applied here. For charged residue naming, "ASP", "GLU", "HIS", "ARG",
 *  "LYS" are used. Corresponding residues in neutral is labled as "AS0", "GL0", "HI0", "AR0", "LY0".
 *
 */

void CEnergy:: charge(bool& bNeutral)
{
    int i;
    float fSiz, fCrg;
    string key;
    stringstream sstmp;
    
#ifdef PRINT_CHARGED_PQR
    
    ofstream output_chrgPQR;
    
    output_chrgPQR.open("new.pqr");
    
#endif
    
    int atom_serial_num = 1;
    for(i=0; i<newPDB.size();++i)
    {
        key = key_to_hashmapQR(newPDB[i]);
        
        if (sizmap.find(key)!=sizmap.end()) {
            fSiz = sizmap.find(key)->second;
        }
        else {
#ifdef MPI_PARALLEL
            if(id==0) {
#endif
                
            cout << "The radii  for atom " << newPDB[i].atom_name << " in residue " << newPDB[i].res_name << " "
            << newPDB[i].res_num << " is not found. The radii  will be set to 0.0 " << endl;
                
#ifdef MPI_PARALLEL
            }
#endif
            fSiz = 0.0000;
        }
        
        if(!bNeutral) {
            
            
            if(vecCrgStat[i]) {
                if (crgmap.find(key)!=crgmap.end()) {
                    fCrg = crgmap.find(key)->second;
                }
                else {
#ifdef MPI_PARALLEL
                    if(id==0) {
#endif
                        
                    cout << "The charge for atom " << newPDB[i].atom_name << " in residue " << newPDB[i].res_name << " "
                    <<  newPDB[i].res_num << " is not found. The charge will be set to 0.0000 "  << endl;
                   
#ifdef MPI_PARALLEL
                    }
#endif
                    fCrg = 0.0000;
                }
            }
            
            else {
                fCrg = 0.0000;
            }
            
        }
        
        if(bNeutral) {
            
            if(vecCrgStat[i]) {
                
                if (newPDB[i].res_name == "ASP") {
                    key = "AS0 " + newPDB[i].atom_name;
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue AS0 " <<  newPDB[i].res_num << " is not found. The charge will be set to 0.0 " << endl;
                        
                        fCrg = 0.0000;
                    }
                }
                
                else if (newPDB[i].res_name == "LYS") {
                    key = "LY0 " + newPDB[i].atom_name;
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue LY0 " <<  newPDB[i].res_num << " is not found. The charge will be set to 0.0 " << endl;
                        
                        fCrg = 0.0000;
                    }
                }
                
                else if (newPDB[i].res_name == "HIS") {
                    key = "HI0 " + newPDB[i].atom_name;
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue HI0 " <<  newPDB[i].res_num << " is not found. The charge will be set to 0.0 " << endl;
                        
                        fCrg = 0.0000;
                    }
                }
                
                else if (newPDB[i].res_name == "ARG") {
                    key = "AR0 " + newPDB[i].atom_name;
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue AR0 " <<  newPDB[i].res_num << " is not found. The charge will be set to 0.0 " << endl;
                        
                        fCrg = 0.0000;
                    }
                }
                
                else if (newPDB[i].res_name == "GLU") {
                    key = "GL0 " + newPDB[i].atom_name;
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue GL0 " <<  newPDB[i].res_num << " is not found. The charge will be set to 0.0 " << endl;
                        
                        fCrg = 0.0000;
                    }
                }
                
                else if (newPDB[i].res_name == "A") {
                    key = "A0 " + newPDB[i].atom_name;
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue A0 is not found. The charge will be set to 0.0 "
                        << endl;
                        fCrg = 0.0000;
                    }
                }
				
                else if (newPDB[i].res_name == "DA") {
                    key = "DA0 " + newPDB[i].atom_name;
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue DA0 is not found. The charge will be set to 0.0 "
                        << endl;
                        fCrg = 0.0000;
                    }
                }
                
                else if (newPDB[i].res_name == "C") {
                    key = "C0 " + newPDB[i].atom_name;
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue C0 is not found. The charge will be set to 0.0 "
                        << endl;
                        fCrg = 0.0000;
                    }
                }
                
                else if (newPDB[i].res_name == "DC") {
                    key = "DC0 " + newPDB[i].atom_name;
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue DC0 is not found. The charge will be set to 0.0 "
                        << endl;
                        fCrg = 0.0000;
                    }
                }
				
                else {
                    key = key_to_hashmapQR(newPDB[i]);
                    if (crgmap.find(key) != crgmap.end()) {
                        fCrg = crgmap.find(key)->second;
                    }
                    else {
                        cout << "The charge for atom " << newPDB[i].atom_name << " in residue " << newPDB[i].res_name << " is not found. The charge will be set to 0.0 "
                        << endl;
                        fCrg = 0.0000;
                    }
                }
                
            }
            
            else
                fCrg = 0.000;
        }
        
        
        sstmp << left << setw(6) << "ATOM"
        << right << setw(5) << atom_serial_num << " "
        << left  << setw(4) << newPDB[i].atom_name << " "
        << left  << setw(3) << newPDB[i].res_name << " "
        << right << setw(1) << newPDB[i].chain_id
        << right << setw(4) << newPDB[i].res_num << "    "
        << right << setw(8) << fixed << setprecision(3) << newPDB[i].coord.X
        << right << setw(8) << fixed << setprecision(3) << newPDB[i].coord.Y
        << right << setw(8) << fixed << setprecision(3) << newPDB[i].coord.Z
        << right << setw(8) << fixed << setprecision(4) << fCrg
        << right << setw(7) << fixed << setprecision(4) << fSiz;
        
        PQR_to_DelPhi[i] = sstmp.str();
        sstmp.str("");
        
#ifdef PRINT_CHARGED_PQR
        
        output_chrgPQR << PQR_to_DelPhi[i] << endl;
        
#endif
        
        atom_serial_num++;
    }
    
    if(bHETATMinPQR) {
        
        int iN = newPDB.size();
        
        for(i=0; i<strHETATM2.size(); i++) {
            PQR_to_DelPhi[iN+i] = strHETATM2[i];
            
#ifdef PRINT_CHARGED_PQR
            
            output_chrgPQR << PQR_to_DelPhi[iN+i] << endl;
            
#endif
            
        }
    }
    
#ifdef PRINT_CHARGED_PQR
    
    output_chrgPQR.close();
    
#endif
    
}





/*
 *  followed with setbase() and setfocus. call delphi class to run the focusing.
 */

bool CEnergy::runFocus(shared_ptr<SPrime> param, SGrid<double>& center1)
{
    bool   bCenterReset = 0;
    int    len1,halflen1, len0, halflen0;
    int    gsize0 = param->igrid1, gsize1 = 65;
    double scale0 = param->scale1, scale1;
    double dist_l, dist_r;
    
    
    if (param->scale == 4)  return false;
    
    param->scale *= 2;
    scale1 = param->scale;
    
    if ((gsize0/scale0) <= (gsize1/scale1)) {
        param->scale *= 2;
        scale1 = param->scale;
        cout << " Scale = 2 The new grid will be out of box. Therefore, the scale will be set to 4.0 for next focussing run. " << endl;
        if ((gsize0/scale0) <= (gsize1/scale1) || scale1 > 4.0) {
            return false;
        }
    }
    
    center0.nX = param->oldmid1.nX; center0.nY = param->oldmid1.nY; center0.nZ = param->oldmid1.nZ;
    
    halflen0 = (gsize0 - 1) / (2 * scale0);
    
    double x_min0 = center0.nX - halflen0;
    double x_max0 = center0.nX + halflen0;
    
    double y_min0 = center0.nY - halflen0;
    double y_max0 = center0.nY + halflen0;
    
    double z_min0 = center0.nZ - halflen0;
    double z_max0 = center0.nZ + halflen0;
    
#ifdef PRIME_DEBUG
    cout << "X0 : " << x_min0 << " " << x_max0 << endl;
    cout << "Y0 : " << y_min0 << " " << y_max0 << endl;
    cout << "Z0 : " << z_min0 << " " << z_max0 << endl;
#endif
    
    len1 = (gsize1 - 1) / scale1; halflen1 = len1 / 2;
    
    double x_min1 = center1.nX - halflen1;
    double x_max1 = center1.nX + halflen1;
    
    double y_min1 = center1.nY - halflen1;
    double y_max1 = center1.nY + halflen1;
    
    double z_min1 = center1.nZ - halflen1;
    double z_max1 = center1.nZ + halflen1;
    
#ifdef PRIME_DEBUG
    cout << "X1 : " << x_min1 << " " << x_max1 << endl;
    cout << "Y1 : " << y_min1 << " " << y_max1 << endl;
    cout << "Z1 : " << z_min1 << " " << z_max1 << endl;
#endif
    
    
    // Check X-boundary
    dist_l = x_min0 - x_min1;   dist_r = x_max0 - x_max1;
    if (dist_l >= 0) {center1.nX += dist_l + 0.01; bCenterReset = true;}
    if (dist_r <= 0) {center1.nX += dist_r - 0.01; bCenterReset = true;}
    
    // Check Y-boundary
    dist_l = y_min0 - y_min1;   dist_r = y_max0 - y_max1;
    if (dist_l >= 0) {center1.nY += dist_l + 0.01; bCenterReset = true;}
    if (dist_r <= 0) {center1.nY += dist_r - 0.01; bCenterReset = true;}
    
    // Check Z-boundary
    dist_l = z_min0 - z_min1;   dist_r = z_max0 - z_max1;
    if (dist_l >= 0) {center1.nZ += dist_l + 0.01; bCenterReset = true;}
    if (dist_r <= 0) {center1.nZ += dist_r - 0.01; bCenterReset = true;}
    
    
    if (bCenterReset)
        cout << "Reset the Acenter. New Acenter : " << center1.nX << " " << center1.nY << " " << center1.nZ << endl;
    
    return true;
}





/*
 *  set up the first step in a three focusing run, with scale = 1
 */

void CEnergy::setbase(shared_ptr<SPrime> param, SGrid<double>& center1)
{
    param->gsize       = 0;
    param->scale       = 1.0;
    param->perfil      = 70.0;
    param->bndcon      = 2;
    
    param->vecStrPDB   = PQR_to_DelPhi;
    param->vecStrFRCIn = PQR_to_DelPhi;
    param->strFRCOut   = "log00.frc";
    
    param->bAcent      = true;
    param->center[0]   = center1.nX;
    param->center[1]   = center1.nY;
    param->center[2]   = center1.nZ;
    
    if(iGaussian) {
        param->iGaussian = 1;
        param->fSigma    = fSigma;
    }
    
    else {
        param->iGaussian = 0;
        param->fSigma    = 0.0;
    }
    
}




/*
 *  generate residue side-chain in free state for desolvation energy calculation
 */
void CEnergy::resSolv(shared_ptr<SPrime> param)
{
    param->bndcon      = 2;
    param->strFRCOut   = "log0.frc";
    param->vecStrPDB   = singleRes_to_DelPhi;
    param->vecStrFRCIn = singleRes_to_DelPhi;
    
}





/*
 *  set up the second and third steps in a three focusing run, with scale = 2 and 4
 */

void CEnergy::setfocus(shared_ptr<SPrime> param)
{
    param->gsize  = 65;
    param->bndcon = 3;
    
    param->strFRCOut = frcnam;
    param->vecStrPDB = PQR_to_DelPhi;
    param->vecStrFRCIn = PQR_to_DelPhi;
    
    if(iGaussian) {
        param->iGaussian = 1;
        param->fSigma    = fSigma;
    }
    else {
        param->iGaussian = 0;
        param->fSigma    = 0.0;
    }
}



// The following HomoDielec function is just a template that shows how to set delphi parameter for homogeneous boundary senario



/*
 *  evoke the homogeneous dielectric model in delphicpp.
 *  since gaussian flag is set to be true by default. this option needs some modification before it works
 *  only for testing purpose, otherwise do not use this function.
 */

void CEnergy :: HomoDielec(shared_ptr<SPrime> param)
{
    param->gsize = 59;
    param->scale = 1.0;
    param->perfil = 70.0;
    param->pdbfile = pdb_input;
    param->strFRCOut = "log.frc";
    param->pdbformat = CommPDB;
    param->vecStrPDB = PQR_to_DelPhi;
    param->vecStrFRCIn = PQR_to_DelPhi;
    
    param->indi = 2;
    param->exdi = 80;
    param->bndcon = 2;
    param->prbrad = 1.4;
    param->ionrad = 2.0;
    param->fMaxc = 0.0001;
    param->bCommFrc = true;
}


// The following GaussDielec function is just a template that shows how to set delphi parameter for gaussian boundary senario


/*
 *  evoke the gaussian dielectric model in delphicpp.
 */

void CEnergy :: GaussDielec(shared_ptr<SPrime> param)
{
    param->iGaussian = 1;
    param->fSigma = 0.95;
    param->fSrfcut = 40.0;
    
    param->scale = 1.0;
    param->perfil = 70.0;
    param->pdbfile = pdb_input;
    param->strFRCOut = "log.frc";
    param->pdbformat = CommPDB;
    param->vecStrPDB = PQR_to_DelPhi;
    param->vecStrFRCIn = PQR_to_DelPhi;
    
    param->indi = 2;
    param->exdi = 80;
    param->bndcon = 2;
    param->prbrad = 1.4;
    param->ionrad = 2.0;
    param->fMaxc = 0.0001;
    param->bCommFrc = true;
    
}





/*
 *  update potential map after each step during a focusing run. only potentials within the newer box will be updated.
 */

void CEnergy :: updateGridPotential(shared_ptr<SPrime> param)

{
    unordered_set<string> setResInBox;
    
    resInBox.resize(newPDB.size());
    
    string key_to_set;
    
    double halflen = (param->gsize - 1) / (2 * param->scale);
    
    double x_min = center1.nX - halflen;
    double x_max = center1.nX + halflen;
    
    double y_min = center1.nY - halflen;
    double y_max = center1.nY + halflen;
    
    double z_min = center1.nZ - halflen;
    double z_max = center1.nZ + halflen;
    
    for (int i=0;i<newPDB.size();i++){
        
        if (newPDB[i].coord.X > x_max || newPDB[i].coord.X < x_min ||
            newPDB[i].coord.Y > y_max || newPDB[i].coord.Y < y_min ||
            newPDB[i].coord.Z > z_max || newPDB[i].coord.Z < z_min)
        {
            string tmp = newPDB[i].res_name + " " + to_string(newPDB[i].res_num);
            
            setResInBox.insert(tmp);
        }
    }
    
    for (int i=0;i<newPDB.size();i++) {
        
        key_to_set = newPDB[i].res_name + " " + to_string(newPDB[i].res_num);
        
        auto itr = setResInBox.find(key_to_set);
        
        if ( itr != setResInBox.end() ) {
            resInBox[i] = false;
        }
        else {
            resInBox[i] = true;
        }
        //        cout << newPDB[i].res_name << " " << newPDB[i].res_num <<" " << newPDB[i].atom_name <<  " " << resInBox[i] << endl;
    }
    
#ifdef PRIME_DEBUG
    
    cout << " center X, Y, Z : " << center1.nX << " " << center1.nY << " " << center1.nZ << endl;
    cout << " Scale = " << param->scale << "   GridSize = " << param->gsize << endl;
    
#endif
    
    for(int itr=0;itr<vecGridPotential.size()-vecCrgHETATM.size();itr++) {
        if (resInBox[itr]) {
            vecGridPotential[itr] = param->vecGridPot[itr];
        }
    }
    
}


