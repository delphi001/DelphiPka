//
//  network.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "network.h"
#include <cmath>
#include <iomanip>



void CNetwork :: importcenter(){
    
    int i;
    
    
    for (i=0;i<newPDB.size();i++) {
        if (newPDB[i].bIonizable)
        {
            if      (newPDB[i].res_name == "ARG" && newPDB[i].atom_name == "CD") {
                vecPDB.push_back(newPDB[i]);
            }
            else if (newPDB[i].res_name == "LYS" && newPDB[i].atom_name == "NZ") {      // originally CD
                vecPDB.push_back(newPDB[i]);
            }
            else if (newPDB[i].res_name == "ASP" && newPDB[i].atom_name == "OD1") {      // originally CG
                vecPDB.push_back(newPDB[i]);
            }
            else if (newPDB[i].res_name == "GLU" && newPDB[i].atom_name == "OE1") {     // originally CD
                vecPDB.push_back(newPDB[i]);
            }
            else if (newPDB[i].res_name == "HIS" && newPDB[i].atom_name == "NE2") {      // originally CG
                vecPDB.push_back(newPDB[i]);
            }
            else if (newPDB[i].res_name == "A"   && newPDB[i].atom_name == "C5" ) {
                vecPDB.push_back(newPDB[i]);
            }
            else if (newPDB[i].res_name == "C"   && newPDB[i].atom_name == "C4" ) {
                vecPDB.push_back(newPDB[i]);
            }
            else if (newPDB[i].res_name == "DA"   && newPDB[i].atom_name == "C5" ) {
                vecPDB.push_back(newPDB[i]);
            }
            else if (newPDB[i].res_name == "DC"   && newPDB[i].atom_name == "C4" ) {
                vecPDB.push_back(newPDB[i]);
            }
        }
    }
}


float CNetwork :: dist2(PDBFORM a, PDBFORM b)
{
    float X = a.coord.X - b.coord.X;
    float Y = a.coord.Y - b.coord.Y;
    float Z = a.coord.Z - b.coord.Z;
    
    return sqrt( X*X + Y*Y + Z*Z );
    
}


void CNetwork :: grouping()
{
    int i, j, k;
    float dist;
    
    vector<int> vectmp;
    vector<string> vectmp1;
    
    vecCluster.resize(vecPDB.size());
    vecClusterChID.resize(vecPDB.size());
    
    for(i=0; i<vecPDB.size(); i++) {
        
        vectmp.push_back(vecPDB[i].res_num);
        vectmp1.push_back(vecPDB[i].chain_id);
        
        for(j=0; j<vecPDB.size(); j++) {
            if (i != j) {
                dist = dist2(vecPDB[i], vecPDB[j]);
                
                if ( dist <= clusterThreshold ) {
                    
                    vectmp.push_back(vecPDB[j].res_num);
                    vectmp1.push_back(vecPDB[j].chain_id);
                }
                
            }
        }
        
        vecCluster[i] = vectmp;
        vecClusterChID[i] = vectmp1;
        
        vectmp.clear();
        vectmp1.clear();
    }
    
    
}

void CNetwork :: run()
{
    
    importcenter();
    
    grouping();
    
    int i,j;
    
#ifdef PRIME_DEBUG
    
    for (i=0; i<vecCluster.size(); i++) {
        cout << " Cluster " << setw(2) << right << i+1 << " : ";
        for (j=0; j<vecCluster[i].size(); j++) {
            cout << vecClusterChID[i][j];
            cout << setfill('0') << setw(4) << vecCluster[i][j] << "  ";cout << setfill(' ');
        }
        cout << endl;
    }
#endif
 
#ifdef MPI_PARALLEL
    
    int id = MPI::COMM_WORLD.Get_rank();
    
    if (0==id) {
        
#endif
        
        ofstream ClusterOUT;
        ClusterOUT.open("clusters.txt");
        
        for (i=0; i<vecCluster.size(); i++) {
            ClusterOUT << " Cluster " << setw(2) << right << i+1 << " : ";
            for (j=0; j<vecCluster[i].size(); j++) {
                ClusterOUT << vecClusterChID[i][j];
                ClusterOUT << setfill('0') << setw(4) << vecCluster[i][j] << "  "; ClusterOUT << setfill(' ');
            }
            ClusterOUT << endl;
        }
        
        cout << " Clusters information is saved in clusters.txt  file ." << endl;
        cout << endl;

#ifdef MPI_PARALLEL
    }
#endif
}

