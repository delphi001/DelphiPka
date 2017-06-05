//
//  clustering.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "clustering.h"
#include <cmath>
#include <bitset>
#include <iomanip>

void CClustering :: importcenter(){
    
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
        }
    }
}


float CClustering :: dist2(PDBFORM a, PDBFORM b)
{
    float X = a.coord.X - b.coord.X;
    float Y = a.coord.Y - b.coord.Y;
    float Z = a.coord.Z - b.coord.Z;
    
    return X*X + Y*Y + Z*Z;
    
}

float CClustering :: randf(float m)
{
    return m * rand() / (RAND_MAX - 1.);
}



int CClustering :: nearest(PDBFORM pt, vector<PDBFORM> cent, int n_cluster, float *d2)
{
    int i, min_index;
    float d, min_dist;
    
    min_index = pt.iGroup;
    min_dist = HUGE_VAL;
    
    for (i=0; i<n_cluster; i++) {
        d = dist2(cent[i], pt);
        
        if (min_dist > d) {
            
            min_dist  = d;
            min_index = i;
        }
    }
    
    if(d2) *d2 = min_dist;
    
    
    return min_index;
}



void CClustering :: kmeanpp(vector<PDBFORM>& pts, int len, vector<PDBFORM> cent, int n_cent)
{
    int i, j, n_cluster;
    float sum;
    vector<float> d(len);
    
    i = rand() % len;
    cent[0] = pts[ len - i - 1 ];
    
    for (n_cluster=1; n_cluster < n_cent; n_cluster++) {
        sum = 0.0;
        for (j=0; j<len; j++) {
            nearest(pts[j], cent, n_cluster, &d[j]);
            sum += d[j];
            
        }
        sum = randf(sum);
        
        for (j=len-1; j>=0; j--) {
            sum -= d[j];
            if (sum > 1.0e-13)  continue;
            cent[n_cluster] = pts[j];
            break;
        }
    }
    
    for (j=0; j<len; j++) {
        pts[j].iGroup = nearest(pts[j], cent, n_cluster, 0);
    }
    
    
}


vector<PDBFORM> CClustering :: centdetect(vector<PDBFORM>& pts, int len, int n_cluster)
{
    int i, j, min_i;
    int changed;
    vector<PDBFORM> cent(n_cluster);
    
    kmeanpp(pts, len, cent, n_cluster);
    
    while (changed > (len >> 10))
    {
        for (i=0; i<n_cluster; i++) {
            cent[i].iGroup = 0;
            cent[i].coord.X = 0.0;
            cent[i].coord.Y = 0.0;
            cent[i].coord.Z = 0.0;
        }
        
        for (j=0; j<len; j++) {
            
            cent[pts[j].iGroup].iGroup++;
            cent[pts[j].iGroup].coord.X += pts[j].coord.X;
            cent[pts[j].iGroup].coord.Y += pts[j].coord.Y;
            cent[pts[j].iGroup].coord.Z += pts[j].coord.Z;
            
        }
        
        for (i=0; i<n_cluster; i++) {
            
            cent[i].coord.X /= cent[i].iGroup;
            cent[i].coord.Y /= cent[i].iGroup;
            cent[i].coord.Z /= cent[i].iGroup;
            
        }
        
        changed = 0;
        
        for (j=0; j<len; j++) {
            
            min_i = nearest(pts[j], cent, n_cluster, 0);
            
            if (min_i != pts[j].iGroup) {
                changed++;
                pts[j].iGroup = min_i;
            }
            
        }
    }
    
    for (i=0; i<n_cluster; i++) {
        cent[i].iGroup = i;
    }
    
    
    
    return cent;
}



void CClustering :: checkcluster(vector<PDBFORM>& pts, int len, int n_cluster)
{
    int i, j, k, group = 0;
    float d0, min_d0;
    float d1, min_d1;
    
    for (j=0; j<len; j++) {

        min_d0 = HUGE_VAL;
        min_d1 = HUGE_VAL;
        
        for (k=0; k<len; k++) {
            if(j!=k) {
                if(pts[j].iGroup == pts[k].iGroup) {
                    d0 = dist2(pts[j], pts[k]);
                    if(d0 < min_d0) min_d0 = d0;
                }
                
                if(pts[j].iGroup != pts[k].iGroup) {
                    d1 = dist2(pts[j], pts[k]);
                    if(d1 < min_d1) { min_d1 = d1; group = pts[k].iGroup; }
                }
            }
        }
        
        if (min_d0 > min_d1)
            pts[j].iGroup = group;
    }
}



void CClustering :: run()
{    
    int i, j, iMaxNum, ntmp;

    vector<int> vectmp;
    vector<string> vectmp1;
    
    cout << " Start grouping the molecule into clusters... " << endl;
    
    if(bClusterAuto)    { cout << " User input cluster number is  AUTO . " << endl; cout << endl; }
    else                { cout << " User input cluster number is " << n_cluster << " . " << endl; cout << endl; }
    
    
    importcenter();
    
    
    if(!bClusterAuto) {
        
        vector<PDBFORM> c = centdetect(vecPDB, vecPDB.size(), n_cluster);
        
        checkcluster(vecPDB, vecPDB.size(), n_cluster);
    
    }
    
    if(bClusterAuto) {
        
        n_cluster = ( vecPDB.size() / 16 ) + 1 ;
        
        iMaxNum = 16;
        
        while (iMaxNum > 15) {
            
            
            vector<PDBFORM> c = centdetect(vecPDB, vecPDB.size(), n_cluster);
            
            checkcluster(vecPDB, vecPDB.size(), n_cluster);
            
            iMaxNum = 0;
            
            for (i=0; i<n_cluster; i++) {
                ntmp = 0;
                for (j=0; j<vecPDB.size(); j++) {
                    if(vecPDB[j].iGroup == i)   ntmp++;
                }
                if(ntmp > iMaxNum)  iMaxNum = ntmp;
            }
            
            n_cluster++;
        
        }
        
    }
    
    if(bClusterAuto)    n_cluster--;
    
    for (i=0; i<n_cluster; i++) {
        for (j=0; j<vecPDB.size(); j++) {
            if(vecPDB[j].iGroup != i)  continue;
            
            vectmp.push_back(vecPDB[j].res_num);
            vectmp1.push_back(vecPDB[j].chain_id);
        }
        if (!vectmp.empty()) {
            vecCluster.push_back(vectmp);
            vecClusterChID.push_back(vectmp1);
        }
        vectmp.clear();
        vectmp1.clear();
        
    }
    
    
    
#ifdef PRIME_DEBUG

    cout << " The Cluster Number is " << vecCluster.size() << " .     " << endl;
 
    ofstream clusterOut;
    clusterOut.open("clusters.txt");
    
    for (i=0; i<vecCluster.size(); i++) {
        clusterOut << " Cluster " << setw(2) << right << i+1 << " : ";
        for (j=0; j<vecCluster[i].size(); j++) {
            clusterOut << vecClusterChID[i][j];
            clusterOut << setfill('0') << setw(4) << vecCluster[i][j] << "  ";clusterOut << setfill(' ');
        }
        clusterOut << endl;
    }
    
    
    cout << " Clusters information is saved in  clusters.txt  file ." << endl;

#endif
    

#ifdef PRIME_DEBUG
    
    int atom_serial_num = 1;
    for(auto itr : vecPDB) {
        
        cout << left << setw(6) << "ATOM"
            << right << setw(5) << atom_serial_num << " "
            << left  << setw(4) << itr.atom_name << " "
            << left  << setw(3) << itr.res_name << " "
            << right << setw(1) << itr.chain_id << " "
            << right << setw(4) << itr.res_num << "   "
            << right << setw(8) << fixed << setprecision(3) << itr.coord.X
            << right << setw(8) << fixed << setprecision(3) << itr.coord.Y
            << right << setw(8) << fixed << setprecision(3) << itr.coord.Z << "    "
            << right << setw(1) << itr.iGroup    // the last column is the group number in the clusters.
            << endl;
        atom_serial_num++;
        
    }
    
#endif
    
}