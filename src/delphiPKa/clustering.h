//
//  clustering.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*  This class is used to cluster the system,
 *  which includes kmean (original, do not use any more) and network partition (current).
 *
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include <iostream>
#include <fstream>
#include "prime_environment.h"
#include "global_type.h"
#include "data_store.h"
using namespace std;

class CClustering {
private:
    
    /* From Data_STORE */
    
    const vector<PDBFORM>& newPDB;
    
    int& n_cluster;
    bool& bClusterAuto;
    vector<vector<int> >& vecCluster;
    vector<vector<string> >& vecClusterChID;
    
    /* END */
    
    vector<PDBFORM> vecPDB;
    
    float randf(float m); // generate random float number
    
    float dist2(PDBFORM a, PDBFORM b); // calculate the distance of two atoms
    
    int nearest(PDBFORM pt, vector<PDBFORM> cent, int n_cluster, float *d2); // find the nearest centroid
    
    vector<PDBFORM> centdetect(vector<PDBFORM>& pts, int len, int n_cent); // find the center of the current cluster
    
    void kmeanpp(vector<PDBFORM>& pts, int len, vector<PDBFORM> cent, int n_cent); // k-Mean algorithm
    
    void checkcluster(vector<PDBFORM>& pts, int len, int n_cluster); // check if the cluster reached convergence
    
    void importcenter();
    
public:
    CClustering(shared_ptr<DATA_STORE> pData) :
    
    newPDB(pData->newPDB),
    n_cluster(pData->n_cluster),
    bClusterAuto(pData->bClusterAuto),
    vecCluster(pData->vecCluster),
    vecClusterChID(pData->vecClusterChID)
    
    {};
    
    void run();
    
    ~CClustering() {};
    
};

#endif // CLUSTERING_H_
