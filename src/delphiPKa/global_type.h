//
//  global_type.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*  
 *  This class is used to create some structs and general string trimming
 */

#ifndef GLOBAL_TYPE_H_
#define GLOBAL_TYPE_H_

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "prime_environment.h"
using namespace std;



struct nVECTOR
{
    float X, Y, Z;
};


struct TOPOLOGY {
    string resname;
    string atom;
    string orbital;
    string bondatom1, bondatom2, bondatom3, bondatom4;
    string conf;
};



struct PDBFORM {
    string res_name;
    string atom_name;
    string chain_id;
    int res_num;
    nVECTOR  coord;
    string conf;
    float fCharge;
    float fRidus;
    int iGroup;
    bool bIonizable;
    bool bCharged;
};

struct Energy {
    int resnum1;
    int resnum2;
    float fEnergy;
};

struct respair {
    int res1;
    int res2;
};

typedef map<string, TOPOLOGY> mapTopology;

typedef map<string, PDBFORM> mapPDBFORM;



//------------------- overload operator + -------------------//

nVECTOR operator+ (const nVECTOR& Vector1, const nVECTOR& Vector2);



nVECTOR operator+ (const nVECTOR& Vector, const float& Value);


nVECTOR operator+ (const float& Value, const nVECTOR& Vector);




//------------------- overload operator - -------------------//

nVECTOR operator- (const nVECTOR& Vector1, const nVECTOR& Vector2);


nVECTOR operator- (const nVECTOR& Vector, const float& Value);



nVECTOR operator- (const float& Value, const nVECTOR& Vector);



nVECTOR  operator- (const nVECTOR& Vector);



//------------------- overload operator * -------------------//

nVECTOR operator* (const nVECTOR& Vector1, const nVECTOR& Vector2);


nVECTOR operator* (const nVECTOR& Vector, const float& Value);


nVECTOR operator* (const float& Value, const nVECTOR& Vector);


//------------------- overload operator / -------------------//


nVECTOR operator/ (const nVECTOR& Vector, const float& Value);



//------------------- special operation on nVECTOR -------------------//

float optDot(const nVECTOR& Vector1, const nVECTOR& Vector2);



nVECTOR optCross(const nVECTOR& Vector1, const nVECTOR& Vector2);


nVECTOR optNorm(const nVECTOR& Vector);


nVECTOR optSqrt(const nVECTOR& Vector);


nVECTOR optABS(const nVECTOR& Vector);


//------------------------------------------------------------------//


extern ostream& operator<< (ostream& os, TOPOLOGY& topology_data);
extern ostream& operator<< (ostream& os, PDBFORM& pdbform_data);
extern ostream& operator<< (ostream& os, nVECTOR& nVector_data);




void remove_all_whitespace_string(string& str);

void remove_lead_whitespace_string(string& str);

void remove_tail_whitespace_string(string& str);

void remove_leadtail_whitespace_string(string& str);

string key_to_hashmapQR(const PDBFORM& inPDB);



#endif //GLOBAL_TYPE_H_
