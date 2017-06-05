//
//  global_type.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "global_type.h"

ostream& operator<< (ostream& os, nVECTOR& nVector_data)
{
    os << "(" << nVector_data.X << ", " << nVector_data.Y << ", " << nVector_data.Z << ")";
    
    return os;
}


ostream& operator<< (ostream& os, TOPOLOGY& topology_data)
{
    os << topology_data.resname << " " << topology_data.atom << " " << topology_data.orbital << " " << topology_data.bondatom1 << " " << topology_data.bondatom2 << " " << topology_data.bondatom3 << " " << topology_data.bondatom4;
    return os;
}


ostream& operator<< (ostream& os, PDBFORM& pdbform_data)
{
    os << pdbform_data.atom_name << " " << pdbform_data.res_name << " " << pdbform_data.chain_id << " " << pdbform_data.res_num << " " << pdbform_data.coord.X << " " << pdbform_data.coord.Y << " " << pdbform_data.coord.Z;
    return os;
}



nVECTOR operator+ (const nVECTOR& Vector1, const nVECTOR& Vector2)
{
    nVECTOR  Vector3;
    Vector3.X = Vector1.X + Vector2.X;
    Vector3.Y = Vector1.Y + Vector2.Y;
    Vector3.Z = Vector1.Z + Vector2.Z;
    return Vector3;
};


nVECTOR operator+ (const nVECTOR& Vector, const float& Value)
{
    nVECTOR ValuePlusVECTOR;
    ValuePlusVECTOR.X = Vector.X + Value;
    ValuePlusVECTOR.Y = Vector.Y + Value;
    ValuePlusVECTOR.Z = Vector.Z + Value;
    return ValuePlusVECTOR;
};

nVECTOR operator+ (const float& Value, const nVECTOR& Vector)
{
    nVECTOR ValuePlusVECTOR;
    ValuePlusVECTOR.X = Value + Vector.X;
    ValuePlusVECTOR.Y = Value + Vector.Y;
    ValuePlusVECTOR.Z = Value + Vector.Z;
    return ValuePlusVECTOR;
};



//------------------- overload operator - -------------------//

nVECTOR operator- (const nVECTOR& Vector1, const nVECTOR& Vector2)
{
    nVECTOR Vector3;
    Vector3.X = Vector1.X - Vector2.X;
    Vector3.Y = Vector1.Y - Vector2.Y;
    Vector3.Z = Vector1.Z - Vector2.Z;
    return Vector3;
};

nVECTOR operator- (const nVECTOR& Vector, const float& Value)
{
    nVECTOR VECTORMinusValue;
    VECTORMinusValue.X = Vector.X - Value;
    VECTORMinusValue.Y = Vector.Y - Value;
    VECTORMinusValue.Z = Vector.Z - Value;
    return VECTORMinusValue;
};


nVECTOR operator- (const float& Value, const nVECTOR& Vector)
{
    nVECTOR ValueMinusVECTOR;
    ValueMinusVECTOR.X = Value - Vector.X;
    ValueMinusVECTOR.Y = Value - Vector.Y;
    ValueMinusVECTOR.Z = Value - Vector.Z;
    return ValueMinusVECTOR;
};


nVECTOR  operator- (const nVECTOR& Vector)
{
    nVECTOR Vector2;
    Vector2.X = - Vector.X;
    Vector2.Y = - Vector.Y;
    Vector2.Z = - Vector.Z;
    return Vector2;
};


//------------------- overload operator * -------------------//

nVECTOR operator* (const nVECTOR& Vector1, const nVECTOR& Vector2)
{
    nVECTOR Vector3;
    Vector3.X = Vector1.X * Vector2.X;
    Vector3.Y = Vector1.Y * Vector2.Y;
    Vector3.Z = Vector1.Z * Vector2.Z;
    return Vector3;
};

nVECTOR operator* (const nVECTOR& Vector, const float& Value)
{
    nVECTOR VECTORMultiplyValue;
    VECTORMultiplyValue.X = Vector.X * Value;
    VECTORMultiplyValue.Y = Vector.Y * Value;
    VECTORMultiplyValue.Z = Vector.Z * Value;
    return VECTORMultiplyValue;
};

nVECTOR operator* (const float& Value, const nVECTOR& Vector)
{
    nVECTOR ValueMultiplyVECTOR;
    ValueMultiplyVECTOR.X = Value * Vector.X;
    ValueMultiplyVECTOR.Y = Value * Vector.Y;
    ValueMultiplyVECTOR.Z = Value * Vector.Z;
    return ValueMultiplyVECTOR;
};

//------------------- overload operator / -------------------//


nVECTOR operator/ (const nVECTOR& Vector, const float& Value)
{
    nVECTOR VECTORDivideValue;
    VECTORDivideValue.X = Vector.X / Value;
    VECTORDivideValue.Y = Vector.Y / Value;
    VECTORDivideValue.Z = Vector.Z / Value;
    return VECTORDivideValue;
};


//------------------- special operation on nVECTOR -------------------//

float optDot(const nVECTOR& Vector1, const nVECTOR& Vector2)
{
    float DotValue;
    DotValue = Vector1.X * Vector2.X + Vector1.Y * Vector2.Y + Vector1.Z * Vector2.Z;
    return DotValue;
    
}


nVECTOR optCross(const nVECTOR& Vector1, const nVECTOR& Vector2)
{
    nVECTOR Vector3;
    Vector3.X = Vector1.Y * Vector2.Z - Vector1.Z * Vector2.Y;
    Vector3.Y = Vector1.Z * Vector2.X - Vector1.X * Vector2.Z;
    Vector3.Z = Vector1.X * Vector2.Y - Vector1.Y * Vector2.X;
    return Vector3;
};

nVECTOR optNorm(const nVECTOR& Vector)
{
    nVECTOR Vector2;
    float dist = sqrt(Vector.X * Vector.X + Vector.Y * Vector.Y + Vector.Z * Vector.Z);
    if(dist < 1.0e-20) {
        Vector2.X = 0.0;
        Vector2.Y = 0.0;
        Vector2.Z = 0.0;
    }
    else {
        Vector2.X = Vector.X / dist;
        Vector2.Y = Vector.Y / dist;
        Vector2.Z = Vector.Z / dist;
    }
    
    return Vector2;
};



nVECTOR optSqrt(const nVECTOR& Vector)
{
    nVECTOR Vector2;
    Vector2.X = sqrt(Vector.X);
    Vector2.Y = sqrt(Vector.Y);
    Vector2.Z = sqrt(Vector.Z);
    return Vector2;
};


nVECTOR optABS(const nVECTOR& Vector)
{
    nVECTOR Vector2;
    Vector2.X = Vector.X > 0 ? Vector.X : -Vector.X;
    Vector2.Y = Vector.Y > 0 ? Vector.Y : -Vector.Y;
    Vector2.Z = Vector.Z > 0 ? Vector.Z : -Vector.Z;
    return Vector2;
    
};




void remove_all_whitespace_string(string& str) {
    str.erase(std::remove(str.begin(), str.end(), ' ' ), str.end());
}

void remove_lead_whitespace_string(string& str) {
    str.erase(0,str.find_first_not_of(' '));
}

void remove_tail_whitespace_string(string& str) {
    str.erase(str.find_last_not_of(' ')+1,str.length()-str.find_last_not_of(' ')-1);
}

void remove_leadtail_whitespace_string(string& str) {
    str.erase(0,str.find_first_not_of(' '));
    str.erase(str.find_last_not_of(' ')+1,str.length()-str.find_last_not_of(' ')-1);
}

string key_to_hashmapQR(const PDBFORM& inPDB) {
    string strtmp = inPDB.res_name + " " + inPDB.atom_name;
    return strtmp;
}

