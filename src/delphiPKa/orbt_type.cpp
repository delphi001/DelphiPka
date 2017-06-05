//
//  orbt_type.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "orbt_type.h"
using namespace std;



void sp3_paint1(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, const nVECTOR& Vector3, nVECTOR& Vector4, const float& BondLenth_04)
{
    
    nVECTOR normVector01, normVector02, normVector03, normVector04;
    normVector01 = optNorm(Vector1 - Vector0);
    normVector02 = optNorm(Vector2 - Vector0);
    normVector03 = optNorm(Vector3 - Vector0);
    normVector04 = optNorm( - (normVector01 + normVector02 + normVector03));
    Vector4 = Vector0 + normVector04 * BondLenth_04;
    
}



void sp3_paint2(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, nVECTOR& Vector3, nVECTOR& Vector4, const float& BondLenth_03, const float& BondAngle_304)
{
    
    float semiAngle = BondAngle_304 / 2.;
    nVECTOR normVector01, normVector02, normVector102, normBisectVector304;
    normVector01 = optNorm(Vector1 - Vector0);
    normVector02 = optNorm(Vector2 - Vector0);
    normBisectVector304 = optNorm( - (normVector01 + normVector02));
    normVector102 = optNorm( optCross(normVector01, normVector02));
    Vector3 = Vector0 + ((normBisectVector304 * cos(semiAngle) + normVector102 * sin(semiAngle)) * BondLenth_03);
    Vector4 = Vector0 + ((normBisectVector304 * cos(semiAngle) - normVector102 * sin(semiAngle)) * BondLenth_03);
}



void sp3_paint3(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, nVECTOR& Vector3, nVECTOR& Vector4, nVECTOR& Vector5, const float& BondLenth_03, const float& BondAngle_103, const float& TorsionAngle_3012)
{
    
    float theta, phi;
    nVECTOR normVector10, normVector12, normVector012, normVector01_012;
    theta = PI - BondAngle_103;
    
    normVector10 = optNorm(Vector0 - Vector1);
    normVector12 = optNorm(Vector2 - Vector1);
    normVector012 = optNorm(optCross(normVector10, normVector12));
    normVector01_012 = optNorm(optCross(normVector10, normVector012));
    
    phi = TorsionAngle_3012 - PI/2.0;
    Vector3 = Vector0 + (normVector012*cos(phi)*sin(theta) + normVector01_012*sin(phi)*sin(theta) + normVector10*cos(theta))*BondLenth_03;
    
    
    phi = TorsionAngle_3012 + PI/6.0;
    Vector4 = Vector0 + (normVector012*cos(phi)*sin(theta) + normVector01_012*sin(phi)*sin(theta) + normVector10*cos(theta))*BondLenth_03;
    
    
    phi = TorsionAngle_3012 + PI*5./6.;
    Vector5 = Vector0 + (normVector012*cos(phi)*sin(theta) + normVector01_012*sin(phi)*sin(theta) + normVector10*cos(theta))*BondLenth_03;
    
}




void sp2_paint1(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, nVECTOR& Vector3, const float& BondLenth_03)
{
    nVECTOR normVector01, normVector02, normVector03;
    
    normVector01 = optNorm(Vector1 - Vector0);
    normVector02 = optNorm(Vector2 - Vector0);
    normVector03 = optNorm( - (normVector01 + normVector02));
    
    Vector3 = Vector0 + normVector03 * BondLenth_03;
    
}



void sp2_paint2(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, nVECTOR& Vector3, nVECTOR& Vector4, const float& BondLenth_03, const float& BondAngle_103, const float& TorsionAngle_3012)
{
    nVECTOR normVector10, normVector12, normVector012, normVector01_012;
    float theta, phi;
    
    normVector10 = optNorm(Vector0 - Vector1);
    normVector12 = optNorm(Vector2 - Vector1);
    
    normVector012 = optNorm(optCross(normVector10, normVector12));
    normVector01_012 = optNorm(optCross(normVector10, normVector012));
    
    theta = PI - BondAngle_103;
    
    phi = TorsionAngle_3012 - PI/2.0;
    Vector3 = Vector0 + (normVector012*cos(phi)*sin(theta) + normVector01_012*sin(phi)*sin(theta) + normVector10*cos(theta)) * BondLenth_03;
    
    phi = TorsionAngle_3012 + PI/2.0;
    Vector4 = Vector0 + (normVector012*cos(phi)*sin(theta) + normVector01_012*sin(phi)*sin(theta) + normVector10*cos(theta)) * BondLenth_03;
    
}