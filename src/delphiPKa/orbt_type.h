//
//  orbt_type.h
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

/*
 *  This class is used to realize various atomic orbital types for protonation (adding position of hydrogens)
 */


#ifndef ORBT_TYPE_
#define ORBT_TYPE_

#include <iostream>
#include "prime_environment.h"
#include "global_type.h"
using namespace std;


void sp3_paint1(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, const nVECTOR& Vector3, nVECTOR& Vector4, const float& BondLenth_04);


void sp3_paint2(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, nVECTOR& Vector3, nVECTOR& Vector4, const float& BondLenth_03, const float& BondAngle_304);


void sp3_paint3(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, nVECTOR& Vector3, nVECTOR& Vector4, nVECTOR& Vector5, const float& BondLenth_03, const float& BondAngle_103, const float& TorsionAngle_3012);


void sp2_paint1(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, nVECTOR& Vector3, const float& BondLenth_03);


void sp2_paint2(const nVECTOR& Vector0, const nVECTOR& Vector1, const nVECTOR& Vector2, nVECTOR& Vector3, nVECTOR& Vector4, const float& BondLenth_03, const float& BondAngle_103, const float& TorsionAngle_3012);




#endif // ORBT_TYPE_
