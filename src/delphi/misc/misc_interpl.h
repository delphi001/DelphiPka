/**
 * @file misc_interpl.h
 * @brief interpolate grid value
 *
 * @author Chuan Li, chuanli@clemson.edu
 */

#ifndef INTERPL_H_
#define INTERPL_H_

#include <vector>
#include <math.h> /* floor */

#include "../interface/environment.h"
#include "../interface/exceptions.h"
#include "misc_grid.h"

using namespace std;

/**
 * function to interpolate the value at a given pt in fmap based on its neighboring values
 */
extern delphi_real interpl(const delphi_integer& ieExtrema, delphi_real *** fMap, const SGrid<delphi_real>& gPoint);

/**
 * function to interpolate the value at a given pt in fmap based on the bool types of its neighbors
 */
extern delphi_real boolinterpl(const delphi_integer& igMaxGrid, const vector<bool>& prgbMap, const SGrid<delphi_real>& gPoint);

#endif // INTERPL_H_
