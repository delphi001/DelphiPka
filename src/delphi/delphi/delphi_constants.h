/**
 * @file delphi_constants.h
 * @brief globally defined constants used for delphicpp
 *
 * @author Chuan Li, chuanli@clemson.edu
 */

#ifndef DELPHI_CONSTANTS_H_
#define DELPHI_CONSTANTS_H_

#include <vector>
#include <string>
#include "../interface/environment.h"
#include "../misc/misc_grid.h"

using namespace std;

#ifdef MX
const delphi_real fPi    =  3.141592653589793; ///< pi
const delphi_real f4Pi   =  1.256637061435917e+1; ///< 4*pi
const delphi_real fSixth =  1.666666666666667e-1; ///< 1/6
const delphi_real fZero  =  1.0e-15; ///< practical zero for comparing delphi_real numbers
#endif

#ifdef SP
const delphi_real fPi    =  3.141593;
const delphi_real f4Pi   =  1.256637e+1;
const delphi_real fSixth =  1.666667e-1;
const delphi_real fZero  =  1.0e-6;
#endif

#ifdef DP
const delphi_real fPi    =  3.141592653589793;
const delphi_real f4Pi   =  1.256637061435917e+1;
const delphi_real fSixth =  1.666666666666667e-1;
const delphi_real fZero  =  1.0e-15;
#endif

#ifdef LD
const delphi_real fPi    =  3.1415926535897932384626433832795;
const delphi_real f4Pi   =  1.2566370614359172953850573533118e+1;
const delphi_real fSixth =  1.6666666666666666666666666666667e-1;
const delphi_real fZero  =  1.0e-15;
#endif

const int    iStatementNum = 51; ///< number of statements increased from 45 to 51(Gaussian MEMPOT)
const int    iFunctionNum_FullName  = 11; ///< number of functions in full name
const int    iFunctionNum_ShortName = 5; ///< number of functions in short name

const int    iDMax = 50;
const double dAbsoluteZero = -273.15; /// temperature of absolute zero
const double dAtomicUnitCrg = 1.602176487e-19; ///< e, Coulomb C
const double dBoltzmannConstant = 1.3806504e-23; ///< k, Joule per Kelvin JK^(-1)
const double dVacuumPermittivity = 8.8541878176e-12; ///< e0, farads per meter Fm^(-1)
const double dEPK = 167100.9162872952; ///< e^2/(4*pi*e0*k*1.0e-10)

const char rgcFileMap[] = "qdiffxas: qdiffxs4 with an improved surfacing routine";

/**
 * struct used to pass atom info from delphicpp to mcce
 */
typedef struct {
   delphi_real x;
   delphi_real y;
   delphi_real z;
   delphi_real crg;
   delphi_real size;
} SAtomInfo;

/**
 * struct used to pass calculated values between delphicpp and mcce
 */
typedef struct {
   //------- IN
   delphi_integer gsize;
   delphi_real scale;
   float grids_per_ang;
   int del_runs;
   int n_retry; // default = 0, reset delphi failure counter for this conformer
   string pdbfile;
   int pdbformat;
   delphi_real indi;
   delphi_real exdi;
   delphi_real ionrad;
   delphi_real prbrad;
   delphi_real salt;
   int  bndcon;
   delphi_real center[3];
   //string frcfile;
   //string phifile;
   string uniqID;

   //-------- INTERMEDIATE
   vector<delphi_real> phimap;
   delphi_real scale1;
   SGrid<delphi_real> oldmid1;
   delphi_integer igrid1;

   //-------- OUT
   delphi_real ergs;
   vector<delphi_real> phiv;
   char del_err;
} SMCCE;


/**
 * struct used to pass calculated values between delphicpp and PrimePKA
 */

typedef struct {
    //-------- Input
    
    int ntimes;
    int iGaussian;
    delphi_real fSigma;
    delphi_real fSrfcut;
    delphi_integer gsize;
    delphi_real scale;
    delphi_real perfil;
    delphi_real fMaxc;
    string pdbfile;
    int pdbformat;
    int  bndcon;
    bool bFrcPDBitself;
    bool bCommFrc;
    bool bAcent;
    bool bSolvEng;
    delphi_real indi;
    delphi_real exdi;
    delphi_real ionrad;
    delphi_real prbrad;
    delphi_real salt;
    delphi_real center[3];
    string uniqID;
    vector<string> vecStrPDB;
    vector<string> vecStrFRCIn;
    vector<string> strAtomDes;
    vector<delphi_real> vecGridPot;
    vector<delphi_real> vecSiteCrg;
    vector<delphi_real> vecAtomPot;
    string strFRCIn;
    string strFRCOut;
    SGrid<delphi_real> oldmid;

    //-------- INTERMEDIATE
    vector<delphi_real> phimap;
    delphi_real scale1;
    SGrid<delphi_real> oldmid1;
    delphi_integer igrid1;


    //-------- Output
    delphi_real ergs;
    delphi_real ergg;

} SPrime;


#endif // DELPHI_CONSTANTS_H_

