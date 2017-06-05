/*
 * delphi_datamarshal_setDefault.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: chuan
 */

#include "delphi_datamarshal.h"

void CDelphiDataMarshal::setDefault()
{
   //----------------------- set by Statments ------------------------//
   bAutoConverge        = true;         // iautocon AUTOC
   iBndyType            = 2;            // ibctyp   BNDCON
   fPercentageFill      = 80.0;         // perfil   PERFIL
   bFixedRelaxParam     = false;        // icheb    CHEBIT
   bOutGraspSurf        = false;        // isrf     CLCSRF
   iConvergeFract       = 1;            // icon2    CONFRA
   iIterateInterval     = 10;           // icon1    CONINT
   bExitUniformDielect  = false;        // iexun    EXITUN
   fExDielec            = 80;           // repsout  EXDI
   bCrgInterplateType   = false;        // isph     FCRG
   fGridConverge        = 0.0;          // gten     GRDCON
   iGrid                = 0;            // igrid    GSIZE
   fInDielec            = 2.0;          // repsin   INDI
   vctfSalt.push_back(0.0);             // conc     SALT,SALT2
   vctfSalt.push_back(0.0);
   fIonRadius           = 2.0;          // exrad    IONRAD
   iLinIterateNum       = 0;            // nlit     LINIT
   bLogGraph            = false;        // igraph   LOGGRP
   bLogPotential        = false;        // ipoten   LOGPOT
   fMaxc                = 0.0;          // res2     MAXC
   iNonIterateNum       = 0;            // nnit     NONIT
   vctbPeriodicBndy.push_back(false);   // iper     PBX,PBY,PBZ
   vctbPeriodicBndy.push_back(false);
   vctbPeriodicBndy.push_back(false);
   vctbPeriodicBndy.push_back(false);
   vctbPeriodicBndy.push_back(false);
   vctbPeriodicBndy.push_back(false);
   bOutCrgDensity       = false;        // iconc    PHICON
   vctfProbeRadius.push_back(1.4);      // radprb   PRBRAD,RADPR2
   vctfProbeRadius.push_back(-1.0);
   fSpectralRadius      = 0.9975;       // uspec    RELFAC
   fRelaxParam          = 1.0;          // relpar   RELPAR
   fRmsc                = 0.0;          // res1     RMSC
   fScale               = 2.0;          // scale    SCALE
   bSolvePB             = true;         // isolv    SOLVPB
   vctiValence1.push_back(1);           // ival     VAL+1,VAL-1
   vctiValence1.push_back(1);
   vctiValence2.push_back(0);           // ival2    VAL+2,VAL-2
   vctiValence2.push_back(0);
   fPotentialUpperBond  = 0.5;          // atompotdist  ATPODS
   fTemper              = 297.3342119;  // temperature  TEMPER
   gfPotentialDrop.nX   = 0.0;          // vdrop
   gfPotentialDrop.nY   = 0.0;
   gfPotentialDrop.nZ   = 0.0;
   bSpectralRadius      = false;        // iuspec
   bManualRelaxParam    = false;        // imanual
   iPhiInType           = 0;
   //-------------------Lin Li: Gaussian & MEMPOT -----------------//
   fCutoff              = 1.0;
   fSigma               = 1.0;
   iInhomo              = 0;
   fSrfcut              = 20.0;
   iGaussian            = 0;
   fRadipz              = -1.0;

   //-------------------------- io file names ------------------------//
   //strParamFile       = "fort.10";    // prmnam
   strSizeFile          = "fort.11";    // siznam
   strCrgFile           = "fort.12";    // crgnam
   strPdbFile           = "fort.13";    // pdbnam
   strPhiFile           = "fort.14";    // phinam
   strFrciFile          = "fort.15";    // frcinam
   strFrcFile           = "fort.16";    // frcnam
   strEpsFile           = "fort.17";    // epsnam
   strPhiiFile          = "fort.18";    // phiinam
   strModifiedPdbFile   = "fort.19";    // mpdbnam
   strUnformatPdbFile   = "fort.20";    // updbnam
   strUnformatFrcFile   = "fort.21";    // ufrcnam
   strGraspFile         = "grasp.srf";  // srfnam
   strEnergyFile        = "energy.dat"; // nrgnam
   strScrgFile          = "scrg.dat";   // scrgnam

   //----------------------- set by functions ------------------------//
   // set by CENTER or CENT function:
   gfOffCenter.nX       = 0.0;          // offset
   gfOffCenter.nY       = 0.0;
   gfOffCenter.nZ       = 0.0;

   // set by ACENTER or ACENT function
   gfAcent.nX           = 0.0;          // acent
   gfAcent.nY           = 0.0;
   gfAcent.nZ           = 0.0;
   bIsAcent             = false;        // iacent

   // set by READ or IN function
   iPdbFormatIn         = 10;           // pdbfrm
   bPdbUnformatIn       = false;        // ipdbrd

   // set by WRITE or OUT function
   bPhimapOut           = false;        // phiwrt
   iPhiFormatOut        = 0;            // phifrm
   bBiosystemOut        = false;        // ibios
   bBemSrfOut           = false;        // ibem
   bSiteOut             = false;        // isite
   iFrcFormatOut        = 0;            // frcfrm
   bEpsOut              = false;        // epswrt
   bModPdbOut           = false;        // iatout
   iModPdbFormatOut     = 0;            // mpdbfrm
   bUnformatPdbOut      = false;        // ipdbwrt
   bUnformatFrcOut      = false;        // ifrcwrt
   bEngOut              = false;        // inrgwrt
   bGridCrgOut          = false;        // iwgcrg
   bHsurf2DatOut        = false;        // iacs
   bDbOut               = false;        // idbwrt
   bSurfEngOut          = false;        // isen
   bSurfCrgOut          = false;        // isch
   iSurfCrgFormatOut    = 0;            // scrgfrm

   // set by ENERGY function
   bGridEng             = false;        // logg
   bSolvEng             = false;        // logs
   bAnalySurfEng        = false;        // logas
   bAnalyEng            = false;        // loga
   bIonsEng             = false;        // logions
   bCoulombEng          = false;        // logc

   // set by SITE function: all MUST be initialized to to false
   bAtomInSite          = false;        // isita
   bCrgInSite           = false;        // isitq
   bGridPotentialInSite = false;        // isitp
   bAtomPotentialInSite = false;        // isitap
   bDebyeFractionInSite = false;        // isitdeb
   bFieldInSite         = false;        // isitf
   bReactPotentialInSite= false;        // isitr
   bCoulombPotentialInSite = false;     // isitc
   bAtomCoordInSite     = false;        // isitx
   bSaltInSite          = false;        // isiti
   bTotalPotentialInSite= false;        // isitt
   bReactForceInSite    = false;        // isitrf
   bCoulombForceInSite  = false;        // isitcf
   bMDInSite            = false;        // isitmd
   bSurfCrgInSite       = false;        // isitsf
   bTotalForceInSite    = false;        // isittf
   bPotentialInSite     = false;        // isitpot
   bReactFieldInFRC     = false;        // irea
   bPDB2FRCInSite       = false;        // iself
   bCommFRCIn           = false;        // isitcomm    Update: Dec 19, 2014 by Lin Wang

   // set by BUFZ function
   eiBuffz.nMin.nX      = 0;            // bufz
   eiBuffz.nMin.nY      = 0;
   eiBuffz.nMin.nZ      = 0;
   eiBuffz.nMax.nX      = 0;
   eiBuffz.nMax.nY      = 0;
   eiBuffz.nMax.nZ      = 0;
   bIsBuffz             = false;        // ibufz

   // set by SURFACE function
   iTypeSurf            = -1;           // iTypeSurf

   //------------------------------ DelPhi ---------------------------//
   fDebyeLength         = 0.0;          // deblen
   fEpsOut              = 80.0;         // epsout
   gfCoordinateRange.nX = 0.0;          // cran
   gfCoordinateRange.nY = 0.0;
   gfCoordinateRange.nZ = 0.0;
   gfGeometricCenter.nX = 0.0;          // pmid
   gfGeometricCenter.nY = 0.0;
   gfGeometricCenter.nZ = 0.0;
   gfBoxCenter.nX       = 0.0;          // oldmid
   gfBoxCenter.nY       = 0.0;
   gfBoxCenter.nZ       = 0.0;
   fIonStrength         = 0.0;          // rionst
   fTaylorCoeff1        = 0.0;          // chi1
   fTaylorCoeff2        = 0.0;          // chi2
   fTaylorCoeff3        = 0.0;          // chi3
   fTaylorCoeff4        = 0.0;          // chi4
   fTaylorCoeff5        = 0.0;          // chi5
   bNonlinearEng        = false;        // lognl
   fEPKT                = 0.0;          // epkt
   fEpsIn               = 2.0;          // epsin
   bFrcUnformatIn       = false;        // ifrcrd
   iDirectEpsMap        = 1;            // iDirectEpsMap
   iMoleculeNum         = 0;            // numbmol
   fMaxRadius           = 0.01;         // rdmx
   bUniformDielec       = true;         // uniformdiel

   //-------------------------------- IO -----------------------------//
   iResidueNum          = 0;            // resnummax
   iMediaNum            = 1;            // nmedia
   iObjectNum           = 1;            // nobject
   iAtomNum             = 0;            // natom
   bOnlyMolecule        = true;         // ionlymol

   //------------------------------ Surface --------------------------//
   iCrgGridNum          = 0;            // nqass
   fNetCrg              = 0.0;          // qnet
   fMinusCrg            = 0.0;          // qmin
   fPlusCrg             = 0.0;          // qplus
   gfPlusCrgCenter.nX   = 0.0;          // cqplus
   gfPlusCrgCenter.nY   = 0.0;
   gfPlusCrgCenter.nZ   = 0.0;
   gfMinusCrgCenter.nX  = 0.0;          // cqmin
   gfMinusCrgCenter.nY  = 0.0;
   gfMinusCrgCenter.nZ  = 0.0;
   gfMinCoordinate.nX   = 6000.0;       // cmin
   gfMinCoordinate.nY   = 6000.0;
   gfMinCoordinate.nZ   = 6000.0;
   gfMaxCoordinate.nX   =-6000.0;       // cmax
   gfMaxCoordinate.nY   =-6000.0;
   gfMaxCoordinate.nZ   =-6000.0;
   iBndyGridNum         = 0;            // ibnum
   iCrg2GridNum         = 0;            // nqgrd

   //------------------------------ Solver ---------------------------//
   iDielecBndySum       = 0;            // icount2b
   iCrgedGridSum        = 0;            // icount1b
   iCrgBdyGrid          = 0;            // ibc

   //------------------------------ Energy ---------------------------//
   fEngGrid             = 0.0;          // test_ergg
   fEngCoul             = 0.0;          // test_ergc
   fEngCorrect          = 0.0;          // test_ergs
   fEngReact            = 0.0;          // test_ergr
   fEngIons             = 0.0;          // test_ergions

   //---------------------------- statements -------------------------//
   rgstrStatement_ShortForm[0]    = "UNUSED";
   rgstrStatement_ShortForm[1]    = "GSIZE";
   rgstrStatement_ShortForm[2]    = "SCALE";
   rgstrStatement_ShortForm[3]    = "PERFIL";
   rgstrStatement_ShortForm[4]    = "INDI";
   rgstrStatement_ShortForm[5]    = "EXDI";
   rgstrStatement_ShortForm[6]    = "PRBRAD";
   rgstrStatement_ShortForm[7]    = "IONRAD";
   rgstrStatement_ShortForm[8]    = "SALT";
   rgstrStatement_ShortForm[9]    = "BNDCON";
   rgstrStatement_ShortForm[10]   = "LINIT";
   rgstrStatement_ShortForm[11]   = "NONIT";
   //rgstrStatement_ShortForm[12] = "MEMDAT"; // OBSOLELE. REMOVED FROM THE LIST.
   rgstrStatement_ShortForm[12]   = "UNUSED";
   rgstrStatement_ShortForm[13]   = "FCRG";
   rgstrStatement_ShortForm[14]   = "LOGPOT";
   rgstrStatement_ShortForm[15]   = "LOGGRP";
   rgstrStatement_ShortForm[16]   = "CONINT";
   rgstrStatement_ShortForm[17]   = "CONFRA";
   rgstrStatement_ShortForm[18]   = "PBX";
   rgstrStatement_ShortForm[19]   = "PBY";
   rgstrStatement_ShortForm[20]   = "PBZ";
   rgstrStatement_ShortForm[21]   = "AUTOC";
   rgstrStatement_ShortForm[22]   = "EXITUN";
   rgstrStatement_ShortForm[23]   = "GRDCON";
   rgstrStatement_ShortForm[24]   = "RELFAC";
   rgstrStatement_ShortForm[25]   = "CHEBIT";
   rgstrStatement_ShortForm[26]   = "SOLVPB";
   rgstrStatement_ShortForm[27]   = "CLCSRF";
   rgstrStatement_ShortForm[28]   = "PHICON";
   //rgstrStatement_ShortForm[29] = "RADPOL"; // ONLY FOR OBJECTS. REMOVED FROM THE LIST.
   rgstrStatement_ShortForm[29]   = "UNUSED";
   rgstrStatement_ShortForm[30]   = "RELPAR";
   rgstrStatement_ShortForm[31]   = "SALT2";
   rgstrStatement_ShortForm[32]   = "RADPR2";
   rgstrStatement_ShortForm[33]   = "VAL+1";
   rgstrStatement_ShortForm[34]   = "VAL-1";
   rgstrStatement_ShortForm[35]   = "VAL+2";
   rgstrStatement_ShortForm[36]   = "VAL-2";
   rgstrStatement_ShortForm[37]   = "RMSC";
   rgstrStatement_ShortForm[38]   = "MAXC";
   //rgstrStatement_ShortForm[39] = "NORMC"; // UNUSED. REMOVED FROM THE LIST.
   rgstrStatement_ShortForm[39]   = "UNUSED";
   rgstrStatement_ShortForm[40]   = "VDROPX";
   rgstrStatement_ShortForm[41]   = "VDROPY";
   rgstrStatement_ShortForm[42]   = "VDROPZ";
   rgstrStatement_ShortForm[43]   = "ATPODS";
   rgstrStatement_ShortForm[44]   = "TEMPER";
   //-------------------Lin Li: Gaussian & MEMPOT -----------------//
   rgstrStatement_ShortForm[45]   = "CUTOFF";
   rgstrStatement_ShortForm[46]   = "SIGMA";
   rgstrStatement_ShortForm[47]   = "INHOMO";
   rgstrStatement_ShortForm[48]   = "SRFCUT";
   rgstrStatement_ShortForm[49]   = "GAUSSIAN";
   rgstrStatement_ShortForm[50]   = "RADIPZ";



   rgstrStatement_2lAbbre[0]      = "UNUSED";
   rgstrStatement_2lAbbre[1]      = "GS";
   rgstrStatement_2lAbbre[2]      = "SC";
   rgstrStatement_2lAbbre[3]      = "PF";
   rgstrStatement_2lAbbre[4]      = "ID";
   rgstrStatement_2lAbbre[5]      = "ED";
   rgstrStatement_2lAbbre[6]      = "PR";
   rgstrStatement_2lAbbre[7]      = "IR";
   rgstrStatement_2lAbbre[8]      = "IS";
   rgstrStatement_2lAbbre[9]      = "BC";
   rgstrStatement_2lAbbre[10]     = "LI";
   rgstrStatement_2lAbbre[11]     = "NI";
   //rgstrStatement_2lAbbre[12]   = "MD"; // OBSOLELE. REMOVED FROM THE LIST.
   rgstrStatement_2lAbbre[12]     = "UNUSED";
   rgstrStatement_2lAbbre[13]     = "FC";
   rgstrStatement_2lAbbre[14]     = "LP";
   rgstrStatement_2lAbbre[15]     = "LG";
   rgstrStatement_2lAbbre[16]     = "CI";
   rgstrStatement_2lAbbre[17]     = "CF";
   rgstrStatement_2lAbbre[18]     = "PX";
   rgstrStatement_2lAbbre[19]     = "PY";
   rgstrStatement_2lAbbre[20]     = "PZ";
   rgstrStatement_2lAbbre[21]     = "AC";
   rgstrStatement_2lAbbre[22]     = "XU";
   rgstrStatement_2lAbbre[23]     = "GC";
   rgstrStatement_2lAbbre[24]     = "RF";
   rgstrStatement_2lAbbre[25]     = "CB";
   rgstrStatement_2lAbbre[26]     = "SP";
   rgstrStatement_2lAbbre[27]     = "CS";
   rgstrStatement_2lAbbre[28]     = "PC"; // new for PHICON
   rgstrStatement_2lAbbre[29]     = "RL";
   rgstrStatement_2lAbbre[30]     = "RR";
   rgstrStatement_2lAbbre[31]     = "S2";
   rgstrStatement_2lAbbre[32]     = "R2";
   rgstrStatement_2lAbbre[33]     = "+1";
   rgstrStatement_2lAbbre[34]     = "-1";
   rgstrStatement_2lAbbre[35]     = "+2";
   rgstrStatement_2lAbbre[36]     = "-2";
   rgstrStatement_2lAbbre[37]     = "MC";
   rgstrStatement_2lAbbre[38]     = "XC";
   rgstrStatement_2lAbbre[39]     = "NC";
   rgstrStatement_2lAbbre[40]     = "VX";
   rgstrStatement_2lAbbre[41]     = "VY";
   rgstrStatement_2lAbbre[42]     = "VZ";
   rgstrStatement_2lAbbre[43]     = "AD";
   rgstrStatement_2lAbbre[44]     = "TE";
  //-------------------Lin Li: Gaussian & MEMPOT -----------------//
   rgstrStatement_2lAbbre[45]     = "CT";
   rgstrStatement_2lAbbre[46]     = "SG";
   rgstrStatement_2lAbbre[47]     = "IH";
   rgstrStatement_2lAbbre[48]     = "SF";
   rgstrStatement_2lAbbre[49]     = "GN";
   rgstrStatement_2lAbbre[50]     = "RZ";

   //------------------------------ functions ------------------------//
   rgstrFunction_FullForm[0]      = "UNUSED";
   rgstrFunction_FullForm[1]      = "CENTER";
   rgstrFunction_FullForm[2]      = "ACENTER";
   rgstrFunction_FullForm[3]      = "READ";
   rgstrFunction_FullForm[4]      = "WRITE";
   rgstrFunction_FullForm[5]      = "ENERGY";
   rgstrFunction_FullForm[6]      = "SITE";
   rgstrFunction_FullForm[7]      = "BUFFZ";
   rgstrFunction_FullForm[8]      = "QPREF";
   rgstrFunction_FullForm[9]      = "INSOBJ";
   rgstrFunction_FullForm[10]     = "SURFACE";
   //rgstrFunction_FullForm[11]   = "SOLVER";

   rgstrFunction_ShortForm[0]     = "UNUSED";
   rgstrFunction_ShortForm[1]     = "CENT";
   rgstrFunction_ShortForm[2]     = "ACENT";
   rgstrFunction_ShortForm[3]     = "IN";
   rgstrFunction_ShortForm[4]     = "OUT";

   //------------------------------ LOCAL ----------------------------//
   strCentFile          = "fort.27";    // centnam (not to be mapped)
   fMaxDimension        = 0.0;          // rmaxdim (not to be mapped)

}
