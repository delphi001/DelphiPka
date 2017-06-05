/*
 * delphi_data_setMap.cpp
 *
 *  Created on: Feb 09, 2014
 *      Author: chuan
 */

#include "delphi_data.h"

void CDelphiData::setMap()
{
   //--------------------- uniform parameters ---------------------//
   myData["biomodel"]    = pddm->strBioModel;
   myData["solver"]      = pddm->strNumSolver;
   //--------------------- set by Statements ----------------------//
   myData["iautocon"]    = pddm->bAutoConverge;
   myData["ibctyp"]      = pddm->iBndyType;
   myData["perfil"]      = pddm->fPercentageFill;
   myData["icheb"]       = pddm->bFixedRelaxParam;
   myData["isrf"]        = pddm->bOutGraspSurf;
   myData["icon2"]       = pddm->iConvergeFract;
   myData["icon1"]       = pddm->iIterateInterval;
   myData["iexun"]       = pddm->bExitUniformDielect;
   myData["repsout"]     = pddm->fExDielec;
   myData["isph"]        = pddm->bCrgInterplateType;
   myData["gten"]        = pddm->fGridConverge;
   myData["igrid"]       = pddm->iGrid;
   myData["repsin"]      = pddm->fInDielec;
   myData["conc"]        = pddm->vctfSalt; // std::vector
   myData["exrad"]       = pddm->fIonRadius;
   myData["nlit"]        = pddm->iLinIterateNum;
   myData["igraph"]      = pddm->bLogGraph;
   myData["ipoten"]      = pddm->bLogPotential;
   myData["res2"]        = pddm->fMaxc;
   myData["nnit"]        = pddm->iNonIterateNum;
   myData["iper"]        = pddm->vctbPeriodicBndy; // std::vector
   myData["iconc"]       = pddm->bOutCrgDensity;
   myData["radprb"]      = pddm->vctfProbeRadius; // std::vector
   myData["uspec"]       = pddm->fSpectralRadius;
   myData["relpar"]      = pddm->fRelaxParam;
   myData["res1"]        = pddm->fRmsc;
   myData["scale"]       = pddm->fScale;
   myData["isolv"]       = pddm->bSolvePB;
   myData["ival"]        = pddm->vctiValence1; // std::vector
   myData["ival2"]       = pddm->vctiValence2; // std::vector
   myData["atompotdist"] = pddm->fPotentialUpperBond;
   myData["temperature"] = pddm->fTemper;
   myData["vdrop"]       = pddm->gfPotentialDrop; // SGrid<delphi_real>
   myData["iuspec"]      = pddm->bSpectralRadius;
   myData["imanual"]     = pddm->bManualRelaxParam;
   myData["phiintype"]   = pddm->iPhiInType; //for focusing

   //---------------------Gaussian & MEMPOT ----------------------//
   //Lin Li: for Gaussian & MEMPOT options
   myData["cutoff"]      = pddm->fCutoff;
   myData["sigma"]       = pddm->fSigma;
   myData["inhomo"]      = pddm->iInhomo;
   myData["srfcut"]      = pddm->fSrfcut;
   myData["gaussian"]    = pddm->iGaussian;
   myData["gepsmp"]      = pddm->vctgfGepsMap;
   myData["gepsmp2"]     = pddm->vctgfGepsMap2;
   myData["ergsgaussian"]= pddm->fErgsgaussian;
   myData["radipz"]      = pddm->fRadipz;


   //-------------------------- io file names ------------------------//
   //myData["prmnam"]    = pddm->strParamFile[0]; // (not to be mapped)
   myData["siznam"]      = pddm->strSizeFile;
   myData["crgnam"]      = pddm->strCrgFile;
   myData["pdbnam"]      = pddm->strPdbFile;
   myData["phinam"]      = pddm->strPhiFile;
   myData["frcinam"]     = pddm->strFrciFile;
   myData["frcnam"]      = pddm->strFrcFile;
   myData["epsnam"]      = pddm->strEpsFile;
   myData["phiinam"]     = pddm->strPhiiFile;
   myData["mpdbnam"]     = pddm->strModifiedPdbFile;
   myData["updbnam"]     = pddm->strUnformatPdbFile;
   myData["ufrcnam"]     = pddm->strUnformatFrcFile;
   myData["srfnam"]      = pddm->strGraspFile;
   myData["nrgnam"]      = pddm->strEnergyFile;
   myData["scrgnam"]     = pddm->strScrgFile;
   //myData["centnam"]   = pddm->strCentFile; // renamed to be fort.27 (not to be mapped)
   //----------------------- set by functions ------------------------//
   /*
    * set by CENTER or CENT function:
    */
   myData["offset"]      = pddm->gfOffCenter; // SGrid<delphi_real>

   /*
    * set by ACENTER or ACENT function
    */
   myData["acent"]       = pddm->gfAcent; // SGrid<delphi_real>
   myData["iacent"]      = pddm->bIsAcent;

   /*
    * set by READ or IN function
    */
   myData["pdbfrm"]      = pddm->iPdbFormatIn;
   myData["ipdbrd"]      = pddm->bPdbUnformatIn;

   /*
    * set by WRITE or OUT function
    */
   myData["phiwrt"]      = pddm->bPhimapOut;
   myData["phifrm"]      = pddm->iPhiFormatOut;
   myData["ibios"]       = pddm->bBiosystemOut;
   myData["ibem"]        = pddm->bBemSrfOut;
   myData["isite"]       = pddm->bSiteOut;
   myData["frcfrm"]      = pddm->iFrcFormatOut;
   myData["epswrt"]      = pddm->bEpsOut;
   myData["iatout"]      = pddm->bModPdbOut;
   myData["mpdbfrm"]     = pddm->iModPdbFormatOut;
   myData["ipdbwrt"]     = pddm->bUnformatPdbOut;
   myData["ifrcwrt"]     = pddm->bUnformatFrcOut;
   myData["inrgwrt"]     = pddm->bEngOut;
   myData["iwgcrg"]      = pddm->bGridCrgOut;
   myData["iacs"]        = pddm->bHsurf2DatOut;
   myData["idbwrt"]      = pddm->bDbOut;
   myData["isen"]        = pddm->bSurfEngOut;
   myData["isch"]        = pddm->bSurfCrgOut;
   myData["scrgfrm"]     = pddm->iSurfCrgFormatOut;

   /*
    * set by ENERGY function
    */
   myData["logg"]        = pddm->bGridEng;
   myData["logs"]        = pddm->bSolvEng;
   myData["logas"]       = pddm->bAnalySurfEng;
   myData["loga"]        = pddm->bAnalyEng;
   myData["logions"]     = pddm->bIonsEng;
   myData["logc"]        = pddm->bCoulombEng;

   /*
    * set by SITE function: all MUST be initialized to to false
    */
   myData["isita"]       = pddm->bAtomInSite;
   myData["isitq"]       = pddm->bCrgInSite;
   myData["isitp"]       = pddm->bGridPotentialInSite;
   myData["isitap"]      = pddm->bAtomPotentialInSite;
   myData["isitdeb"]     = pddm->bDebyeFractionInSite;
   myData["isitf"]       = pddm->bFieldInSite;
   myData["isitr"]       = pddm->bReactPotentialInSite;
   myData["isitc"]       = pddm->bCoulombPotentialInSite;
   myData["isitx"]       = pddm->bAtomCoordInSite;
   myData["isiti"]       = pddm->bSaltInSite;
   myData["isitt"]       = pddm->bTotalPotentialInSite;
   myData["isitrf"]      = pddm->bReactForceInSite;
   myData["isitcf"]      = pddm->bCoulombForceInSite;
   myData["isitmd"]      = pddm->bMDInSite;
   myData["isitsf"]      = pddm->bSurfCrgInSite;
   myData["isittf"]      = pddm->bTotalForceInSite;
   myData["isitpot"]     = pddm->bPotentialInSite;
   myData["irea"]        = pddm->bReactFieldInFRC;
   myData["iself"]       = pddm->bPDB2FRCInSite;

   /*
    * set by BUFFZ function
    */
   myData["buffz"]       = pddm->eiBuffz;
   myData["ibufz"]       = pddm->bIsBuffz;

   /*
    * set by SURFACE function
    */
   myData["isurftype"]   = pddm->iTypeSurf;
   //----------------------- set by DelPhi ------------------------//
   myData["deblen"]      = pddm->fDebyeLength;
   myData["epsout"]      = pddm->fEpsOut;
   myData["cran"]        = pddm->gfCoordinateRange;
   myData["pmid"]        = pddm->gfGeometricCenter;
   myData["oldmid"]      = pddm->gfBoxCenter;
   myData["rionst"]      = pddm->fIonStrength;
   myData["chi1"]        = pddm->fTaylorCoeff1;
   myData["chi2"]        = pddm->fTaylorCoeff2;
   myData["chi3"]        = pddm->fTaylorCoeff3;
   myData["chi4"]        = pddm->fTaylorCoeff4;
   myData["chi5"]        = pddm->fTaylorCoeff5;
   myData["lognl"]       = pddm->bNonlinearEng;
   myData["epkt"]        = pddm->fEPKT;
   myData["epsin"]       = pddm->fEpsIn;
   myData["ifrcrd"]      = pddm->bFrcUnformatIn;
   myData["idirectalg"]  = pddm->iDirectEpsMap;
   myData["numbmol"]     = pddm->iMoleculeNum;
   myData["rdmx"]        = pddm->fMaxRadius;
   myData["uniformdiel"] = pddm->bUniformDielec;
   myData["limobject"]   = pddm->vctefExtrema;     // std::vector< SExtrema<delphi_real> >
   myData["xn1"]         = pddm->vctgfAtomCoordA;  // std::vector< SGrid<delphi_real> >
   myData["xn2"]         = pddm->vctgfAtomCoordG;  // std::vector< SGrid<delphi_real> >
   //----------------------- set by IO class ------------------------//
   myData["resnummax"]   = pddm->iResidueNum;
   myData["nmedia"]      = pddm->iMediaNum;
   myData["medeps"]      = pddm->vctfMediaEps;     // std::vector<delphi_real>
   myData["nobject"]     = pddm->iObjectNum;
   myData["dataobject"]  = pddm->vctstrObject;     // std::vector<string>
   myData["natom"]       = pddm->iAtomNum;
   myData["delphipdb"]   = pddm->vctapAtomPdb;     // std::vector<CAtomPdb>
   myData["iatmmed"]     = pddm->vctiAtomMediaNum; // std::vector<delphi_integer>
   myData["ionlymol"]    = pddm->bOnlyMolecule;
   //myData["ndistr"]    = pddm->iCrgDistrNum;
   //myData["datadistr"] = pddm->prgstrCrgDistr;
   //------------------- set by Surface class ---------------------//
   myData["nqass"]       = pddm->iCrgGridNum;
   myData["qnet"]        = pddm->fNetCrg;
   myData["qmin"]        = pddm->fMinusCrg;
   myData["qplus"]       = pddm->fPlusCrg;
   myData["cqplus"]      = pddm->gfPlusCrgCenter;
   myData["cqmin"]       = pddm->gfMinusCrgCenter;
   myData["cmin"]        = pddm->gfMinCoordinate;
   myData["cmax"]        = pddm->gfMaxCoordinate;
   myData["ibnum"]       = pddm->iBndyGridNum;
   myData["iepsmp"]      = pddm->vctgiEpsMap;      // std::vector< SGrid<delphi_integer> >
   myData["idebmap"]     = pddm->vctbDielecMap;    // std::vector<bool>
   myData["ibgrd"]       = pddm->vctgiBndyGrid;    // std::vector< SGrid<delphi_integer> >
   myData["nqgrd"]       = pddm->iCrg2GridNum;
   myData["chrgv2"]      = pddm->vctgvfCrg2Grid;   // std::vector< SGridValue<delphi_real> >
   myData["nqgrdtonqass"]= pddm->vctiCrg2GridMap;  // std::vector<delphi_integer>
   myData["atmcrg"]      = pddm->vctgvfAtomCrg;    // std::vector< SGridValue<delphi_real> >
   myData["chgpos"]      = pddm->vctgfCrgPoseA;    // std::vector< SGrid<delphi_real> >
   myData["scspos"]      = pddm->vctgfSurfCrgA;    // std::vector< SGrid<delphi_real> >
   myData["crgatn"]      = pddm->vctiCrgAt;        // std::vector<delphi_integer>
   myData["atsurf"]      = pddm->vctiAtSurf;       // std::vector<delphi_integer>
   myData["atndx"]       = pddm->vctiAtNdx;        // std::vector<delphi_integer>
   myData["scsnor"]      = pddm->vctgfSurfCrgE;    // std::vector< SGrid<delphi_real> >
   myData["atmeps"]      = pddm->vctfAtomEps;      // std::vector<delphi_real>
   //------------------- set by Solver class ---------------------//
   myData["icount2b"]    = pddm->iDielecBndySum;
   myData["icount1b"]    = pddm->iCrgedGridSum;
   myData["gchrg"]       = pddm->vctfGridCrg;      // std::vector<delphi_real>
   myData["gchrgp"]      = pddm->vctgiGridCrgPose; // std::vector< SGrid<delphi_integer> >
   myData["ibc"]         = pddm->iCrgBdyGrid;
   myData["cgbp"]        = pddm->vctdgvCrgBndyGrid;// std::vector<SDoubleGridValue>
   myData["phimap"]      = pddm->vctfPhiMap;       // std::vector<delphi_real>
   myData["phimap_pre"]  = pddm->vctfPhiMap_Pre;       // previous phimap for focusing

   //------------------- set by Energy class ---------------------//
   myData["schrg"]       = pddm->vctfSurfCrgE;     // std::vector<delphi_real>
   myData["ergg"]        = pddm->fEngGrid;
   myData["ergc"]        = pddm->fEngCoul;
   myData["ergs"]        = pddm->fEngCorrect;
   myData["ergr"]        = pddm->fEngReact;
   myData["ergions"]     = pddm->fEngIons;
   //------------------- New var by Interface  ---------------------//
#ifdef PRIME
   myData["vecfrcin"]    = pddm->strCommFRCIn;    // update: Dec 19, 2014 by Lin Wang
#endif
   myData["isitcomm"]    = pddm->bCommFRCIn;      // update: Dec 19, 2014 by Lin Wang

}



