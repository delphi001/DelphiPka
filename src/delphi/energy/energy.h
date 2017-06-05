/*
 * energy.h
 *
 *  Created on: Feb 23, 2014
 *      Author: Lin Wang, lwang3@g.clemson.edu
 *
 *	This file declares the public and private variables in CDelphiEnergy Class.
 *
 */


#ifndef ENERGY_H_
#define ENERGY_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include "../interface/interface_abstractmodule.h"
#include "../misc/misc_timer.h"
#include "../misc/misc_grid.h"
#include "energy_exceptions.h"

using namespace std;

class CDelphiEnergy:virtual public IAbstractModule
{
   private:

      shared_ptr<CTimer> pTimer;

      const int& iBndyType;
      const string& strEnergyFile;
      const string& strScrgFile;
      const bool& bEngOut;
      const bool& bSurfEngOut;
      const bool& bSurfCrgOut;
      const int& iSurfCrgFormatOut;
      const bool& bGridEng;
      const bool& bSolvEng;
      const bool& bAnalySurfEng;
      const bool& bAnalyEng;
      const bool& bIonsEng;
      const bool& bCoulombEng;
      const bool& bPotentiallnSite;
      const bool& bReactFieldlnFRC;
      const bool& bBuffz;
      const int& iCrgGridNum;
      const delphi_real& fEpsOut;
      const bool& bNonlinearEng;
      const int& iMediaNum;
      const int& iGaussian; //for Gaussian
      const int& inhomo;

      const delphi_real& prgrMediaEps;
      const int& prgiAtomMediaNum;
      const int& iTotalBdyGridNum;
      const int& iCrgBdyGrid;
      const vector< SGrid<int> >& prgigBndyGrid;
      const vector< SGridValue<delphi_real> >& prggvAtomicCrg;
      const int& iGrid;
      const int& iCrgedGridB;
      const delphi_real& fEpsin;
   //	const bool& DEVELOPER;    // Developer is not from datacontainer now. It changed to the flag.
      const vector<delphi_real>& prgfPhimap;
      const delphi_real& fScale;
      const SGrid<delphi_real>& fgPotentialDrop;
      const vector< SGrid<delphi_real> >& prgfgCrgPoseA;
      const vector< SGrid<int> >& prgigGridCrgPose;
      const vector<delphi_real>& prgfGridCrg;
      const delphi_real& fIonStrength;
      const SGrid<delphi_real>& fgBoxCenter;
      const vector<bool>& prgbDielecMap;
      const vector<delphi_real>& prgfAtomEps;
      const SExtrema<int>& ieBuffz;

      const int& iObjectNum;
      const int& iAtomNum;
      const int& iDielecBndyOdd;
      const delphi_real& fEPKT;
      const vector<CAtomPdb>& prgapAtomPdb;
      vector< SDoubleGridValue >& prgdgvCrgBndyGrid;
      const vector< SGrid<delphi_real> >& prgfgSurfCrgA;
      const vector<delphi_integer>& prgiCrgAt;
      const vector<delphi_integer>& atsurf;


      const delphi_real& fTaylorCoeff1;
      const delphi_real& fTaylorCoeff2;
      const delphi_real& fTaylorCoeff3;
      const delphi_real& fTaylorCoeff4;
      const delphi_real& fTaylorCoeff5;

      bool debug_energy;
      vector < SGridValue<delphi_real> > sout;
      SGrid<int> lim_min, lim_max;

      vector <delphi_real>& schrg;
      delphi_real& ergg;
      delphi_real& ergc;
      delphi_real& ergs;
      delphi_real& ergr;
      delphi_real& ergions;
      delphi_real& ergsgaussian;

      //***************//

      void energy_clb(delphi_real& ergc);
      void energy_clbmedia(delphi_real& ergc);
      void energy_clbnonl(delphi_real& ergc, delphi_real& ergest, int& igridout);
      void energy_clbtotal(delphi_real& ergest, delphi_real& ergc);
      void energy_react(delphi_real& ergs, delphi_real& ergas, int& iisitpot);
      void energy_nonl(delphi_real& ergnl, int& igridout);


   public:
      CDelphiEnergy(shared_ptr<IDataContainer> pdc,shared_ptr<CTimer> pt):
         IAbstractModule(pdc),
         pTimer(pt),

         iBndyType(pdc->getKey_constRef<int>("ibctyp")),
         strEnergyFile(pdc->getKey_constRef<string>("nrgnam")),
         strScrgFile(pdc->getKey_constRef<string>("scrgnam")),
         bEngOut(pdc->getKey_constRef<bool>("inrgwrt")),
         bSurfEngOut(pdc->getKey_constRef<bool>("isen")),
         bSurfCrgOut(pdc->getKey_constRef<bool>("isch")),
         iSurfCrgFormatOut(pdc->getKey_constRef<int>("scrgfrm")),
         bGridEng(pdc->getKey_constRef<bool>("logg")),
         bSolvEng(pdc->getKey_constRef<bool>("logs")),
         bAnalySurfEng(pdc->getKey_constRef<bool>("logas")),
         bAnalyEng(pdc->getKey_constRef<bool>("loga")),
         bIonsEng(pdc->getKey_constRef<bool>("logions")),
         bCoulombEng(pdc->getKey_constRef<bool>("logc")),
         bPotentiallnSite(pdc->getKey_constRef<bool>("isitpot")),
         bReactFieldlnFRC(pdc->getKey_constRef<bool>("irea")),
         bBuffz(pdc->getKey_constRef<bool>("ibufz")),
         iCrgGridNum(pdc->getKey_constRef<delphi_integer>("nqass")),
         fEpsOut(pdc->getKey_constRef<delphi_real>("epsout")),
         bNonlinearEng(pdc->getKey_constRef<bool>("lognl")),
         iMediaNum(pdc->getKey_constRef<delphi_integer>("nmedia")),
         iGaussian(pdc->getKey_constRef<int>("gaussian")),
         inhomo(pdc->getKey_constRef<int>("inhomo")),
         prgrMediaEps(pdc->getKey_constRef<delphi_real>("medeps")),
         prgiAtomMediaNum(pdc->getKey_constRef<delphi_integer>("iatmmed")),
         iTotalBdyGridNum(pdc->getKey_constRef<delphi_integer>("ibnum")),
         iCrgBdyGrid(pdc->getKey_constRef<delphi_integer>("ibc")),
         prgigBndyGrid(pdc->getKey_constRef< vector< SGrid<int> > >("ibgrd")),
         iGrid(pdc->getKey_constRef<delphi_integer>("igrid")),
         iCrgedGridB(pdc->getKey_constRef<delphi_integer>("icount1b")),
         fEpsin(pdc->getKey_constRef<delphi_real>("epsin")),
   //		DEVELOPER(pdc->getKey_constRef<bool>("ideveloper")),
         prgfPhimap(pdc->getKey_constRef< vector<delphi_real> >("phimap")),  // phimap read function to convert to 3d array. //
         fScale(pdc->getKey_constRef<delphi_real>("scale")),
         fgPotentialDrop(pdc->getKey_constRef< SGrid<delphi_real> >("vdrop")),
         prgfgCrgPoseA(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("chgpos")),
         prgfGridCrg(pdc->getKey_Ref< vector<delphi_real> >("gchrg")),
         prgigGridCrgPose(pdc->getKey_constRef< vector< SGrid<int> > >("gchrgp")),
         prggvAtomicCrg(pdc->getKey_constRef< vector< SGridValue<delphi_real> > >("atmcrg")),
         fIonStrength(pdc->getKey_constRef<delphi_real>("rionst")),
         fgBoxCenter(pdc->getKey_constRef< SGrid<delphi_real> >("oldmid")),
         prgbDielecMap(pdc->getKey_constRef< vector<bool> >("idebmap")),
         prgfAtomEps(pdc->getKey_constRef< vector<delphi_real> >("atmeps")),
         ieBuffz(pdc->getKey_constRef< SExtrema<int> >("buffz")),

         iObjectNum(pdc->getKey_constRef<int>("nobject")),
         iAtomNum(pdc->getKey_constRef<int>("natom")),
         iDielecBndyOdd(pdc->getKey_constRef<int>("icount2b")),
         fEPKT(pdc->getKey_constRef<delphi_real>("epkt")),
         prgdgvCrgBndyGrid(pdc->getKey_Ref< vector<SDoubleGridValue> >("cgbp")),
         prgapAtomPdb(pdc->getKey_Ref< vector<CAtomPdb> >("delphipdb")),
         prgfgSurfCrgA(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("scspos")),
         prgiCrgAt(pdc->getKey_constRef< vector<delphi_integer> >("crgatn")),
         atsurf(pdc->getKey_constRef< vector<delphi_integer> >("atsurf")),
         fTaylorCoeff1(pdc->getKey_constRef<delphi_real>("chi1")),
         fTaylorCoeff2(pdc->getKey_constRef<delphi_real>("chi2")),
         fTaylorCoeff3(pdc->getKey_constRef<delphi_real>("chi3")),
         fTaylorCoeff4(pdc->getKey_constRef<delphi_real>("chi4")),
         fTaylorCoeff5(pdc->getKey_constRef<delphi_real>("chi5")),
         schrg(pdc->getKey_Ref< vector<delphi_real> >("schrg")),
         ergg(pdc->getKey_Ref<delphi_real>("ergg")),
         ergc(pdc->getKey_Ref<delphi_real>("ergc")),
         ergs(pdc->getKey_Ref<delphi_real>("ergs")),
         ergr(pdc->getKey_Ref<delphi_real>("ergr")),
         ergions(pdc->getKey_Ref<delphi_real>("ergions")),
         ergsgaussian(pdc->getKey_Ref<delphi_real>("ergsgaussian"))
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                CDelphiEnergy is constructed                  *\n";
         cout << "****************************************************************\n";
#endif
      };

       ~CDelphiEnergy()
       {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                 CDelphiEnergy is destroyed                   *\n";
         cout << "****************************************************************\n";
#endif
       };

       virtual void validateInput(){};

       virtual void run();


};



#endif /* ENERGY_H_ */
