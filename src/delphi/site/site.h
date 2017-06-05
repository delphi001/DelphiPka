/*
 * site.h
 *
 *  Created on: Feb 17, 2014
 *      Author: chuan
 */

#ifndef SITE_H_
#define SITE_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <memory>
#include <cstdio>

#include "../interface/interface_datacontainer.h"
#include "../misc/misc_timer.h"
#include "../misc/misc_interpl.h"
#include "../delphi/delphi_constants.h"
#include "../io/io.h"
#include "site_exceptions.h"

using namespace std;

class CSite:public CIO
{
   private:
      /*********************************************************************************************
       *                                                                                           *
       *              references to the variables obtained from the data container                 *
       *                                                                                           *
       ********************************************************************************************/
      shared_ptr<CTimer> pTimer;
      //----- set by statements
      const delphi_integer& iGrid;                             // igrid
      const delphi_real&    fPercentageFill;                   // perfil
      const delphi_real&    fExDielec;                         // repsout
      const delphi_real&    fIonRadius;                        // exrad
      const vector<delphi_real>& rgfProbeRadius;               // radprb
      const int&     iLinIterateNum;                    // nlit
      const int&     iNonIterateNum;                    // nnit
      const int&     iBndyType;                         // ibctyp
      const delphi_real&    fPotentialUpperBond;               // atompotdist
      const bool&    bOutCrgDensity;                    // iconc
      const float&        radipz ;                      //radipz
      //----- io file names
      const string&  strFrciFile;                       // frcinam
      const string&  strFrcFile;                        // frcnam
      const string&  strPhiFile;                        // phinam
      //----- set by functions
      const bool&    bAtomInSite;                       // isita
      const bool&    bSaltInSite;                       // isiti
      const bool&    bMDInSite;                         // isitmd
      const bool&    bPotentialInSite;                  // isitpot
      const int&     iPhiFormatOut;                     // phifrm
      const bool&    bBiosystemOut;                     // ibios
      //----- set by DelPhi
      const delphi_integer& iNatom;
      const delphi_real&    fIonStrength;                      // rionst
      const SGrid<delphi_real>& fgBoxCenter;                   // oldmid
      const vector< SGrid<delphi_real> >& prgfgAtomCoordA;     // xn1(natom)
      const delphi_real&    fTaylorCoeff1;                     // chi1
      const delphi_real&    fTaylorCoeff2;                     // chi2
      const delphi_real&    fTaylorCoeff3;                     // chi3
      const delphi_real&    fTaylorCoeff4;                     // chi4
      const delphi_real&    fTaylorCoeff5;                     // chi5
      const vector < SGrid<delphi_real> >& xn2;                //xn2
      //----- set by Surface class
      const delphi_integer& iBndyGridNum;                      // ibnum
      const delphi_integer& iCrgGridNum;                       // nqass
      const vector< SGridValue<delphi_real> >& prggvAtomicCrg; // atmcrg(nqass)
      const vector< SGrid<delphi_real> >& prgfgCrgPoseA;       // chgpos(ibnum)
      const vector< SGrid<delphi_real> >& prgfgSurfCrgA;       // scspos(ibnum)
      const vector<delphi_integer>& prgiCrgAt;                 // crgatn(nqass)
      const vector<delphi_integer>& prgiAtSurf;                // atsurf(ibnum)
      const vector< SGrid<delphi_real> >& prgfgSurfCrgE;       // scsnor(ibnum)
      const vector<bool>& prgbDielecMap;                // idebmap(igrid,igrid,igrid)
      const vector< SGrid<delphi_integer> >& prgigBndyGrid;    // ibgrd(ibnum)
      const vector<delphi_real>& prgfAtomEps;                  // atmeps(nqass)
      const vector<delphi_integer>& prgiAtNdx;                 // atndx(ibnum)
      //----- set by Solver class
      const delphi_integer& iDielecBndySum;                    // icount2b
      //----- set by Energy class
      const vector<delphi_real>& prgfSurfCrgE;                 // schrg(ibnum)
      //++++++++++++++++ reference to read-and-write variables from data container +++++++++++++++//
      bool&          bAtomCoordInSite;                  // isitx
      bool&          bCrgInSite;                        // isitq
      bool&          bFieldInSite;                      // isitf
      bool&          bGridPotentialInSite;              // isitp
      bool&          bReactPotentialInSite;             // isitr
      bool&          bCoulombPotentialInSite;           // isitc
      bool&          bAtomPotentialInSite;              // isitap
      bool&          bDebyeFractionInSite;              // isitdeb
      bool&          bSurfCrgInSite;                    // isitsf
      bool&          bTotalForceInSite;                 // isittf
      bool&          bReactForceInSite;                 // isitrf
      bool&          bTotalPotentialInSite;             // isitt
      bool&          bCoulombForceInSite;               // isitcf
      bool&          bPDB2FRCInSite;                    // iself
      int&           iFrcFormatOut;                     // frcfrm
      delphi_real&          fScale;                            // scale
      vector<delphi_real>&  prgfPhiMap;                        // phimap
#ifdef PRIME
    
      const vector<string>& strCommFRCIn;

#endif
      const bool&    bCommFRCIn;
      /*********************************************************************************************
       *                                                                                           *
       *                            variables defined in this class                                *
       *                                                                                           *
       ********************************************************************************************/
      delphi_real *** phimap;

      vector< SGrid<delphi_real> > rforceeps1();

      vector< SGrid<delphi_real> > rforce();

      //delphi_real debinterpl(const SGrid<delphi_real> & gPoint);

      delphi_real tops(const SGrid<delphi_real>& xxo,const SGrid<delphi_real>& xxu,const delphi_real& crg,const delphi_real& eps,const int& flag);

      void phicon();

      void writePotential_insight(vector<delphi_real>& phimapIn);

      void writePotential_grasp(vector<delphi_real>& phimapIn);

      void writePotential_ccp4(vector<delphi_real>& phimapIn);

      void writePotential_fromPrevious(vector<delphi_real>& phimapIn);

      void writePotential_cube();

      void writePotential_delphi();

      void writePhiMap(const int& formatflag,vector<delphi_real>& phimap4,ofstream& ofFileStream);

      void expand(const int& mgrid, vector<delphi_real>& phimapIn);



   public:

#ifdef MCCE
      vector<delphi_real> mcce_phiv; // to save phiv produced in site_writeSite()
#endif

#ifdef PRIME
      vector<delphi_real> prime_grdphiv;
      vector<string>      prime_atomdes;
      vector<delphi_real> prime_crhgv;
#endif

      CSite(shared_ptr<IDataContainer> pdc,shared_ptr<CTimer> pt):
         CIO(pdc->getKey_Val<delphi_real>("repsin"),pdc->getKey_Val<delphi_real>("epkt")),
         pTimer(pt),
         //----- set by statements
         iGrid(pdc->getKey_constRef<delphi_integer>("igrid")), // modified in setFocusBndy
         fPercentageFill(pdc->getKey_constRef<delphi_real>("perfil")),
         fExDielec(pdc->getKey_constRef<delphi_real>("repsout")),
         fIonRadius(pdc->getKey_constRef<delphi_real>("exrad")),
         rgfProbeRadius(pdc->getKey_constRef< vector<delphi_real> >("radprb")),
         iLinIterateNum(pdc->getKey_constRef<int>("nlit")),
         iNonIterateNum(pdc->getKey_constRef<int>("nnit")),
         iBndyType(pdc->getKey_constRef<int>("ibctyp")),
         fPotentialUpperBond(pdc->getKey_constRef<delphi_real>("atompotdist")),
         bOutCrgDensity(pdc->getKey_constRef<bool>("iconc")),
         radipz (pdc->getKey_constRef<float>("radipz")),
         //----- io file names
         strFrciFile(pdc->getKey_constRef<string>("frcinam")),
         strFrcFile(pdc->getKey_constRef<string>("frcnam")),
         strPhiFile(pdc->getKey_constRef<string>("phinam")),
         //----- set by functions
         bAtomInSite(pdc->getKey_constRef<bool>("isita")),
         bSaltInSite(pdc->getKey_constRef<bool>("isiti")),
         bMDInSite(pdc->getKey_constRef<bool>("isitmd")),
         bPotentialInSite(pdc->getKey_constRef<bool>("isitpot")),
         iPhiFormatOut(pdc->getKey_constRef<int>("phifrm")),
         bBiosystemOut(pdc->getKey_constRef<bool>("ibios")),
         //----- set by DelPhi
         iNatom (pdc->getKey_constRef<delphi_integer>("natom")),
         fIonStrength(pdc->getKey_constRef<delphi_real>("rionst")),
         fgBoxCenter(pdc->getKey_constRef< SGrid<delphi_real> >("oldmid")),
         prgfgAtomCoordA(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("xn1")),
         fTaylorCoeff1(pdc->getKey_constRef<delphi_real>("chi1")),
         fTaylorCoeff2(pdc->getKey_constRef<delphi_real>("chi2")),
         fTaylorCoeff3(pdc->getKey_constRef<delphi_real>("chi3")),
         fTaylorCoeff4(pdc->getKey_constRef<delphi_real>("chi4")),
         fTaylorCoeff5(pdc->getKey_constRef<delphi_real>("chi5")),
         xn2(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("xn2")),

         //----- set by Surface class
         iBndyGridNum(pdc->getKey_constRef<delphi_integer>("ibnum")),
         iCrgGridNum(pdc->getKey_constRef<delphi_integer>("nqass")),
         prggvAtomicCrg(pdc->getKey_constRef< vector< SGridValue<delphi_real> > >("atmcrg")),
         prgfgCrgPoseA(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("chgpos")),
         prgfgSurfCrgA(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("scspos")),
         prgiCrgAt(pdc->getKey_constRef< vector<delphi_integer> >("crgatn")),
         prgiAtSurf(pdc->getKey_constRef< vector<delphi_integer> >("atsurf")),
         prgfgSurfCrgE(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("scsnor")),
         prgbDielecMap(pdc->getKey_constRef< vector<bool> >("idebmap")),
         prgigBndyGrid(pdc->getKey_constRef< vector< SGrid<delphi_integer> > >("ibgrd")),
         prgfAtomEps(pdc->getKey_constRef< vector<delphi_real> >("atmeps")),
         prgiAtNdx(pdc->getKey_constRef< vector<delphi_integer> >("atndx")),

         //----- set by Solver class
         iDielecBndySum(pdc->getKey_constRef<delphi_integer>("icount2b")),
         //----- set by Energy class
         prgfSurfCrgE(pdc->getKey_constRef< vector<delphi_real> >("schrg")),
         //++++++++++++++ reference to read-and-write variables from data container ++++++++++++++//
         bAtomCoordInSite(pdc->getKey_Ref<bool>("isitx")),
         bCrgInSite(pdc->getKey_Ref<bool>("isitq")),
         bFieldInSite(pdc->getKey_Ref<bool>("isitf")),
         bGridPotentialInSite(pdc->getKey_Ref<bool>("isitp")),
         bReactPotentialInSite(pdc->getKey_Ref<bool>("isitr")),
         bCoulombPotentialInSite(pdc->getKey_Ref<bool>("isitc")),
         bAtomPotentialInSite(pdc->getKey_Ref<bool>("isitap")),
         bDebyeFractionInSite(pdc->getKey_Ref<bool>("isitdeb")),
         bSurfCrgInSite(pdc->getKey_Ref<bool>("isitsf")),
         bTotalForceInSite(pdc->getKey_Ref<bool>("isittf")),
         bReactForceInSite(pdc->getKey_Ref<bool>("isitrf")),
         bTotalPotentialInSite(pdc->getKey_Ref<bool>("isitt")),
         bCoulombForceInSite(pdc->getKey_Ref<bool>("isitcf")),
         bPDB2FRCInSite(pdc->getKey_Ref<bool>("iself")),
         iFrcFormatOut(pdc->getKey_Ref<int>("frcfrm")),
         fScale(pdc->getKey_Ref<delphi_real>("scale")),
#ifdef PRIME
         strCommFRCIn(pdc->getKey_constRef< vector<string> >("vecfrcin")),
#endif
         bCommFRCIn(pdc->getKey_constRef<bool>("isitcomm")),
         prgfPhiMap(pdc->getKey_Ref< vector<delphi_real> >("phimap"))

      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                    CSite is constructed                      *\n";
         cout << "****************************************************************\n";
#endif

         //----- variables inherited from CIO class
         iResidueNum  = pdc->getKey_Val<delphi_integer>("resnummax");
         vctfMediaEps = pdc->getKey_Val< vector<delphi_real> >("medeps");
         iMediaNum    = pdc->getKey_Val<delphi_integer>("nmedia");
         iAtomNum     = pdc->getKey_Val<delphi_integer>("natom");
         vctapAtomPdb = pdc->getKey_Val< vector<CAtomPdb> >("delphipdb");

         //----- local variables
         phimap = pdc->getKey_Ptr<delphi_real>("phimap",iGrid,iGrid,iGrid); // pointer to 3D phimap

         //cout << "fDielec = " << fDielec << endl;
      };

      /**
       * destructor
       */
      ~CSite()
      {
         /*
          * delete delphi_real *** phimap without deleting underneath prgfPhiMap in data container
          */
         for(int i = 0; i != iGrid; ++i)
         {
            //for(int j = 0; j != iGrid; ++j)
            //{
            //    delete[] phimap[i][j];
            //}
            delete[] phimap[i];
         }
         delete[] phimap;

#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                    CSite is destroyed                        *\n";
         cout << "****************************************************************\n";
#endif
      };

      void writeSite(const int& iisitsf);

      void writePhi();

      void wirtePAnalysis(); //Lin Li: radipz, panalysis
    
#ifdef PRIME
      void clearIO();
#endif
    
};



#endif /* SITE_H_ */
