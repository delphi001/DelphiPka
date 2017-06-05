/**
 * @file delphi_data.h
 * @brief class CDelphiData
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * This file declares the class CDelphiData, which inherits from the interface class IDataContainer.
 * Class CDelphiData implements a map of global variables passing to other classes.
 */

#ifndef CDELPHIDATA_H_
#define CDELPHIDATA_H_

#include <iostream>
#include <string>
#include <vector>
#include <memory> // std::shared_ptr

#include "../misc/misc_grid.h"
#include "../misc/misc_timer.h"
#include "../io/io_datatype.h"
#include "../interface/environment.h"
#include "../interface/interface_datacontainer.h"
#include "delphi_datamarshal.h"

using namespace std;

//-----------------------------------------------------------------------//

class CDelphiData:virtual public IDataContainer
{	
   private:
      /**
       * a shared pointer pointing to an object of class CDelphiDataMarshal
       */
      shared_ptr<CDelphiDataMarshal> pddm;

      /**
       * member function to print a flash of program information
       */
      void flash() const;

      /**
       * member function to implement the virtual function setMap() declared in IDataContainer
       */
      virtual void setMap();

	public:     
      /**
       * constructor to generate regular stand-alone executable delphicpp.
       *
       * @note This constructor depends on and is paired with one constructor, CDelphiDataMarshal(argc,argv,pTimer),
       * of class CDelphiDataMarshal
       *
       * @param[in] argc   Number of parameters in the command line
       * @param[in] argv[] paraemters in the command line
       * @param[in] pTimer pointer to an object of class CTimer to report execution time
       */
      CDelphiData(int argc,char* argv[],shared_ptr<CTimer> pTimer):pddm(new CDelphiDataMarshal(argc,argv,pTimer))
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                 CDelphiData is constructed                   *\n";
         cout << "****************************************************************\n";
#endif

         flash();
      
         pTimer->start(); // start to record program execution time

         pddm->read(pddm->strParamFile);
         
         pddm->updateParameters();
                 
         setMap();
                 
#ifdef DEBUG_DELPHI_MAP
         showMap("delphicpp_datacontainer.dat");
#endif   

         pddm.reset(); // destroy the CDelphiDataMarshal object
      };
      
      /**
       * constructor to allow delphicpp to be compiled with mcce in order to avoid intensive IO operations.
       *
       * @note This constructor depends on and is paired with one constructor, CDelphiDataMarshal(mcce_data,pTimer),
       * of class CDelphiDataMarshal
       *
       * @param[in] mcce_data A pointer to the interface struct SMCCE
       * @param[in] pTimer pointer to an object of class CTimer to report execution time
       */
      CDelphiData(SMCCE * mcce_data,shared_ptr<CTimer> pTimer):pddm(new CDelphiDataMarshal(pTimer))
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                 CDelphiData is constructed                   *\n";
         cout << "****************************************************************\n";
#endif

         flash();

         pTimer->start(); // start to record program execution time

         /*
          * Get values from struct mcce_data instead of reading the parameter file.
          * The same constrains described in getStatement and getFunction must be applied here as well
          * in order to maintain consistency.
          */

         // equivalent to gsize = value in the parameter file
         pddm->iGrid = mcce_data->gsize;
         if (5 > pddm->iGrid || 2000 < pddm->iGrid) COutOfRange_GSIZE warning(pddm->iGrid);

         // equivalent to scale = value in the parameter file
         pddm->fScale = (mcce_data->grids_per_ang+0.01*(float)mcce_data->n_retry)/pow(2,mcce_data->del_runs-1);
         if (0.0 >= pddm->fScale || 40.0 <= pddm->fScale) COutofRange_SCALE warning(pddm->fScale);

         // equivalent to in(modpdb,file=\"apbs.pqr\",format=pqr) in the parameter file
         pddm->strPdbFile = mcce_data->pdbfile;
         pddm->iPdbFormatIn = mcce_data->pdbformat;

         // equivalent to indi = value in the parameter file
         pddm->fInDielec = mcce_data->indi;
         if (1000.0 < pddm->fInDielec) COutOfRange_INDI warning(pddm->fInDielec);

         // equivalent to exdi = value in the parameter file
         pddm->fExDielec = mcce_data->exdi;
         if (0.0 >= pddm->fExDielec || 1000.0 < pddm->fExDielec) COutOfRange_EXDI warning(pddm->fExDielec);

         // equivalent to ionrad = value in the parameter file
         pddm->fIonRadius = mcce_data->ionrad;
         if (0.0 > pddm->fIonRadius || 10.0 < pddm->fIonRadius) COutOfRange_IONRAD warning(pddm->fIonRadius);

         // equivalent to prbrad = value in the parameter file
         pddm->vctfProbeRadius[0] = mcce_data->prbrad;
         if (0.0 > pddm->vctfProbeRadius[0] || 10.0 <= pddm->vctfProbeRadius[0]) COutOfRange_PRBRAD warning(pddm->vctfProbeRadius[0]);

         // equivalent to salt = value in the parameter file
         pddm->vctfSalt[0] = mcce_data->salt;
         if (0.0 > pddm->vctfSalt[0] || 10.0 < pddm->vctfSalt[0]) COutOfRange_SALT warning(pddm->vctfSalt[0]);

         // equivalent to bndcon = value in the parameter file
         pddm->iBndyType = mcce_data->bndcon;

         /*
         if (3 == mcce_data->bndcon)
         {
            // equivalent to in(phi,file="run.phi") in the parameter file
            pddm->vctfPhiMap  = mcce_data->phimap;
            pddm->fScale      = mcce_data->scale;
            pddm->gfBoxCenter = mcce_data->oldmid;
            pddm->iGrid       = mcce_data->igrid;
         }
         */
          

         // equivalent to center(x,y,z) in the parameter file
         pddm->gfOffCenter.nX = mcce_data->center[0]; pddm->gfOffCenter.nY = mcce_data->center[1]; pddm->gfOffCenter.nZ = mcce_data->center[2];

         // equivalent to out(frc,file="filename") in the parameter file
         pddm->bSiteOut = true; //pddm->strFrcFile = mcce_data->frcfile;

         // equivalent to out(phi,file="filename") in the parameter file
         //pddm->bPhimapOut = true; pddm->strPhiFile = mcce_data->phifile;

         // equivalent to site(a,c,p) in the parameter file
         pddm->bAtomInSite = true; pddm->bCoulombPotentialInSite = true; pddm->bGridPotentialInSite = true;

         // equivalent to energy(g,s) in the parameter file
         pddm->bGridEng = true; pddm->bSolvEng = true;

         // equivalent to relaxationfactor=0.8 in the parameter file
         pddm->bSpectralRadius = true; pddm->fSpectralRadius = 0.8;

         pddm->updateParameters();

         setMap();

         pddm.reset(); // destroy the CDelphiDataMarshal object
      };

    
    CDelphiData(shared_ptr<SPrime> prime_data,shared_ptr<CTimer> pTimer):pddm(new CDelphiDataMarshal(pTimer))
    {
#ifdef PRIME
        
        flash();
        
        pTimer->start();
        
        pddm->iGaussian = prime_data->iGaussian;
        
        pddm->fSigma = prime_data->fSigma;
        
        pddm->fSrfcut = prime_data->fSrfcut;
        
        pddm->fMaxc = prime_data->fMaxc;
        
        
        // equivalent to gsize = value in the parameter file
        pddm->iGrid = prime_data->gsize;
        if (5 > pddm->iGrid || 2000 < pddm->iGrid) COutOfRange_GSIZE warning(pddm->iGrid);
        
        // equivalent to scale = value in the parameter file
        pddm->fScale = prime_data->scale;
        if (0.0 >= pddm->fScale || 40.0 <= pddm->fScale) COutofRange_SCALE warning(pddm->fScale);
        
        pddm->fPercentageFill = prime_data->perfil;
//        if (0.0 >= pddm->fPercentageFill || 100.0 <= pddm->fPercentageFill) COutOfRange_PERFIL warning(pddm->fPercentageFill);
        
        // equivalent to in(modpdb,file=\"apbs.pqr\",format=pqr) in the parameter file
        pddm->strPdbFile = prime_data->pdbfile;
        pddm->iPdbFormatIn = prime_data->pdbformat;
        
        // equivalent to indi = value in the parameter file
        pddm->fInDielec = prime_data->indi;
        if (1000.0 < pddm->fInDielec) COutOfRange_INDI warning(pddm->fInDielec);
        
        // equivalent to exdi = value in the parameter file
        pddm->fExDielec = prime_data->exdi;
        if (0.0 >= pddm->fExDielec || 1000.0 < pddm->fExDielec) COutOfRange_EXDI warning(pddm->fExDielec);
        
        // equivalent to ionrad = value in the parameter file
        pddm->fIonRadius = prime_data->ionrad;
        if (0.0 > pddm->fIonRadius || 10.0 < pddm->fIonRadius) COutOfRange_IONRAD warning(pddm->fIonRadius);
        
        // equivalent to prbrad = value in the parameter file
        pddm->vctfProbeRadius[0] = prime_data->prbrad;
        if (0.0 > pddm->vctfProbeRadius[0] || 10.0 <= pddm->vctfProbeRadius[0]) COutOfRange_PRBRAD warning(pddm->vctfProbeRadius[0]);
        
        // equivalent to salt = value in the parameter file
        pddm->vctfSalt[0] = prime_data->salt;
        if (0.0 > pddm->vctfSalt[0] || 10.0 < pddm->vctfSalt[0]) COutOfRange_SALT warning(pddm->vctfSalt[0]);
        
        // equivalent to bndcon = value in the parameter file
        pddm->iBndyType = prime_data->bndcon;
        
/*
         if (3 == prime_data->bndcon)
         {
         // equivalent to in(phi,file="run.phi") in the parameter file
         pddm->vctfPhiMap  = prime_data->phimap;
         pddm->fScale      = prime_data->scale1;
         pddm->gfBoxCenter = prime_data->oldmid1;
         pddm->iGrid       = prime_data->igrid1;

         }
*/
        
        pddm->bIsAcent = prime_data->bAcent;
        // equivalent to center(x,y,z) in the parameter file
        pddm->gfAcent.nX = prime_data->center[0]; pddm->gfAcent.nY = prime_data->center[1]; pddm->gfAcent.nZ = prime_data->center[2];
        
        // equivalent to out(frc,file="filename") in the parameter file
        pddm->bSiteOut = true;
        
        pddm->strFrciFile = prime_data->strFRCIn;
        
        pddm->strFrcFile = prime_data->strFRCOut;
        
        
        
        // equivalent to out(phi,file="filename") in the parameter file
        //pddm->bPhimapOut = true; pddm->strPhiFile = prime_data->phifile;
        
        // equivalent to site(a,c,p) in the parameter file
        pddm->bAtomInSite             = true;
        pddm->bCoulombPotentialInSite = false;
        pddm->bGridPotentialInSite    = true;
        pddm->bAtomPotentialInSite    = false;
        pddm->bCrgInSite              = true;
        
        // equivalent to energy(g,s,c) in the parameter file
        pddm->bGridEng = true; pddm->bSolvEng = prime_data->bSolvEng; pddm->bCoulombEng = false;
        
        // equivalent to relaxationfactor=0.8 in the parameter file
        pddm->bSpectralRadius = true; pddm->fSpectralRadius = 0.8;
        
        pddm->strCommPDB     = prime_data->vecStrPDB;
        pddm->bPDB2FRCInSite = false;  //prime_data->bFrcPDBitself;
        
        pddm->bCommFRCIn     = prime_data->bCommFrc;
        pddm->strCommFRCIn   = prime_data->vecStrFRCIn;
        
        pddm->updateParameters();
        
        setMap();
        
        pddm.reset(); // destroy the CDelphiDataMarshal object
#endif
    };
    
    /**
     * destructor
     */
    ~CDelphiData()
    {
#ifdef DEBUG_OBJECT
        cout << endl;
        cout << "****************************************************************\n";
        cout << "*                  CDelphiData is destroyed                    *\n";
        cout << "****************************************************************\n";
#endif
    }

      /**
       * member function to implement the virtual function showMap() declared in IDataContainer
       *
       * @note used for debug purpose only
       *
       * @param[in] strMapFile Name of the txt file containing the values of variables in CDelphiData
       */
      virtual void showMap(const string& strMapFile);

      /**
       * member function to implement the virtual function reset() declared in IDataContainer
       *
       * @note used for debug purpose only
       *
       * @param[in] strF95File Name of the txt file containing the values of variables generated by delphi95
       */
      virtual void reset(const string& strF95File);
};




#endif // CDELPHIDATA_H_
