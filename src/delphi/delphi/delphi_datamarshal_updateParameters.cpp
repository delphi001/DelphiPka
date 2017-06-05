/*
 * delphi_datamarshal_updateParameters.cpp
 *
 *  Created on: Mar 18, 2014
 *      Author: chuan
 */

#include "delphi_datamarshal.h"

//-----------------------------------------------------------------------//
void CDelphiDataMarshal::updateParameters()
{
   // ++++++ update parameters after reading the parameter file ++++++//
   //*****************************************************************//
   //                                                                 //
   //        perform updates after reading the parameter file         //
   //                                                                 //
   //*****************************************************************//

   if (-1 == vctfProbeRadius[1])
      vctfProbeRadius[1] = vctfProbeRadius[0];

   /*
    * parameters set at the end of subroutine rdprm
    */
   if (0 > fExDielec || 0 > fInDielec)
   {
      CMinusDielec warning(fExDielec, fInDielec);

      fExDielec = abs(fExDielec);
      fInDielec = abs(fInDielec);
   }

   delphi_real fZ1Plus  = vctiValence1[0]; // valence 1
   delphi_real fZ1Minus = vctiValence1[1];
   delphi_real fZ2Plus  = vctiValence2[0]; // valence 2
   delphi_real fZ2Minus = vctiValence2[1];

   /*
    * concentration of positive ion
    */
   delphi_real fZ1PlusConcentrate = vctfSalt[0]*fZ1Minus;
   delphi_real fZ2PlusConcentrate = vctfSalt[1]*fZ2Minus;

   fIonStrength = (fZ1PlusConcentrate*fZ1Plus*(fZ1Plus+fZ1Minus) + fZ2PlusConcentrate*fZ2Plus*(fZ2Plus+fZ2Minus))/2.0;

   /*
    * coefficients in Taylor series of the charge concentration apart from n! (order >=1)
    * (chuan 2012Apr24) Correct coefficients in Taylor series. NOT in compact form just for clean math formula
    */
   fTaylorCoeff1 = -2.0*fIonStrength;

   fTaylorCoeff2 =  ( fZ1PlusConcentrate*fZ1Plus *pow(fZ1Plus, 2) -
                      fZ1PlusConcentrate*fZ1Minus*pow(fZ1Minus,2) +
                      fZ2PlusConcentrate*fZ2Plus *pow(fZ2Plus, 2) -
                      fZ2PlusConcentrate*fZ2Minus*pow(fZ2Minus,2) )/2.0;

   fTaylorCoeff3 = -( fZ1PlusConcentrate*fZ1Plus *pow(fZ1Plus, 3) +
                      fZ1PlusConcentrate*fZ1Minus*pow(fZ1Minus,3) +
                      fZ2PlusConcentrate*fZ2Plus *pow(fZ2Plus, 3) +
                      fZ2PlusConcentrate*fZ2Minus*pow(fZ2Minus,3) )/6.0;

   fTaylorCoeff4 =  ( fZ1PlusConcentrate*fZ1Plus *pow(fZ1Plus, 4) -
                      fZ1PlusConcentrate*fZ1Minus*pow(fZ1Minus,4) +
                      fZ2PlusConcentrate*fZ2Plus *pow(fZ2Plus, 4) -
                      fZ2PlusConcentrate*fZ2Minus*pow(fZ2Minus,4) )/24.0;

   fTaylorCoeff5 = -( fZ1PlusConcentrate*fZ1Plus *pow(fZ1Plus, 5) +
                      fZ1PlusConcentrate*fZ1Minus*pow(fZ1Minus,5) +
                      fZ2PlusConcentrate*fZ2Plus *pow(fZ2Plus, 5) +
                      fZ2PlusConcentrate*fZ2Minus*pow(fZ2Minus,5) )/120.0;

   /*
    * convert ionic strength to debye length
    */
   if (fZero < fIonStrength)
   {
      delphi_real fDebyeFactor = 0.01990076478*sqrt(fTemper*fExDielec);
      fDebyeLength = fDebyeFactor/sqrt(fIonStrength);

      if (0 < iNonIterateNum) bNonlinearEng = true;
   }
   else
   {
      bIonsEng = false;
      fDebyeLength = 1.0e6;
   }

   /*
    * epkt assignment as a function of temperature
    */
   fEPKT = dEPK/fTemper;

   /*
    * set epsin and epsout (= epkt adjusted dielectrics such that all distances are in angstroms, charges in e)
    */
   fEpsIn  = fInDielec/fEPKT;
   fEpsOut = fExDielec/fEPKT;

   /*
    * test for unformatted pdb and frc files
    * the same tests have been implemented in class CIO. However, here users may make mistakes in parameter file to
    * require wrong format of the files. So check again to reset format flags.
    */
   string strASCI = "1234567890 .-+#,$asdfghjklzxcvbnmqwertyuiopASDFGHJKLZXCVBNMQWERTYUIOP)(}{][/";

   ifstream ifFileHandle;
   char cTestChar[80];

   if (!bPdbUnformatIn) // PDB file
   {
      ifFileHandle.open(strPdbFile.c_str());

      if (!ifFileHandle.is_open()) throw CUnknownFile(strPdbFile);

      ifFileHandle.read(cTestChar,80);

      int iCount = 0;

      for (int i = 0; i < 80; i++)
      {
         if (string::npos == strASCI.find(cTestChar[i]) ) iCount += 1;
      }

      if (10 < iCount) // unformatted PDB
      {
         bPdbUnformatIn = true;
         CToUnformattedFile warning(strPdbFile,strASCI);
      }

      ifFileHandle.close();
   }

   if (bUnformatFrcOut) // FRC file
   {
      ifFileHandle.open(strFrcFile.c_str());

      if (!ifFileHandle.is_open()) throw CUnknownFile(strFrcFile);

      ifFileHandle.read(cTestChar,80);

      int iCount = 0;

      for (int i = 0; i < 80; i++)
      {
         if (string::npos == strASCI.find(cTestChar[i]) ) iCount += 1;
      }

      if (10 < iCount) // unformatted FRC
      {
         CToUnformattedFile warning(strFrcFile,strASCI);
         bFrcUnformatIn = true;
      }

      ifFileHandle.close();
   }

   //*****************************************************************//
   //                                                                 //
   //               read size, charge, PDB files                      //
   //                                                                 //
   //*****************************************************************//
   unique_ptr<CIO> pIO(new CIO(fInDielec,fEPKT)); // smart unique_ptr

   pIO->setDelphiAtom(bSolvePB,bSurfCrgInSite,strSizeFile,strCrgFile,strPdbFile,iPdbFormatIn,bPdbUnformatIn,strCommPDB);

   iMediaNum        = pIO->iMediaNum;        // nmedia
   iObjectNum       = pIO->iObjectNum;       // nobject
   iAtomNum         = pIO->iAtomNum;         // natom
   iResidueNum      = pIO->iResidueNum;      // resnummax
   bOnlyMolecule    = pIO->bOnlyMolecule;    // ionlymol
   vctapAtomPdb     = pIO->vctapAtomPdb;     // delphipdb(natom)
   vctfMediaEps     = pIO->vctfMediaEps;     // medeps(0:nmediamax)
   vctstrObject     = pIO->vctstrObject;     // dataobject(nobjectmax,2)
   vctiAtomMediaNum = pIO->vctiAtomMediaNum; // iatmmed(Natom+Nobjectmax)

   /*
    * write unformatted PDB file
    */
   if (bUnformatPdbOut) pIO->writeUnformatPdb(strUnformatPdbFile);

   /*
    * write mod/pqr/mod4/pqr4-type PDB file
    */
   if (bModPdbOut)
      pIO->writeModifiedPdb(strModifiedPdbFile,iModPdbFormatOut);

   fEpsIn = pIO->vctfMediaEps[1];

   if (1 < iMediaNum) iDirectEpsMap = 1;

   //--------------------------- extrmobjects ------------------------//

   /*
    * extrmobjects: find extrema of each object and  according to them limobject contains extreme
    * values of each object for a molecule it has extreme but without radii
    */
   iMoleculeNum = 0;    // numbmol
   fMaxRadius   = 0.01; // rdmx

   SExtrema<delphi_real> tmpExtrema;
   for (unsigned int i = 0; i < vctstrObject.size(); i=i+2)
   {
      string strLine = vctstrObject[2*i]; // vctstrObject[0][nobject]

      if (0 != strLine.compare(0,4,"is a")) throw CIsAnObjectType(strLine);

      iMoleculeNum += 1;

#ifdef VERBOSE
      cout << "Object number " << i << " is a molecule\n";
#endif

      //cout << "## iAtomNum: " << iAtomNum << endl;
      if (0 == iAtomNum) throw CNoAtomsInMolecule(i);

      SGrid<delphi_real> gMinCoord = vctapAtomPdb[0].getPose();
      SGrid<delphi_real> gMaxCoord = vctapAtomPdb[0].getPose();
      fMaxRadius            = vctapAtomPdb[0].getRadius();
      for (delphi_integer j = 1; j < iAtomNum; j++)
      {
         gMinCoord = optMin<delphi_real>(gMinCoord,vctapAtomPdb[j].getPose());
         gMaxCoord = optMax<delphi_real>(gMaxCoord,vctapAtomPdb[j].getPose());
         fMaxRadius = max(fMaxRadius,vctapAtomPdb[j].getRadius());
      }

      tmpExtrema.nMin = gMinCoord; tmpExtrema.nMax = gMaxCoord;
      vctefExtrema.push_back(tmpExtrema);
   }

   //------------------------------ extrm ---------------------------//

   /*
    * find extrema and calculate scale according to them and to the percent box fill
    */
   gfMinCoordinate.nX = 6000.0; gfMinCoordinate.nY = 6000.0; gfMinCoordinate.nZ = 6000.0;

   gfMaxCoordinate.nX =-6000.0; gfMaxCoordinate.nY =-6000.0; gfMaxCoordinate.nZ =-6000.0;

   for (delphi_integer i = 0; i < iAtomNum; i++)
   {
      gfMinCoordinate = optMin<delphi_real>(gfMinCoordinate,vctapAtomPdb[i].getPose()-vctapAtomPdb[i].getRadius());
      gfMaxCoordinate = optMax<delphi_real>(gfMaxCoordinate,vctapAtomPdb[i].getPose()+vctapAtomPdb[i].getRadius());
   }

   gfGeometricCenter = (gfMinCoordinate + gfMaxCoordinate)/2.0;

   //------------------------------- off ----------------------------//
   if (bIsAcent)
      gfBoxCenter = gfAcent;
   else
   {
      gfBoxCenter = gfGeometricCenter - gfOffCenter/fScale;

      if (fZero > abs(gfOffCenter.nX-999.0) || fZero > abs(gfOffCenter.nX-777.0))
      {
         if (fZero > abs(gfOffCenter.nX-999.0))
         {
            cout << "modifying midpoints using frc input file \n";
            gfBoxCenter = pIO->readFrcFile(strFrciFile,gfOffCenter,fScale);
         }
         else
         {
            cout << "modifying midpoints using fort.27 \n";
            gfBoxCenter = pIO->readFrcFile(strCentFile,gfOffCenter,fScale);
         }
      } // ---------- end of if (fZero > abs(gfOffCenter-999.0) || fZero > abs(gfOffCenter-777.0))
   } // ---------- end of if (bIsAcent)

   gfCoordinateRange = gfMaxCoordinate - gfMinCoordinate;

   {
      SGrid<delphi_real> gVec1 = 2.0*optABS<delphi_real>(gfMaxCoordinate-gfBoxCenter);
      delphi_real fMaxVal1 = optMax<delphi_real>(gVec1);

      SGrid<delphi_real> gVec2 = 2.0*optABS<delphi_real>(gfMinCoordinate-gfBoxCenter);
      delphi_real fMaxVal2 = optMax<delphi_real>(gVec2);

      fMaxDimension = (fMaxVal1 > fMaxVal2) ? fMaxVal1 : fMaxVal2;
   }

   if (0 == iGrid)
   {
      if (fZero > abs(fScale-10000.0) )          fScale = 2.0;
      if (fZero > abs(fPercentageFill-10000.0) ) fPercentageFill = 80.0;
      iGrid = fScale*100.0/fPercentageFill*fMaxDimension;
   }
   else if (fZero > abs(fScale-10000.0) )
   {
      if (fZero > abs(fPercentageFill-10000.0) )
      {
         fScale          = 2.0;
         fPercentageFill = 100.0*fMaxDimension*fScale/(iGrid-1);
      }
      else
         fScale = (iGrid-1)*fPercentageFill/(100.0*fMaxDimension);
   }
   else
      fPercentageFill = 100.0*fMaxDimension*fScale/(iGrid-1);

   if (0 == iGrid%2) iGrid += 1;

   vctfMediaEps[0] = fEpsOut;

   if (bPdbUnformatIn && bUnformatPdbOut)
      CReadWriteUnformatPdb warning(bUnformatPdbOut);

   if (bFrcUnformatIn && bUnformatFrcOut)
      CReadWriteUnformatFrc warning(bUnformatPdbOut);


   SGrid<delphi_real> gLeftBndy  = gfBoxCenter - (1.0/fScale)*(iGrid+1)*0.5;
   SGrid<delphi_real> gRightBndy = gfBoxCenter + (1.0/fScale)*(iGrid+1)*0.5;
   if (optORLT(gfMinCoordinate,gLeftBndy) || optORGT(gfMaxCoordinate,gRightBndy))
      CSystemOutsideBox warning;

   /*
    * convert atom coordinates from angstroms to grid units
    * Added allocation of xn1 array in order to deal with different names in real (atpos) and dummy (xn1)
    * arguments in following subroutines
    */
   for (delphi_integer i = 0; i < iAtomNum; i++)
   {
      vctgfAtomCoordA.push_back( vctapAtomPdb[i].getPose() );
      vctgfAtomCoordG.push_back( (vctapAtomPdb[i].getPose()-gfBoxCenter)*fScale+(delphi_real)((iGrid+1)/2) );
   }

   /*
    * verify if dielectric is uniform
    */
   bUniformDielec = true;
   for (delphi_integer i = 0; i < iMediaNum; i++)
   {
      if (vctfMediaEps[i] != vctfMediaEps[i+1]) bUniformDielec = false;
   }

   /*
    * now pass uniformdiel to epsmak make the epsmap, and also a listing of boundary elements, and the
    * second epsmap used for the molecular surface scaling
    */
   //prgiEpsMap.assign(iGrid*iGrid*iGrid,0);
   //vctbDielecMap.assign(iGrid*iGrid*iGrid,false);

   /*
    * new updates in c++ for more rigid check of input parameter values
    */
   if (false == bAutoConverge && 0 == iLinIterateNum && 0 == iNonIterateNum)
      throw CBadAutoConvergence(bAutoConverge);

   showParameters();

   pIO.reset();

}

