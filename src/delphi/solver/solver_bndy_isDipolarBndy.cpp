/*
 * solver_bndy_isDipolarBndy.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

bool CDelphiFastSOR::isDipolarBndy(delphi_real *** phimap)
{
   int ix,iy,iz;
   delphi_real fEpsTemp;//when Gaussian or Convolution is on, the epsvalue need to be determined.
   delphi_real fDist,fTempVal1,fTempVal2;
   SGrid<delphi_real> fgXYZ,fgTempGrid;

   fNetCrg = fMinusCrg + fPlusCrg;
   fEpsTemp=fEpsOut;
   if (( iGaussian==1 || iConvolute != 0 ) && inhomo==1) fEpsTemp=fEpsIn;

   //cout << "Lin Li: fEpsOut: " << setprecision(10) << fEpsOut << endl;
   //cout << "Lin Li: fEpsIn: " << fEpsIn << endl;
   //cout << "Lin Li: fEpsTemp: " << fEpsTemp << endl;
   //cout << "Lin Li: igaussian: " << iGaussian << "   inhomo: " << inhomo << endl;


   for (iz = 0; iz < iGrid; iz += 1)
   {
      for (iy = 0; iy < iGrid; iy += 1)
      {
          for (ix = 0; ix < iGrid; ix += iGrid-1)
          {
             fgXYZ.nX = (delphi_real)(ix+1); fgXYZ.nY = (delphi_real)(iy+1); fgXYZ.nZ = (delphi_real)(iz+1);

             fgTempGrid = fgPlusCrgCenter - fgXYZ;
             fDist = sqrt(optDot<delphi_real>(fgTempGrid,fgTempGrid))/fScale;
             fTempVal1 = fPlusCrg*exp(-fDist/fDebyeLength)/(fDist*fEpsTemp);

             fgTempGrid = fgMinusCrgCenter - fgXYZ;
             fDist = sqrt(optDot<delphi_real>(fgTempGrid,fgTempGrid))/fScale;
             fTempVal2 = fMinusCrg*exp(-fDist/fDebyeLength)/(fDist*fEpsTemp);

             phimap[iz][iy][ix] = fTempVal1 + fTempVal2;
          }
      }
   }

   for (iz = 0; iz < iGrid; iz += 1)
   {
      for (iy = 0; iy < iGrid; iy += iGrid-1)
      {
          for (ix = 0; ix < iGrid; ix += 1)
          {
             fgXYZ.nX = (delphi_real)(ix+1); fgXYZ.nY = (delphi_real)(iy+1); fgXYZ.nZ = (delphi_real)(iz+1);

             fgTempGrid = fgPlusCrgCenter - fgXYZ;
             fDist = sqrt(optDot<delphi_real>(fgTempGrid,fgTempGrid))/fScale;
             fTempVal1 = fPlusCrg*exp(-fDist/fDebyeLength)/(fDist*fEpsTemp);

             fgTempGrid = fgMinusCrgCenter - fgXYZ;
             fDist = sqrt(optDot<delphi_real>(fgTempGrid,fgTempGrid))/fScale;
             fTempVal2 = fMinusCrg*exp(-fDist/fDebyeLength)/(fDist*fEpsTemp);

             phimap[iz][iy][ix] = fTempVal1 + fTempVal2;
          }
      }
   }

   for (iz = 0; iz < iGrid; iz += iGrid-1)
   {
      for (iy = 0; iy < iGrid; iy += 1)
      {
         for (ix = 0; ix < iGrid; ix += 1)
         {
            fgXYZ.nX = (delphi_real)(ix+1); fgXYZ.nY = (delphi_real)(iy+1); fgXYZ.nZ = (delphi_real)(iz+1);

            fgTempGrid = fgPlusCrgCenter - fgXYZ;
            fDist = sqrt(optDot<delphi_real>(fgTempGrid,fgTempGrid))/fScale;
            fTempVal1 = fPlusCrg*exp(-fDist/fDebyeLength)/(fDist*fEpsTemp);

            fgTempGrid = fgMinusCrgCenter - fgXYZ;
            fDist = sqrt(optDot<delphi_real>(fgTempGrid,fgTempGrid))/fScale;
            fTempVal2 = fMinusCrg*exp(-fDist/fDebyeLength)/(fDist*fEpsTemp);

            phimap[iz][iy][ix] = fTempVal1 + fTempVal2;
         }
      }
   }

   return true;
}
