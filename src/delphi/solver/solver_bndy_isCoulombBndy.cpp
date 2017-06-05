/*
 * solver_bndy_isCoulombBndy.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

bool CDelphiFastSOR::isCoulombBndy(delphi_real *** phimap)
{
   int ix,iy,iz,ic;
   delphi_real fDist,fTempVal;
   SGrid<delphi_real> fgXYZ,fgTempGrid;

   fNetCrg = fMinusCrg + fPlusCrg;

   for (iz = 0; iz < iGrid; iz += 1)
   {
      for (iy = 0; iy < iGrid; iy += 1)
      {
          for (ix = 0; ix < iGrid; ix += iGrid-1)
          {
             fgXYZ.nX = (delphi_real)(ix+1); fgXYZ.nY = (delphi_real)(iy+1); fgXYZ.nZ = (delphi_real)(iz+1);

             for (ic = 0; ic < iCrgGridNum; ic += 1)
             {
                fgTempGrid = prggvCrgedAtom[ic].nGrid - fgXYZ;
                fDist      = sqrt(optDot<delphi_real>(fgTempGrid,fgTempGrid))/fScale;
                fTempVal   = prggvCrgedAtom[ic].nValue*exp(-fDist/fDebyeLength)/(fDist*fEpsOut);
                phimap[iz][iy][ix] += fTempVal;
             }
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

             for (ic = 0; ic < iCrgGridNum; ic += 1)
             {
                fgTempGrid = prggvCrgedAtom[ic].nGrid - fgXYZ;
                fDist      = sqrt(optDot<delphi_real>(fgTempGrid,fgTempGrid))/fScale;
                fTempVal   = prggvCrgedAtom[ic].nValue*exp(-fDist/fDebyeLength)/(fDist*fEpsOut);
                phimap[iz][iy][ix] += fTempVal;
             }
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

            for (ic = 0; ic < iCrgGridNum; ic += 1)
            {
               fgTempGrid = prggvCrgedAtom[ic].nGrid - fgXYZ;
               fDist      = sqrt(optDot<delphi_real>(fgTempGrid,fgTempGrid))/fScale;
               fTempVal   = prggvCrgedAtom[ic].nValue*exp(-fDist/fDebyeLength)/(fDist*fEpsOut);
               phimap[iz][iy][ix] += fTempVal;
            }
         }
      }
   }

   return true;
}
