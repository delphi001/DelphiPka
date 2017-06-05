/*
 * solver_fastSOR_relfac.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

delphi_real CDelphiFastSOR::calculateRelaxFactor()
{
   delphi_real fRelaxFactor = 0.0;

   int ix,iy,iz;
   delphi_integer iw;
   delphi_real temp,temp2,temp3;

   //----- setup sn array for lowest eigenstate
   vector<delphi_real> sn1(iGrid,0.0);

   for(ix = 1; ix < iGrid-1; ix++)
   {
      temp = fPi*ix/(iGrid-1);
      sn1[ix] = sqrt(2.0)*sin(temp)/sqrt((delphi_real)(iGrid-1));
   }

   vector<delphi_real> sn2(sn1),sn3(sn1);

   delphi_real recipr = 1.0/sqrt((delphi_real)iGrid);
   if (rgbPeriodicBndy[0]) sn1.assign(iGrid,recipr);
   if (rgbPeriodicBndy[1]) sn2.assign(iGrid,recipr);
   if (rgbPeriodicBndy[0]) sn3.assign(iGrid,recipr);

   //----- map sn arrays to prgfPhiMap
   iw = 0;
   for (iz = 0; iz < iGrid; iz++)
   {
      temp3 = sn3[iz];
      for (iy = 0; iy < iGrid; iy++)
      {
         temp2 = temp3*sn2[iy];
         for (ix = 0; ix < iGrid; ix++)
         {
            prgfPhiMap[iw] = temp2*sn1[ix];
            iw++;
         }
      }
   }

   //----- setup periodic boundaries, start and stop vectors etc. for odd/even loops
   initOddEvenItr(0); // forWhom = 0

   //----- iterate once over odd and even points
   itrOddPoints(0,10001); // forWhom = 0

   itrEvenPoints(0,10002); // forWhom = 0

   //----- caculate the estimated spectral radius
   fSpec = 0.0;

   for (unsigned int i = 0; i != phimap1.size(); i++)
      fSpec += phimap1[i]*phimap2[i];

   fSpec = 2.0*fSpec;

   //----- following needed as spec exceeds 1.0 occasionally in focussing calculations (SS May 8, 1998)
   //if (1.0 < fSpec) fSpec = 0.995;

#ifdef VERBOSE
   cout << "\n gauss-seidel spectral radius is " << fSpec << endl;
#endif

#ifdef DEBUG_DELPHI_SOLVER_RELFAC
   {
      string strTestFile = "test_relfac.dat";
      ofstream ofTestStream(strTestFile.c_str());
      ofTestStream << boolalpha;
      ofTestStream << fixed << setprecision(7);

      for (iw = 0; iw < iHalfGridNum; iw++)
         ofTestStream << "phimap1(" << setw(6) << right << iw+1 << ") = " << setw(11) << phimap1[iw] << endl;

      for (iw = 0; iw < iHalfGridNum; iw++)
         ofTestStream << "phimap2(" << setw(6) << right << iw+1 << ") = " << setw(11) << phimap2[iw] << endl;

      for (iw = 0; iw < iGrid*iGrid*iGrid; iw++)
         ofTestStream << "phimap3(" << setw(6) << right << iw+1 << ") = " << setw(11) << prgfPhiMap[iw] << endl;



   }
#endif // DEBUG_DELPHI_SOLVER


   prgfPhiMap.assign(iGrid*iGrid*iGrid,0.0); // restore prgfPhiMap

   return fRelaxFactor;
}
