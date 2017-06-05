/*
 * solver_fastSOR_initOddEvenItr.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

void CDelphiFastSOR::initOddEvenItr(const int& forWhom)
{
   int ix,iy,iz;
   delphi_integer ip,iq;
   delphi_integer iadd1,iadd2,itemp1,itemp2,star,fin;

   /*
    * setup periodic boundaries
    */
   if (rgbPeriodicBndy[0])
   {
      for (iz = 1; iz < iGrid-1; iz++)
      {
         iadd1 = iz*iGrid*iGrid;
         for (iy = 1; iy < iGrid-1; iy++)
         {
            iadd2 = (iadd1+iy*iGrid+2)/2;
            ibndx.push_back(iadd2);
         }
      }

      idif1x = (iGrid-2)/2; idif2x = idif1x+1; inc1xa = 1; inc1xb = 0; inc2xa = 0; inc2xb = 1;
   }

   if (rgbPeriodicBndy[1])
   {
      for (iz = 1; iz < iGrid-1; iz++)
      {
         iadd1 = iz*iGrid*iGrid;
         for (ix = 1; ix < iGrid-1; ix++)
         {
            iadd2 = (iadd1+ix+2)/2;
            ibndy.push_back(iadd2);
         }
      }

      idif1y = iGrid*(iGrid-2)/2; idif2y = idif1y+1; inc1ya = iGrid/2+1; inc1yb = inc1ya-1; inc2ya = inc1yb; inc2yb = inc1ya;
   }

   if (rgbPeriodicBndy[2])
   {
      for (ix = 1; ix < iGrid-1; ix++)
      {
         iadd1 = ix+2;
         for (iy = 1; iy < iGrid-1; iy++)
         {
            iadd2 = (iadd1+iy*iGrid)/2;
            ibndz.push_back(iadd2);
         }
      }

      idif1z = iGrid*iGrid*(iGrid-2)/2; idif2z = idif1z+1; inc1za = iGrid*iGrid/2+1; inc1zb = inc1za; inc2za = inc1zb; inc2zb = inc1za;
   }

   /*
    * setup start and stop vectors
    */
   sta1.assign(iGrid,0); sta2.assign(iGrid,0); fi1.assign(iGrid,0); fi2.assign(iGrid,0);

   sta1[1] = (iGrid*iGrid+iGrid+4)/2; sta2[1] = sta1[1]-1; fi1[1] = iGrid*iGrid-(iGrid+1)/2; fi2[1] = fi1[1];

   itemp1 = iGrid+2; itemp2 = iGrid*iGrid-iGrid-2;
   for (delphi_integer i = 2; i < iGrid-1; i++)
   {
      sta1[i] = fi1[i-1]+itemp1; sta2[i] = fi2[i-1]+itemp1; fi1[i] = sta1[i-1]+itemp2; fi2[i] = sta2[i-1]+itemp2;
   }

   lat1 = (iGrid-1)/2; lat2 = (iGrid+1)/2; long1 = (iGrid*iGrid-1)/2; long2 = (iGrid*iGrid+1)/2;

   /*
    * split prgfPhiMap to odd and even vectors
    */
    //cout << "ip,iq: " << ip << " " << iq <<endl;
   for (ip = 0; ip < iHalfGridNum; ip++)
   {
      iq = ip*2;
        //cout << "ip,iq: " << ip << " " << iq <<endl;
      if (prgfPhiMap.size() > iq)
         phimap1.push_back(prgfPhiMap[iq]);
      else
         phimap1.push_back(0.0);

      if (prgfPhiMap.size() > iq+1)
         phimap2.push_back(prgfPhiMap[iq+1]);
      else
         phimap2.push_back(0.0);
   }

    //for (int ii=1;ii<=phimap2.size();ii++)
    //for (int ii=1;ii<=700;ii++)
    //    {
    //        cout << "Aphimap2: " << setw(6) << right << ii << " " << setw(11) << right << phimap2[ii-1] << endl;
    //    }
   /*
    * setup vectors for restoring x boundary values
    */
   star = (iGrid+1)/2; iq = iGrid*(iGrid+1)/2-iGrid+1; fin = (iGrid*(iGrid-1)-2)/2;
   for(ip = 0; ip < fin-star+1; ip++)
   {
      iq += iGrid;
      bndx1.push_back(phimap1[iq-1]);
      bndx2.push_back(phimap1[iq+((iGrid+1)/2-1)-1]);
   }

   star = (iGrid+2)/2; iq = iGrid*(iGrid+2)/2-iGrid+1; fin = (iGrid*(iGrid-1)-1)/2;
   for(ip = 0; ip < fin-star+1; ip++)
   {
      iq += iGrid;
      bndx3.push_back(phimap2[iq-1]);
      bndx4.push_back(phimap2[iq+((iGrid+1)/2-1)-1]);
   }

   if (0 == forWhom || (1 == forWhom && bFixedRelaxParam))
   {
      om2 = 1.0; sixth = fSixth;
   }
   else
   {
      om2=2.0/(1.0+sqrt(1.0-fSpec));

      for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
         *it = (*it)*om2;

      for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
         *it = (*it)*om2;

      for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
         *it = (*it)*om2;

      for (iy = 0; iy < iDielecBndyOdd; iy++)
         for (ix = 0; ix < 6; ix++){
            //cout << "iy,ix1: " << iy << " " << ix  << endl;
            //cout << "iy,ix2: " << iy << " " << ix << " " << prgfBndyDielec[iy][ix] << " " << iDielecBndyOdd << endl;
            prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix]*om2;

         }


      sixth = fSixth*om2;
   }

   om1 = 1.0 - om2;

}
