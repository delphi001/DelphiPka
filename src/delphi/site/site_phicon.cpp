/*
 * site_phicon.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: chuan
 *
 *  Description:
 *     convert potentials to concentrations
 */

#include "site.h"

void CSite::phicon()
{
   delphi_real tmp,temp,phi,phisq;
   delphi_integer ix,iy,iz,iw;

   tmp = abs(fTaylorCoeff2*fTaylorCoeff4);

   if (0.0 < fIonStrength)
   {
      cout << "\nconverting potentials to \n";
      cout << "net charge concentrations...\n\n";

      /*
       * if nonlinear equation is used then use the exponential form, otherwise use the linear form
       * use same number of terms in expansion of sinh as for iteration in itit.f
       */
      {
         CUntestedPhicon warning;
      }

      if (0 != iNonIterateNum)
      {
         if (fZero > tmp)
         {
            for (iz = 0; iz < iGrid; iz++)
            {
               for (iy = 0; iy < iGrid; iy++)
               {
                  for (ix = 0; ix < iGrid; ix++)
                  {
                     //----- use first three terms of sinh
                     iw = iz*iGrid*iGrid+iy*iGrid+ix;
                     if (prgbDielecMap[iw])
                     {
                        phi   = phimap[iz][iy][ix];
                        phisq = phi*phi;
                        //----- Horner scheme for charge and osmotic term
                        temp  = phisq*fTaylorCoeff5 + fTaylorCoeff3;
                        temp  = temp*phisq          + fTaylorCoeff1;
                        phimap[iz][iy][ix] = temp*phi;
                     }
                     else
                        phimap[iz][iy][ix] = 0.0;
                  }
               }
            }
         }
         else
         {
            //----- asymmetric salt
            for (iz = 0; iz < iGrid; iz++)
            {
               for (iy = 0; iy < iGrid; iy++)
               {
                  for (ix = 0; ix < iGrid; ix++)
                  {
                     iw = iz*iGrid*iGrid+iy*iGrid+ix;
                     if (prgbDielecMap[iw])
                     {
                        phi  = phimap[iz][iy][ix];
                        //----- Horner scheme for charge and osmotic term
                        temp = phi*fTaylorCoeff5 + fTaylorCoeff4;
                        temp = phi*temp          + fTaylorCoeff3;
                        temp = phi*temp          + fTaylorCoeff2;
                        temp = phi*temp          + fTaylorCoeff1;
                        phimap[iz][iy][ix] = temp*phi;
                     }
                     else
                        phimap[iz][iy][ix] = 0.0;
                  }
               }
            }
         }
      }
      else
      {
         for (iz = 0; iz < iGrid; iz++)
         {
            for (iy = 0; iy < iGrid; iy++)
            {
               for (ix = 0; ix < iGrid; ix++)
               {
                  iw = iz*iGrid*iGrid+iy*iGrid+ix;
                  if (prgbDielecMap[iw])
                  {
                     phi = phimap[iz][iy][ix];
                     phimap[iz][iy][ix] = fTaylorCoeff1*phi;
                  }
                  else
                     phimap[iz][iy][ix] = 0.0;
               }
            }
         }
      }
   }

   if (fZero > abs(fIonStrength)) CNoPotential2CrgConcentrate warning;
}
