/*
 * site_writePhi.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: chuan
 */

#include "site.h"

void CSite::writePhi()
{
   vector<delphi_real> phimap3,phimap4;

   if (2 != iPhiFormatOut)
      phimap4.assign(prgfPhiMap.begin(),prgfPhiMap.end()); // need a copy of phimap

   if (bOutCrgDensity)
   {
      phimap3.assign(prgfPhiMap.begin(),prgfPhiMap.end());
      phicon();
   }

   if (bBiosystemOut)
      writePotential_insight(phimap4);
   else if (2 == iPhiFormatOut)
      writePotential_grasp(phimap4);
   else if (3 == iPhiFormatOut)
      writePotential_ccp4(phimap4);
   else if (4 == iPhiFormatOut)
      writePotential_fromPrevious(phimap4);
   else if (5 == iPhiFormatOut)
      writePotential_cube();
   else
      writePotential_delphi();

   if (2 != iPhiFormatOut)
      vector<delphi_real>().swap(phimap4);

   if (bOutCrgDensity)
   {
      prgfPhiMap.assign(phimap3.begin(),phimap3.end());
      vector<delphi_real>().swap(phimap3);
   }
}
