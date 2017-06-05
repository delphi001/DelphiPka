/*
 * site_rforce.cpp
 *
 *  Created on: Feb 22, 2014
 *      Author: chuan
 */

#include "site.h"

vector< SGrid<delphi_real> > CSite::rforce()
{
   vector< SGrid<delphi_real> > afield;

   SGrid<delphi_real> sfield = {0.0,0.0,0.0},vtemp;
   delphi_real sc,dist,temp;
   delphi_integer iat,i,j;

   afield.assign(iBndyGridNum,sfield);

   for (i = 0; i < iBndyGridNum; i++)
   {
      sfield.nX = 0.0; sfield.nY = 0.0; sfield.nZ = 0.0;
      sc = prgfSurfCrgE[i];

      for (j = 0; j < iCrgGridNum; j++) // calculate total field on this surface element due to ALL charges
      {
         vtemp  = prgfgSurfCrgA[i] - prgfgCrgPoseA[j];
         dist   = optDot(vtemp,vtemp);
         temp   = prggvAtomicCrg[j].nValue/(dist*sqrt(dist));
         sfield = sfield + temp*vtemp;
         iat    = prgiCrgAt[j];
         if (0 >  iat) continue;
         if (0 == iat) throw CCrgatnError();
         afield[iat-1] = afield[iat-1] - (sc*temp)*vtemp;
      }

      sfield = sfield*prgfSurfCrgE[i];
      j = prgiAtSurf[i];
      afield[j-1] = afield[j-1] + sfield;
   }

   return afield;
}
