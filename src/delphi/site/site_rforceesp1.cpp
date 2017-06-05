/*
 * site_rforceesp1.cpp
 *
 *  Created on: Feb 22, 2014
 *      Author: chuan
 */

#include "site.h"

vector< SGrid<delphi_real> > CSite::rforceeps1()
{
   vector< SGrid<delphi_real> > afield;

   SGrid<delphi_real> sfield = {0.0,0.0,0.0},xyz,vtemp;
   SGrid<int> ixyz;
   delphi_real fact,fact1,sc,vvtemp,deltas,sigmap,realds,rmyds,trullo,dist,temp;
   delphi_integer p,i,j,cont,iat;

   afield.assign(iBndyGridNum,sfield);
   fact = -2.0*fPi*80.0/79.0;
   realds = 0.0; rmyds = 0.0; cont = 0;

   for (p = 0; p < iBndyGridNum; p++)
   {
      sfield.nX = 0.0; sfield.nY = 0.0; sfield.nZ = 0.0;

      //----- rounds to the nearest whole number
      ixyz.nX = floor((prgfgSurfCrgA[p]-fgBoxCenter).nX*fScale+0.5);
      ixyz.nY = floor((prgfgSurfCrgA[p]-fgBoxCenter).nY*fScale+0.5);
      ixyz.nZ = floor((prgfgSurfCrgA[p]-fgBoxCenter).nZ*fScale+0.5);

      xyz.nX = (delphi_real)ixyz.nX/fScale+fgBoxCenter.nX;
      xyz.nY = (delphi_real)ixyz.nY/fScale+fgBoxCenter.nY;
      xyz.nZ = (delphi_real)ixyz.nZ/fScale+fgBoxCenter.nZ;

      sc = 0.5*prgfSurfCrgE[p];
      vtemp = prgfgSurfCrgA[p]-xyz; vvtemp = optDot(vtemp,prgfgSurfCrgE[p]);
      fact1 = 0.8/(fScale*fScale)-vvtemp*vvtemp;
      deltas = fPi*fact1; sigmap = prgfSurfCrgE[p]/deltas;
      realds += prgfSurfCrgE[p]/0.0393; rmyds += deltas;
      trullo = fact*sigmap*sigmap*deltas;

      if (400 == (p+1))
      {
         cout << "normale " << prgfgSurfCrgE[p].nX << " " << prgfgSurfCrgE[p].nY << " " << prgfgSurfCrgE[p].nZ << endl;
         cout << "P " << vtemp << endl;
         cout << prgfgSurfCrgA[p].nX << " " << prgfgSurfCrgA[p].nY << " " << prgfgSurfCrgA[p].nZ << endl;
         cout << "area " << prgfSurfCrgE[p]/0.0393 << " " << prgfSurfCrgE[p]*fScale*fScale/0.0393 << endl;
      }

      j = prgiAtSurf[p];

      for (iat = 0; iat < iAtomNum; iat++)
      {
         if ((j-1) == iat)
         {
            cont++; afield[iat].nX += trullo;
         }
      }

      for (i = 0; i < iCrgGridNum; i++) // calculate total field on this surface element due to ALL charges
      {
         vtemp  = prgfgSurfCrgA[p] - prgfgCrgPoseA[i];
         dist   = optDot(vtemp,vtemp);
         temp   = prggvAtomicCrg[j].nValue/(dist*sqrt(dist));
         sfield = sfield + temp*vtemp;

         iat    = prgiCrgAt[i];
         if (0 >  iat) continue;
         if (0 == iat) throw CCrgatnError();
         if (j == iat) afield[iat-1] = afield[iat-1] - (sc*temp)*vtemp+fact1*prgfgSurfCrgE[p];
      }

      sfield = sfield*sc;
   }

   cout << "supcalc " << realds/(16.0*fPi) << " mia " << rmyds/(16.0*fPi) << endl;

   return afield;
}
