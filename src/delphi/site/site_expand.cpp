/*
 * site_expand.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: chuan
 *
 *  Description:
 *     expands igrid**3 grid to mgrid**3
 *     grid (65**3 default) to give compatibilty with previous phimap and epsmap formats using trilinear interpolation
 */

#include "site.h"

void CSite::expand(const int& mgrid, vector<delphi_real>& phimapIn)
{
   delphi_integer iw;
   SGrid<delphi_real> gc;
   delphi_real rscale,phiv;

   rscale = (iGrid-1.0)/(mgrid-1.0);

   /*
    * do high end first to prevent overwriting
    * find small grid values and interpolate into big grid
    */
   if (iGrid == mgrid)
   {
      phimapIn.assign(prgfPhiMap.begin(),prgfPhiMap.end());
   }
   else
   {
      if (0 == phimapIn.size()) phimapIn.assign(mgrid*mgrid*mgrid,0.0); // 65X65X65

      if (!bBiosystemOut)
      {
         for (int iz = mgrid; iz >0; iz--)
         {
            gc.nZ = (iz-1)*rscale+1.0;
            for (int iy = mgrid; iy >0; iy--)
            {
               gc.nY = (iy-1)*rscale+1.0;
               for (int ix = mgrid; ix >0; ix--)
               {
                  gc.nX = (ix-1)*rscale+1.0;
                  phiv = interpl(iGrid,phimap,gc);
                  iw = (iz-1)*mgrid*mgrid+(iy-1)*mgrid+(ix-1);
                  phimapIn[iw] = phiv;
               }
            }
         }
      }
   }

   fScale = fScale/rscale;
   cout << "new scale is " << fScale << " grids/ang" << endl;
}
