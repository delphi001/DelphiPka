/*
 * site_tops.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: chuan
 *
 *  Description:
 *     sc     = float = scale at which we need potential (removed)
 *     xo     = float(3) = position of charge
 *     crg    = float = value of charge
 *     phi    = float (0:64,0:64,0:64) = potential map for lookup table
 *     pot    = float = potential we need
 *     trgelf = float(3) = electric field at target point (used and removed)
 *     flag   = delphi_integer = 1 ==> calc. only potential
 *                      = 2 ==> calc. only field
 *                      = 3 ==> calc. both potential & field
 */

#include "site.h"

delphi_real CSite::tops(const SGrid<delphi_real>& xxo,const SGrid<delphi_real>& xxu,const delphi_real& crg,const delphi_real& eps,const int& flag)
{
   delphi_real pot = 0.0;

   delphi_real phi[65][65][65];
   delphi_real xo[3],axo[3],bxo[3],crgridpos[3][8],faxo[3],fbxo[3],mfo[8],crgrid[8];
   delphi_real xu[3],axu[3],bxu[3],trgridpos[3][8],faxu[3],fbxu[3],mfu[8];
   delphi_real v[8],e[3],dummy,elf[3][8];
   int  vec[3],si,a[3],b[3];
   int i,j,k,l;

   xo[0] = xxo.nX; xo[1] = xxo.nY; xo[2] = xxo.nZ;
   xu[0] = xxu.nX; xu[1] = xxu.nY; xu[2] = xxu.nZ;

   string strFileName = "lkphi.dat";
   ifstream ifFileStream;
   ifFileStream.open(strFileName.c_str());
   if (!ifFileStream.is_open()) throw CUnknownGridEngFile(strFileName);

   cout << "reading data file for analytic grid energies \n";

   for (k = 0; k < 65; k++)
      for (j = 0; j < 65; j++)
         for (i = 0; i < 65; i++)
            ifFileStream.read(reinterpret_cast<char*>(&phi[i][j][k]),sizeof(delphi_real));

   ifFileStream.close();

   //---------------------------------------------------------------------------------------------//

   for (i = 0; i < 3; i++)
   {
      axo[i] = (delphi_real)floor(xo[i]);
      if (0.0 > xo[i]) bxo[i] = axo[i] - 1.0;
      else             bxo[i] = axo[i] + 1.0;
   }

   for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
         crgridpos[j][i] = axo[j];

   for (i = 4; i < 8; i++)
      crgridpos[2][i] = bxo[2];

   for (i = 0; i < 2; i++)
      crgridpos[i][6] = bxo[i];

   crgridpos[0][1] = bxo[0]; crgridpos[1][3] = bxo[1];
   crgridpos[0][2] = bxo[0]; crgridpos[1][2] = bxo[1];

   for (i = 4; i < 8; i++)
      for (j = 0; j < 2; j++)
         crgridpos[j][i] = crgridpos[j][i-4];

   for (i = 0; i < 3; i++)
   {
      faxo[i] = abs(xo[i]-axo[i]); fbxo[i] = abs(xo[i]-bxo[i]);
   }

   mfo[0] = fbxo[0]*fbxo[1]*fbxo[2]; mfo[1] = faxo[0]*fbxo[1]*fbxo[2]; mfo[2] = faxo[0]*faxo[1]*fbxo[2];
   mfo[3] = fbxo[0]*faxo[1]*fbxo[2]; mfo[4] = fbxo[0]*fbxo[1]*faxo[2]; mfo[5] = faxo[0]*fbxo[1]*faxo[2];
   mfo[6] = faxo[0]*faxo[1]*faxo[2]; mfo[7] = fbxo[0]*faxo[1]*faxo[2];

   for (i = 0; i < 8; i++)
      crgrid[i] = crg*mfo[i];

   /*
    * crgrid(8) now contains the trilinearly interpolated charge at the eight vertices surrounding the charge
    */
   for (i = 0; i < 3; i++)
   {
      axu[i] = (delphi_real)floor(xu[i]);
      if (0.0 > xu[i]) bxu[i] = axu[i] - 1.0;
      else             bxu[i] = axu[i] + 1.0;
   }

   for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
         trgridpos[j][i] = axu[j];

   for (i = 4; i < 8; i++)
      trgridpos[2][i] = bxu[2];

   for (i = 0; i < 2; i++)
      trgridpos[i][6] = bxu[i];

   trgridpos[0][1] = bxu[0]; trgridpos[1][3] = bxu[1];
   trgridpos[0][2] = bxu[0]; trgridpos[1][2] = bxu[1];

   for (i = 4; i < 8; i++)
      for (j = 0; j < 2; j++)
         trgridpos[j][i] = trgridpos[j][i-4];

   for (i = 0; i < 3; i++)
   {
      faxu[i] = abs(xu[i]-axu[i]); fbxu[i] = abs(xu[i]-bxu[i]);
   }

   mfu[0] = fbxu[0]*fbxu[1]*fbxu[2]; mfu[1] = faxu[0]*fbxu[1]*fbxu[2]; mfu[2] = faxu[0]*faxu[1]*fbxu[2];
   mfu[3] = fbxu[0]*faxu[1]*fbxu[2]; mfu[4] = fbxu[0]*fbxu[1]*faxu[2]; mfu[5] = faxu[0]*fbxu[1]*faxu[2];
   mfu[6] = faxu[0]*faxu[1]*faxu[2]; mfu[7] = fbxu[0]*faxu[1]*faxu[2];

   /*
    * mfu contains the trilinear interpolation proportionality constants for interpolating the potential
    * onto the target point once we know the potentials of the neighboring eight vertices to the target point
    */

   if (1 == flag || 3 == flag) { for (i = 0; i < 8; i++) v[i] = 0.0;}

   for (i = 0; i < 8; i++)
   {
      if (1 < flag) { for (j = 0; j < 3; j++) e[j] = 0.0; }

      for (j = 0; j < 8; j++)
      {
         if (fZero < abs(crgrid[j]))
         {
            for (k = 0; k < 3; k++) vec[k] = trgridpos[k][i] - crgridpos[k][j];

            if (2 != flag)
            {
               dummy = phi[(int)abs(vec[0])][(int)abs(vec[1])][(int)abs(vec[2])];
               dummy = dummy*crgrid[j];
               v[i] += dummy/eps*fScale;
            }

            if (1 < flag)
            {
               for (k = 0; k < 3; k++)
               {
                  si = 1;
                  if (0 > vec[k]) si = -1;
                  for (l = 0; l < 3; l++) { a[l] = (int)abs(vec[l]); b[l] = (int)abs(vec[l]); }
                  if (64 > a[k]) a[k] += 1;
                  if (0  < b[k]) b[k] -= 1;
                  else           b[k] += 1;
                  dummy = (phi[a[0]][a[1]][a[2]]-phi[b[0]][b[1]][b[2]])/eps*fScale;
                  e[k] -= si*crgrid[j]*dummy;
               }
            }
         }
      }

      if (1 < flag)
      {
         for (j = 0; j < 3; j++) elf[j][i] = e[j];
      }
   }

   pot = 0.0;
   for (i = 0; i < 3; i++) pot += v[i]*mfu[i];

   return pot;
}
