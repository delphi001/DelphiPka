/*
 * site_writePotential_grasp.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: chuan
 */

#include "site.h"

void CSite::writePotential_grasp(vector<delphi_real>& phimapIn)
{
   cout << "\nwriting potential map in GRASP format\n\n";
   cout << "potential map written to file " << strPhiFile << endl << endl;

   int  mgrid = 65;
   char nxtlbl[10],head[] = "now starting phimap ",tail[] = " end of phimap ";
   delphi_real scalesingle = fScale;
   SGrid<delphi_real> oldmidsingle = fgBoxCenter;

   expand(mgrid,phimapIn);

   ofstream ofPhiStream(strPhiFile.c_str(),ios::binary); int formatflag = 0; //unformatted

   if (bOutCrgDensity && (0 != fIonStrength))
      strcpy(nxtlbl,"concentrat");
   else
      strcpy(nxtlbl,"potential ");

   ofPhiStream.write(head,sizeof(head));
   ofPhiStream.write(nxtlbl,10);
   ofPhiStream.write(rgcFileMap,sizeof(rgcFileMap));
   writePhiMap(formatflag,phimapIn,ofPhiStream);
   ofPhiStream.write(tail,sizeof(tail));
   ofPhiStream.write(reinterpret_cast<char*>(&scalesingle),sizeof(delphi_real));
   ofPhiStream.write(reinterpret_cast<char*>(&oldmidsingle),sizeof(SGrid<delphi_real>));

   ofPhiStream.close();

   vector<delphi_real>().swap(phimapIn);
}
