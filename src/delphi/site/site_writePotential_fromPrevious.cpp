/*
 * site_writePotential_fromPrevious.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: chuan
 */

#include "site.h"

void CSite::writePotential_fromPrevious(vector<delphi_real>& phimapIn)
{
   cout << "\nreading previous .phi file\n\n";

   char uplbl[20],nxtlbl[10],toplbl[53],botlbl[16];
   delphi_real fData,scalesingle1;
   SGrid<delphi_real> oldmidsingle1;
   delphi_integer igrid1,iw;

   string strPreviousFile = "previous.phi";
   ifstream ifPreviousStream;
   ifPreviousStream.open(strPreviousFile.c_str());
   if (!ifPreviousStream.is_open()) throw CUnknownPreviousPhiFile(strPreviousFile);

   ifPreviousStream.read(uplbl, 20); // uplbl
   ifPreviousStream.read(nxtlbl,10); // nxtlbl
   ifPreviousStream.read(toplbl,53); // toplbl

   //----- read previous phimap stored in unformatted file
   for (vector<delphi_real>::iterator it = phimapIn.begin(); it != phimapIn.end(); ++it)
   {
      ifPreviousStream.read(reinterpret_cast<char*>(&fData),sizeof(delphi_real));
      *it = fData;
   }

   ifPreviousStream.read(botlbl,16); // botlbl
   ifPreviousStream.read(reinterpret_cast<char*>(&scalesingle1),sizeof(delphi_real));
   ifPreviousStream.read(reinterpret_cast<char*>(&(oldmidsingle1.nX)),sizeof(delphi_real));
   ifPreviousStream.read(reinterpret_cast<char*>(&(oldmidsingle1.nY)),sizeof(delphi_real));
   ifPreviousStream.read(reinterpret_cast<char*>(&(oldmidsingle1.nZ)),sizeof(delphi_real));
   ifPreviousStream.read(reinterpret_cast<char*>(&igrid1),sizeof(delphi_integer));

   ifPreviousStream.close();

   if ( fZero < abs(fScale-scalesingle1) || fZero < abs(oldmidsingle1.nX-fgBoxCenter.nX) ||  fZero < abs(oldmidsingle1.nY-fgBoxCenter.nY) ||
        fZero < abs(oldmidsingle1.nZ-fgBoxCenter.nZ) ||  fZero < abs(iGrid-igrid1) ) throw CUnmatchPotentialMap();

   /*
    * subroutine submaps(mgrid,phimap,phimap4,opt) is removed due to the fact that realsiz cannot be 4
    */
   iw = 0;
   for (vector<delphi_real>::iterator it = phimapIn.begin(); it != phimapIn.end(); ++it)
   {
      *it = prgfPhiMap[iw] - *it;
      iw++;
   }

   cout << "potential map written to file " << strPhiFile << endl << endl;

   char head[] = "now starting phimap ",tail[] = " end of phimap ";
   delphi_real scalesingle = fScale;
   SGrid<delphi_real> oldmidsingle = fgBoxCenter;
   delphi_integer igrid = iGrid;

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
   ofPhiStream.write(reinterpret_cast<char*>(&igrid),sizeof(delphi_integer));

   ofPhiStream.close();
}
