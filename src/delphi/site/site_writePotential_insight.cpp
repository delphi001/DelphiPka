/*
 * site_writePotential_insight.cpp
 *
 *  Created on: Feb 26, 2014
 *      Author: chuan
 */

#include "site.h"

void CSite::writePotential_insight(vector<delphi_real>& phimap4)
{
   cout << "potential map written in INSIGHT format to file " << strPhiFile << endl << endl;

   int  ivary,nbyte,intdat,intx,inty,intz;
   delphi_real xang,yang,zang,xmax,range,extent;
   SGrid<delphi_real> xyzstart,xyzend;

   ivary = 0;         nbyte = 4;         intdat = 0;
   xang  = 90.0;      yang  = 90.0;      zang   = 90.0;
   intx  = iGrid - 1; inty  = iGrid - 1; intz   = iGrid - 1;

   xmax = optMax(fgBoxCenter); range = (iGrid-1.0)/(2.0*fScale); extent = range + xmax;
   xyzstart = (fgBoxCenter-range)/extent; xyzend = (fgBoxCenter-range)/extent;

   ofstream ofPhiStream(strPhiFile.c_str()); int formatflag = 1; //formatted output

   ofPhiStream << rgcFileMap << endl;
   ofPhiStream << ivary << " " << nbyte << " " << intdat << " " << extent << " " << extent << " " << extent << " "
               << xang << " " << yang << "" << zang << " " << xyzstart.nX << " " << xyzend.nX << " "
               << xyzstart.nY << " " << xyzend.nY << " " << xyzstart.nZ << " " << xyzend.nZ << " "
               << intx << " " << inty << " " << intz << endl;
   writePhiMap(formatflag,phimap4,ofPhiStream);

   ofPhiStream.close();
}
