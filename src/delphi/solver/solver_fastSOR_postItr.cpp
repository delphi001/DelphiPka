/*
 * solver_fastSOR_postItr.cpp
 *
 *  Created on: Feb 12, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

const string strEmpty60 = "                                                           ";
const string strLine60  = "-----------------------------------------------------------";

void CDelphiFastSOR::conplt(const vector<delphi_real>& array,const string& title,const int& iclr,const int& iscl,const int& imk,
                            const int& iplt,const char symb,const int& ixmin,const int& ixmax,vector<string>& iplot,delphi_real& ymin,delphi_real& ymax)
{
   delphi_integer iyup,iylw,iyran;
   delphi_real yminl,ymaxl;

   if (1 == iclr) { ymin = 1.0e10; ymax = 0.0; iplot.assign(nyran,strEmpty60); } //----- clear plot

   if (1 == iscl)
   {
      //----- scale plot and find max, min of data
      for (vector<delphi_real>::const_iterator it = array.cbegin(); it != array.cend(); ++it)
      {
         if (0.0 < *it)
         {
            ymin = min(ymin,*it);
            ymax = max(ymax,*it);
         }
      }

      //----- find y plot range in log scale
      yminl = log10(ymin); ymaxl = log10(ymax);
      iyup  = (delphi_integer)(1.0 + ymaxl); iylw = (delphi_integer)(yminl - 1.0);
      iyran = iyup - iylw;
   }

   if (1 == imk)
   {
      //----- make plot and stick x values in the appropriate bins after clipping
      yminl = log10(ymin); ymaxl = log10(ymax);
      iyup  = (delphi_integer)(1.0 + ymaxl); iylw = (delphi_integer)(yminl - 1.0);
      iyran = iyup - iylw;

      delphi_real temp,temp1;
      int ibin;
      for (vector<delphi_real>::const_iterator it = array.cbegin(); it != array.cend(); ++it)
      {
         if (ymin <= *it && *it <= ymax)
         {
            temp  = log10(*it);
            temp1 = (temp - iylw)/iyran;
            ibin  = (int)(temp1*(nyran - 1) + 1);
            iplot[ibin-1].assign(nxran,symb);
         }
      }
   }

   if (1 == iplt) //----- draw out plot
   {
      cout << endl;
      cout << title << endl;
      cout << " " << scientific << setprecision(2) << setw(9) << ymax << " |-" << strLine60 << "-| \n";

      cout << " " << scientific << setprecision(2) << setw(9) << ymin << " |-" << strLine60 << "-| \n";
      cout << "           " << " | " << strEmpty60 << " | \n";
      cout << "      " << setw(5) << ixmin << strEmpty60 << setw(5) << ixmax << endl;
   }
}


void CDelphiFastSOR::postItr(const vector<delphi_real>& rmaxl,const vector<delphi_real>& rmsl)
{
   /*
    * remap into phimap
    */
   for (delphi_integer iy = 0; iy < (iGrid*iGrid*iGrid-1)/2; iy++)
   {
      delphi_integer ix = iy*2;
      prgfPhiMap[ix]   = phimap1[iy];
      prgfPhiMap[ix+1] = phimap2[iy];
   }
   prgfPhiMap[iGrid*iGrid*iGrid-1] = phimap1[iHalfGridNum-1];

#ifdef VERBOSE
   cout << endl;
   cout << "finished qdiffx linear iterations "; pTimer->showTime(); cout << endl;
   cout << "total time elapsed so far: "; pTimer->showElapse(); cout << endl;
   //cout << "mean,max change (kT/e)   : " << rmsch2 << " " << rmxch2 << endl;
#endif

   /*
    * plot convergence history
    */
   if (bLogGraph)
   {
      int iclr = 1, iscl = 1, imk = 0, iplt = 0;
      char symb = 'M';
      string title = "    linear iteration convergence history   ";
      vector<string> iplot(nyran,strEmpty60);
      delphi_real ymin,ymax;

      conplt(rmaxl,title,iclr,iscl,imk,iplt,symb,1,iLinIterateNum,iplot,ymin,ymax);

      iclr = 0;
      conplt(rmsl,title,iclr,iscl,imk,iplt,symb,1,iLinIterateNum,iplot,ymin,ymax);

      iscl = 0; imk = 1;
      conplt(rmaxl,title,iclr,iscl,imk,iplt,symb,1,iLinIterateNum,iplot,ymin,ymax);

      symb = 'A'; iplt = 1;
      conplt(rmsl,title,iclr,iscl,imk,iplt,symb,1,iLinIterateNum,iplot,ymin,ymax);
   }

   /*
    * give some intermediate output of phi
    */
   if (bLogPotential)
   {
      int m,n,nn,ii;
      const delphi_real *** phimap = pdc->getKey_constPtr<delphi_real>("phimap",iGrid,iGrid,iGrid); // pointer to 3D phimap
      int midg = (iGrid+1)/2;
      for (m = 1; m <= 5; m++)
      {
         n = (iGrid-1)/4; nn = (m-1)*n+1;
         cout << "phi " << "  " << right << nn << "   " << right << midg << endl;
         for (ii = 1; ii <= iGrid; ii++)
            cout << right << phimap[ii-1][midg-1][nn-1] << "   ";
         cout << endl << endl;
      }
   }

}
