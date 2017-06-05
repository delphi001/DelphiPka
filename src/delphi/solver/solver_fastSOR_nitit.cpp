/*
 * solver_fastSOR_nitit.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

void CDelphiFastSOR::nitit(const delphi_real& qfact)
{
   delphi_integer ix,iy,iz;
   vector<delphi_real> rmsl(nxran,0.0),rmaxl(nxran,0.0);
   bool ichangeom,istop,inewfirst,inew;
   delphi_real relparprev,factor,fraction,derprec,der;
   int itr,itnum,ires,icountplus;
   delphi_real rmsch,rmxch;
   delphi_real conv[3] = {0.0,0.0,0.0};
   string nlstr;

   if (0 == iConvergeFract) { iIterateInterval = 10; iConvergeFract = 1; }

   if (iLinIterateNum < iIterateInterval) iIterateInterval = iLinIterateNum;

   debmap1.assign(iHalfGridNum,0.0); debmap2.assign(iHalfGridNum,0.0);
   for (ix = 0; ix < iHalfGridNum; ix++)
   {
      iy = ix*2;
      if (prgbDielecMap[iy])   debmap1[ix] = 1.0;
      if (prgbDielecMap[iy+1]) debmap2[ix] = 1.0;
   }

   initOddEvenItr(2); // forWhom = 2

   cout << " linear rel. parameter = " << om2 << endl;

   if (bManualRelaxParam)
   {
      ichangeom = true;
      cout << " non linear fixed rel. parameter = " << fRelaxParam << endl;
   }
   else
   {
      ichangeom = false;
      cout << " non linear initial rel. parameter = " << fRelaxParam << endl;
      cout << " q factor " << qfact << endl;
   }

   cout << "\n\n  rms-change     max change       #iterations" << endl;

#ifdef DEBUG_DELPHI_SOLVER
   {
      string strTestFile = "test_nitit.dat";
      ofstream ofTestStream(strTestFile.c_str());
      ofTestStream << boolalpha;
      ofTestStream << fixed << setprecision(7);

      ix = 1;
      for (vector<delphi_real>::iterator it = debmap1.begin(); it != debmap1.end(); ++it)
      {
         ofTestStream << "debmap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
         ix++;
      }

      ix = 1;
      for (vector<delphi_real>::iterator it = debmap2.begin(); it != debmap2.end(); ++it)
      {
         ofTestStream << "debmap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
         ix++;
      }

      ix = 1;
      for (vector<delphi_real>::iterator it = bndx1.begin(); it != bndx1.end(); ++it)
      {
         ofTestStream << "bndx1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
         ix++;
      }

      ix = 1;
      for (vector<delphi_real>::iterator it = bndx2.begin(); it != bndx2.end(); ++it)
      {
         ofTestStream << "bndx2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
         ix++;
      }

      ix = 1;
      for (vector<delphi_real>::iterator it = bndx3.begin(); it != bndx3.end(); ++it)
      {
         ofTestStream << "bndx3[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
         ix++;
      }

      ix = 1;
      for (vector<delphi_real>::iterator it = bndx4.begin(); it != bndx4.end(); ++it)
      {
         ofTestStream << "bndx4[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
         ix++;
      }

      ix = 1;
      for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
      {
         ofTestStream << "phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
         ix++;
      }

      ix = 1;
      for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
      {
         ofTestStream << "phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
         ix++;
      }

      ofTestStream.close();
   }

   bool testflag = true;
#endif // DEBUG_DELPHI_SOLVER

   itr = 1; itnum = 0; ires = 0; istop = true; inewfirst = false; inew = true;
   factor = 1.0; nlstr = "                  "; icountplus = 0; relparprev = om2; fraction=0.0;

   do
   {
      rmsch = 0.0; rmxch = 0.0;

      /*
       * iterate over odd points
       */
      itrOddPoints(2,itr); // forWhom = 2

#ifdef DEBUG_DELPHI_SOLVER
      if (false)
      {
         testflag = false;

         string strTestFile = "test_nitit.dat";
         ofstream ofTestStream(strTestFile.c_str());
         ofTestStream << boolalpha;
         ofTestStream << fixed << setprecision(7);

         ix = 1;
         for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
         {
            ofTestStream << "phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
         }

         ix = 1;
         for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
         {
            ofTestStream << "phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
         }

         ofTestStream.close();
      }
#endif // DEBUG_DELPHI_SOLVER

      /*
       * Next update phimap2 using the new phimap1
       */
      itrEvenPoints(2,itr); // forWhom = 2

#ifdef DEBUG_DELPHI_SOLVER
      if (false)
      {
         testflag = false;

         string strTestFile = "test_nitit.dat";
         ofstream ofTestStream(strTestFile.c_str());
         ofTestStream << boolalpha;
         ofTestStream << fixed << setprecision(7);

         ix = 1;
         for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
         {
            ofTestStream << "phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
         }

         ix = 1;
         for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
         {
            ofTestStream << "phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
         }

         ofTestStream.close();
      }
#endif // DEBUG_DELPHI_SOLVER

      /*
       * we also save time by only checking convergence every 10 iterations, rather than every single iteration.
       * store phi2 in phi3 to compare against next iteration
       */
      if (iIterateInterval-1 == itr%iIterateInterval) // itr = 9,19,29,...
      {
         for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract) prgfPhiMap[ix] = phimap2[ix];
      }

      /*
       * check to see if accuracy is sufficient
       */
      if (0 == itr%iIterateInterval)
      {
         delphi_real rnorm2=0.0, temp2;

         for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract)
         {
            temp2   = prgfPhiMap[ix] - phimap2[ix];
            rnorm2 += temp2*temp2;
            rmxch   = max(rmxch,abs(temp2));
         }

         conv[2] = conv[1]; conv[1] = conv[0]; conv[0] = rmxch;

         rmsch = sqrt((delphi_real)iConvergeFract*rnorm2/((iGrid-2)*(iGrid-2)*(iGrid-2)));
         //rnormch = sqrt(rnorm2)/fRelaxParam;

         if ((fRmsc > rmsch || fMaxc > rmxch) && (22 < itnum)) ires = 1;

         if (0 == itnum) cout << scientific << rmsch << "  "  << rmxch << "  at  " << setw(5) << left << itr << " iterations\n";

         istop = !(0.22 < rmxch);

         if (bLogGraph)
         {
            int ibin;
            for (int j = itr-9; j <= itr; j++)
            {
               ibin = (j-1)*(60-1)/(iLinIterateNum-1)+1;
               rmsl[ibin-1]  = rmsch;
               rmaxl[ibin-1] = rmxch;
            }
         }

         //----- optimization of the first step
         inewfirst = (itr > iLinIterateNum-iLinIterateNum%iIterateInterval) && istop && (0 == itnum);

         if (!bManualRelaxParam && inewfirst && (3.2 < qfact))
         {
            factor    = exp(-qfact*2.1)+fZero;
            ichangeom = true;
         }

         //----- nonlinear part
         inew = inew || inewfirst;

         if (0 != itnum)
         {
            cout << scientific << rmsch << " " << rmxch << " " << itnum << " it. " << nlstr << endl;

            if (!bManualRelaxParam)
            {
               derprec = der; der = (conv[0]-conv[1])/conv[1];

               if (fZero > rmxch)
               {
                  factor = 1.2; ichangeom = true;
               }
               else
               {
                  if (0.0 < der && !inew)
                  {
                     icountplus++;
                     factor = factor*pow(1.0-der,0.99);
                     if (0.55 < der) { ichangeom = true; if (1.0 <= der) factor = 1.0e-5; }
                     if (0.35 < der && 0.1 < conv[0]) { ichangeom = true; factor = pow(factor*0.05/conv[0],4.0); }
                  }

                  if ((0.0 < der && 0.1 < conv[0]) && inew && !inewfirst)
                  {
                     ichangeom = true;
                     factor = min(factor*0.05/conv[0],factor*pow(1.0-der,0.86));
                  }

                  if (0.0 >= der)
                  {
                     icountplus = 0;
                     factor = 1.0;
                     if (24 < itnum && itnum < 24+0.75*(iNonIterateNum-24) && 0.03 > rmxch && 0.0 >= derprec)
                     {
                        if (-0.2 < der && -0.2 < derprec)
                        {
                           cout << "Trying to speed up the convergence process\n";
                           factor = 1.1;
                           ichangeom = true;
                           if (0.2 > fRelaxParam && -0.05 < der && -0.05 < derprec) factor = 1-45.226*(fRelaxParam-0.2);
                        }
                     }
                  }

                  if (2 <= icountplus) ichangeom = true;
               }

               inewfirst = false;
            }
         }
      } //----- end of convergence check

      itr++;

      if ((iLinIterateNum >= itr || !istop) && fZero > abs(ires) ) continue;

      if (0 < iNonIterateNum && 0 == itnum)
      {
         cout << "\nnow for the non-linear iterations\n";
         cout << "\n  rms-change     max change         #iterations\n";
      }

      iIterateInterval = 10; //----- icon1 = how many blocks each iteration convergence occurs
      iLinIterateNum   = 10; //----- nlit  = How many iterations in the block
      itnum++;

      if (iNonIterateNum < itnum || (1 == ires && !inew)) break;

      itr = 1; inew = ichangeom;
      if (ichangeom)
      {
         delphi_real omcomp;

         fRelaxParam = fRelaxParam*factor;

         if (1.0e-4 > fRelaxParam)
         {
            cout << "estimation " << fRelaxParam << " 1E-4 preferred\n";
            fRelaxParam = 1.0e-4;
         }

         factor = 1.0;
         cout << "                 New relaxation parameter = " << fRelaxParam << endl;
         ichangeom = false;
         icountplus = 0;
         omcomp = fRelaxParam/relparprev;
         relparprev = fRelaxParam;
         om1 = 1.0 - fRelaxParam;

         for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
            *it = (*it)*omcomp;

         for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
            *it = (*it)*omcomp;

         for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
            *it = (*it)*omcomp;

         for (delphi_integer iy = 0; iy < iDielecBndyOdd; iy++)
            for (delphi_integer ix = 0; ix < 6; ix++)
               prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix]*omcomp;

         sixth = sixth*omcomp;
      }

      fraction += 0.05;
      if (1.0 < fraction) {fraction = 1.0; nlstr = "full non-linearity";}

      {
         delphi_real temp1,temp2,temp3,temp4;
         delphi_real fac1 = fraction*fDebFct/(2.0*fIonStrength*fEpsOut);

         for (ix = 0; ix < iHalfGridNum; ix++)
         {
            temp1 = phimap1[ix]*debmap1[ix]; temp3 = temp1*temp1;
            temp2 = phimap2[ix]*debmap2[ix]; temp4 = temp2*temp2;
            qmap1[ix] = fac1*temp3*(fTaylorCoeff2+temp1*(fTaylorCoeff3+temp1*(fTaylorCoeff4+temp1*fTaylorCoeff5)));
            qmap2[ix] = fac1*temp4*(fTaylorCoeff2+temp2*(fTaylorCoeff3+temp2*(fTaylorCoeff4+temp2*fTaylorCoeff5)));

            //----- threshold setting to increase stability of the convergence
            if (2500.0 < temp3) qmap1[ix] = fac1*2500.0*(fTaylorCoeff2+50.0*(fTaylorCoeff3+50.0*(fTaylorCoeff4+50.0*fTaylorCoeff5)));
            if (2500.0 < temp4) qmap1[ix] = fac1*2500.0*(fTaylorCoeff2+50.0*(fTaylorCoeff3+50.0*(fTaylorCoeff4+50.0*fTaylorCoeff5)));
         }
      }
   } while(true);

   if (0.05 > fRelaxParam) {CSmallRelaxParam waring;}

   postItr(rmaxl,rmsl);

#ifdef DEBUG_DELPHI_SOLVER
   if (false)
   {
      string strTestFile = "test_nitit.dat";
      ofstream ofTestStream(strTestFile.c_str());
      ofTestStream << boolalpha;
      ofTestStream << fixed << setprecision(7);

      const delphi_real *** phimap = pdc->getKey_constPtr<delphi_real>("phimap",iGrid,iGrid,iGrid); // const pointer to 3D phimap
      for (iz = 0; iz < iGrid; iz++)
      {
         for (iy = 0; iy < iGrid; iy++)
         {
            for (ix = 0; ix < iGrid; ix++)
            {
               ofTestStream << "phimap[" << setw(6) << right << iz+1 << "," << setw(6) << right << iy+1 << ","
                            << setw(6) << right << ix+1 << "] = " << setw(11) << right << phimap[iz][iy][ix] << endl;
            }
         }
      }

      ofTestStream.close();
   }
#endif // DEBUG_DELPHI_SOLVER

}
