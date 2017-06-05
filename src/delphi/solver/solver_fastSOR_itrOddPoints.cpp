/*
 * solver_fastSOR_itrOddPoints.cpp
 *
 *  Created on: Feb 7, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

#ifdef PARALLEL_OMP
#include <omp.h>
#endif

void CDelphiFastSOR::itrOddPoints(const int& forWhom, const int& flag)
{
   delphi_integer n,ix,iy,iz;
   delphi_integer star,fin;
   delphi_real temp1,temp2,temp3,temp4;
   delphi_integer itemp1,itemp2,itemp3,itemp4;

    //cout << "### oddpoints phimap1: " << flag << endl;
#ifdef PARALLEL_OMP
   int omp_num_threads,omp_thread_id;

   /*
    * set number of threads = number of processors
    */
   //omp_set_num_threads(2);
   omp_set_num_threads(omp_get_num_procs());

   #pragma omp parallel default(shared) private(omp_thread_id,n,ix,iy,star,fin,temp1,temp2,temp3)
   {
      delphi_integer omp_index;

      omp_thread_id = omp_get_thread_num();

      if (0 == omp_thread_id) omp_num_threads = omp_get_num_threads();

      //cout << "thread " << omp_thread_id << " of " << omp_num_threads << " is alive\n";
#endif

      /* the following loops are about four times faster than the original loop over all grid points for
       * several reasons, the biggest being that we are only solving laplace's equation (unless salt is present),
       * which numerically much simpler, hence faster. we put all we leave out, back in below, ending up with
       * an equivalent calculation, but much faster.
       */
      if (fZero < abs(fIonStrength))  //----- the main loop is as below:
      {
#ifdef PARALLEL_OMP
         #pragma omp for schedule(auto)
#endif
         for (n = 1; n < iGrid-1; n++)
         {
            star = sta1[n]; fin = fi1[n];
            for (ix = star; ix <= fin; ix++)
            {
               temp1 = phimap2[ix-1]         + phimap2[(ix-1)-1];
               temp2 = phimap2[(ix-1)+lat1]  + phimap2[(ix-1)-lat2];
               temp3 = phimap2[(ix-1)+long1] + phimap2[(ix-1)-long2];
               phimap1[ix-1] = phimap1[ix-1]*om1 + (qmap1[ix-1]+temp1+temp2+temp3)*prgfSaltMap1[ix-1];
            }
         }
      }
      else //----- if there is no salt then the main loop is executed without sf saving about 15% in execution time
      {
#ifdef PARALLEL_OMP
         #pragma omp for schedule(auto)
#endif
         for (n = 1; n < iGrid-1; n++)
         {
            star = sta1[n]; fin = fi1[n];
            for (ix = star; ix <= fin; ix++)
            {
               temp1 = phimap2[ix-1]         + phimap2[(ix-1)-1];
               temp2 = phimap2[(ix-1)+lat1]  + phimap2[(ix-1)-lat2];
               temp3 = phimap2[(ix-1)+long1] + phimap2[(ix-1)-long2];
               phimap1[ix-1] = phimap1[ix-1]*om1 + (temp1+temp2+temp3)*sixth;
                //cout << "phimap1: " << right << setw(10) << flag << setw(10) << ix << setw(20) << setprecision(5) << fixed << phimap1[ix-1] << endl;
               //if(flag==1)cout << "1phimap1: " << right << setw(10) << flag << setw(10) << ix << setw(20) << setprecision(5) << fixed << phimap1[ix-1] << endl;

               //if( flag==2 && ix==498 )
               //cout << "phimap1: " << right << setw(10) << flag << setw(8) << ix << setw(20) << setprecision(5) << fixed << phimap1[ix-1]
               //     << " " << om1 << " " << temp1  << " " << temp2 << " " << temp3 << " " << sixth <<endl;

            }
         }
      }

#ifdef PARALLEL_OMP
      //#pragma omp barrier
#endif

      /*
       * first we add back the dielectric boundary points, by recalculating them individually. note this is still
       * vectorised by means of a gathering load by the compiler.
       */
#ifdef PARALLEL_OMP
      #pragma omp for schedule(auto)
#endif
      for (n = 0; n < iDielecBndyEven; n++)
      {
         ix = prgiBndyDielecIndex[n];
         temp1 = phimap2[(ix-1)-1]*prgfBndyDielec[n][0]     + phimap2[ix-1]*prgfBndyDielec[n][1];
         temp2 = phimap2[(ix-1)-lat2]*prgfBndyDielec[n][2]  + phimap2[(ix-1)+lat1]*prgfBndyDielec[n][3];
         temp3 = phimap2[(ix-1)-long2]*prgfBndyDielec[n][4] + phimap2[(ix-1)+long1]*prgfBndyDielec[n][5];
         phimap1[ix-1] += temp1 + temp2 + temp3;
        /*
        if(flag==1)cout << "2phimap1: " << right << setw(10) << flag << setw(10) << ix << setw(10) << setprecision(5) << fixed << phimap1[ix-1]
            <<setw(10) << phimap2[(ix-1)-long2] <<setw(10) << prgfBndyDielec[n][4]
            <<setw(10) << phimap2[(ix-1)+long1] <<setw(10) <<prgfBndyDielec[n][5]
            <<setw(10) << (ix-1)-long2 <<setw(10) << (ix-1)+long1
            << endl;
        */
        //if( flag==1 && ix==498 )
          //  cout << "phimap1: " << right << setw(10) << ix << setw(20) << setprecision(5) << fixed << phimap1[ix-1]
            //    << " " << temp1  << " " << temp2 << " " << temp3 <<endl;
      }


      /*
       * Now reset boundary values altered in above loops.
       */
#ifdef PARALLEL_OMP
      star = (iGrid+1)/2; fin = (iGrid*(iGrid-1)-2)/2; omp_index = iGrid*(iGrid+1)/2-iGrid+1;//iy = iGrid*(iGrid+1)/2-iGrid+1;
      #pragma omp for schedule(auto)
      for (n = 0; n < fin-star+1; n++)
      {
         iy = omp_index+(n+1)*iGrid;
         phimap1[iy-1] = bndx1[n];
         phimap1[iy+((iGrid+1)/2-1)-1] = bndx2[n];
      }
#else
      star = (iGrid+1)/2; fin = (iGrid*(iGrid-1)-2)/2; iy = iGrid*(iGrid+1)/2-iGrid+1;
      for (n = 0; n < fin-star+1; n++)
      {
         iy = iy+iGrid; phimap1[iy-1] = bndx1[n]; phimap1[iy+((iGrid+1)/2-1)-1] = bndx2[n];
      }
#endif

      /*
       * next we add back an adjustment to all the charged grid points due to the charge assigned. the compiler
       * directive just reassures the vector compiler that all is well as far as recurrence is concerned, i.e. it
       * would think there is a recurrence below, where as in fact there is none.
       */
      if (0 != forWhom)
      {
#ifdef PARALLEL_OMP
         #pragma omp for schedule(auto)
#endif
         for (n = 0; n < iCrgedGridEven; n++)
         {
            ix = prgiCrgPose[n]; phimap1[ix-1] += prgfCrgValA[n];
            //if(flag==1)cout << "3phimap1: " << right << setw(10) << flag << setw(10) << ix << setw(20) << setprecision(5) << fixed << phimap1[ix-1] << endl;

         }
      }

#ifdef PARALLEL_OMP
   } // end of #pragma omp parallel
#endif

   /*
    * if periodic boundary condition option, force periodicity using wrap around update of boundary values:
    *    2nd slice-->last
    *    last-1 slice-->first
    */
   if (rgbPeriodicBndy[2]) //----- z periodicity
   {
      for (iz = 0; iz < (iGrid-2)*(iGrid-2); iz += 2)
      {
         temp1 = ibndz[iz];      itemp1 = (delphi_integer)temp1;
         temp2 = temp1 + idif1z; itemp2 = (delphi_integer)temp2;
         temp3 = temp2 + inc1za; itemp3 = (delphi_integer)temp3;
         temp4 = temp1 + inc1zb; itemp4 = (delphi_integer)temp4;
         phimap1[itemp1-1] = phimap2[itemp2-1];
         phimap1[itemp3-1] = phimap2[itemp4-1];
      }
   }

   if (rgbPeriodicBndy[1]) //----- y periodicity
   {
      for (iy = 0; iy < (iGrid-2)*(iGrid-2); iy += 2)
      {
         temp1 = ibndy[iy];      itemp1 = (delphi_integer)temp1;
         temp2 = temp1 + idif1y; itemp2 = (delphi_integer)temp2;
         temp3 = temp2 + inc1ya; itemp3 = (delphi_integer)temp3;
         temp4 = temp1 + inc1yb; itemp4 = (delphi_integer)temp4;
         phimap1[itemp1-1] = phimap2[itemp2-1];
         phimap1[itemp3-1] = phimap2[itemp4-1];
      }
   }

   if (rgbPeriodicBndy[0]) //----- x periodicity
   {
      for (ix = 0; ix < (iGrid-2)*(iGrid-2); ix += 2)
      {
         temp1 = ibndx[ix];      itemp1 = (delphi_integer)temp1;
         temp2 = temp1 + idif1x; itemp2 = (delphi_integer)temp2;
         temp3 = temp2 + inc1xa; itemp3 = (delphi_integer)temp3;
         temp4 = temp1 + inc1xb; itemp4 = (delphi_integer)temp4;
         phimap1[itemp1-1] = phimap2[itemp2-1];
         phimap1[itemp3-1] = phimap2[itemp4-1];
      }
   }
}
