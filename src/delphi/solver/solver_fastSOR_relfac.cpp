/*
 * solver_fastSOR_relfac.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"
#include <cmath>
delphi_real CDelphiFastSOR::calculateRelaxFactor()
{
   delphi_real fRelaxFactor = 0.0;
   delphi_integer iw;
   delphi_real temp,temp1,temp2,temp3,onesixth;
   int i,k,n,ix,iy,iz,icgrid,ihgd,ihgd2;
   int star,fin;
   int itemp1,itemp2,lat1,lat2;
   int long1,long2;

    delphi_integer& icount2a = iDielecBndyEven;
    delphi_integer& icount2b = iDielecBndyOdd;
    vector<delphi_integer>& idpos = prgiBndyDielecIndex;
    vector< vector<delphi_real> >& db = prgfBndyDielec;


   vector < int > sta1,sta2,fi1,fi2;


   sta1.assign(iGrid,0);
   sta2.assign(iGrid,0);
   fi1.assign(iGrid,0);
   fi2.assign(iGrid,0);

   //solver_pdc->showMap("test_delphicpp_1.dat");


   icgrid=pow(iGrid,3); ihgd=(iGrid+1)/2;
   //phimap1.assign(icgrid/2,0.0);
   //phimap2.assign(icgrid/2,0.0);
   phimap1.assign((icgrid+1)/2,0.0);
   phimap2.assign((icgrid+1)/2,0.0);


   //cout << "Lin relfac 1: phimap1.size(): " << phimap1.size() << endl;
   //cout << "Lin relfac 1: phimap2.size(): " << phimap2.size() << endl;
   phimap3.assign(iGrid*iGrid*iGrid,0.0);
   //prgfPhiMap.assign(iGrid*iGrid*iGrid,0.0);


   // b++++++++++++debug++++++++++++++++++++++++
   //
   // set up start and stop vectors
   onesixth=1./6.;
   sta1[1]=(iGrid*iGrid + iGrid + 4)/2; sta2[1]=sta1[1]-1;
   fi1[1]=pow(iGrid,2) - (iGrid+1)/2;   fi2[1]=fi1[1];
   itemp1=iGrid + 2;  itemp2=iGrid*iGrid -iGrid -2;
   for( i=2;i < iGrid-1;i++)
   {
      sta1[i]=fi1[i-1] + itemp1;  sta2[i]=fi2[i-1] + itemp1;
      fi1[i]=sta1[i-1] + itemp2;  fi2[i]=sta2[i-1] + itemp2;
   }

   //c also
   lat1= (iGrid-1)/2; lat2= (iGrid+1)/2;
   long1= (iGrid*iGrid - 1)/2; long2= (iGrid*iGrid + 1)/2;


   //----- setup sn array for lowest eigenstate
   vector<delphi_real> sn1(iGrid,0.0);

   for(ix = 1; ix < iGrid-1; ix++)
   {
      temp = fPi*ix/(iGrid-1);
      sn1[ix] = sqrt(2.0)*sin(temp)/sqrt((delphi_real)(iGrid-1));
      //cout << setprecision(3) << "ix: " << ix << " sn1: " << sn1[ix] << endl;
   }

   vector<delphi_real> sn2(sn1),sn3(sn1);

   delphi_real recipr = 1.0/sqrt((delphi_real)iGrid);
   if (rgbPeriodicBndy[0]) sn1.assign(iGrid,recipr);
   if (rgbPeriodicBndy[1]) sn2.assign(iGrid,recipr);
   if (rgbPeriodicBndy[0]) sn3.assign(iGrid,recipr);
/*
      //debug++++++++++++++ Lin Li
      ofstream ofMapStream1("test_delphicpp_fi_sta.txt");
      for (i = 0; i < iGrid; i++)
          {
             ofMapStream1 << setw(10) << right << i+1 << ": "
                         << setw(10) << sta1[i] << setw(10) << sta2[i] << setw(10) << fi1[i] << setw(10) << fi2[i] << endl;

          }
      ofMapStream1 << "lat1, lat2, long1, long2:" << right << setw(10) << lat1 << setw(10) << lat2 << setw(10) << long1 << setw(10) << long2 << endl;
      for (i = 0; i < iGrid; i++)
          {
             ofMapStream1 << "i,sn1,sn2,sn3:" << setw(10) << i+1 << fixed<< setprecision(7) << setw(20) << sn1[i]<< setw(20) << sn2[i] << setw(20) << sn3[i] << endl;

          }
      ofMapStream1 << "rionst:" << fixed<< setprecision(7) << setw(20) << fIonStrength << endl; 
      ofMapStream1 << "iper:" << fixed << setw(10) << rgbPeriodicBndy[0] << setw(10) << rgbPeriodicBndy[1] << setw(10) << rgbPeriodicBndy[2] << endl; 
      ofMapStream1 << "sixth:" << fixed<< setprecision(7) << setw(20) << onesixth << endl; 
      ofMapStream1 << "icount2a:"  << setw(20) << icount2a << endl; 
      ofMapStream1 << "icount2b:"  << setw(20) << icount2b << endl; 

      ofMapStream1 << "idpos:" << endl; 
      for(k=1;k<=icount2a;k++){
           ofMapStream1 << setw(10) << k << setw(10) << idpos[k-1] << endl; 
      }

      ofMapStream1 << "db:" << endl; 
      for(k=1;k<=icount2a;k++){
           //ofMapStream1 << setw(10) << k << setw(20) << setprecision(7) << fixed << db[0][k-1] << setw(20) << db[1][k-1] << setw(20) << db[2][k-1] << setw(20) << db[3][k-1] << setw(20) << db[4][k-1] << setw(20) << db[5][k-1] << setw(20) << endl; 
           ofMapStream1 << setw(10) << k << setw(20) << setprecision(7) << fixed << db[k-1][0] << setw(20) << db[k-1][1] << setw(20) << db[k-1][2] << setw(20) << db[k-1][3] << setw(20) << db[k-1][4] << setw(20) << db[k-1][5] << setw(20) << endl; 
 
      }

      ofMapStream1.close();
*/


   //----- map sn arrays to prgfPhiMap
   iw = 0;
   for (iz = 0; iz < iGrid; iz++)
   {
      temp3 = sn3[iz];
      for (iy = 0; iy < iGrid; iy++)
      {
         temp2 = temp3*sn2[iy];
         for (ix = 0; ix < iGrid; ix++)
         {
            prgfPhiMap[iw] = temp2*sn1[ix];
            phimap3[iw] = temp2*sn1[ix]; 

        //cout << setprecision(16)<< "iw: " << iw << " prgfPhiMap[iw]: " << prgfPhiMap[iw]  << endl; //LinLi test
            iw++;
         }
      }
   }

      temp=0.0;
      //for( ix=2; ix<=icgrid-1;ix+=2){
      for( ix=1; ix<=icgrid-2;ix+=2){
         iy=ix/2; 
         phimap2[iy] = phimap3[ix];
         temp = temp + phimap3[ix] * phimap3[ix];
         //cout << "iy,ix,phimap2,phimap3: " << setw(6) << fixed << setprecision(7) << iy << " " << setw(6) << ix << " " << phimap2[iy] << " " << phimap3[ix] << endl;
       }


/*
      ofstream ofMapStream3("test_delphicpp_phimap3.txt");
      ofMapStream3 << "            phimap3 : " << endl;
         for (unsigned int k = 0; k < iGrid; k++)
         {
            for (unsigned int j = 0; j < iGrid; j++)
            {
               for (unsigned int i = 0; i < iGrid; i++)
               {
                  ofMapStream3 << setw(3) << right << i << "," << setw(3) << right << j << "," << setw(3) << right << k << ": "
                              << setw(20) << fixed << setprecision(7)  << phimap3[k*iGrid*iGrid+j*iGrid+i] << endl;
               }
            }
         }
      ofMapStream3.close();
      ofstream ofMapStream2("test_delphicpp_phimap2.txt");
      ofMapStream2 << "            phimap2 : " << endl;
      for (unsigned int i = 0; i < icgrid/2; i++)
          {
             ofMapStream2 << setw(10) << right << i << ": " << setw(10) << fixed << setprecision(7)  << phimap2[i] << endl;
          }

      ofMapStream2.close();
*/
   //solver_pdc->showMap("test_delphicpp_2.dat");
   //----- setup periodic boundaries, start and stop vectors etc. for odd/even loops



if(fIonStrength>0.0) {
    for(n=1;n<=iGrid-2;n++){
        star=sta1[n];
        fin=fi1[n];
        for(ix=star-1;ix<=fin-1;ix++){
            temp1 = phimap2[ix] + phimap2[ix-1];
            temp2 = phimap2[ix+lat1] + phimap2[ix-lat2];
            temp3 = phimap2[ix+long1] + phimap2[ix-long2];
            phimap1[ix] = (temp1+temp2+temp3)*prgfSaltMap1[ix];
            //cout << ix << phimap1[ix] << prgfSaltMap1[ix] << endl;
        }// do
    }// do
    //!c otherwise the main loop is as below:
}
else{
    for(n=1;n<=iGrid-2;n++){
        star=sta1[n];
        fin=fi1[n];
        for(ix=star-1;ix<=fin-1;ix++){
            temp1 = phimap2[ix] + phimap2[ix-1];
            temp2 = phimap2[ix+lat1] + phimap2[ix-lat2];
            temp3 = phimap2[ix+long1] + phimap2[ix-long2];
            phimap1[ix] = (temp1+temp2+temp3)*onesixth;
        }// do
    }// do
}// if

 
/*
// ++++++++++++++++debug Lin Li ++++++++++++
      ofstream ofMapStream4("test_delphicpp_phimap2_2.txt");
      ofMapStream4 << "            phimap2 : " << endl;
      for (i = 0; i < icgrid/2; i++)
          {
             ofMapStream4 << setw(10) << right << i << ": " << setw(10) << fixed << setprecision(7)  << phimap2[i] << endl;
          }

      ofMapStream4.close();

// ++++++++++++++++debug Lin Li ++++++++++++
      ofstream ofMapStream5("test_delphicpp_phimap1_1.txt");
      ofMapStream5 << "            phimap1 : " << endl;
      for (i = 0; i < icgrid/2; i++)
          {
             ofMapStream5 << setw(10) << right << i << ": " << setw(10) << fixed << setprecision(7)  << phimap1[i] << endl;
          }

      ofMapStream5.close();
*/



//!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (rgbPeriodicBndy[0]) {
/*
    //!c calculating first slice
    ix=1+iGrid*lat2;
    if(fIonStrength>0.0) {
        for(n=1;n<=(iGrid-3)*lat2+1;n++){
            temp1 = phimap2[ix]+phimap2[ix-1+lat1];
            temp2 = phimap2[ix+lat1]+phimap2[ix-lat2];
            rtemp2=phimap2[ix+lat1+lat1]+phimap2[ix-lat2+lat1];
            temp3 = phimap2[ix+long1]+phimap2[ix-long2];
            rtemp3=phimap2[ix+long1+lat1]+phimap2[ix-long2+lat1];
            phimap1[ix]=(temp1+.5*(temp2+temp3+rtemp2+rtemp3))*(prgfSaltMap1[ix]+prgfSaltMap1[ix+lat1])*.5;
            //c now updating last slice
            phimap1[ix+lat1]=phimap1[ix];
            ix=ix+iGrid;
        }// do
    }
    else{
        for(n=1;n<=(iGrid-3)*lat2+1;n++){
            temp1 = phimap2[ix]+phimap2[ix-1+lat1];
            temp2 = phimap2[ix+lat1] + phimap2[ix-lat2];
            temp3 = phimap2[ix+long1] + phimap2[ix-long2];
            phimap1[ix]=(temp1+temp2+temp3)*sixth;
            //!c now updating last slice
            phimap1[ix+lat1]=phimap1[ix];
            ix=ix+iGrid;
        }// do
    }// if
*/
}// if
//!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!C$DIR NO_RECURRENCE
for(k=1;k<=icount2a;k++){
    ix=idpos[k-1];
    temp1=phimap2[ix-1-1]*db[k-1][0]+phimap2[ix-1]*db[k-1][1];
    temp2=phimap2[ix-lat2-1]*db[k-1][2]+phimap2[ix+lat1-1]*db[k-1][3];
    temp3=phimap2[ix-long2-1]*db[k-1][4]+phimap2[ix+long1-1]*db[k-1][5];
    phimap1[ix-1]= phimap1[ix-1] + temp1+temp2+temp3;
        //cout << "k,ix,...:" << setw(6) << k << setw(6) << ix << fixed << setprecision(7) << setw(11) << temp1 << setw(11) << temp2 << setw(11) << temp3 << endl;

}// do
//!c }// if


/*
// ++++++++++++++++debug Lin Li ++++++++++++
      ofstream ofMapStream6("test_delphicpp_phimap1_2.txt");
      ofMapStream6 << "            phimap1 : " << endl;
      for (i = 0; i < icgrid/2; i++)
          {
             ofMapStream6 << setw(10) << right << i << ": " << setw(10) << fixed << setprecision(7)  << phimap1[i] << endl;
          }

      ofMapStream6.close();

*/

//!c Now reset boundary values altered in above loops.
star=iGrid*(iGrid+1)/2;
fin=iGrid*(iGrid*(iGrid-1)-2)/2;
//!C$DIR NO_RECURRENCE
for(ix=star-1;ix<=fin-1;ix+=iGrid){
    phimap1[ix+1]=0.0;
    phimap1[ix+ihgd]=0.0;
}// do
/* temp is not used, so commented out-- Lin Li
temp=0.0;
for(ix=0;ix<=((icgrid-1)/2)-1 ;ix++){
    temp=temp + phimap1[ix]*phimap3[2*ix-1];
}// do
*/
//!c if periodic boundary condition option
//!c force periodicity using wrap around update of boundary values:
//!c 2nd slice-->last
//!c last-1 slice-->first
//!c z periodicity

/*
// ++++++++++++++++debug Lin Li ++++++++++++
      ofstream ofMapStream7("test_delphicpp_phimap1_3.txt");
      ofMapStream7 << "            phimap1 : " << endl;
      for (i = 0; i < icgrid/2; i++)
          {
             ofMapStream7 << setw(10) << right << i << ": " << setw(10) << fixed << setprecision(7)  << phimap1[i] << endl;
          }

      ofMapStream7.close();
*/


//========================================================================!
//!c Next update phimap3 using the new phimap1
if(fIonStrength>0.0) {
    for(n=1;n<=iGrid-2;n++){
        star=sta2[n];
        fin=fi2[n];
        for(ix=star-1;ix<=fin-1;ix++){
            temp1 = phimap1[ix] + phimap1[ix+1];
            temp2 = phimap1[ix+lat2] + phimap1[ix-lat1];
            temp3 = phimap1[ix+long2] + phimap1[ix-long1];
            phimap3[ix] =(temp1+temp2+temp3)*prgfSaltMap2[ix];
            //cout << "n,ix...:" << setw(6) << n << setw(6) << ix << setprecision(7) << fixed << setw(15) << temp1 << setw(15) << temp2 << setw(15) << temp3 << setw(15) << prgfSaltMap2[ix] << endl;
        }// do
    }// do
}
else{
    for(n=1;n<=iGrid-2;n++){
        star=sta2[n];
        fin=fi2[n];
        for(ix=star-1;ix<=fin-1;ix++){
            temp1 = phimap1[ix] + phimap1[ix+1];
            temp2 = phimap1[ix+lat2] + phimap1[ix-lat1];
            temp3 = phimap1[ix+long2] + phimap1[ix-long1];
            phimap3[ix] = (temp1+temp2+temp3)*onesixth;
        }// do
    }// do
}// if

/*
      ofstream ofMapStream8("test_delphicpp_phimap3_2.txt");
      ofMapStream8 << "            phimap3 : " << endl;
      for (i = 0; i < icgrid; i++)
          {
             ofMapStream8 << setw(10) << right << i << ": " << setw(10) << fixed << setprecision(7)  << phimap3[i] << endl;
          }

      ofMapStream8.close();
*/


//!c b+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (rgbPeriodicBndy[0]) {
/*
	//!c calculating first slice
	ix=long2;
	if(fIonStrength>0.0) {
		for(n=1;n<=long1-iGrid-1;n++){
			ix=ix+iGrid;
			temp1 = phimap1[ix+lat1] + phimap1[ix+1];
			temp2 = phimap1[ix+lat2] + phimap1[ix-lat1];
			rtemp2 = phimap1[ix+lat2+lat1] + phimap1[ix];
			temp3 = phimap1[ix+long2] + phimap1[ix-long1];
			rtemp3 = phimap1[ix+long2+lat1] + phimap1[ix-long1+lat1];
			phimap3[ix]=(temp1+.5*(temp2+temp3+rtemp2+rtemp3))*(prgfSaltMap2[ix]+prgfSaltMap2[ix+lat1])*.5;
			//!c now updating last slice
			phimap3[ix+lat1]=phimap3[ix];
		}// do
	}
	else{
		for(n=1;n<=long1-iGrid-1;n++){
			ix=ix+iGrid;
			temp1 = phimap1[ix+lat1] + phimap1[ix+1];
			temp2 = phimap1[ix+lat2] + phimap1[ix-lat1];
			temp3 = phimap1[ix+long2] + phimap1[ix-long1];
			phimap3[ix]=(temp1+temp2+temp3)*sixth;
			//!c now updating last slice
			phimap3[ix+lat1]=phimap3[ix];
		}// do
	}// if
*/
}// if
//!c e+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for(k=icount2a+1;k<=icount2b;k++){
	ix=idpos[k-1];
	temp1=phimap1[ix-1]*db[k-1][0]+phimap1[ix]*db[k-1][1];
	temp2=phimap1[ix-1-lat1]*db[k-1][2]+phimap1[ix-1+lat2]*db[k-1][3];
	temp3=phimap1[ix-1-long1]*db[k-1][4]+phimap1[ix-1+long2]*db[k-1][5];
	phimap3[ix-1]=phimap3[ix-1] + temp1+temp2+temp3;
        //cout << "k,ix,...:" << setw(6) << k << setw(6) << ix << fixed << setprecision(7) << setw(11) << temp1 << setw(11) << temp2 << setw(11) << temp3 << endl;
}// do
//!c reset boundary condition



star=(iGrid+2)/2;
iy=(iGrid*(iGrid+2)/2) - iGrid +1;
fin=(iGrid*(iGrid-1)-1)/2;
ihgd2=ihgd-1;
//!C$DIR NO_RECURRENCE
for(ix=star-1;ix<=fin-1;ix++){
	iy=iy+iGrid;
	phimap3[iy-1]=0.0;
	phimap3[iy-1+ihgd2]=0.0;
}// do

/*
      ofstream ofMapStream9("test_delphicpp_phimap3_3.txt");
      ofMapStream9 << "            phimap3 : " << endl;
      for (i = 0; i < icgrid; i++)
          {
             ofMapStream9 << setw(10) << right << i << ": " << setw(10) << fixed << setprecision(7)  << phimap3[i] << endl;
          }

      ofMapStream9.close();
*/

/*
   initOddEvenItr(0); // forWhom = 0

   //----- iterate once over odd and even points
   itrOddPoints(0,10001); // forWhom = 0

   itrEvenPoints(0,10002); // forWhom = 0
*/


   //solver_pdc->showMap("test_delphicpp_3.dat");
   //----- caculate the estimated spectral radius
   temp = 0.0;

   for (i = 0; i <(icgrid-1)/2; i++)
      temp += phimap3[i]*phimap2[i];
   fSpec = 2.0*temp;
   //----- following needed as spec exceeds 1.0 occasionally in focussing calculations (SS May 8, 1998)
   if (1.0 < fSpec) fSpec = 0.995;



#ifdef VERBOSE
   //cout << "\n gauss-seidel spectral radius is " << setw(15) << setprecision(10) << fSpec << endl;
   cout << "\n gauss-seidel spectral radius is " << setw(6) << setprecision(3) << fSpec << endl;
#endif

#ifdef DEBUG_DELPHI_SOLVER_RELFAC
   {
      string strTestFile = "test_relfac.dat";
      ofstream ofTestStream(strTestFile.c_str());
      ofTestStream << boolalpha;
      ofTestStream << fixed << setprecision(7);

      for (iw = 0; iw < iHalfGridNum; iw++)
         ofTestStream << "phimap1(" << setw(6) << right << iw+1 << ") = " << setw(11) << phimap1[iw] << endl;

      for (iw = 0; iw < iHalfGridNum; iw++)
         ofTestStream << "phimap2(" << setw(6) << right << iw+1 << ") = " << setw(11) << phimap2[iw] << endl;

      for (iw = 0; iw < iGrid*iGrid*iGrid; iw++)
         ofTestStream << "phimap3(" << setw(6) << right << iw+1 << ") = " << setw(11) << prgfPhiMap[iw] << endl;

   }
#endif // DEBUG_DELPHI_SOLVER
   //phimap1.clear();
   //phimap2.clear();
   //phimap3.clear();
   bndx1.clear();
   bndx2.clear();
   bndx3.clear();
   bndx4.clear();

   fRelaxFactor=fSpec;

   prgfPhiMap.assign(iGrid*iGrid*iGrid,0.0); // restore prgfPhiMap

   return fRelaxFactor;
}

