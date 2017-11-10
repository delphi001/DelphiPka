/*
 * misc_interpl.cpp
 *
 *  Created on: Apr 14, 2013
 *      Author: chuan
 *
 * interpolates the potential at any point inside a cubical volume using the potential
 * values at the 8 vertices by means of a trilinear function:
 *
 *	   W = A1.XYZ + A2.XY + A3.XZ + A4.YZ + A5.X +A6.Y + A7.Z + A8
 *
 * where Ai coefficients are linear combinations of Wi at the cube corners
 */

#include "misc_interpl.h"

//-----------------------------------------------------------------------//
class CPointOutBox : public CWarning
{
   public:
      CPointOutBox(const SGrid<delphi_real> & gPoint, const delphi_integer& igMaxGrid)
      {
#ifdef VERBOSE
         cerr << "PAY ATTENSION, POINT OUT OF THE CUBE!! \n";
         cerr << "VALUES: " << gPoint.nX << " " << gPoint.nY << " " << gPoint.nZ << " GRID: " << igMaxGrid << endl;
         cerr << "\t VALUE AT THIS POINT IS SET = 0." << endl;
#endif		 
      }
};

//-----------------------------------------------------------------------//
delphi_real interpl(const delphi_integer& igMaxGrid, delphi_real *** prgfMap, const SGrid<delphi_real>& gPoint)
{
   delphi_real fInterpl = 0.0; // returned interpolation value

   SGrid<delphi_real> fgMaxGrid = {(delphi_real)(igMaxGrid),(delphi_real)(igMaxGrid),(delphi_real)(igMaxGrid)};


   /*
    * return 0.0 if outside grid
    */
    
   if (optORLT(gPoint,1.0) || optORGT(gPoint,fgMaxGrid))
   {
#ifndef PRIME
       
      CPointOutBox waring(gPoint,igMaxGrid);
       
#endif
      return 0.0;
   }

    
   /*
    * find lower left bottom grid point
    */
   SGrid<delphi_integer> iFloor = {(delphi_integer)floor(gPoint.nX),(delphi_integer)floor(gPoint.nY),(delphi_integer)floor(gPoint.nZ)}; //nx,ny,nz
   
   SGrid<delphi_integer> iCeiling = iFloor+1; //nx1,ny1,nz1
   if (iCeiling.nX > igMaxGrid) iCeiling.nX = iFloor.nX;
   if (iCeiling.nY > igMaxGrid) iCeiling.nY = iFloor.nY;
   if (iCeiling.nZ > igMaxGrid) iCeiling.nZ = iFloor.nZ;
   
   /*
    * calculate cube coordinates of point
    */
   SGrid<delphi_real> fFloor = {(delphi_real)iFloor.nX,(delphi_real)iFloor.nY,(delphi_real)iFloor.nZ};
   fFloor = gPoint - fFloor; // xgr,ygr,zgr

   //cout << "Lin Li test: fFllor: " << fFloor << endl;
   /*
    * calculate coefficients of trilinear function
    * notice that
    * 1. 3D map prgfMap is prgfMap[iz][iy][ix]
    * 2. C++ index starts with 0 (therefore - 1)
    */
   delphi_real a8 = prgfMap[iFloor.nZ-1][iFloor.nY-1][iFloor.nX-1],
        a7 = prgfMap[iCeiling.nZ-1][iFloor.nY-1][iFloor.nX-1]-a8,
        a6 = prgfMap[iFloor.nZ-1][iCeiling.nY-1][iFloor.nX-1]-a8,
        a5 = prgfMap[iFloor.nZ-1][iFloor.nY-1][iCeiling.nX-1]-a8,
        a4 = prgfMap[iCeiling.nZ-1][iCeiling.nY-1][iFloor.nX-1]-a8-a7-a6,
        a3 = prgfMap[iCeiling.nZ-1][iFloor.nY-1][iCeiling.nX-1]-a8-a7-a5,
        a2 = prgfMap[iFloor.nZ-1][iCeiling.nY-1][iCeiling.nX-1]-a8-a6-a5,
        a1 = prgfMap[iCeiling.nZ-1][iCeiling.nY-1][iCeiling.nX-1]-a8-a7-a6-a5-a4-a3-a2;
        
   // determine value of phi     
   fInterpl = a1*fFloor.nX*fFloor.nY*fFloor.nZ + a2*fFloor.nX*fFloor.nY + a3*fFloor.nX*fFloor.nZ + a4*fFloor.nY*fFloor.nZ +
              a5*fFloor.nX + a6*fFloor.nY + a7*fFloor.nZ + a8;

   return fInterpl;
}

//---cubic interpolation: --------------------------------------------------------------------//
delphi_real cubicInterpl(delphi_real p[4], delphi_real x)
{
//      for cubic interpolation 
//      Approach 1:(delivered by Lin Li)
//      a=p[3]/6-p[2]/2+p[1]/2-p[0]/6
//      b=p[0]/2+p[2]/2-p[1]
//      c=-p[3]/6+p[2]-p[1]/2-p[0]/3
//      d=p[1]
//
	//return x*x*x*(p[3]-p[0]+3*p[1]-3*p[2])/6.0 +x*x*(p[0]+p[2]-2*p[1])/2.0+ x*(6*p[2]-p[3]-3*p[1]-2*p[0])/6.0+p[1];
//
//
//      for cubic interpolation 
//      Approach 2:(delivered by Pauli http://www.paulinternet.nl/?page=bicubic)
//      I think it's not as accurate as my method --- Lin Li
//	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));

//      for quadratic interpolation 
//      Approach 1:(delivered by Lin Li)
//      a=p[3]/2-p[2]+p[1]/2
//      b=2p[2]-3/2*p[1]-p[3]/2
//      c=p[1]

        return x*x*(p[3]-2*p[2]+p[1])/2 + x*(4*p[2]-3*p[1]-p[3])/2 + p[1];
}

delphi_real bicubicInterpl (delphi_real p[4][4], delphi_real x, delphi_real y) 
{
	delphi_real arr[4];
	arr[0] = cubicInterpl(p[0], y);
	arr[1] = cubicInterpl(p[1], y);
	arr[2] = cubicInterpl(p[2], y);
	arr[3] = cubicInterpl(p[3], y);
	return cubicInterpl(arr, x);
}

delphi_real tricubicInterpl(const delphi_integer& igMaxGrid, delphi_real *** prgfMap, const SGrid<delphi_real>& gPoint)
{
   delphi_real fInterpl = 0.0; // returned interpolation value
   int i,j,k;
   //compared to trilinear interplation, tricubic interpolation need extra grids.
   SGrid<delphi_real> fgMaxGrid = {(delphi_real)(igMaxGrid-1),(delphi_real)(igMaxGrid-1),(delphi_real)(igMaxGrid-1)};

   /*
    * return 0.0 if outside grid
    */
    
   if (optORLT(gPoint,2.0) || optORGT(gPoint,fgMaxGrid))
   {
#ifndef PRIME
       
      CPointOutBox waring(gPoint,igMaxGrid);
       
#endif
      return 0.0;
   }

    
   /*
    * find lower left bottom grid point
    */
   SGrid<delphi_integer> iFloor = {(delphi_integer)floor(gPoint.nX),(delphi_integer)floor(gPoint.nY),(delphi_integer)floor(gPoint.nZ)}; //nx,ny,nz

   //seems like we can remove the code below:   
   SGrid<delphi_integer> iCeiling = iFloor+1; //nx1,ny1,nz1
   if (iCeiling.nX > igMaxGrid) iCeiling.nX = iFloor.nX;
   if (iCeiling.nY > igMaxGrid) iCeiling.nY = iFloor.nY;
   if (iCeiling.nZ > igMaxGrid) iCeiling.nZ = iFloor.nZ;
   
   /*
    * calculate cube coordinates of point
    */
   SGrid<delphi_real> fFloor = {(delphi_real)iFloor.nX,(delphi_real)iFloor.nY,(delphi_real)iFloor.nZ};
   fFloor = gPoint - fFloor; // xgr,ygr,zgr

   delphi_real p[4][4][4];
//   delphi_real p[4][4][4]={{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32}}; //p[z][y][x]
//   cout << "p[0][0][1],p[0][0][2] " << p[0][0][1] << p[0][0][2] << endl;
//   cout << "p[0][1][1],p[0][1][2] " << p[0][1][1] << p[0][1][2] << endl;
//   cout << "p[1][1][1],p[1][1][2] " << p[1][1][1] << p[1][1][2] << endl;
   for(i=iFloor.nX-1;i<=iFloor.nX+2;i++)
   {
      for(j=iFloor.nY-1;j<=iFloor.nY+2;j++)
      {
         for(k=iFloor.nZ-1;k<=iFloor.nZ+2;k++)
         {
             p[i-iFloor.nX+1][j-iFloor.nY+1][k-iFloor.nZ+1]=prgfMap[k-1][j-1][i-1];//notice that the x,y,z are reverted
             //p[i-iFloor.nX+1][j-iFloor.nY+1][k-iFloor.nZ+1]=prgfMap[iFloor.nZ-1][iFloor.nY-1][iFloor.nX-1];//notice that the x,y,z are reverted
             //p[i-iFloor.nZ+1][j-iFloor.nY+1][k-iFloor.nX+1]=prgfMap[iFloor.nZ-1][iFloor.nY-1][iFloor.nX-1];//notice that the x,y,z are reverted
         }
      }
   }
   delphi_real arr[4];
   arr[0] = bicubicInterpl(p[0], fFloor.nY, fFloor.nZ);
   arr[1] = bicubicInterpl(p[1], fFloor.nY, fFloor.nZ);
   arr[2] = bicubicInterpl(p[2], fFloor.nY, fFloor.nZ);
   arr[3] = bicubicInterpl(p[3], fFloor.nY, fFloor.nZ);
   /*cout << "p[0][0][0]: " << p[0][0][0] << endl; 
   cout << "p[0][0][1]: " << p[0][0][1] << endl; 
   cout << "p[0][0][2]: " << p[0][0][2] << endl; 
   cout << "p[0][0][3]: " << p[0][0][3] << endl; 
   cout << "p[0][1][0]: " << p[0][1][0] << endl; 
   cout << "p[0][1][1]: " << p[0][1][1] << endl; 
   cout << "p[0][1][2]: " << p[0][1][2] << endl; 
   cout << "p[0][1][3]: " << p[0][1][3] << endl; 
   cout << "arr[0]: " << arr[0] << endl; 
   cout << "arr[1]: " << arr[1] << endl; 
   cout << "arr[2]: " << arr[2] << endl; 
   cout << "arr[3]: " << arr[3] << endl; 
    */
   fInterpl=cubicInterpl(arr, fFloor.nX);
   //cout << "Lin Li: fFloor: " << fFloor << endl;

   return fInterpl;
}



//---bool interpl--------------------------------------------------------------------//
delphi_real boolinterpl(const delphi_integer& igMaxGrid, const vector<bool>& prgbMap, const SGrid<delphi_real>& gPoint)
{
   delphi_real fInterpl = 0.0; // return interploation value

   SGrid<delphi_real> fgMaxGrid = {(delphi_real)(igMaxGrid),(delphi_real)(igMaxGrid),(delphi_real)(igMaxGrid)};

   // return 0.0 if outside grid
   if (optORLT(gPoint,0.0) || optORGT(gPoint,fgMaxGrid))
   {
      CPointOutBox waring(gPoint,igMaxGrid);
      return fInterpl = 0.0;
   }

   // find lower left bottom grid point
   SGrid<delphi_integer> iFloor = {(delphi_integer)floor(gPoint.nX),(delphi_integer)floor(gPoint.nY),(delphi_integer)floor(gPoint.nZ)}; //nx,ny,nz

   SGrid<delphi_integer> iCeiling = iFloor+1; //nx1,ny1,nz1
   if (iCeiling.nX > igMaxGrid) iCeiling.nX = iFloor.nX;
   if (iCeiling.nY > igMaxGrid) iCeiling.nY = iFloor.nY;
   if (iCeiling.nZ > igMaxGrid) iCeiling.nZ = iFloor.nZ;

   // calculate cube coordinates of point
   SGrid<delphi_real> fFloor = {(delphi_real)iFloor.nX,(delphi_real)iFloor.nY,(delphi_real)iFloor.nZ};
   fFloor = gPoint - fFloor; // xgr,ygr,zgr

   // convert logical idebmap to delphi_real
   delphi_real debm[8];
   delphi_integer nx  = iFloor.nX,   ny  = iFloor.nY,   nz  = iFloor.nZ;
   delphi_integer nx1 = iCeiling.nX, ny1 = iCeiling.nY, nz1 = iCeiling.nZ;

   delphi_integer nw;

   nw = nz*igMaxGrid*igMaxGrid+ny*igMaxGrid+nx;
   if (prgbMap[nw]) debm[0] = 1.0;

   nw = nz1*igMaxGrid*igMaxGrid+ny*igMaxGrid+nx;
   if (prgbMap[nw]) debm[1] = 1.0;

   nw = nz*igMaxGrid*igMaxGrid+ny1*igMaxGrid+nx;
   if (prgbMap[nw]) debm[2] = 1.0;

   nw = nz*igMaxGrid*igMaxGrid+ny*igMaxGrid+nx1;
   if (prgbMap[nw]) debm[3] = 1.0;

   nw = nz1*igMaxGrid*igMaxGrid+ny1*igMaxGrid+nx;
   if (prgbMap[nw]) debm[4] = 1.0;

   nw = nz1*igMaxGrid*igMaxGrid+ny*igMaxGrid+nx1;
   if (prgbMap[nw]) debm[5] = 1.0;

   nw = nz*igMaxGrid*igMaxGrid+ny1*igMaxGrid+nx1;
   if (prgbMap[nw]) debm[6] = 1.0;

   nw = nz1*igMaxGrid*igMaxGrid+ny1*igMaxGrid+nx1;
   if (prgbMap[nw]) debm[7] = 1.0;

   // calculate coefficients of trilinear function
   delphi_real a8 = debm[0],
        a7 = debm[1]-a8,
        a6 = debm[2]-a8,
        a5 = debm[3]-a8,
        a4 = debm[4]-a8-a7-a6,
        a3 = debm[5]-a8-a7-a5,
        a2 = debm[6]-a8-a6-a5,
        a1 = debm[7]-a8-a7-a6-a5-a4-a3-a2;

   // determine value of phi
   fInterpl = a1*fFloor.nX*fFloor.nY*fFloor.nZ + a2*fFloor.nX*fFloor.nY + a3*fFloor.nX*fFloor.nZ + a4*fFloor.nY*fFloor.nZ +
              a5*fFloor.nX + a6*fFloor.nY + a7*fFloor.nZ + a8;

   return fInterpl;
}

