/*
 * site_wirtePAnalysis.cpp
 *
 *  Created on: 01-04, 2015
 *      Author: Lin Li
 */

#include "site.h"

void CSite::wirtePAnalysis()
{
    cout << "#######in writePAnalysis ######" << endl;
    //LinLi:2013 March: for analyzing the average potential against z-axis.
    int i,j,k,m,xbox,ybox,zbox,lbox;
    float rmid,distance;
    delphi_real psumarr[10000];
    delphi_integer npoint[10000];
    //---------------------------------------------------------------------------
    rmid=float((iGrid+1)/2);
    distance=(radipz*fScale)*(radipz*fScale);

    for(i=0; i<10000; i++)
    {
        psumarr[i]=0.;
        npoint[i]=0;
    }
    lbox=int(radipz*fScale);

    cout << "LinLi,radipz,scale,distance,lbox: " << radipz<< " " << fScale<< " " << distance<< " " << lbox << endl;

    ofstream pzfile;
    pzfile.open ("pz.txt");


    //open(111,file="pz.txt",form="formatted")
    for(k = 1; k<=iGrid; k++)
    {
        for(i = 1; i<=iGrid; i++)
        {
            for(j = 1; j<=iGrid; j++)
            {
                psumarr[k]=psumarr[k]+prgfPhiMap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)];
                npoint[k]=npoint[k]+1;
                //k_j_i to i_j_k
                //cout <<"LinLi,i,j,k,phimap: " << i << " " << j << " " << k << " " << prgfPhiMap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)] << endl;
            }// end do
        }//end do
        //cout << k << " " << psumarr[k] << endl;
    }//end do


    for (m=0; m<iNatom; m++) // LinLi: 0 or 1??
    {
        //cout << xn2[m].nX << " " << xn2[m].nY << " " << xn2[m].nZ << endl;
        for (i=int(xn2[m].nX)-lbox; i<=int(xn2[m].nX)+lbox+1; i++)
        {
            for (j=int(xn2[m].nY)-lbox; j<=int(xn2[m].nY)+lbox+1; j++)
            {
                for (k=int(xn2[m].nZ)-lbox; k<=int(xn2[m].nZ)+lbox+1; k++)
                {

                    if((i-xn2[m].nX)*(i-xn2[m].nX)+(j-xn2[m].nY)*(j-xn2[m].nY)+(k-xn2[m].nZ)*(k-xn2[m].nZ) < distance
                    && (prgfPhiMap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)]>0.1
                    ||prgfPhiMap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)]<-0.00001))
                    {
                        psumarr[k]=psumarr[k]-prgfPhiMap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)];
                        npoint[k]=npoint[k]-1;
                        prgfPhiMap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)]=0.0;

                        //k_j_i to i_j_k
                        //cout << i << " " << j << " " << k << " " << prgfPhiMap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)] << endl;
                    }//endif
                }//end do
            }//end do
        }//end do
    }//end do

    for(k=1; k<=iGrid; k++)
    {
        //write(111,'(a3,i5,f8.3,a8,i10,a8,f12.4,a12,f12.4,a4)'),'z:'&
        //&,k,(k-rmid)/scale+oldmid.nZ,'n:',npoint(k),'Pz:'&
        //&,psumarr[k]/npoint(k),'kt/e  or',psumarr[k]/npoint(k)*25.85,'mv'

        pzfile  << "z: "
                << setw(5) << k << fixed
                << setw(8) << setprecision(3) << (k-rmid)/fScale+fgBoxCenter.nZ
                << "   n: "
                << setw(10) << npoint[k]
                << "   Pz: "
                << setw(12) << setprecision(4) << psumarr[k]/npoint[k]
                << " kt/e or  "
                << setw(12) << setprecision(4) << psumarr[k]/npoint[k]*25.85
                << " mv"
                << endl;
    }

    pzfile.close();
}
