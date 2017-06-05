#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"


#include "../space/space.h"

using namespace std;
//#####################################################################
// ATTENTION! This file is part of epsmakmod module.
// Do not compile separately!
//#####################################################################

//---------------------------------------------------------------------
void CDelphiSpace::cubedata(delphi_real fac, delphi_real cbln)
{



//2011-05-27 Declarations added due to IMPLICIT NONE
    delphi_real off;
    SGrid <delphi_real> xyzp;
    off=0.1;
    //cout << fac << " " << cbln << endl;

    if(debug_space) cout << "### in cubedata: fac,cbln: " << fac << " " << cbln << endl;
    //fac=1.0; // fac=1.0 was passed as parameter, now removed passing process and set fac=1.0 below:
    xyzo=cMin-((fac*cbln)+off);
    xyzp=cMax+((fac*cbln)+off);
    //cout << "cMin << cMax << off: " << cMin << cMax << off << endl;

    lmncb=optCast<delphi_integer,delphi_real> ((xyzp-xyzo)/cbln);
    //cout << "xyzp,xyzo,cbln,lmncb: " << xyzp << xyzo << cbln << lmncb << endl;
    lcb=lmncb.nX;
    mcb=lmncb.nY;
    ncb=lmncb.nZ;
    if(debug_space) cout << "lmncb" << lmncb << endl;

    cbai=1./cbln;

}// void cubedata;



//---------------------------------------------------------------------
// 2011-05-27 Other parameters transfered via modules qlog and pointers
// rda changed to sDelPhiPDB[].radius

void CDelphiSpace::cube()
{

//here rda equals rad3
//2011-05-27 Arrays declared in pointers and allocated in calling
//void
    //delphi_integer cbn1[0:lcb,0:mcb,0:ncb],cbn2[0:lcb,0:mcb,0:ncb];
    //delphi_integer ***cbn1,***cbn2;


    //cbn1=get_pt3d<delphi_integer>(lcb+1,mcb+1,ncb+1);
    //cbn2=get_pt3d<delphi_integer>(lcb+1,mcb+1,ncb+1);

    get_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
    if(debug_space) cout << "### alocated cbn1 ###" << endl;
    get_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);

    delphi_integer newatm,itmp,kind;


//2011-05-27 Variables are accessible via qlog and pointers modules
// delphi_integer iNObject,numbmol
// delphi_integer cbal[1]

/**
 * creating a set of fictious atoms occupying little cubes
*/
    //SGrid <delphi_integer> icbn[iNatom+(iNObject-numbmol)*(lcb+1)*(mcb+1)*(ncb+1)];
    SGrid <delphi_integer> icbn[iNatom+(iNObject-numbmol)*(lcb+1)*(mcb+1)*(ncb+1)+1];

//iatmobj connects the fictious atoms to the objects
    delphi_integer iatmobj[(iNObject-numbmol)*(lcb+1)*(mcb+1)*(ncb+1)];

    string strtmp;

    SGrid <delphi_real> xq,xyz;
    SGrid <delphi_integer> ixyz;

    delphi_real shift;
    //delphi_real dist;
    delphi_real cost;

//2011-05-27 Declarations added due to IMPLICIT NONE
    delphi_real prbrd,cbln;
    delphi_integer i,j,k,ii,ix,iy,iz,jx,jy,jz,icum;


    if(debug_space) cout << "### in cube: ###" << endl;


    for(i=0;i<=iNatom+(iNObject-numbmol)*(lcb+1)*(mcb+1)*(ncb+1);i++){
        icbn[i]=sgrid_temp_int;
    }

    //cout << "Lin Li 1: icbn[818]:" << icbn[818] << endl;

    prbrd=fRadPrb[1];

    cbln=1./cbai;


//2011-05-27 Changed to array operations
    //cbn1=1;
    //cbn2=0;

    for (i=0;i<=lcb;i++){
        for (j=0;j<=mcb;j++){
            for (k=0;k<=ncb;k++){
                cbn1[i][j][k]=1;
                cbn2[i][j][k]=0;
            }
        }
    }


    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {

/**
 * xyzo variable is declared in pointers module, cbai in qlog
*/
            xyz=(xn1[i]-xyzo)*cbai;
            ixyz=optCast<delphi_integer,delphi_real>(xyz);
            if(debug_space) {


                if(optORLT(ixyz,1)) cout <<"ix << iy << iz: " << ixyz << endl;
                if(optORGE(ixyz,lmncb)) cout <<"ix << iy << iz: " << ixyz << endl;
            }
            ix=ixyz.nX;
            iy=ixyz.nY;
            iz=ixyz.nZ;

            for(jz=iz-1; jz<=iz+1; jz++)
            {
                for(jy=iy-1; jy<=iy+1; jy++)
                {
                    for(jx=ix-1; jx<=ix+1; jx++)
                    {
                        cbn2[jx][jy][jz]=cbn2[jx][jy][jz]+1;
                    }// do
                }// do
            }// do

            icbn[i]=ixyz;
            //cout << "i,icbn: "<< i << " " << icbn[i] << endl;
        }// if
    }// do

    if(debug_space) cout << "cbln: " << cbln << endl;
    newatm=iNatom;
    cost=cbln*.87+1./fScale;
    shift=cost+prbrd;


    if(debug_space) cout << "dataobject_v.size(): " << dataobject_v.size() << endl;
    for (i=0;i<=iNObject-1;i++){
        if(debug_space) cout << "i,dataobject_v[i]: " << i << " " <<  dataobject_v[i] << endl;

    }

   // cout << "Lin Li 2: icbn[818]:" << icbn[818] << endl;

/**
 * icbn will contain also coord center of fictious atoms
*/
    for(ii=1; ii<=iNObject; ii++)
    {
        //strtmp=dataobject[ii-1][0];

        strtmp=dataobject_v[(ii-1)*2];

        //read(strtmp(16:18),*)kind;
        //kind=strtmp.substr(15,3);
        kind = atoi(strtmp.substr(15,3).c_str());
        if (strtmp.substr(0,4)!="is a" && kind!=2)
        {
            if(debug_space) cout << "### Warning: strtmp.substr(0,4)!= is a && kind!=2 ###" << endl;
            itmp=ii+iNatom;

            for(iz=0; iz<=ncb; iz++)
            {
                for(iy=0; iy<=mcb; iy++)
                {
                    for(ix=0; ix<=lcb; ix++)
                    {
                        ixyz=int_coord(ix,iy,iz);
                        xq=((optCast<delphi_real,delphi_integer>(ixyz)+0.5)*cbln)+xyzo;


                        //if ((sLimObject[ii].nMin.vandle.(xq+shift))&& (sLimObject[ii].nMax.vandge.(xq-shift)))
                        if (optANDLE(sLimObject[ii].nMin,(xq+shift))&& optANDGE(sLimObject[ii].nMax,(xq-shift)))
                        {
                            //call distobj(xq,dist,dxyz,ii,prbrd,true);

                            //if (dist<=cost)
                            if (false) // never go to this if statement
                            {
                                newatm=newatm+1;
                                cbn2[ix][iy][iz]=cbn2[ix][iy][iz]+1;
                                icbn[newatm]=ixyz;
                                iatmobj[newatm-iNatom]=itmp;
                            }// if
                        }// if
                    }// do
                }// do
            }// do
        }// if
    }// do


    icum=1;
    for(iz=0; iz<=ncb; iz++)
    {
        for(iy=0; iy<=mcb; iy++)
        {
            for(ix=0; ix<=lcb; ix++)
            {
                if (cbn2[ix][iy][iz]>0)
                {
                    cbn1[ix][iy][iz]=icum;
                    icum=icum+cbn2[ix][iy][iz];
                    //cout << "iz,iy,ix,icum: " << iz << " " << iy<< " " << ix<< " " << icum << endl;
                }// if
            }// do
        }// do
    }// do


    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
        }// if
    }// do

    //cout << "Lin Li 1: cbal[3955] :" << cbal[3955] << endl;
    for(i=iNatom+1; i<=newatm; i++) //This will not be excuted
    {
        cout << "Warning: newatm > iNatom" << endl;
        ix=icbn[i].nX;
        iy=icbn[i].nY;
        iz=icbn[i].nZ;
        cbal[cbn1[ix][iy][iz]]=iatmobj[i-iNatom];
        cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
    }// do


    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            //cout << "i,ix,iy,iz,cbn1[ix][iy][iz],icbn[i]: " << " " << i << " " << ix << " " << iy << " " << iz << " " << cbn1[ix][iy][iz] << " " << icbn[i]<< endl;

        }// if
    }// do

//1,0,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            //cout << "i,ix,iy,iz,cbn1[ix][iy][iz],icbn[i]: " << " " << i << " " << ix << " " << iy << " " << iz << " " << cbn1[ix][iy][iz] << " " << icbn[i]<< endl;
        }// if
    }// do

//0,-1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            iy=iy-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            icbn[i].nY=iy;
        }// if
    }// do



//0,1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }// if
    }// do

//0,0,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy-1;
            iz=iz-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
            icbn[i].nZ=iz;
        }// if
    }// do

//0,0,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iz=iz+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nZ=iz;
        }// if
    }// do

//nn=2
//1,0,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do

//-1,0,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do

//0,1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+1;
            iy=iy+1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            icbn[i].nY=iy;
        }// if
    }// do
    //cout << "Lin Li 63: cbal[3955] :" << cbal[3955] << endl;
//0,-1,1
    for(i=1; i<=iNatom; i++)
    {

        if (sDelPhiPDB[i].radius>0.0)
        {

            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy-2;
            cbal[cbn1[ix][iy][iz]]=i;
            //if(cbn1[ix][iy][iz]==3955) cout << "Lin Li: cbn1[ix][iy][iz]:" << cbn1[ix][iy][iz] << " " << ix << " " << iy << " " << iz << endl;
            //if(cbn1[ix][iy][iz]==3955) cout << "Lin Li: icbn[i],i1:" << icbn[i] << " " << i << endl;
            //if(i==818) cout << "Lin Li: icbn[i],i:"                 << icbn[i] << " " << i << endl;

            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }// if
    }// do
    //cout << "Lin Li 65: cbal[3955] :" << cbal[3955] << endl;
//-1,-1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            iz=iz-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            icbn[i].nZ=iz;
        }// if
    }// do

//1,-1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do

//1,1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }// if
    }// do

//-1,1,0
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do

//-1,0,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iz=iz-1;
            iy=iy-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
            icbn[i].nZ=iz;
        }// if
    }// do

//1,0,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do

//0,1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            iy=iy+1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
            icbn[i].nY=iy;
        }// if
    }// do

//0,-1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }// if
    }// do

//nn=3
//-1,-1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-1;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do

//1,-1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do

//1,1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }// if
    }// do

//-1,1,-1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do

//-1,1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iz=iz+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nZ=iz;
        }// if
    }// do

//1,1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix+2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do

//1,-1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            iy=iy-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nY=iy;
        }// if
    }// do

//-1,-1,1
    for(i=1; i<=iNatom; i++)
    {
        if (sDelPhiPDB[i].radius>0.0)
        {
            ix=icbn[i].nX;
            iy=icbn[i].nY;
            iz=icbn[i].nZ;
            ix=ix-2;
            cbal[cbn1[ix][iy][iz]]=i;
            cbn1[ix][iy][iz]=cbn1[ix][iy][iz]+1;
            icbn[i].nX=ix;
        }// if
    }// do


//reset cbn1
    icum=1;
    for(iz=0; iz<=ncb; iz++)
    {
        for(iy=0; iy<=mcb; iy++)
        {
            for(ix=0; ix<=lcb; ix++)
            {
                if (cbn2[ix][iy][iz]>0)
                {
                    cbn1[ix][iy][iz]=icum;
                    icum=icum+cbn2[ix][iy][iz];
                    cbn2[ix][iy][iz]=icum-1;
                    //cout << "iz,iy,ix,icum,cbn2[ix][iy][iz]: " << iz << " " << iy << " " << ix << " " << icum << " " << cbn2[ix][iy][iz] << endl;
                }// if
            }// do
        }// do
    }// do


    icum=icum-1;
    if(debug_space) cout << "icum: " << icum << endl;


/*
    for (i=0;i<=lcb;i++){
        for (j=0;j<=mcb;j++){
            for (k=0;k<=ncb;k++){
                cbn1_v[i*mcb+j*ncb+k]=cbn1[i][j][k];
                cbn2_v[i*mcb+j*ncb+k]=cbn2[i][j][k];

            }
        }
    }
*/
}// void cube;


