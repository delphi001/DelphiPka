//#define DEBUG_DELPHI_SPACE

#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"

#include "../space/space.h"
#ifdef PARALLEL_OMP
#include <omp.h>
#endif

using namespace std;
//#####################################################################
//ATTENTION! This file is part of epsmakmod module.
// Do not compile separately!
//
// 2011-05-12 Parameters transfered via modules
//#####################################################################

void CDelphiSpace::setout()
{

//2011-05-12 Some local variables

    SGrid <delphi_real>  sqtemp[31],rad2aavtemp[31];
    SGrid <delphi_real> * sq=sqtemp+15;
    SGrid <delphi_real> * rad2aav=rad2aavtemp+15;
    SGrid <delphi_integer>* ioff;   ioff = NULL;

    bool itobig,itest2;
    string strtmp,strtmp1;
//2011-05-12 Non-standard int variable, thus necessary
    delphi_integer epsdim;
//2011-05-12 Non-standard float variable, thus necessary
    //delphi_real modul,modul2, mod2,modx,mody;
//2011-05-12 Non-standard type of variable, thus necessary
    //SGrid <delphi_real> xa,xb,xc,xd,xp;
    SGrid <delphi_real> ddist,dxyz,ddxyz,xn;
    //SGrid <delphi_real> tmpvect1,tmpvect2;
    SGrid <delphi_real> rad2av,fxn,vtemp;
    SGrid <delphi_integer> ismin,ismax,idist,ixyz,itest,ixn,i123;
/**
 * here fRadPrb is not zero only if one wants to map the ext}//ed
 *surf. that has to be done with fRadPrb(1) if solute-solvent
 *interface is concerned imedia = medium number for a object
 *a non-zero entry in iEpsMap indicates an atom # plus 1 (to
 *properly treat re-entrant mid-points later (15 Aug 93)
*/
//2011-05-27 Declarations added due to IMPLICIT NONE
    //delphi_integer iboxt;
    delphi_integer iac,ibox,igrdc,i,iv,ix,iy,iz;
    delphi_integer limmax,lim;
    delphi_real dis2min2,dis2min1,distsq,dist,rad,rad2;
    //delphi_real rad4,rad5;
    delphi_real rad2a,radmax2,fRMid,radp2,radtest,radprobe;
    delphi_real radp,temp;


    if(debug_space) cout << "############ beginning of setout: ##############" << endl;

#ifdef VERBOSE
    cout << "Starting creating Van der Waals Epsilon Map: " << endl;
#endif
    radprobe=0; // this radprobe seems not useful.
    epsdim=iNatom+iNObject+2;
    //iboxt=0;
    radmax2=0.0;
    fRMid=float((iGrid+1)/2);
    itest2=false;


#ifdef DEBUG_DELPHI_SPACE
    cout << "pdb radius 0 & iNatom-1: " << sDelPhiPDB[0].radius << " " << sDelPhiPDB[iNatom-1].radius << endl;
#endif // DEBUG_DELPHI_SPACE

    if (iNatom>0)
    {
//2011-05-12 Label 691 removed
        for(ix=1; ix<=iNatom; ix++)
        {
//2011-05-13 Chaged to derive-type array sDelPhiPDB (pointers module)
            radmax2=max(radmax2,sDelPhiPDB[ix].radius);

        }// do

#ifdef DEBUG_DELPHI_SPACE
        cout << "radmax2: " << radmax2 << endl;
#endif // DEBUG_DELPHI_SPACE

/**
 * this is probably the best way to do it,dep}//ing on which
 *surf. is desired
*/
        temp=max(radprobe,fExternRadius);

#ifdef DEBUG_DELPHI_SPACE
        cout << "radprobe: " << radprobe << " fExternRadius: " << fExternRadius << endl;
#endif // DEBUG_DELPHI_SPACE


        radmax2=fScale*(radmax2+temp);

        lim=1+radmax2;
        limmax = 12;
        itobig=false;
        if(lim>limmax) itobig=true;
        igrdc=pow((2*lim+1),3);
        ioff = new SGrid <int> [igrdc];

        if (!itobig)
        {
//2011-05-12 removed goto 7878 statement
            radtest= pow( (radmax2 + 0.5*sqrt(3.0)),2 );
            ibox=-1;

//2011-05-12 Strange statement. May allocate or may not
//allocate array that used later in the program
//irrespectively of itobig value, thus moved array
//allocation before if condition


            if(debug_space) cout << "radmax2: " << radmax2 << " radtest: " << radtest << " lim: " << lim << endl;
#ifdef PARALLEL_OMP
        #pragma omp for schedule(auto)
#endif

            for(ix=-lim; ix<=lim; ix++)
            {
                for(iy=-lim; iy<=lim; iy++)
                {
                    for(iz=-lim; iz<=lim; iz++)
                    {
//2011-05-13 Replaced by faster operation
//2011-05-13 Using operations on coord and
//int_coord type variables defined in module
//operators_on_coordinates
                        idist=int_coord(ix,iy,iz);
                        dist=float ( optDot(idist,idist) );
                        ddist = dist + float(0.25) + optCast <delphi_real,delphi_integer> (idist);

                        if ((dist<radtest)|| optORLT( ddist,radtest ))
                        {
                            ibox++;
                            ioff[ibox]=idist;
                            //cout << ibox <<": " << ioff[ibox] << idist << endl;

                        }// if
                    }// do
                }// do
            }// do
        }// if
    }// if
    if(debug_space) cout << ibox <<": " << ioff[ibox] << idist << endl;

    /*

    //set interiors in OBJECTS
        for(ii=1; ii<=iNObject; ii++)
        {
            ipore=false;
            strtmp=dataobject[ii][1];
            cout <<" " << endl;

            if (strtmp(1:4)!='is a')
            {
                strtmp1=dataobject[ii][2];
                read(strtmp(8:10),'(I3)')objecttype;
                read(strtmp(12:14),'(I3)')imedia;

    //completing iAtomMed with imedia data
                iAtomMed[iNatom+ii]=imedia;
                if(bVerbose)cout << "(setout) object number" << ii << " in medium" << imedia << endl;

    //check if object sits within limits of box and calculating
    //ismin and ismax accordingly
                exrd=fExternRadius*fScale;
                temp=radprobe*fScale+exrd;

    //2011-05-13 Using operations on coord and int_coord type
    //variables defined in module operators_on_coordinates
                ismin=optCast <delphi_integer,delphi_real> (sLimGridUnit[ii]%min-temp-1.);
                ismin=min(ismin,iGrid);
                ismin=max(ismin,1);

                ismax=optCast <delphi_integer,delphi_real> (sLimGridUnit[ii]%max+temp+1.);
                ismax=min(ismax,iGrid);
                ismax=max(ismax,1);

    //2011-05-13 Changed multiple IF's to switch
                switch (objecttype)
                {
                case 1 //if (objecttype==1) {
    //dealing with a sphere
                        read(strtmp(16:80),*)kind,xb,radius;
                    ipore=(kind==2);
                    radius=radius*fScale;
                    radp=radius+exrd;
                    rad2=radius*radius;
                    radp2=radp*radp;

    //2011-05-13 Using operations on coord and int_coord
    //type variables defined in module
    //operators_on_coordinates
                    xb=(xb-cOldMid)*fScale+fRMid;

                    for(ix=ismin.nX; ix<=ismax.nX; ix++)
                    {
                        for(iy=ismin.nY; iy<=ismax.nY; iy++)
                        {
                            for(iz=ismin.nZ; iz<=ismax.nZ; iz++)
                            {
                                ixyz=int_coord(ix,iy,iz);
                                xyz2=(optCast <delphi_real,delphi_integer> (ixyz)-xb)*(float(ixyz)-xb);
                                if(sum(xyz2)<=radp2) bDebMap[ix][iy][iz]=ipore;
                                xyz2.nX=xyz2.nX+0.25;

                                ddist=sum(xyz2)+optCast <delphi_real,delphi_integer> (ixyz)-xb;

                                if(ddist.nX<=rad2) iEpsMap[ix][iy][iz].nX=iNatom+1+ii+imedia*epsdim;
                                if(ddist.nY<=rad2) iEpsMap[ix][iy][iz].nY=iNatom+1+ii+imedia*epsdim;
                                if(ddist.nZ<=rad2) iEpsMap[ix][iy][iz].nZ=iNatom+1+ii+imedia*epsdim;
                            }// do
                        }// do
                    }// do
                    break;
                case 2 //if (objecttype==2) {
    //dealing with a cylinder
                        read(strtmp(16:80),*)kind,xa,xb,radius;
                    ipore=(kind==2);
                    read(strtmp1,'(5f8.3)')vectz,modul,modul2;

    //here we work in grid units
                    radius=radius*fScale;
                    rad2=radius*radius;
                    radp=radius+exrd;
                    radp2=radp*radp;
                    modul=modul*fScale;
                    modul2=modul*modul;
                    tmp=exrd*modul;

    //2011-05-13 Using operations on coord and int_coord
    //type variables defined in module
    //operators_on_coordinates
                    xb=(xb-cOldMid)*fScale + fRMid;
                    vectz=vectz*fScale;

                    for(ix=ismin.nX; ix<=ismax.nX; ix++)
                    {
                        for(iy=ismin.nY; iy<=ismax.nY; iy++)
                        {
                            for(iz=ismin.nZ; iz<=ismax.nZ; iz++)
                            {
    //vectz=A-B
                                modul=|A-B|;
                                tmpvect1=P-B;

    //tmp1=(A-B)(P-B)
    //mod2=(P-B)**2
                                modul2=(A-B)**2,;
    //tmp=fExternRadius*|A-B|.
    //2011-05-13 Using operations on coord and
    //int_coord type variables defined in module
    //operators_on_coordinates
                                ixyz=int_coord(ix,iy,iz);
                                tmpvect1=optCast <delphi_real,delphi_integer> (ixyz)-xb;
                                tmp1=tmpvect1.dot.vectz;

                                if((tmp1>=-tmp)&&(tmp1<=modul2+tmp))
                                {
                                    mod2=tmpvect1.dot.tmpvect1;
                                    if((mod2-(tmp1/modul)**2)<=radp2) bDebMap[ix][iy][iz]=ipore;
                                }// if

    //x-offset
                                tmpvect1.nX=tmpvect1.nX+.5;
                                tmp1=tmpvect1.dot.vectz;
                                if ((tmp1>=0.0)&&(tmp1<=modul2))
                                {
                                    mod2=tmpvect1.dot.tmpvect1;
                                    if ((mod2-(tmp1/modul)**2)<=rad2)
                                    {
                                        iEpsMap[ix][iy][iz].nX=iNatom+1+ii+imedia*epsdim;
                                    }// if
                                }// if

    //y-offset
                                tmpvect1.nX=tmpvect1.nX-.5;
                                tmpvect1.nY=tmpvect1.nY+.5;
                                tmp1=tmpvect1.dot.vectz;
                                if ((tmp1>=0.0)&&(tmp1<=modul2))
                                {
                                    mod2=tmpvect1.dot.tmpvect1;
                                    if ((mod2-(tmp1/modul)**2)<=rad2)
                                    {
                                        iEpsMap[ix][iy][iz].nY=iNatom+1+ii+imedia*epsdim;
                                    }// if
                                }// if

    //z-offset
                                tmpvect1.nY=tmpvect1.nY-.5;
                                tmpvect1.nZ=tmpvect1.nZ+.5;
                                tmp1=tmpvect1.dot.vectz;
                                if ((tmp1>=0.0)&&(tmp1<=modul2))
                                {
                                    mod2=tmpvect1.dot.tmpvect1;
                                    if ((mod2-(tmp1/modul)**2)<=rad2)
                                    {
                                        iEpsMap[ix][iy][iz].nZ=iNatom+1+ii+imedia*epsdim;
                                    }// if
                                }// if
                            }// do
                        }// do
                    }// do
                    break;
                case 3 //if (objecttype==3) {
    //dealing with a cone
                        read(strtmp(16:80),*)kind,xa,xb,alpha;
                    ipore=(kind==2);

    //conversion degrees --> radiants
                    alpha=alpha*3.14159265359/180.;
                    tan2=tan(alpha*.5)**2;
                    read(strtmp1,'(5f8.3)')vectz,modul,modul2;

    //here we work in grid units
                    modul=modul*fScale;
                    modul2=modul*modul;

    //2011-05-13 Using operations on coord and int_coord
    //type variables defined in module
    //operators_on_coordinates
                    xb=(xb-cOldMid)*fScale + fRMid;
                    vectz=vectz*fScale;

                    for(ix=ismin.nX; ix<=ismax.nX; ix++)
                    {
                        for(iy=ismin.nY; iy<=ismax.nY; iy++)
                        {
                            for(iz=ismin.nZ; iz<=ismax.nZ; iz++)
                            {
                                ixyz=int_coord(ix,iy,iz);
                                tmpvect1=(optCast <delphi_real,delphi_integer> (ixyz)-fRMid)/fScale+cOldMid;

    //2011-05-13 Parameters transfered via module
    //architecture
                                call distobj(tmpvect1,dist,vecdist,ii,fExternRadius,true);

                                if (dist<=0.0) bDebMap[ix][iy][iz]=ipore;

    //vectz=A-B
                                tmpvect1=P-B;
                                tmp1=(A-B)(P-B);
    //mod2=(P-B)**2
                                modul2=(A-B)**2.;

    //x-offset
    //2011-05-13 Using operations on coord and
    //int_coord type variables defined in module
    //operators_on_coordinates
                                tmpvect1=optCast <delphi_real,delphi_integer> (ixyz)-xb;
                                tmpvect1.nX=tmpvect1.nX+0.5;
                                tmp1=tmpvect1.dot.vectz;

                                if ((tmp1>=0.0)&&(tmp1<=modul2))
                                {
                                    mod2=(tmpvect1.dot.tmpvect1)*modul2;
                                    tmp2=(1+tan2)*tmp1*tmp1;
                                    if (mod2<=tmp2)
                                    {
                                        iEpsMap[ix][iy][iz].nX=iNatom+1+ii+imedia*epsdim;
                                    }// if
                                }// if

    //y-offset
                                tmpvect1.nX=tmpvect1.nX-0.5;
                                tmpvect1.nY=tmpvect1.nY+0.5;
                                tmp1=tmpvect1.dot.vectz;

                                if ((tmp1>=0.0)&&(tmp1<=modul2))
                                {
                                    mod2=(tmpvect1.dot.tmpvect1)*modul2;
                                    tmp2=(1+tan2)*tmp1*tmp1;
                                    if (mod2<=tmp2)
                                    {
                                        iEpsMap[ix][iy][iz].nY=iNatom+1+ii+imedia*epsdim;
                                    }// if
                                }// if

    //z-offset
                                tmpvect1.nY=tmpvect1.nY-0.5;
                                tmpvect1.nZ=tmpvect1.nZ+0.5;
                                tmp1=tmpvect1.dot.vectz;
                                if ((tmp1>=0.0)&&(tmp1<=modul2))
                                {
                                    mod2=(tmpvect1.dot.tmpvect1)*modul2;
                                    tmp2=(1+tan2)*tmp1*tmp1;
                                    if (mod2<=tmp2)
                                    {
                                        iEpsMap[ix][iy][iz].nZ=iNatom+1+ii+imedia*epsdim;
                                    }// if
                                }// if

                            }// do
                        }// do
                    }// do
                    break;
                case 4 //if (objecttype==4) {
    //dealing with a parallelepiped
                        read(strtmp(16:80),*)kind,xa,xb,xc,xd;
                    ipore=(kind==2);
                    read(strtmp1,'(12f8.3)')vectz,modul,vecty,mody,vectx,modx;

    //conversion to axial symmetry points
    //2011-05-13 Using operations on coord and int_coord
    //type variables defined in module
    //operators_on_coordinates
                    xb=0.5*(xc+xb);
                    xb=(xb-cOldMid)*fScale + fRMid;

                    modul=modul*fScale;
                    modul2=modul*modul;
                    vectz=vectz*fScale;

                    modx=modx*fScale;
                    tmp1=modx*modx/2.;
                    vectx=vectx*fScale;

                    mody=mody*fScale;
                    tmp2=mody*mody/2.;
                    vecty=vecty*fScale;

                    for(ix=ismin.nX; ix<=ismax.nX; ix++)
                    {
                        for(iy=ismin.nY; iy<=ismax.nY; iy++)
                        {
                            for(iz=ismin.nZ; iz<=ismax.nZ; iz++)
                            {
    //vectz=A-B
                                vectx=C-D;
                                vecty=(C+D)/2-B;

    //tmp1=|C-D|/2
    //modul2=(A-B)**2, tmp2=mody/2
    //dot=(P-B)(..-..)

                                ixyz=int_coord(ix,iy,iz);
                                xp=(optCast <delphi_real,delphi_integer> (ixyz)-fRMid)/fScale + cOldMid;

    //2011-05-13 Parameters transfered via module
    //architecture
                                call distobj(xp,dist,vecdist,ii,fExternRadius,true);

                                if (dist<=0.0) bDebMap[ix][iy][iz]=ipore;

    //now xp=P-B

    //x-offset
                                xp=optCast <delphi_real,delphi_integer> (ixyz)-xb;
                                xp.nX=xp.nX+0.5;
                                dot=vectz.dot.xp;

                                if ((dot>=0.0)&&(dot<=modul2))
                                {
                                    dot=vectx.dot.xp;
                                    if (abs(dot)<=tmp1)
                                    {
                                        dot=vecty.dot.xp;
                                        if (abs(dot)<=tmp2)
                                        {
                                            iEpsMap[ix][iy][iz].nX=iNatom+1+ii+imedia*epsdim;
                                        }// if
                                    }// if
                                }// if

    //y-offset
                                xp.nX=xp.nX-0.5;
                                xp.nY=xp.nY+0.5;
                                dot=vectz.dot.xp;

                                if ((dot>=0.0)&&(dot<=modul2))
                                {
                                    dot=vectx.dot.xp;
                                    if (abs(dot)<=tmp1)
                                    {
                                        dot=vecty.dot.xp;
                                        if (abs(dot)<=tmp2)
                                        {
                                            iEpsMap[ix][iy][iz].nY=iNatom+1+ii+imedia*epsdim;
                                        }// if
                                    }// if
                                }// if

    //z-offset
                                xp.nY=xp.nY-0.5;
                                xp.nZ=xp.nZ+0.5;
                                dot=vectz.dot.xp;

                                if ((dot>=0.0)&&(dot<=modul2))
                                {
                                    dot=vectx.dot.xp;
                                    if (abs(dot)<=tmp1)
                                    {
                                        dot=vecty.dot.xp;
                                        if (abs(dot)<=tmp2)
                                        {
                                            iEpsMap[ix][iy][iz].nZ=iNatom+1+ii+imedia*epsdim;
                                        }// if
                                    }// if
                                }// if

                            }// do
                        }// do
                    }// do
                }// select;
            }// if
        }// do //end of setting in OBJECTS;


    */
/**
 * set interiors in MOLECULES
*/
    if(itest2||itobig) cout <<"setout method 1 " << itest2 << " " << itobig << endl;

    //DoATOMS:
    for( iv=1; iv<=iNatom; iv++)
    {
/**
 * restore values
*/
        rad= sDelPhiPDB[iv].radius;

        xn=xn2[iv];
        //cout << "iv,xn1[iv],xn2[iv]: " <<iv<< xn1[iv+1] << xn2[iv+1] << endl;

//2011-05-13 Removed GOTO statement
        if (rad<1.e-6)
        {
            continue;
        }// if

//fScale radius to grid
        rad=rad*fScale;
        //rad5=pow( (rad+0.5),2);
        radp=rad+fExternRadius*fScale;
        rad=rad+radprobe*fScale;
        //rad4=pow( (rad+0.5),2);
        rad2=rad*rad;
        radp2=radp*radp;

/**
 * set dielectric map
 *check if sphere sits within limits of box
*/
        itest2=false;

        ismin=optCast <delphi_integer,delphi_real> (xn-radmax2-1.0);
        ismax=optCast <delphi_integer,delphi_real> (xn+radmax2+1.0);
        itest=ismin;
        ismin=optMin(ismin,iGrid);
        ismin=optMax(ismin,1);
        if(itest!=ismin) itest2=true;
        itest=ismax;
        ismax=optMin(ismax,iGrid);
        ismax=optMax(ismax,1);
        if(itest != ismax) itest2=true;

        if (itest2||itobig)   //slow method;
        {
//2011-05-13 Seems redundant statement
            if(debug_space) cout << "### slow method:"  << endl;
            if(debug_space) cout << "itest, itest2,itobig: " << itest << " " << itest2 << " " << itobig << endl;
            rad2a = rad2 - 0.25;
#ifdef KJI
#ifdef PARALLEL_OMP
        #pragma omp for schedule(auto)
#endif


            for(iz=ismin.nZ; iz<=ismax.nZ; iz++)
            {
                for(iy=ismin.nY; iy<=ismax.nY; iy++)
                {
                    for(ix=ismin.nX; ix<=ismax.nX; ix++)
                    {

                        ixyz=int_coord(ix,iy,iz);
                        dxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
                        distsq=optDot(dxyz,dxyz);
                        dxyz=dxyz+distsq;

                        if (dxyz.nX<rad2a)
                        {
                            iepsmp[iz][iy][ix].nX=iv+1+iAtomMed[iv]*epsdim;
                        }// if

                        if (dxyz.nY<rad2a)
                        {
                            iepsmp[iz][iy][ix].nY=iv+1+iAtomMed[iv]*epsdim;
                        }// if

                        if (dxyz.nZ<rad2a)
                        {
                            iepsmp[iz][iy][ix].nZ=iv+1+iAtomMed[iv]*epsdim;
                        }// if

                        if(distsq<radp2) idebmap[ix][iy][iz] =false;
                    }// do
                }// do
            }// do
#endif // KJI

            for(iz=ismin.nZ; iz<=ismax.nZ; iz++)
            {
                for(iy=ismin.nY; iy<=ismax.nY; iy++)
                {
                    for(ix=ismin.nX; ix<=ismax.nX; ix++)
                    {
                        ixyz=int_coord(ix,iy,iz);
                        dxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
                        distsq=optDot(dxyz,dxyz);
                        dxyz=dxyz+distsq;

                        if (dxyz.nX<rad2a)
                        {
                            iepsmp[iz][iy][ix].nX=iv+1+iAtomMed[iv]*epsdim;
                        }// if

                        if (dxyz.nY<rad2a)
                        {
                            iepsmp[iz][iy][ix].nY=iv+1+iAtomMed[iv]*epsdim;
                        }// if

                        if (dxyz.nZ<rad2a)
                        {
                            iepsmp[iz][iy][ix].nZ=iv+1+iAtomMed[iv]*epsdim;
                        }// if

                        if(distsq<radp2) idebmap[iz][iy][ix] =false;
                    }// do
                }// do
            }// do


        } //if
        else  /**faster method;*/
        {
//IT HAS PROBLEMS!!!! Walter (be careful before using also in multidilectric case!!!&&!isitmd
            //cout << "####faster method:" << endl;
            rad2a=rad2-0.25;

            ixn=optRound(xn);

            fxn=optCast <delphi_real,delphi_integer> (ixn)-xn;
            rad2av=rad2a-fxn;

            for(ix=-lim; ix<=lim; ix++)
            {
                vtemp= double(ix)+fxn;
                //sq[ix]=vtemp*vtemp;
                //rad2aav[ix]=rad2a-vtemp;
                sqtemp[ix+15]=vtemp*vtemp;
                rad2aavtemp[ix+15]=rad2a-vtemp;




            }// do

//adjust inter-atom, different epsilon bgps+++04/2004 Walter

            if (iNMedia>1&&bOnlyMol)
            {
                if(debug_space)cout << "#### multi media:" << endl;
                for(i=0; i<=ibox; i++)
                {

                    i123=ioff[i];
                    ixyz=ixn+i123;
                    ix=ixyz.nX;
                    iy=ixyz.nY;
                    iz=ixyz.nZ;
                    //distsq=sq[i123.nX].nX+sq[i123.nY].nY+sq[i123.nZ].nZ;
		    distsq = sqtemp[i123.nX+15].nX +sqtemp[i123.nY+15].nY + sqtemp[i123.nZ+15].nZ;
                    //if (distsq<rad2aav[i123.nX].nX)
                    if (distsq<rad2aavtemp[i123.nX+15].nX)
		    {
                        //iac=(iEpsMap[ix][iy][iz].nX % epsdim)-1;
                        iac=(iepsmp[ix][iy][iz].nX % epsdim)-1;

                        if (iac==-1||iac>iNatom)
                        {
//occhio! non so cosa fa con i pori!!
                            //iEpsMap[ix][iy][iz].nX=iv+1+iAtomMed[iv]*epsdim;
                            iepsmp[ix][iy][iz].nX=iv+1+iAtomMed[iv]*epsdim;
                        }
                        else
                        {
//2011-05-14 Using operations on coord and int_coord type variables defined in module operators_on_coordinates
                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
                            ddxyz.nX= ddxyz.nX+0.5;
                            dis2min1=optDot(ddxyz,ddxyz)-rad2;

                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn2[iac];
                            ddxyz.nX= ddxyz.nX+0.5;
                            dis2min2=optDot(ddxyz,ddxyz)-pow((sDelPhiPDB[iac].radius*fScale),2);

                            if (dis2min2>dis2min1) iac=iv;
                            //iEpsMap[ix][iy][iz].nX=iac+1+iAtomMed[iac]*epsdim;
                            iepsmp[ix][iy][iz].nX=iac+1+iAtomMed[iac]*epsdim;
                        }// if
                    }// if

                    //if (distsq<rad2aav[i123.nY].nY)
		    if (distsq<rad2aavtemp[i123.nY+15].nY)
                    {
                        //iac=(iEpsMap[ix][iy][iz].nY % epsdim)-1;
                        iac=(iepsmp[ix][iy][iz].nY % epsdim)-1;
                        if (iac==-1||iac>iNatom)
                        {
//occhio! non so cosa fa con i pori!!
                            //iEpsMap[ix][iy][iz].nY=iv+1+iAtomMed[iv]*epsdim;
                            iepsmp[ix][iy][iz].nY=iv+1+iAtomMed[iv]*epsdim;
                        }
                        else
                        {
//2011-05-14 Using operations on coord and int_coord type variables defined in module operators_on_coordinates
                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
                            ddxyz.nY= ddxyz.nY+0.5;
                            dis2min1=optDot(ddxyz,ddxyz)-rad2;

                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn2[iac];
                            ddxyz.nY= ddxyz.nY+0.5;
                            dis2min2=optDot(ddxyz,ddxyz)-pow((sDelPhiPDB[iac].radius*fScale),2);

                            if (dis2min2>dis2min1) iac=iv;
                            //iEpsMap[ix][iy][iz].nY=iac+1+iAtomMed[iac]*epsdim;
                            iepsmp[ix][iy][iz].nY=iac+1+iAtomMed[iac]*epsdim;
                        }// if
                    }// if

                    //if (distsq<rad2aav[i123.nZ].nZ)
		    if (distsq<rad2aavtemp[i123.nZ+15].nZ)
                    {
                        //iac=(iEpsMap[ix][iy][iz].nZ%epsdim)-1;
                        iac=(iepsmp[ix][iy][iz].nZ%epsdim)-1;
                        if (iac==-1||iac>iNatom)
                        {
//occhio! non so cosa fa con i pori!!
                            //iEpsMap[ix][iy][iz].nZ=iv+1+iAtomMed[iv]*epsdim;
                            iepsmp[ix][iy][iz].nZ=iv+1+iAtomMed[iv]*epsdim;
                        }
                        else
                        {

                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
                            ddxyz.nZ= ddxyz.nZ+0.5;
                            dis2min1=optDot(ddxyz,ddxyz)-rad2;

                            ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn2[iac];
                            ddxyz.nZ=ddxyz.nZ+0.5;
                            dis2min2=optDot(ddxyz,ddxyz)-pow((sDelPhiPDB[iac].radius*fScale),2);

                            if (dis2min2>dis2min1) iac=iv;
                            //iEpsMap[ix][iy][iz].nZ=iac+1+iAtomMed[iac]*epsdim;
                            iepsmp[ix][iy][iz].nZ=iac+1+iAtomMed[iac]*epsdim;
                        }// if
                    }// if

                    //if(distsq<radp2) bDebMap[ix][iy][iz]=false;
                    if(distsq<radp2) idebmap[ix][iy][iz]=false;
                }// do
            }
            else
            {
                if(debug_space) cout << "### fast method + single media:" << endl;

                //cout << "ix,rad2aav[ix]" << ix <<" " << rad2aav[ix] << endl;
		//cout << "#### ibox: " << ibox << endl;

	//cout << "### omp start: " << endl;
    	//omp_set_num_threads(4);
        //#pragma omp parallel for schedule(static) private(xn,rad,radp,rad2,radp2,itest2,ismin,ismax,itest,rad2a,ixn,fxn,rad2av,vtemp,rad2aavtemp,i123,ixyz,ix,iy,iz,distsq,sqtemp)
        //#pragma omp parallel for schedule(static) private(xn,rad,radp,rad2,radp2,itest2,ismin,ismax,itest,rad2a,ixn,fxn,rad2av,vtemp,rad2aavtemp,i123,ixyz,ix,iy,iz,distsq,sqtemp) shared(idebmap)

#ifdef PARALLEL_OMP
        #pragma omp for schedule(auto)
        //#pragma omp parallel for schedule(static) private(xn,rad,radp,rad2,radp2,itest2,ismin,ismax,itest,rad2a,ixn,fxn,rad2av,vtemp,rad2aavtemp,i123,ixyz,ix,iy,iz,distsq,sqtemp)
#endif

                for(i=0; i<=ibox; i++)
                {

                    i123=ioff[i];
                    ixyz=ixn+i123;
                    ix=ixyz.nX;
                    iy=ixyz.nY;
                    iz=ixyz.nZ;
                    distsq = sq[i123.nX].nX +sq[i123.nY].nY + sq[i123.nZ].nZ;
                    //cout << "i123,ixn,ix,iy,iz,distsq: " << i123 << " " <<ixn << " " <<ix << " " <<iy << " " <<iz << " " <<distsq << endl;
                    if (distsq<rad2aav[i123.nX].nX)
                    {
                        //Lin Li: all indexes -1, because of C++ arrays' indexes start from 0;
                        //Lin Li: because iv start from 0, so +2:
                        //iEpsMap[ix-1][iy-1][iz-1].nX=iv+2+iAtomMed[iv]*epsdim;
                        iepsmp[iz][iy][ix].nX=iv+1+iAtomMed[iv]*epsdim;

                        //cout << "ix,iy,iz,iEpsMap[ix-1][iy-1][iz-1].nX: " << ix<< " " << iy << " " << iz << " " <<iEpsMap[ix-1][iy-1][iz-1].nX << endl;
                        //cout << "ix,iy,iz,iv,iAtomMed[iv],epsdim " << ix<< " " << iy << " " << iz << " " << iv << " " << iAtomMed[iv]<< " " << epsdim << endl;
                    }// if

                    if (distsq<rad2aav[i123.nY].nY)
                    {
                        //iEpsMap[ix-1][iy-1][iz-1].nY=iv+2+iAtomMed[iv]*epsdim;
                        iepsmp[iz][iy][ix].nY=iv+1+iAtomMed[iv]*epsdim;
                    }// if

                    if (distsq<rad2aav[i123.nZ].nZ)
                    {
                        //iEpsMap[ix-1][iy-1][iz-1].nZ=iv+2+iAtomMed[iv]*epsdim;
                        iepsmp[iz][iy][ix].nZ=iv+1+iAtomMed[iv]*epsdim;
                    }// if
#ifdef IJK
                    if (distsq<radp2) idebmap[ix][iy][iz]=false;
#endif // IJK
                    if (distsq<radp2) idebmap[iz][iy][ix]=false;


                }// do


            }// if
        }// if
    }// do DoATOMS //end do of atoms;

    if (ioff != NULL)    delete [] ioff;

#ifdef DEBUG_DELPHI_SPACE
    ofstream epsfile;
    ofstream debfile;
    epsfile.open ("epsmap_win.txt");
    debfile.open ("debmap_win.txt");
    for (ix=1; ix<=iGrid; ix++)
    {
        for (iy=1; iy<=iGrid; iy++)
        {
            for (iz=1; iz<=iGrid; iz++)
            {
                epsfile << setw(5) << iEpsMap[ix-1][iy-1][iz-1].nX
                        << setw(5) << iEpsMap[ix-1][iy-1][iz-1].nY
                        << setw(5) << iEpsMap[ix-1][iy-1][iz-1].nZ
                        << endl;
                debfile << bDebMap[ix-1][iy-1][iz-1] << endl;
                /*
                if(bDebMap[ix-1][iy-1][iz-1]){
                   debfile << "T" << " " << bDebMap[ix-1][iy-1][iz-1] <<endl;

                }
                else{
                   debfile << "F" << " " << bDebMap[ix-1][iy-1][iz-1] <<endl;
                }*/

            }
        }
    }

    epsfile.close();
    debfile.close();

#endif //DEBUG_DELPHI_SPACE

    if(false)
    {
        // testing iepsmp and idebmap:
        ofstream epsfile;
        ofstream debfile;
        epsfile.open ("epsmap_win_bastar.txt");
        debfile.open ("debmap_win_bastar.txt");
        for (ix=1; ix<=iGrid; ix++)
        {
            for (iy=1; iy<=iGrid; iy++)
            {
                for (iz=1; iz<=iGrid; iz++)
                {
                    epsfile << setw(5) << iepsmp[iz][iy][ix].nX
                            << setw(5) << iepsmp[iz][iy][ix].nY
                            << setw(5) << iepsmp[iz][iy][ix].nZ
                            << endl;
#ifdef IJK
                    debfile << idebmap[ix][iy][iz] << endl;
#endif // IJK
                    debfile << idebmap[iz][iy][ix] << endl;

                    /*
                    if(bDebMap[ix-1][iy-1][iz-1]){
                       debfile << "T" << " " << bDebMap[ix-1][iy-1][iz-1] <<endl;

                    }
                    else{
                       debfile << "F" << " " << bDebMap[ix-1][iy-1][iz-1] <<endl;
                    }*/

                }
            }
        }

        epsfile.close();
        debfile.close();

    }

/*
#ifdef DEBUG
    //open(52,file='iepsmapnewfirst',form='formatted');
    ofstream file1;
    file1.open("iepsmap.txt");
    //write (52,*)iGrid,iGrid*iGrid;
    file1 << "TEST" << endl;
    for(ix=1; ix<=iGrid; ix++)
    {
        for(iy=1; iy<=iGrid; iy++)
        {
            iz=(iGrid+1)/2;
            file1 << ix << " " << iy << " " << iz << " " << iEpsMap[ix][iy][iz].nZ << endl;
        }// do
    }// do
    //close (52);
    file1.close();
#endif // DEBUG
*/
#ifdef VERBOSE
    cout <<"Ending creating Van der Waals Epsilon Map " << endl;
#endif

    //test_pdc->showMap("test_setout.txt");
    return;

}// void setout;

