//Tips for Gaussian:
//the new statement can be added by modifying following files:
//
//delphi_constants.h : increase iStatementNum
//delphi_datamarshal.h
//delphi_data_setMap
//delphi_data_setDefault
//delphi_data_getstatement

// new variables and arrays increased in Gaussian & MEMPOT:
// cutoff,sigma,srfcut,radipz,inhomo,gaussian,ergsgaussian
// gepsmp(:,:,:),gepsmp2(:,:,:)



#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"
#include "space.h"

using namespace std;


void CDelphiSpace::setGaussian()
{


    //SGrid <float> sq[-150:150],rad2aav[-150:150],vtemp2[-150:150]//Gaussian enlarged the values;

    SGrid <delphi_real>  sqtemp[301],rad2aavtemp[301],vtemp2temp[301];
    delphi_real eps_min,eps_max,eps_diff,bnd_buff;
    SGrid <delphi_real> * sq=sqtemp+150;
    SGrid <delphi_real> * rad2aav=rad2aavtemp+150;
    SGrid <delphi_real> * vtemp2=vtemp2temp+150;
    SGrid <delphi_integer>* ioff;

    ioff=NULL;
    bool itobig,itest2,ipore,bOnlyMolbDebug;
    string strtmp,strtmp1;
// 2011-05-12 Non-standard int variable, thus necessary
    //int epsdim, objecttype,iflag;
    delphi_integer epsdim;

// 2011-05-12 Non-standard float variable, thus necessary
    float modul,modul2, mod2,modx,mody,dentemp;
// 2011-05-12 Non-standard type of variable, thus necessary
    SGrid <delphi_real> xa,xb,xc,xd,xp,ddist,xyz2,dxyz,ddxyz,xn,vecdist;
    SGrid <delphi_real> tmpvect1,tmpvect2,origin;
    SGrid <delphi_real> vectx,vecty,vectz,rad2av,fxn,vtemp;
    SGrid <delphi_integer> ismin,ismax,idist,idist1,ixyz,itest,ixn,i123;
    delphi_real coeff,stepsize;
    delphi_integer longint;


//------------------------------------------------------------------------------------------
// 2011-05-27 Declarations added due to IMPLICIT NONE
    int iboxt,iac,ibox,ii,igrdc,i,j,k,imedia,iv,ix,iy,iz,kind;
    delphi_integer limmax,lim,n;
    //integer,dimension(1:6) ::inwater //mid point;
    //logical,dimension(1:6) ::ifinwater //mid point;
    delphi_real alpha,dis2min2,dis2min1,dot,distsq,exrd,dist,rad,rad2;
    delphi_real rad2a,rad4,radius,radmax2,fRMid,rad5,radp2,radtest,radprobe;
    delphi_real radp,tan2,temp,tmp,tmp1,tmp2,epstemp;
    delphi_real peak,pi,den,distance2,sigmatime,radsq; //sigma=sigma*radius. peak=4/(3*pi**0.5*sigma**3)

    if(debug_space)cout << "######## start Gaussian: #########" << endl;
    eps_diff=1.2;
    //coeff=0.5291772108;
    //stepsize=1.0/fScale;
    //origin=(cOldMid-stepsize*(iGrid-1)/2)/coeff;
    //cout << "origin,cOldMid: " << origin << " "  << endl;
    //cout << "origin,cOldMid: " << cOldMid.nX << " " << cOldMid.nY << " " << cOldMid.nZ << endl;

    if(!(iGaussian==1&&inhomo==0&&logs))
    {
        //goto 1010;
        //if

        pi=3.14159265359;
        sigmatime=3.0;
        for(i=1; i<=iGrid; i++)
        {
            for(j=1; j<=iGrid; j++)
            {
                for(k=1; k<=iGrid; k++)
                {
                    gepsmp2[i][j][k].nX=repsout;
                    gepsmp2[i][j][k].nY=repsout;
                    gepsmp2[i][j][k].nZ=repsout;

                    gepsmp[i][j][k].nX=0.0;
                    gepsmp[i][j][k].nY=0.0;
                    gepsmp[i][j][k].nZ=0.0;
                }
            }
        }

        if(debug_space) cout << "############ beginning of setout: ##############" << endl;
#ifdef VERBOSE
        cout << "Starting creating Van der Waals Epsilon Map: " << endl;
#endif
        radprobe=0; // this radprobe seems not useful.
        epsdim=iNatom+iNObject+2;
        //iboxt=0; //not useful
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


            //radmax2=fScale*(radmax2+temp);
            radmax2=sigmatime*fScale*(radmax2*sigma+temp); //Gaussian: 3 sigma plus temp. sigma=sigma*radius.Now radmax is 3 sigma.
            lim=1+radmax2;
            limmax = 1200; //Gaussian:original value:12

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



        //set interiors in OBJECTS

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
            rad= sigmatime*sigma*rad; // 3 sigma for Gaussian

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
            //rad4=pow( (rad+0.5),2); // not used
            rad2=rad*rad;
            radp2=radp*radp;

            radsq=rad2/(sigmatime*sigmatime); // for Gaussian

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

            if(debug_space) cout << "### itobig: " << itobig  << endl;
            //if (itest2||itobig)   //slow method;
            if (itobig)   //slow method;
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
                if(debug_space) cout << "####faster method:" << endl;
                rad2a=rad2-0.25;

                ixn=optRound(xn);

                fxn=optCast <delphi_real,delphi_integer> (ixn)-xn;
                rad2av=rad2a-fxn;

                for(ix=-lim; ix<=lim; ix++)
                {
                    vtemp= double(ix)+fxn;
                    //sq[ix]=vtemp*vtemp;
                    //rad2aav[ix]=rad2a-vtemp;
                    sqtemp[ix+150]=vtemp*vtemp;
                    rad2aavtemp[ix+150]=rad2a-vtemp;
                    vtemp2[ix]=vtemp;
                    //if(debug_space) cout << "#### ix: " << ix << endl;

                }// do

//adjust inter-atom, different epsilon bgps+++04/2004 Walter
                if(debug_space)cout << "#### iNMedia,bOnlyMol: " << iNMedia << " " << bOnlyMol << endl;
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
                        distsq = sqtemp[i123.nX+150].nX +sqtemp[i123.nY+150].nY + sqtemp[i123.nZ+150].nZ;
                        //if (distsq<rad2aav[i123.nX].nX)
                        if (distsq<rad2aavtemp[i123.nX+150].nX)
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
                        if (distsq<rad2aavtemp[i123.nY+150].nY)
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
                        if (distsq<rad2aavtemp[i123.nZ+150].nZ)
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
                        //cout << "i,ix,iy,iz: " << i << " " << ix << " " << iy << " " << iz << endl;
                        if(ix<iGrid&&ix>1 && iy<iGrid&&iy>1 && iz<iGrid&&iz>1)
                        {
                            distsq = sq[i123.nX].nX +sq[i123.nY].nY + sq[i123.nZ].nZ;
                            //if(debug_space) cout << iv << " " << i << " " << distsq << " " << sq[i123.nX].nX << "+" << sq[i123.nY].nY << "+" << sq[i123.nZ].nZ<< endl;
                        }
                        else // for grid outside the box
                        {
                            //distsq=10000.0;
                            continue;
                        }

//---------------Lin Li: key section for Gaussian:
                        //cout << "i123,ixn,ix,iy,iz,distsq: " << i123 << " " <<ixn << " " <<ix << " " <<iy << " " <<iz << " " <<distsq << endl;
                        if (distsq<rad2aav[i123.nX].nX)
                        {
                            //iepsmp(ix,iy,iz)%i=iv+1+iatmmed(iv)*epsdim
                            //iepsmp[iz][iy][ix].nX=iv+1+iAtomMed[iv]*epsdim;
                            //cout << "distsq<rad2aav[i123.nX].nX: " << distsq << " " << rad2aav[i123.nX].nX << endl;

                            distance2=distsq+0.25+vtemp2[i123.nX].nX;
                            den=exp(-(distance2/(sigma*sigma*radsq)));
                            // cout << ix << " " << iy << " " << iz << " " << distsq << " " << vtemp2[i123.nX].nX << " " << distance2 << " " << den << " "  << endl;


                            //gepsmp[iz][iy][ix].nX=1-(1-gepsmp[iz][iy][ix].nX)*(1-den);
                            gepsmp[ix][iy][iz].nX=1-(1-gepsmp[ix][iy][iz].nX)*(1-den);
                            //if(gepsmp[iz][iy][ix].nX > 0.2) cout << iz << " " << iy << " " << ix << " " << gepsmp[iz][iy][ix].nX << endl;

                            //cout << ix << " " << iy << " " << iz << " " << distsq << " " << vtemp2[i123.nX].nX << " " << distance2 << " " << den << " " << gepsmp[ix][iy][iz].nX << endl;

                        }// if

                        if (distsq<rad2aav[i123.nY].nY)
                        {
                            //iEpsMap[ix-1][iy-1][iz-1].nY=iv+2+iAtomMed[iv]*epsdim;
                            //iepsmp[iz][iy][ix].nY=iv+1+iAtomMed[iv]*epsdim;

                            distance2=distsq+0.25+vtemp2[i123.nY].nY;
                            den=exp(-(distance2/(sigma*sigma*radsq)));
                            gepsmp[ix][iy][iz].nY=1-(1-gepsmp[ix][iy][iz].nY)*(1-den);


                        }// if

                        if (distsq<rad2aav[i123.nZ].nZ)
                        {
                            //iEpsMap[ix-1][iy-1][iz-1].nZ=iv+2+iAtomMed[iv]*epsdim;
                            //iepsmp[iz][iy][ix].nZ=iv+1+iAtomMed[iv]*epsdim;
                            distance2=distsq+0.25+vtemp2[i123.nZ].nZ;
                            den=exp(-(distance2/(sigma*sigma*radsq))); //Gaussian: make density at center is 1.
                            gepsmp[ix][iy][iz].nZ=1-(1-gepsmp[ix][iy][iz].nZ)*(1-den);

                        }// if
#ifdef IJK
                        if (distsq<radp2) idebmap[ix][iy][iz]=false;
#endif // IJK
                        //if (distsq<radp2) idebmap[iz][iy][ix]=false; //Don't need it in Gaussian

                    }// do

                }// if
            }// if
        }// do DoATOMS //end do of atoms;








// cutoff=0.9
// cout << << "cutoff << sigma:" << cutoff << sigma
        if(false)  //cutoff: make the top flat;
        {
            for(ix=1; ix<=iGrid; ix++)
            {
                for(iy=1; iy<=iGrid; iy++)
                {
                    for(iz=1; iz<=iGrid; iz++)
                    {
                        gepsmp[ix][iy][iz].nX=gepsmp[ix][iy][iz].nX/cutoff;
                        gepsmp[ix][iy][iz].nY=gepsmp[ix][iy][iz].nY/cutoff;
                        gepsmp[ix][iy][iz].nZ=gepsmp[ix][iy][iz].nZ/cutoff;
                        if(gepsmp[ix][iy][iz].nX>1.0) gepsmp[ix][iy][iz].nX=1.0;
                        if(gepsmp[ix][iy][iz].nY>1.0) gepsmp[ix][iy][iz].nY=1.0;
                        if(gepsmp[ix][iy][iz].nZ>1.0) gepsmp[ix][iy][iz].nZ=1.0;
                    }//do
                }//do
            }//do
        }//if
//}
    } //if(!(iGaussian==1&&inhomo==0&&logs))1010 continue;


// cout << << "bDebug: srfcut << inhomo:" << srfcut << inhomo
    for(ix=1; ix<=iGrid; ix++)
    {
        for(iy=1; iy<=iGrid; iy++)
        {
            for(iz=1; iz<=iGrid; iz++)
            {
                gepsmp2[ix][iy][iz].nX=gepsmp[ix][iy][iz].nX*repsin+(1-gepsmp[ix][iy][iz].nX)*repsout;
                gepsmp2[ix][iy][iz].nY=gepsmp[ix][iy][iz].nY*repsin+(1-gepsmp[ix][iy][iz].nY)*repsout;
                gepsmp2[ix][iy][iz].nZ=gepsmp[ix][iy][iz].nZ*repsin+(1-gepsmp[ix][iy][iz].nZ)*repsout;
//################### for set epsout in protein larger than in water:##########
                if(gepsmp[ix][iy][iz].nX<0.02)
                {
                    gepsmp2[ix][iy][iz].nX=80.0;
                }//if
                if(gepsmp[ix][iy][iz].nY<0.02)
                {
                    gepsmp2[ix][iy][iz].nY=80.0;
                }//if
                if(gepsmp[ix][iy][iz].nZ<0.02)
                {
                    gepsmp2[ix][iy][iz].nZ=80.0;
                }//if
//################### }// for this epsout in protein different than in water######
            }//do
        }//do
    }//do


    if(inhomo==1)  //reduce epsilon out side protein;
    {
// epstemp=repsin*(1-srfcut)+repsout*srfcut
        epstemp=srfcut;
        if(epstemp<repsin)
        {
            cout << "srfcut is lower than epsin." << endl;
            exit(0);
        }//if
// print* << "epstemp" << epstemp
        for(ix=1; ix<=iGrid; ix++)
        {
            for(iy=1; iy<=iGrid; iy++)
            {
                for(iz=1; iz<=iGrid; iz++)
                {
                    if(gepsmp2[ix][iy][iz].nX>epstemp)
                    {
                        gepsmp2[ix][iy][iz].nX=repsin;
                    }// if
                    if(gepsmp2[ix][iy][iz].nY>epstemp)
                    {
                        gepsmp2[ix][iy][iz].nY=repsin;
                    }// if
                    if(gepsmp2[ix][iy][iz].nZ>epstemp)
                    {
                        gepsmp2[ix][iy][iz].nZ=repsin;
                    }// if
                }//do
            }//do
        }//do
    }//if

    iBoundNum=0;
    longint=0;
    dentemp=0.1;
// ifinwater=true
    //for(i=2; i<=iGrid; i++)
    for(i=2; i<iGrid; i++)
    {
        //for(j=2; j<=iGrid; j++)
        for(j=2; j<iGrid; j++)
        {
            //for(k=2; k<=iGrid; k++)
            for(k=2; k<iGrid; k++)
            {
                if(gepsmp[i][j][k].nX>dentemp||gepsmp[i][j][k].nY>dentemp ||gepsmp[i][j][k].nZ>dentemp||gepsmp[i-1][j][k].nX>dentemp ||gepsmp[i][j-1][k].nY>dentemp||gepsmp[i][j][k-1].nZ>dentemp)
                {

                    //idebmap[i][j][k]=false; //iGaussian change method of generating bDebMap;
                    idebmap[k][j][i]=false; //iGaussian change method of generating bDebMap;

                    // see if the epslon distribution is flat or not
                    eps_min=gepsmp[i][j][k].nX;
                    eps_max=gepsmp[i][j][k].nX;

                    eps_min=min(eps_min,gepsmp[i][j][k].nY);
                    eps_min=min(eps_min,gepsmp[i][j][k].nZ);
                    eps_min=min(eps_min,gepsmp[i-1][j][k].nX);
                    eps_min=min(eps_min,gepsmp[i][j-1][k].nY);
                    eps_min=min(eps_min,gepsmp[i][j][k-1].nZ);

                    eps_max=max(eps_max,gepsmp[i][j][k].nY);
                    eps_max=max(eps_max,gepsmp[i][j][k].nZ);
                    eps_max=max(eps_max,gepsmp[i-1][j][k].nX);
                    eps_max=max(eps_max,gepsmp[i][j-1][k].nY);
                    eps_max=max(eps_max,gepsmp[i][j][k-1].nZ);

                    if(eps_max/eps_min > eps_diff) longint=longint+1;

                }//if
            }//do
        }//do
    }//do

    iBoundNum=longint;
    //allocate(ibgrd[iBoundNum]);
    if(debug_space) cout << "###### iboundnum: " <<iBoundNum << endl;
    if(debug_space)cout << "inhomo: " << inhomo << endl;
    ibgrd_v.assign(iBoundNum, sgrid_temp_int);
    ibgrd=&ibgrd_v[-1];

    n=0;
    //for(i=2; i<=iGrid; i++)
    for(i=2; i<iGrid; i++)
    {
        //for(j=2; j<=iGrid; j++)
        for(j=2; j<iGrid; j++)
        {
            //for(k=2; k<=iGrid; k++)
            for(k=2; k<iGrid; k++)
            {

                if(gepsmp[i][j][k].nX>dentemp||gepsmp[i][j][k].nY>dentemp ||gepsmp[i][j][k].nZ>dentemp||gepsmp[i-1][j][k].nX>dentemp ||gepsmp[i][j-1][k].nY>dentemp||gepsmp[i][j][k-1].nZ>dentemp)


                {
                    eps_min=gepsmp[i][j][k].nX;
                    eps_max=gepsmp[i][j][k].nX;

                    eps_min=min(eps_min,gepsmp[i][j][k].nY);
                    eps_min=min(eps_min,gepsmp[i][j][k].nZ);
                    eps_min=min(eps_min,gepsmp[i-1][j][k].nX);
                    eps_min=min(eps_min,gepsmp[i][j-1][k].nY);
                    eps_min=min(eps_min,gepsmp[i][j][k-1].nZ);

                    eps_max=max(eps_max,gepsmp[i][j][k].nY);
                    eps_max=max(eps_max,gepsmp[i][j][k].nZ);
                    eps_max=max(eps_max,gepsmp[i-1][j][k].nX);
                    eps_max=max(eps_max,gepsmp[i][j-1][k].nY);
                    eps_max=max(eps_max,gepsmp[i][j][k-1].nZ);

                    if(eps_max/eps_min > eps_diff)
                    {

                        n=n+1;
                        ibgrd[n].nX=i;
                        ibgrd[n].nY=j;
                        ibgrd[n].nZ=k;
                    }

                }//if
            }//do
        }//do
    }//do

// Array deallocation is made by ordinary F95 statement
    //if(allocated(ioff)) deallocate(ioff);

    if (ioff != NULL)
    {
        delete [] ioff;
        ioff=NULL;
    }
    //cout << "#### ioff: "<<ioff << endl;
    /*
    //!c b+++++bDebug+++++++++++++++++++++++++++++++
        if (false)   //Gaussian;
        {
            open(52,file='linlieps_k',form='formatted');
            write (52,*)iGrid,iGrid*iGrid;
            for(ix=1; ix<=iGrid; ix++)
            {
                for(iy=1; iy<=iGrid; iy++)
                {
                    for(iz=1; iz<=iGrid; iz++)
                    {
                        write (52,*)ix,iy,iz,iepsmp[ix][iy][iz].nZ;
                    }// do
                }// do
            }// do
            close (52);

        }// if
    */


    if(debug_space)cout << "####test" << endl;
    //if (true)   //bDebug: output density and eps;
    if(debug_space)
    {
        //open(14,file="cube_density", form='formatted');
        coeff=0.5291772108;
        stepsize=1.0/fScale;
        origin=(cOldMid-stepsize*(iGrid-1)/2)/coeff;

        // -------density file: --------
        ofstream densfile;
        densfile.open ("cube_density.txt");

        densfile << "qdiffxs4 with an improved surfacing routine" << endl;
        densfile << "Gaussian cube format phimap" << endl;

        densfile << fixed << setprecision(6);
        densfile << setw(5) << 1 << setw(12) << origin.nX << setw(12) << origin.nY << setw(12) << origin.nZ << endl;
        densfile << setw(5) << iGrid << setw(12) << stepsize/coeff << setw(12) << 0.0 << setw(12) << 0.0 << endl;
        densfile << setw(5) << iGrid << setw(12) << 0.0 << setw(12) << stepsize/coeff << setw(12) << 0.0 << endl;
        densfile << setw(5) << iGrid << setw(12) << 0.0 << setw(12) << 0.0 << setw(12) << stepsize/coeff << endl;
        densfile << setw(5) << 1 << setw(12) << 0.0 << setw(12) << 0.0 << setw(12) << 0.0 << setw(12) << 0.0 << endl;

        densfile << setprecision(5);
        for(ix=1; ix<=iGrid; ix++)
        {
            for(iy=1; iy<=iGrid; iy++)
            {
                for(iz=1; iz<=iGrid; iz++)
                {
                    densfile << setw(13) << scientific << gepsmp[ix][iy][iz].nZ << " " ;
                    if(iz%6 == 0) densfile << endl;
                }
                densfile << endl;
            }
        }
        densfile.close();

        //-------epsilon file: --------

        ofstream epsfile;
        epsfile.open ("cube_eps.txt");

        epsfile << "qdiffxs4 with an improved surfacing routine" << endl;
        epsfile << "Gaussian cube format phimap" << endl;

        epsfile << fixed << setprecision(6);
        epsfile << setw(5) << 1 << setw(12) << origin.nX << setw(12) << origin.nY << setw(12) << origin.nZ << endl;
        epsfile << setw(5) << iGrid << setw(12) << stepsize/coeff << setw(12) << 0.0 << setw(12) << 0.0 << endl;
        epsfile << setw(5) << iGrid << setw(12) << 0.0 << setw(12) << stepsize/coeff << setw(12) << 0.0 << endl;
        epsfile << setw(5) << iGrid << setw(12) << 0.0 << setw(12) << 0.0 << setw(12) << stepsize/coeff << endl;
        epsfile << setw(5) << 1 << setw(12) << 0.0 << setw(12) << 0.0 << setw(12) << 0.0 << setw(12) << 0.0 << endl;

        epsfile << setprecision(5);
        for(ix=1; ix<=iGrid; ix++)
        {
            for(iy=1; iy<=iGrid; iy++)
            {
                for(iz=1; iz<=iGrid; iz++)
                {
                    epsfile << setw(13) << scientific << gepsmp2[ix][iy][iz].nZ << " " ;
                    if(iz%6 == 0) epsfile << endl;
                }
                epsfile << endl;
            }
        }
        epsfile.close();
    }

    /*
            open(14,file="cube_eps", form='formatted');

            write(14,*)'qdiffxs4 with an improved surfacing routine';
            write(14,*) 'Gaussian cube format phimap';
            coeff=0.5291772108;
            stepsize=1.0/fScale;
            origin=(cOldMid-stepsize*(iGrid-1)/2)/coeff;
    // do i=1,3
    // origin(i)=(cOldMid(i)-stepsize*(iGrid-1)/2)/coeff
    // }//do
            write(14,'(i5,3f12.6)') 1, origin;
            write(14,'(i5,3f12.6)') iGrid, stepsize/coeff,0.0,0.0;
            write(14,'(i5,3f12.6)') iGrid, 0.0,stepsize/coeff,0.0;
            write(14,'(i5,3f12.6)') iGrid, 0.0,0.0,stepsize/coeff;
            write(14,'(i5,4f12.6)') 1,0.0,0.0,0.0,0.0;

            for(ix=1; ix<=iGrid; ix++)
            {
                for(iy=1; iy<=iGrid; iy++)
                {
                    write(14,'(6E13.5)')(gepsmp2[ix][iy][iz].nZ,iz=1,iGrid);
                }// do
            }// do
            close(14);
        }//if

    */


//!c e++++++++++++++++++++++++++++++++++++++++++
#ifdef VERBOSE
    cout <<"Ending creating Van der Waals Epsilon Map " << endl;
#endif // VERBOSE


    if(debug_space)cout << "###out of Gaussian#############" << endl;
    return;
}// void setiGaussian;

