#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "../space/space.h"

using namespace std;
//#####################################################################
//ATTENTION! This file is part of epsmakmod module.
// Do not compile separately!
//
//program to reposition the boundary grid points on the molecular surface
//(S. Sridharan May 1994)
//
//2011-05-19 Other parameters to void are transferred via qlog
// and pointers modules
// Parameter iBoundNum is now local because it can take different
// values when sclbp is called from msrf and vwtms
//2011-05-19 Arrays declared in pointers module allocated in vwtms
// void
//#####################################################################

void CDelphiSpace::sclbp() // scspos-> scspos; scsnor-> scsnor; iBoundNum-> iBoundNum
{

    //delphi_integer iab1[0:lcb1,0:mcb1,0:ncb1],iab2[0:lcb1,0:mcb1,0:ncb1];
    //SGrid <delphi_real> scsnor[iBoundNum], scspos[iBoundNum];
    delphi_integer nbra[1001];
    //bool out,outcb[-2:2,-2:2,-2:2];
    bool out;
    bool outcb[5][5][5]; //outcb index should be modified

//2011-05-19 Leaving only variables not declared in qlog and pointers modules
    delphi_integer iaprec,dim1,kind;
    delphi_integer dim,prevmed,med;
    delphi_integer epsdim,imezzo[7],iente[7],ix,iy,iz;
    SGrid <delphi_integer > ixyz;
    delphi_integer iac1,iac2;
    bool lga,lgd,iflag,precedenza,vicinanza,flag;
    string strtmp;
    //delphi_real disev,disodd;
    SGrid <delphi_real> vnor,s123;
    //delphi_real tmp[7];
    delphi_real radpmax,dst;
    SGrid <delphi_real> x1xyz,xg, u123,xxyyzz;
    delphi_real temp;
    SGrid <delphi_real> xq,dxyz,dixyz,dx123, dr;
    //delphi_real tan2,dx,dy,dz;
    delphi_real hgs,ds2min1,ds2min2;
    SGrid <delphi_integer > it,jxyz;
//2011-05-27 Declarations added due to IMPLICIT NONE
    delphi_integer iac,ia,i,j,k,iacl,iii,ii,jjx,jjy,jjz,jx,jy,jz;
    delphi_integer jzi,jyi,jxi,liml,limu,kk,nnbr,ncbp;
    delphi_real x1,cbln,del,dis2,dis,dist,dmn,dcr,ctf,dmn1,dmn2;
    delphi_real cba,ds2,dsr,rmn,rdist,dmx;

    vector <bool> internal;

//2011-05-19 Array allocated by ordinary F95 allocate statement

    if(debug_space) cout << "### in sclbp: ###" << endl;
    //allocate(internal(iNObject));
    internal.assign(iNObject+1,false);
    epsdim=iNatom+iNObject+2;
    iflag=false;
    iac1=0;
    iac2=0;

    iall=0;

/**
 * hgs= half grid spacing
*/
    hgs=1./(2.*fScale);

//2011-05-17 Changed to array operations
    //outcb=true; initiallized above.
    for (i=0; i<=4; i++)
    {
        for (j=0; j<=4; j++)
        {
            for (k=0; k<=4; k++)
            {
                outcb[i][j][k]=true;
            }
        }
    }
    for (i=1; i<=3; i++)
    {
        for (j=1; j<=3; j++)
        {
            for (k=1; k<=3; k++)
            {
                outcb[i][j][k]=false;
            }
        }
    }
    //outcb[-1:1][-1:1][-1:1]=false;

/**
 *convertion from grid to delphi_real coordinates(can also use routine
 *gtoc)
*/
    x1=1.0/fScale;

//2011-05-19 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
    x1xyz=cOldMid-(0.5*x1*(iGrid+1));
    radpmax=max(fRadPrb[1],fRadPrb[2]);

    if (extot==0&&radpmax>0.0&&(iNObject>1||iNatom>1))
    {
//find extrema
//here one should consider the global system (Walter)
        cout <<"Scaling routine in action//" << endl;

//2011-05-19 Using operations on coord and int_coord type variables defined in module operators_on_coordinates
        //cMin= {6000.,6000.,6000.};
        cMin.nX=6000.;
        cMin.nY=6000.;
        cMin.nZ=6000.;

        //cMax= {-6000.,-6000.,-6000.};
        cMax.nX=-6000.;
        cMax.nY=-6000.;
        cMax.nZ=-6000.;


        for(ii=1; ii<=iNObject; ii++)
        {
            cMin=optMin(cMin,sLimObject[ii].nMin);
            cMax=optMax(cMax,sLimObject[ii].nMax);
        }// do

//2011-05-19 All other parameters are transfered via modules
//qlog and pointers
        sas();
    }// if

    del=radpmax;
    del=max(del,1./(2.*fScale));
    cbln=fRMax+del;

    cubedata(2.0,cbln);

    //dim=(lcb+1)*(mcb+1)*(ncb+1); //cbn1,cbn2 are semiglobal
    //allocate(cbn1[dim],cbn2[dim]);


    dim1=27;
    if ((iNObject-numbmol)>0) dim1=max(dim,27);
    //allocate(cbal[dim1*(iNatom+iNObject-numbmol)]);
    cbal.assign(dim1*(iNatom+iNObject-numbmol)+1,0);

//2011-05-19 All parameters transfered via modules qlog and
//pointers
    cube();

    ncbp=0;

//D500:
    //do i=1,iBoundNum;
    cout << "iBoundNum: " << iBoundNum << endl;
    for(i=1; i<=iBoundNum; i++)
    {

//+per trattare molecole con diversa epsilon++01/2002+
        if (iBoundNum!=iBoundNumsurf&&numbmol>1)
        {
//2011-05-19 Converted to int_coord derived type
            cout << "###### sclbp iente iepsmp: " << endl;
            ixyz=optCast <delphi_integer,delphi_real> (scspos[i]);

            ix=ixyz.nX;
            iy=ixyz.nY;
            iz=ixyz.nZ;

            iflag=false;

            iente[1]=iepsmp[ix][iy][iz].nX%epsdim;
            iente[2]=iepsmp[ix][iy][iz].nY%epsdim;
            iente[3]=iepsmp[ix][iy][iz].nZ%epsdim;
            iente[4]=iepsmp[ix-1][iy][iz].nX%epsdim;
            iente[5]=iepsmp[ix][iy-1][iz].nY%epsdim;
            iente[6]=iepsmp[ix][iy][iz-1].nX%epsdim;

            imezzo[1]=iepsmp[ix][iy][iz].nX/epsdim;
            imezzo[2]=iepsmp[ix][iy][iz].nY/epsdim;
            imezzo[3]=iepsmp[ix][iy][iz].nZ/epsdim;
            imezzo[4]=iepsmp[ix-1][iy][iz].nX/epsdim;
            imezzo[5]=iepsmp[ix][iy-1][iz].nY/epsdim;
            imezzo[6]=iepsmp[ix][iy][iz-1].nZ/epsdim;

            //iente[1]=iEpsMap[ix][iy][iz].nX%epsdim;
            //iente[2]=iEpsMap[ix][iy][iz].nY%epsdim;
            //iente[3]=iEpsMap[ix][iy][iz].nZ%epsdim;
            //iente[4]=iEpsMap[ix-1][iy][iz].nX%epsdim;
            //iente[5]=iEpsMap[ix][iy-1][iz].nY%epsdim;
            //iente[6]=iEpsMap[ix][iy][iz-1].nX%epsdim;

            //imezzo[1]=iEpsMap[ix][iy][iz].nX/epsdim;
            //imezzo[2]=iEpsMap[ix][iy][iz].nY/epsdim;
            //imezzo[3]=iEpsMap[ix][iy][iz].nZ/epsdim;
            //imezzo[4]=iEpsMap[ix-1][iy][iz].nX/epsdim;
            //imezzo[5]=iEpsMap[ix][iy-1][iz].nY/epsdim;
            //imezzo[6]=iEpsMap[ix][iy][iz-1].nZ/epsdim;

//guardo se ho due molecole con diversa epsilon nel punto,interno
            if(imezzo[1]!=imezzo[6]&&imezzo[1]*imezzo[6]!=0) iflag=(iente[1]<=iNatom+1&&iente[6]<=iNatom+1);

//iflag sar?vero se il bgp ?interno ma non con un oggetto
            for(ii=2; ii<=6; ii++)
            {
                if(imezzo[ii]!=imezzo[ii-1] && imezzo[ii]*imezzo[ii-1]!=0) iflag=iflag||(iente[ii]<=iNatom+1 && iente[ii-1]<=iNatom+1);
            }// do
        }// if


//2011-05-19 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
        xg=scspos[i]*x1+x1xyz;

/**
 * find the closest surface atom to the gridpoint
*/
        it=optCast <delphi_integer,delphi_real> ((xg-xyzo)*cbai);

        dmn=100.;
        ds2min1=1000.;
        ds2min2=1000.;
        prevmed=0;
        iac=0;
        nnbr=0;
        //lmncb={lcb,mcb,ncb};
        lmncb.nX=lcb;
        lmncb.nY=mcb;
        lmncb.nZ=ncb;

        if (optORLT(it,0)||optORGT(it,lmncb))
        {
/**
 * if the bgp is outside the cube, probably it is due to some object
*/
            for(ii=1; ii<=iNObject; ii++)
            {

                //strtmp=dataobject[ii-1][0];
                strtmp=dataobject_v[(ii-1)*2];

                //read(strtmp(16:18),*)kind;
                //kind=strtmp.substr(15,3);
                kind = atoi(strtmp.substr(15,3).c_str());

                //if (strtmp[1:4]!='is a'&&kind!=2)
                if (strtmp.substr(0,4)!="is a" && kind!=2)

                {
                    if ( optANDLE(xg,(sLimObject[ii].nMax+x1) )&& optANDGT(xg,(sLimObject[ii].nMin-x1)) )
                    {
                        nnbr=nnbr+1;
                        nbra[nnbr]=ii+iNatom;
                        liml=1;
                        limu=0;
                    }// if
                }// if
            }// do

            if(liml!=1||limu!=0) cout <<"bgp close to nothing" << endl;
        }
        else
        {
//2011-05-19 Changed 1d array to 3d array as in cube
//void
            //liml=cbn1[it.nX+1+(lcb+1)*it.nY+(lcb+1)*(mcb+1)*it.nZ];
            //limu=cbn2[it.nX+1+(lcb+1)*it.nY+(lcb+1)*(mcb+1)*it.nZ];

            liml=cbn1[it.nX][it.nY][it.nZ];
            limu=cbn2[it.nX][it.nY][it.nZ];

        }// if

        iaprec=0;

        for(kk=liml; kk<=limu; kk++)
        {
            ia=cbal[kk];

            if (ia<=iNatom)
            {
/**
 *b+aggiunto iflag per salvare comunque in atsurf valore
 *del + vicino (01/02)
 *non sono sicurissimo perche' non ricordo esattmente che
 *fa poi con prevmed...
*/
                if (iflag)
                {
                    dx123=xg-xn1[ia];
                    dis2=optDot(dx123,dx123)-sDelPhiPDB[ia].radius*sDelPhiPDB[ia].radius;

/**
 *dis2, and ds2min are distance**2 from center -
 *radius**2 so they can be <0
*/
                    if (dis2<ds2min1)
                    {
                        iac2=iac1;
                        ds2min2=ds2min1;
                        iac1=ia;
                        ds2min1=dis2;
                    }
                    else if (dis2<=ds2min2)
                    {
                        iac2=ia;
                        ds2min2=dis2;
                    }// if

                }
                else
                {
                    if (ast[ia]==0)
                    {
                        nnbr=nnbr+1;
                        nbra[nnbr]=ia;
                    }// if
                }// if
            }
            else
            {
                if (ia!=iaprec)
                {
                    iaprec=ia;
                    nnbr=nnbr+1;
                    nbra[nnbr]=ia;
                }// if
            }// if
        }// do

        if (iflag)
        {
            atsurf[i]=iac1;

            if (iac1*iac2==0||iac1==iac2)
            {
                cout <<"Problems in Scaling multidielectric Boundary Grid Points" << endl;
                exit(0);
            }// if

            atndx[i]=-1;

/**
 * looks like atndx is used to build Delunay surface, so excluding these bgps
*/

            dx123=xn1[iac2]-xn1[iac1];
            temp=optDot(dx123,dx123);
            temp=0.5*(ds2min2-ds2min1)/temp;
            scspos[i]=xg+(temp*dx123);
            //if(i==190) cout << "Lin Li 1: scspos "<< scspos[i] << endl;
            //scsnor[i]={0.,0.,0.};
            scsnor[i]=sgrid_temp_real;
            //cycle D500;
            continue;
        }
        else
        {
            for(ii=1; ii<=nnbr; ii++)
            {
                ia=nbra[ii];
                med=iAtomMed[ia];
                lgd=(med!=prevmed);

                if (ia>iNatom)
                {
                    iii=ia-iNatom;
                    xq=xg;

/**
 * try to find closest VdW surface, better if it is buried internal is used for the object to which surface the bgp is closer
*/
                    cout << "$$$$$$ warning: distobj is called: $$$$$$$$$$" << endl;
                   // call distobj(xq,dist,dixyz,iii,0.0,false);

                    precedenza=ia>iac&&(iac>iNatom||iac==0);
                    vicinanza=abs(dist)<abs(dmn);
                    lga=(precedenza&&(vicinanza||dist<0.))||(vicinanza&&dmn>0.);

                    if ((dist<dmn&&!lgd)||(lga&&lgd))
                    {
                        dmn=dist;
                        iac=ia;
                        prevmed=med;
                        dr=dixyz*(dist-fRadPrb[1]);
                        vnor=dixyz;
                    }// if

                    internal[iii]=(dist<0.0);
                }
                else
                {
                    dx123=xg-xn1[ia];
                    dis=sqrt(optDot(dx123,dx123) )-sDelPhiPDB[ia].radius;
                    precedenza=ia>iac||iac>iNatom;
                    vicinanza=abs(dis)<abs(dmn);
                    lga=(precedenza&&(vicinanza||dis<0.))||(vicinanza&&dmn>0.);

                    if ((dis<dmn&&!lgd)||(lga&&lgd))
                    {
                        prevmed=med;
                        dmn=dis;
                        iac=ia;
                    }// if
                }// if
            }// do

            atsurf[i]=iac;
            //cout << "i,atsurf: " << i << " " << atsurf[i] << endl;

        }// if


        if (iac==0&&iac1==0)
        {
            cout <<"no close atom or object for boundary pointeger " << i << endl;
            exit(0);
        }// if

/**
 * if iac is an object dr has alredy been calculated and HAS a DIFFERENT value!!!!!!
*/
        if (iac<=iNatom)
        {
            dr=xg-xn1[iac];
        }// if

        dsr=sqrt( optDot(dr,dr));
        out=true;

        if (radpmax>0.0)
        {
//u should have the same value as previous one
            if (iac<=iNatom)
            {
                u123=xn1[iac]+(((r0[iac]*dr)/dsr));
            }
            else
            {
                u123=xg-dr;
            }// if

            it=optCast <delphi_integer,delphi_real> ((u123-xyzo)*cbai);
            nnbr=0;

/**
 * 2011-05-19 Changed 1d to 3d array as in cube void
*/
            //liml=cbn1[it.nX+1+(lcb+1)*it.nY+(lcb+1)*(mcb+1)*it.nZ];
            //limu=cbn2[it.nX+1+(lcb+1)*it.nY+(lcb+1)*(mcb+1)*it.nZ];
            liml=cbn1[it.nX][it.nY][it.nZ];
            limu=cbn2[it.nX][it.nY][it.nZ];

            for(kk=liml; kk<=limu; kk++)
            {
                ia=cbal[kk];
                if (ia<=iNatom)
                {
                    dx123=u123-xn1[ia];
                    ds2=optDot(dx123,dx123);
                    if(ds2<rs2[ia])out=false;
                }
                else
                {
                    if (ia!=iac&&(!internal[ia-iNatom]))
                    {
//I want to know if u is within the shell sorrounding the object
                        xq=u123;

                        cout << "$$$$$$ warning: distobj is called: $$$$$$$$$$" << endl;
                        //call distobj(xq,dist,dixyz,ia-iNatom,0.,true);

                        if (dist>0.0&&dist<fRadPrb[1]-1.e-6) out=false;
                    }// if
                }// if
            }// do
        }// if

        if (out)
        {
            ncbp=ncbp+1;
            if (iac<=iNatom)
            {
                scspos[i]=xn1[iac]+(dr*(sDelPhiPDB[iac].radius/dsr));
                scsnor[i]=dr/dsr;

            }
            else
            {

                scspos[i]=(xg-(fRadPrb[1]*vnor))-dr;
                scsnor[i]=vnor;

            }// if
            //if(i==197) cout << "Lin Li 2: scspos "<< scspos[i] << endl;
            atndx[i]=iac;

            //cout << "i,atndx: " << i << " " << atndx[i] << endl;
        }
        else
        {
            atndx[i]=0;
        }// if


    }// do D500;

    //if(allocated(cbn1)) deallocate(cbn1);
    //if(allocated(cbn2)) deallocate(cbn2);
    //if(allocated(cbal)) deallocate(cbal);
    //release memory in vw2ms





//fScale the re-entrant points with respect to expos if fRadPrb = 0.0 we are done.
    if (radpmax>0.0)
    {
        iall=0;
        cba=1./grdi;

//D700:
        //do i=1,iBoundNum;
        for (i=1; i<=iBoundNum; i++)
        {
//b+++mol. con diversa eps ++++01/02+++++++++++++++
            if(atndx[i]==-1) continue;

            if (atndx[i]==0)
            {
                s123=scspos[i]*x1+x1xyz;

//2011-05-19 mn(x,y,z) and grdi were assigned values in INDVER void now coord type variable mnxyz is declared in pointers
//module and delphi_real grdi declared and thus accessible in qlog module
                xxyyzz=(s123-mnxyz)*grdi;
                jxyz=optCast <delphi_integer,delphi_real> (xxyyzz);
                jx=jxyz.nX;
                jy=jxyz.nY;
                jz=jxyz.nZ;

                dxyz=xxyyzz-optCast<delphi_real,delphi_integer>(jxyz);

                dmn1=min(dxyz.nX,min(dxyz.nY,dxyz.nZ));
                dmx=max(dxyz.nX,max(dxyz.nY,dxyz.nZ));
                dmn2=1.0-dmx;
                dcr=min(dmn1,dmn2);
                ctf=cba*(1+dcr);

                //if(i==197) cout << "Lin Li dmn1,dmn2: " << dmn1<< " " << dmn2<< endl;

                //if(i==197) cout << "Lin Li ctf,cba,dcr: " << ctf<< " " << cba<< " " << dcr << endl;

                ctf=ctf*ctf;
                iacl=0;
                rmn=100.;
                for(jjx=jx-1; jjx<=jx+1; jjx++)
                {
                    for(jjy=jy-1; jjy<=jy+1; jjy++)
                    {
                        for(jjz=jz-1; jjz<=jz+1; jjz++)
                        {
                            for(ii=iab1[jjx][jjy][jjz]; ii<=iab2[jjx][jjy][jjz]; ii++)
                            {
                                iac= icume[ii];
                                dist=optDot( (s123-expos[iac]),(s123-expos[iac]) );

                                if (dist<rmn)
                                {
                                    rmn=dist;
                                    iacl=iac;
                                }// if
                            }// do
                        }// do
                    }// do
                }// do
                //if(i==197) cout << "Lin Li 1: iacl: " << iacl << endl;
                //if(i==197) cout << "Lin Li iacl,rmn,ctf: " << iacl<< " " << rmn<< " " << ctf << endl;

                if (!(iacl>0&&rmn<ctf))
                {
                    //if(i==197) cout << "Lin Li 1: in if1: " << endl;
                    for(jxi=-2; jxi<=2; jxi++)
                    {
                        for(jyi=-2; jyi<=2; jyi++)
                        {
                            for(jzi=-2; jzi<=2; jzi++)
                            {
                                //if (outcb[jxi][jyi][jzi])

                                //if(i==197) cout << "jxi,jyi,jzi:outcb: " <<jxi<< " " << jyi << " " << jzi<< " " << outcb[jxi+2][jyi+2][jzi+2] << endl;
                                if (outcb[jxi+2][jyi+2][jzi+2]) // index of outcb is modified
                                {
                                    //if(i==197) cout << "Lin Li 1: in if2: " << endl;
                                    jjx=jx+jxi;
                                    if (jjx>=0&&jjx<=lcb1)
                                    {
                                        //if(i==197) cout << "Lin Li 1: in if3: " << endl;
                                        jjy=jy+jyi;
                                        if (jjy>=0&&jjy<=mcb1)
                                        {
                                            jjz=jz+jzi;
                                            if (jjz>=0&&jjz<=ncb1)
                                            {
                                                //if(i==197) cout << "Lin Li 1: in if4: " << endl;
                                                ////if(i==197) cout << "iab1:" << iab1[jjx][jjy][jjz] << " " << iab2[jjx][jjy][jjz] << endl;
                                                ////if(i==197) cout << "jjx,jjy,jjz: " << " " << jjx<< " " << jjy << " " << jjz << endl;

                                                for(ii=iab1[jjx][jjy][jjz]; ii<=iab2[jjx][jjy][jjz]; ii++)
                                                {

                                                    iac= icume[ii];
                                                    dist= optDot( (s123-expos[iac]),(s123-expos[iac]) );
                                                    ////if(i==197) cout << "Lin Li:dist,rmn : " << dist << " " << rmn << endl;
                                                    if (dist<rmn)
                                                    {
                                                        //if(i==197) cout << "Lin Li 1: in if5: " << endl;
                                                        rmn=dist;
                                                        iacl=iac;
                                                    }// if
                                                }// do
                                            }// if
                                        }// if
                                    }// if
                                }// if
                            }// do
                        }// do
                    }// do
                    //if(i==197) cout << "Lin Li 2: iacl: " << iacl << endl;

                    if (iacl<=0)
                    {
                        iall=iall+1;
                        for(iac=1; iac<=extot; iac++)
                        {
                            dist=optDot((s123-expos[iac]),(s123-expos[iac]));
                            if (dist<rmn)
                            {
                                rmn=dist;
                                iacl=iac;
                            }// if
                        }// do
                    }// if
                }// if
                //if(i==197) cout << "Lin Li 3: iacl: " << iacl << endl;
                dxyz=s123-expos[iacl];
                rdist=sqrt(optDot(dxyz,dxyz) );

                if (rdist==0)
                {
                    dist=0.0;
                }
                else
                {
//if inside any object  fRadPrb[2]...
                    dst=0.;
                    flag=true;

//D400:
                    //do ii=1,iNObject;

                    for(ii=1; ii<=iNObject; ii++)
                    {


                        //strtmp=dataobject[ii-1][0];
                        strtmp=dataobject_v[(ii-1)*2];
                        //read(strtmp(16:18),*)kind;
                        //kind=strtmp.substr(15,3);
                        kind = atoi(strtmp.substr(15,3).c_str());
                        //if (strtmp[1:4]!='is a'&&kind!=2)
                        if (strtmp.substr(0,4)!="is a" && kind!=2)


                        {
                            xq=s123;
                            cout << "$$$$$$ warning: distobj is called: $$$$$$$$$$" << endl;
                           // call distobj(xq,dist,dixyz,ii,0.,true);

//assuming that if the VdW pointeger is half grid space into an object that means that this belongs to an atom buried in the object
                            if (dst<-hgs)
                            {
                                dist=fRadPrb[2]/rdist;
                                flag=false;
                                break;
                            }// if
                        }// if
                    }// do D400;

                    if(flag) dist=fRadPrb[1]/rdist;
                }// if

                scspos[i]=expos[iacl]+(dxyz*dist);
                //if(i==197) cout << "Lin Li 3: scspos expos[iacl],dxyz,dist,iacl: " << scspos[i] << " " << expos[iacl] << " " <<dxyz << " " <<dist << " " << iacl<< endl;
                if (rdist>1.0e-8)
                {
                    scsnor[i]=(-dxyz)/rdist;
                }
                else
                {
                    cout <<"bdp close to arcp " << i << rdist << endl;
                }// if
            }// if
        }// do D700;
    }// if

    //if (allocated(internal)) deallocate(internal);
    if(internal.size()>0) vector <bool>().swap(internal);

//do i=1,iBoundNum
// write(2,'(i8,3f12.6,i5)') &
// &i][scspos[1][i]][scspos[1][i]][scspos[1][i]][atndx(i)
//}// do
//close (2)
//Varrebbe la pena di capire % of ... cosa significa.
//cout <<"% of boundary points contacting solvent = " << &
// &Int2Float(ncbp)/float(iBoundNum)*100.

}// void sclbp;

