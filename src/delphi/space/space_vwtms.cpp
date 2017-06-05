#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "../space/space.h"

using namespace std;

void CDelphiSpace::VdwToMs()
{
    //bool outcb[-2:2,-2:2,-2:2]; //outcb is not used in this function, thus removed

    string line;
    delphi_integer nbra[1001];
//2011-05-17 bndeps is declared in pointers module
    delphi_integer eps[7],nt;
    delphi_integer itmp[7],dim1,cont;
    delphi_integer iaprec,dim,epsdim,isign;
//    delphi_integer imap[4][6]= {0}; //use same indexes
    delphi_integer imap[5][7]= {0};
    delphi_integer kind,eps2[7];
    bool remov;
    string strtmp;
    SGrid <delphi_real> xq,dxyz,dr123,dx123,u123;
    SGrid <delphi_real> goff[7]= {0.,0.,0.},xg,x1xyz,s123,xxyyzz;
    SGrid <delphi_integer> ixyz,it123,jxyz;
    //type(int_coord), allocatable :: ibgrd_temp(:);
    //SGrid <delphi_integer> * ibgrd_temp;
    //bool out,intb;
    bool exists,flag,cycle_flag=false;
    bool nbe[7]= {false};
//2011-05-27 Declarations added due to IMPLICIT NONE
    delphi_real offf,cbln,cba,del,dis,dmn,dist=0,ds2,off,dsr,r0a;
    delphi_real x1,prbrd12,prbrd22,rsm,prbrd2;
    delphi_integer ia,i,iac,ibgp,ii,iext,iarv,iacv,imedia;
    delphi_integer iiord,ix2,iy2,iz2,ix,iy,iz;
    delphi_integer iii,iord,j,jj,jjj,jx,jy,jz,limu,liml,kk,m,k,n;
    delphi_integer mr,mpr,n2,ndv,ncav,nnn,nnbr;
    //delphi_integer nmt,nmmt,ncms;
    delphi_integer nn,ntt,n1,Nar;
    delphi_integer **** bndeps;  //bndeps is a local variable in C++
    delphi_integer ibmx=1000000; //ibmx is a local variable in C++

    if(iGrid > 300){

	ibmx=5000000;
    }


//2011-05-17 Arrays allocated by ordinary F95 allocate statement
    //allocate(ibnd(ibmx),r0(iNatom),r02(iNatom),rs2(iNatom));
    if(debug_space) cout << "############ start vdwToMs: ##############" << endl;




    if(debug_space) cout << "iAtomMed_v.size(): " << iAtomMed_v.size() << endl;
/*
    for (i=0;i<=iAtomMed_v.size()-1;i++){
        cout << i << " " << iAtomMfed[i+1] << " " << endl;

    }
    cout << endl;
*/


    iarv=0; // Lin Li
    SGrid <delphi_integer> * ibnd = new SGrid <delphi_integer> [ibmx];

    r0.assign(iNatom+1,0.0);
    r02.assign(iNatom+1,0.0);
    rs2.assign(iNatom+1,0.0);
    ast.assign(iNatom+1,0);

    //bndeps = get_pt4d <delphi_integer> (iGrid+1,iGrid+1,iGrid+1,3);
    get_pt4d <delphi_integer> (bndeps,iGrid+1,iGrid+1,iGrid+1,3);



    offf=(iGrid+1.)/2.;
    epsdim=iNatom+iNObject+2;

/**
 * imap maps from midpoint position to iepsmap entry positions
*/

//2011-05-17 Changed to array operations


    imap[1][4]=-1;
    imap[2][5]=-1;
    imap[3][6]=-1;
    imap[4][1]=1;
    imap[4][2]=2;
    imap[4][3]=3;
    imap[4][4]=1;
    imap[4][5]=2;
    imap[4][6]=3;

//2011-05-17 Changed to array operations
    //outcb=true;
    //outcb[-1:1][-1:1][-1:1]=false;
    //nbe=false;
    nbe[1]=true;
    nbe[2]=true;
    nbe[3]=true;
    nbe[4]=true;
    nbe[5]=true;

    //goff=coord(0.,0.,0.)
    off=0.5/fScale;

//2011-05-17 Changed to coord derived variable type
    goff[1].nX=off;
    goff[2].nY=off;
    goff[3].nZ=off;
    goff[4].nX=-off;
    goff[5].nY=-off;
    goff[6].nZ=-off;

    radpmax=max(fRadPrb[1],fRadPrb[2]); //transferred via qlog module;

/**
 *convertion from grid to float coordinates(can also use routine
 *gtoc)
*/
//2011-05-17 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
    x1=1.0/fScale;
    // x1xyz=cOldMid-(0.5*x1*(iGrid+1));
    x1xyz=cOldMid-double(0.5*x1*(iGrid+1));

/**
 * find extrema
*/

    if(debug_space) cout << "cMin,cMax: " << cMin << " " << cMax <<endl;

/**
 * find global extrema
*/
    //cMin= {6000.0,6000.0,6000.0};
    //cMax= {-6000,-6000,-6000};

        //cMin= {6000.,6000.,6000.};
        cMin.nX=6000.;
        cMin.nY=6000.;
        cMin.nZ=6000.;

        //cMax= {-6000.,-6000.,-6000.};
        cMax.nX=-6000.;
        cMax.nY=-6000.;
        cMax.nZ=-6000.;

    for(ii=0; ii<=iNObject-1; ii++)
    {
        cMin=optMin(cMin,sLimObject[ii].nMin);
        cMax=optMax(cMax,sLimObject[ii].nMax);

    }// do

    if(debug_space) cout << "cMin,cMax: " << cMin << " " << cMax <<endl;


/**
 * find vanderwaals boundary
*/
    n=0;
    nn=0;
    //nmt=0;
    //nmmt=0;

/**
 * NB change limits to those of the molecule. set for iepsmp NOT equal to unity
*/
    if(debug_space) cout << "LimEps: " << LimEps.nMin << LimEps.nMax << endl;

    if(debug_space) cout << "HERE are problems: "<< endl;

    //Lin Li: determine cube 2 grids smaller than limeps;
    for(k=LimEps.nMin.nZ+1; k<=LimEps.nMax.nZ-1; k++)
    {
        for(j=LimEps.nMin.nY+1; j<=LimEps.nMax.nY-1; j++)
        {
            for(i=LimEps.nMin.nX+1; i<=LimEps.nMax.nX-1; i++)
            {
/**
 * one distinguishes between internal,external,internal bgp and external bgp
*/
                iext=0;
                ibgp=0;

                //Lin Li: six neighbors: itmp is iatmmed or 0;
                itmp[1]=abs(iepsmp[k][j][i].nX)/epsdim; //right
                itmp[2]=abs(iepsmp[k][j][i].nY)/epsdim; //front
                itmp[3]=abs(iepsmp[k][j][i].nZ)/epsdim; //up
                itmp[4]=abs(iepsmp[k][j][i-1].nX)/epsdim; //left
                itmp[5]=abs(iepsmp[k][j-1][i].nY)/epsdim; //back
                itmp[6]=abs(iepsmp[k-1][j][i].nZ)/epsdim; //down

                //cout << itmp[1] << " "<< itmp[2] << " "<< itmp[3] << " "<< itmp[4] << " "<< itmp[5] << " "<< itmp[6] << endl;
                if(itmp[1]==0) iext=1; //external point
                if(itmp[1]!=itmp[6]) ibgp=1; //bnd point

                for(cont=2; cont<=6; cont++)
                {
                    if(itmp[cont]==0) iext=1; //external point
                    if(itmp[cont]!=itmp[cont-1]) ibgp=1; //bnd point
                }// do
/**
 * assignement of right values to bndeps according to the point nature from now iBoundNum is
 * the total number of internal and external boundary grid points
*/
                if (ibgp>0)
                {

                    n=n+1; //boundary bnd point
                    bndeps[i][j][k][1]=n;
                    bndeps[i][j][k][2]=iext;

                    if (iext>0) nn=nn+1; //bnd point and external point: boundary + surface
                    ibnd[n]=int_coord(i,j,k);

                }
                else
                {


                    bndeps[i][j][k][1]=0;
                    bndeps[i][j][k][2]=0;

                }// if

                {
/*
#ifdef DEBUG
                    nt=0;
//passed from mod to dim, I need the medium
                    if((iepsmp[k][j][i].nX/epsdim)>0)nt=nt+1;
                    if((iepsmp[k][j][i].nY/epsdim)>0)nt=nt+1;
                    if((iepsmp[k][j][i].nZ/epsdim)>0)nt=nt+1;
                    if((iepsmp[k][j][i-1].nX/epsdim)>0)nt=nt+1;
                    if((iepsmp[k][j-1][i].nY/epsdim)>0)nt=nt+1;
                    if((iepsmp[k-1][j][i].nZ/epsdim)>0)nt=nt+1;

                    if (nbe[nt] != ((iext==1)&&(ibgp==1)))
                    {
                        cout <<"PROBLEMS1 " << i << j << k << endl;

//2011-05-17 converts grid to float coordinates
//using operations on coord and int_coord type
//variables defined in module
//operators_on_coordinates
                        itemp=Int2Float(int_coord(i,j,k));
                        rtemp=((itemp-offf)/fScale)+cOldMid;

                        cout <<rtemp << endl;
                        //cout <<iepsmp[k][j][i] << iepsmp[i-1][j][k].nX << iepsmp[i << j-1 << k].nY << iepsmp[i << j << k-1].nZ << endl;
                    }// if
#endif //DEBUG
 */
                };
            }// do
        }// do
    }// do
    if(debug_space) cout << "Problems up to here: "<< endl;

    iBoundNum=n;
    iBoundNumsurf=nn;
    nn=0;
    if(debug_space) cout << "iboundNum,iBoundNumsurf: " << iBoundNum << " " << iBoundNumsurf << endl;


#ifdef VERBOSE
    cout <<"boundary points facing continuum solvent= " << iBoundNumsurf << endl;
    cout <<"total number of boundary points before elab.= " << iBoundNum << endl;
#endif //VERBOSE


    if (iBoundNum>ibmx)
    {
        cout <<"iBoundNum= " << iBoundNum << " is greater than ibmx = " << ibmx << endl;
        cout <<"increase ibmx in vwtms.f" << endl;
        exit(0);
    }// if

    {
        /*
            if (bDebug)
            {
                open(52,file='bgpPrima.log',form='formatted');

        //grid coordinates
                write(52,*)iGrid,iBoundNum;
                for(iiii=1; iiii<=iBoundNum; iiii++)
                {
                    write(52,*)iiii,ibnd(iiii);
                }// do
                close (52);
            }// if //end bDebug;
        */

    };


//2011-05-17 Arrays allocated by ordinary F95 allocate statement


    if (radpmax<1.e-6) //prob radius = 0
    {
        //allocate(ibgrd(iBoundNum));
        //ibgrd_v.assign(iBoundNum+1, {0,0,0});
        ibgrd_v.assign(iBoundNum+1, sgrid_temp_int);
        ibgrd=&ibgrd_v[-1];

//2011-05-17 Changed to array operations, but keeping some assignment in a cycle due to array size mismatch
        //ast=0; already initialized.
        for(i=1; i<=iBoundNum; i++)
        {
            ibgrd[i]=ibnd[i];
        }// do
    }
    else
    {
        //cout << "bOnlyMol,fRadPrb[1] != fRadPrb[2]: " << bOnlyMol << " " <<fRadPrb[1] << " " << fRadPrb[2] << endl;
        if (!bOnlyMol && fRadPrb[1] != fRadPrb[2]) // !bOnlyMol: means object exist, so delete
        {

            /*
            for(i=1; i<=iNatom; i++)
            {

                xq=xn1[i];
                r0a=sDelPhiPDB[i].radius+fRadPrb[1];

                for(ii=1; ii<=iNObject; ii++)
                {
                    strtmp=dataobject[ii][1];
                    read(strtmp(16:18),*)kind;
                    if (strtmp(1:4)!='is a'&&kind!=2)
                    {
                        if (((xq-sDelPhiPDB[i].radius).vandlt.&sLimObject[ii].nMax)&&((xq+&sDelPhiPDB[i].radius).vandgt.sLimObject[ii].nMin))
                        {

                            call distobj(xq,dist,dxyz,ii,0.,true);

            //only full buried atoms have the proberadius
            //changed to fRadPrb[2]
                            if (dist<-sDelPhiPDB[i].radius) r0a=sDelPhiPDB[i].radius+fRadPrb[2];
                        }// if
                    }// if
                }// do

                r0[i]=r0a;
                r02[i]=r0a*r0a;
                rs2[i]=0.99999*r02[i];
            }// do
            */
        }
        else
        {
            for(i=1; i<=iNatom; i++)
            {
                r0a=sDelPhiPDB[i].radius+fRadPrb[1];
                r0[i]=r0a; //atom radius+ probe radius
                r02[i]=r0a*r0a;
                rs2[i]=0.99999*r02[i];
            }// do
        }// if



        /*
                if(iacs) write(6,'(a21)') " opening surface file";

                if(iacs) open(40,file='hsurf2.dat');
        */


/**
 *make a list of accessible points..,expos. all scaling of grid points will be done to thses points..
*/
        prbrd12=fRadPrb[1]*fRadPrb[1]; //fRadPrb[1]: probe radius;fRadPrb[2]=0, no idea.
        prbrd22=fRadPrb[2]*fRadPrb[2];

/**
 * calculate an identity for this conformation
*/
        rsm=0.0;
        if (bOnlyMol)
        {
            for(i=1; i<=iNatom; i++)
            {
                //rsm=rsm+sDelPhiPDB[i].radius*sum(abs(xn1[i]));
                rsm=rsm+sDelPhiPDB[i].radius*optSum(optABS(xn1[i])); //Lin: is it a mistake? rsm? If rms is just for ARCDAT....
                //cout << "i,rsm,radius[i],xn1[i]: "<< i << " " << rsm << " " << sDelPhiPDB[i].radius << " " << xn1[i-1] << " " << endl;

            }// do
        }// if
        if(debug_space) cout << "rsm: " << rsm << endl;
        //cout << "extot: " << extot << endl;
        //inquire(file='ARCDAT',exist=exists);
        flag=true;
        /*
        if (exists&&bOnlyMol) // ARCDAT related, forget it now...
        {

            open(1,file='ARCDAT',form='unformatted',status='old',iostat=ios);
            read(1)natm,radp,nacc,rsm1;

            //maybe it could be improved taking into account objects
            //(Walter 1999)
            if (natm==iNatom&&radp==fRadPrb[1]&& abs(rsm1-rsm)<1.0e-8)
            {
                if(bVerbose)cout <<"reading accessible surface arcs data from file ARCDAT" << endl;
                extot=nacc;
                allocate(expos(extot));

            //2011-05-19 Subroutine arcio is too short and called
            //only from this void, thus no meaning to be
            //separate piece of code
                read(1)ast;
                read(1) expos;
                if(bVerbose) cout <<"no. of arc points read = " << nacc << endl;
                close (1);
                flag=false;

            }
            else
            {
            //close (1);
            }// if

        }// if
        */
        //if (flag) //if no ARCDAT
        if(true)
        {
            //cout << "going to sas..." << endl;
            //cout << "before sas: extot: " << extot << endl;
            sas(); //Lin Li: Is this always called?
            //cout << "after sas: extot: " << extot << endl;

            //writing ARCDAT: //Lin Li: forget this now....
            /*
            if (extot > 0 && bOnlyMol)
            {

                cout <<"writing accessible surface arcs data to ARCDAT" << endl;
                open(1,file='ARCDAT',form='unformatted',iostat=ios);

                if (ios!=0)
                {
                    cout << "error opening ARCDAT file" << endl;
                    exit (0);
                }// if

                write(1)iNatom,fRadPrb[1],extot,rsm;

                if (iNatom==0)
                {
                    write(1)0;
                    write(1) expos;
                }
                else
                {
                    write(1)ast;
                    write(1)expos;
                }// if

                close(1);

            }// if
            */

            /*
            if (bDebug)
            {
                open(52,file='Vertices.txt',form='formatted');
                for(iiii=1; iiii<=extot; iiii++)
                {
                    write (52,*) expos(iiii);
                }// do
                close (52);
            }// if //end bDebug;
            */
        }// if


        del=1./fScale;
        del=max(del,radpmax);
        cbln=fRMax+del;

        cubedata(2.0,cbln);

        dim=(lcb+1)*(mcb+1)*(ncb+1);

//2011-05-17 Array allocation is done by ordinary F95 ALLOCATE
//statement Array allocated as 3D as in the cube void
        //allocate(cbn1_v(dim),cbn2_v(dim));
        //cbn1_v = new int[dim];
        //cbn2_v = new int[dim];

        cbn1_v.assign(dim+1,0);
        cbn2_v.assign(dim+1,0);
        dim1=27;
        if ((iNObject-numbmol)>0) dim1=max(dim,27);
//cout << "### flag0: " << endl;
        //allocate(cbal(dim1*(iNatom+iNObject-numbmol)));
        //cbal = new int[dim1*(iNatom+iNObject-numbmol)];
        cbal.assign(dim1*(iNatom+iNObject-numbmol)+1,0);

        cube();

        //link the accessible points into iab1 and iab2
        indverdata(radpmax);

        cba=1./grdi;
        nnn=(lcb1+1)*(mcb1+1)*(ncb1+1);
        //allocate(iab1(nnn),iab2(nnn));
        //allocate(icume(extot));
        //iab1 = new int[nnn];
        //iab2 = new int[nnn];
        //icume = new int[extot];
        if(debug_space) cout << "nnn: " << nnn << endl;
        //iab1.assign(nnn+1,0); //Lin Li: make them to be pointers
        //iab2.assign(nnn+1,0);
        icume.assign(extot+1,0);

#ifdef VERBOSE
        cout << " grid for indexing accessible points =  " << cba << endl;
#endif // VERBOSE
        //cout << "Lin Li: extot: " << extot << endl;
        indver(extot);

//write out surface data


        /*
                if (iacs)
                {
                    line= ' ';
                    line(1:6)='ATOM ';
                    line(14:14)='O';
                    line(18:20)='SP ';

                    for(i=1; i<=extot; i++)
                    {
                        xg=expos[i];
                        iv=1;

        //2011-05-24 Subroutine watput is small and called only
        //once, thus the code transfered into calling void
                        write(line(7:11),'(i5)')i;
                        write(line(24:26),'(i3)')iv;

        //2011-05-24 In original coordinates from array xo were
        //written into th line.
        //Logic of this piece, however, requires coordinates
        //from xg to be written
                        write(line(31:54),'(3(f8.3))')xg;
                        write(40,'(a80)') line;
                    }// do
                    close (40);
                }// if


        */
//now start the expansion
//m1= the number of boundary points removed from list
        ncav=0;
        n1=1;
        n2=iBoundNum;

/**
 * m= number of new boundary elements..
*/
        mpr=100000;
        ndv=0;

        if(debug_space) cout << "### start while loop: ###" << endl;
//D100:
        while(true)
        {
            m=0;
            mr=0;

            for(i=n1; i<=n2; i++)
            {
                ixyz=ibnd[i];
                ix=ixyz.nX;
                iy=ixyz.nY;
                iz=ixyz.nZ;

//considering both internal and external b.g.p.
                if (bndeps[ix][iy][iz][1]!=0)
                {

/**
 * still has to be considered what is external and what internal!!!!!WWW
 * remov is true if it is an internal midpoint close to an interface where a molecule is present
 * (expansion has to take place also in objects)
*/
                    remov=false;

//tengo il mod perche' deve pr}//ere solo punti in atomi
                    eps[1]=(iepsmp[iz][iy][ix].nX%epsdim);
                    eps[2]=(iepsmp[iz][iy][ix].nY%epsdim);
                    eps[3]=(iepsmp[iz][iy][ix].nZ%epsdim);
                    eps[4]=(iepsmp[iz][iy][ix-1].nX%epsdim);
                    eps[5]=(iepsmp[iz][iy-1][ix].nY%epsdim);
                    eps[6]=(iepsmp[iz-1][iy][ix].nZ%epsdim);

                    remov=((eps[1]>1&&eps[1]<=iNatom+1)||(eps[2]>1&&eps[2]<=iNatom+1));
                    remov=((eps[3]>1&&eps[3]<=iNatom+1)||(eps[4]>1&&eps[4]<=iNatom+1))||remov;
                    remov=((eps[5]>1&&eps[5]<=iNatom+1)||(eps[6]>1&&eps[6]<=iNatom+1))||remov;

//da farsi solo se pores eps2 contiene il mezzo
                    eps2[1]=(iepsmp[iz][iy][ix].nX/epsdim);
                    eps2[2]=(iepsmp[iz][iy][ix].nY/epsdim);
                    eps2[3]=(iepsmp[iz][iy][ix].nZ/epsdim);
                    eps2[4]=(iepsmp[iz][iy][ix-1].nX/epsdim);
                    eps2[5]=(iepsmp[iz][iy-1][ix].nY/epsdim);
                    eps2[6]=(iepsmp[iz-1][iy][ix].nZ/epsdim);

/**cWWW there is still an issue in case there are both molecules and objects: since parent object of
 * reentrant points is only known in sclbp, filling reentrant regions due to molecules in objects might fail
*/
                    remov=remov&&(numbmol>0);

                    xg=(x1*optCast<delphi_real,delphi_integer>(ixyz))+x1xyz;

        //if(i==119) cout<< "flag 2:" << endl;
//D200:


                    for (j=1; j<=6; j++)
                    {
                        //cout << "i,j: " << i << " " <<  j << endl;
//essere in poro ==> eps2=0 and eps >0
            //if(i==119) cout << "j,T or F? " << j << " " << (eps[j]==0||(remov&&eps[j]>iNatom+1)||(eps2[j]==0&&eps[j]>0)) << endl;
                        if (eps[j]==0||(remov&&eps[j]>iNatom+1)||(eps2[j]==0&&eps[j]>0))
                        {
                                    //if(i==119) cout<< "flag 6:" << endl;
                            prbrd2=prbrd22;
                            if (eps[j]==0||eps2[j]==0) prbrd2=prbrd12;

//add midpoint offset to grid point..
                            s123=xg+goff[j];
                            //cout << "xg,goff[j]: " << xg << " " << goff[j] << endl;
//determine if this virgin midpoint is in or out
//2011-05-18 mn(x,y,z) and grdi were assigned values in INDVER void now coord type variable mnxyz is declared
//in pointers module and float grdi declared and thus accessible in qlog module
                            xxyyzz=(s123-mnxyz)*grdi;
                            jxyz=optCast<delphi_integer,delphi_real>(xxyyzz);
                            jx=jxyz.nX;
                            jy=jxyz.nY;
                            jz=jxyz.nZ;
                            //cout << "s123,jxyz: " << s123 << " " << jxyz << endl;
//2011-05-18 Indexes lcb1, mcb1, ncb1 are transfred via qlog module and are set in INDVER void
                            if (optORLE(jxyz,0)||optORGE(jxyz,lmncb1))
                            {
                                cout <<"midpoint out of cube" << endl;
                                //write(6,'(2i5,3f8.3,3i6,3i8)')i,j,xxyyzz,jxyz,lmncb1;
                                cout <<iepsmp[ iz][ iy ][ix ].nX << endl;
                                cout <<iepsmp[ iz][ iy ][ix ].nY << endl;
                                cout <<iepsmp[ iz][ iy ][ix ].nZ << endl;
                                cout <<iepsmp[ iz][ iy ][ix-1 ].nX << endl;
                                cout <<iepsmp[ iz][ iy-1 ][ix ].nY << endl;
                                cout <<iepsmp[ iz-1][ iy ][ix ].nZ << endl;
                            }// if
        //if(i==119) cout<< "flag 7:" << endl;
                            dmn=1000.;
                            iacv=0;


//2011-05-18 Repeating piece of the code now is in the separate file


                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //cout << "dist,prbrd2,dmn,iarv: " << dist << " " << prbrd2 << " " << dmn << " " << iarv << endl;

                            //-1,0,0
                            jx=jx-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,0,0
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

        //if(i==119) cout<< "flag 7.5.1:" << endl;
                            //0,-1,0
                            jx=jx-1;
                            jy=jy-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,1,0
                            jy=jy+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,0,-1
                            jy=jy-1;
                            jz=jz-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,0,1
                            jz=jz+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //nn=2
                            //1,0,1
                            jx=jx+1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,0,1
                            jx=jx-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,1,1
                            jx=jx+1;
                            jy=jy+1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,-1,1
                            jy=jy-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,-1,0
                            jz=jz-1;
                            jx=jx-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,-1,0
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,1,0
                            jy=jy+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,1,0
                            jx=jx-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,0,-1
                            jz=jz-1;
                            jy=jy-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,0,-1
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,1,-1
                            jx=jx-1;
                            jy=jy+1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //0,-1,-1
                            jy=jy-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //nn=3
                            //-1,-1,-1
                            jx=jx-1;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,-1,-1
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,1,-1
                            jy=jy+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,1,-1
                            jx=jx-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,1,1
                            jz=jz+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,1,1
                            jx=jx+2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //1,-1,1
                            jy=jy-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //-1,-1,1
                            jx=jx-2;
                            VdwToMs_piece(cycle_flag, j, jx, jy, jz, prbrd2, s123, eps, eps2, iacv ,dmn);
                            if(cycle_flag) continue;

                            //cout << "dist,prbrd2,dmn,iarv: " << dist << " " << prbrd2 << " " << dmn << " " << iarv << endl;


//it might be in the contact region findthe closest atom surface

                            it123=optCast<delphi_integer,delphi_real>((s123-xyzo)*cbai);
        //if(i==119) cout<< "flag 7.5:" << endl;
                            dmn=100.;
                            iac=0;
                            nnbr=0;
                            lmncb=int_coord(lcb,mcb,ncb);

                            if(optORLT(it123,0)||optORGT(it123,lmncb))
                            {
//if the bgp is outside the cube,probably it is due to some object

                                for(ii=iNObject; ii<=1; ii--)
                                {
                                    //strtmp=dataobject[ii][1];
                                    strtmp=dataobject_v[(ii-1)*2];
                                    //read(strtmp(16:18),*)kind;
                                    kind = atoi(strtmp.substr(15,3).c_str());
                                    //if(strtmp(1:4)!='is a'&&kind!=2)
                                    if (strtmp.substr(0,4)!="is a" && kind!=2)
                                    {
                                        if ( optANDLE(s123,(sLimObject[ii].nMax+x1) ) && optANDGT(s123,(sLimObject[ii].nMin-x1)) )

                                        {
                                            nnbr=nnbr+1;
                                            nbra[nnbr]=ii+iNatom;
                                            liml=0;
                                            limu=0;
                                        }// if
                                    }// if
                                }// do

                                if(liml!=0||limu!=0) cout <<"a bgp close to nothing" << endl;

                            }
                            else
                            {
//2011-05-19 Changed 1d array to 3d array as /in cube void
                                //liml=cbn1_v[it123.nX+1+(lcb+1)*it123.nY+(lcb+1)*(mcb+1)*it123.nZ];
                                //limu=cbn2_v[it123.nX+1+(lcb+1)*it123.nY+(lcb+1)*(mcb+1)*it123.nZ];
                                liml=cbn1[it123.nX][it123.nY][it123.nZ];
                                limu=cbn2[it123.nX][it123.nY][it123.nZ];

                            }// if

                            iaprec=0;
                            //cout << "liml,limu: " << liml << " " << limu << endl;

//DOKK:

                            //if(liml==3911) cout << "liml;limu: " << liml << " " << limu << endl;
                            for( kk=liml; kk<=limu; kk++)
                            {
                                ia=cbal[kk];
                                //if(liml==3911) cout << "kk,ia: " << kk << " " << ia << endl;
                                //if (ia==0) cout << "ia: " << ia << endl;
                                if (ia==0) cout <<"problems with cube" << endl;
                                //if (ia==0) cout << "liml;limu: " << liml << " " << limu << endl;

                                if (ia<=iNatom&&ia>0)
                                {
                                    if (ast[ia]==0)
                                    {
                                        nnbr=nnbr+1;
                                        nbra[nnbr]=ia;
                                    }// if
                                }
                                else
                                {
                                    if (ia!=iaprec&&eps[j]==0)
                                    {
//different from sclbp, I have to consider the object only if the
//midpoint is external O SE C'E' UN PORO E VEDERE IN SEGUITO
                                        iaprec=ia;

//assuming any object is not buried
                                        nnbr=nnbr+1;
                                        nbra[nnbr]=ia;
                                    }// if
                                }// if
                            }// do DOKK;
                            //cout << "nnbr,ia: " << nnbr << " " << ia << endl;



                            //cout << "nbra0-7: " << nbra[0] << " " << nbra[1] << " " << nbra[2] << " " << nbra[3] << " " << nbra[4] << " " << nbra[5] << " " << nbra[6] << " " << nbra[7] << endl;
//DOII:
                            for( ii=1; ii<=nnbr; ii++)
                            {
                                ia=nbra[ii];

                                if (ia>iNatom)
                                {
                                    iii=ia-iNatom;

                                    xq=s123;
                                    cout << "call distobj: ii, iii: " << ii << " " << iii << endl;
/**
 *check if bgp is inside VdW of object
*/
                                    //call distobj(xq,dist,dxyz,iii,0.0,false);

//an object can compete with an atom for a midpoint only if this midpoint is out of the object itself
                                    if (dist>=0.&&dist<dmn)
                                    {
                                        dmn=dist;
                                        iac=ia;
                                        dr123=(-dxyz)*(fRadPrb[1]-dist);
                                    }// if
                                }
                                else
                                {
                                    dx123=s123-xn1[ia];
                                    ds2=optDot(dx123,dx123);
                                    dis=sqrt(ds2)-sDelPhiPDB[ia].radius;

//dis= distance to atom surface
                                    if (dis<dmn)
                                    {
                                        dmn=dis;
                                        iac=ia;
                                    }// if
                                    //cout << "dmn,iac,ia,iNatom,ii: " << dmn << " " << iac << " " << ia << " " << iNatom  << " " << ii << endl;
                                }// if
                            }// do DOII;


                            if (iac==0)
                            {
/*
#ifdef DEBUG

                                cout <<"bgp:" << i << " might be a cavity point" << ix << iy << iz << endl;
                                cout <<"midpoint" << j << " in position [A]" << s1 << s2 << s3 << endl;
                                cout <<"it1:" << it1 << " it2:" << it2 << " it3:" << it3 << endl;

#endif //DEBUG
*/
                                ncav=ncav+1;

//possibly a cavity point
                            }
                            else
                            {
/**
 *check to see if it is in the contact region of that atom or object by
 *projecting it to the atom's acc surface and checking against the acc volumes of nearby atoms
*/
                                if (iac<=iNatom)
                                {
                                    dr123=s123-xn1[iac];
                                    dsr=sqrt(optDot(dr123,dr123));
                                    u123=xn1[iac]+((r0[iac]*dr123)/dsr);
                                }
                                else
                                {
                                    u123=s123-dr123;
                                }// if

                                it123=optCast<delphi_integer,delphi_real>((u123-xyzo)*cbai);

                                //liml=cbn1_v[it123.nX+1+(lcb+1)*it123.nY+(lcb+1)*(mcb+1)*it123.nZ];
                                //limu=cbn2_v[it123.nX+1+(lcb+1)*it123.nY+(lcb+1)*(mcb+1)*it123.nZ];
                                liml=cbn1[it123.nX][it123.nY][it123.nZ];
                                limu=cbn2[it123.nX][it123.nY][it123.nZ];


//DLIM:
                                for( kk=liml; kk<=limu; kk++)
                                {
                                    ia=cbal[kk];
                                    if (ia<=iNatom)
                                    {
                                        dx123=u123-xn1[ia];
                                        ds2=optDot(dx123,dx123);
                                        flag=true;

                                        if (ds2<rs2[ia])
                                        {
                                            flag=false;
                                            break;
                                        }// if
                                    }
                                    else
                                    {
                                        if (ia!=iac&&eps[j]==0)
                                        {
                                            xq=u123;
                                            cout << "call distobj: ia:" << ia << endl;
                                            //call distobj(xq,dist,dxyz,ia-iNatom,fRadPrb[1],true);

                                            flag=true;
                                            if (dist<-1.e-6)
                                            {
                                                flag=false;
                                                break;
                                            }// if

//oriented distance from ext}//ed
//object surface if negative => reentrant region
                                        }// if
                                    }// if
                                }// do DLIM;

/*
 * it is in the contact region. flag the midpoint so it is not checked again iac is atom number...NOT increased by 1
*/

//2011-05-18 To get rid of goto 201
//statements above
                                if (flag)
                                {
                                    eps[j]=-iac;
                                    eps2[j]=-iac;
                                    continue;
                                }// if
                            }// if

                            eps[j]=1; //eps = 1 means cavity or reentrant;

//remap iepsmp
                            if (iac==0)
                            {
/**
 *this is an assumption, still to deeply understand meaning of cavity here and to improve this choice!WWW
*/
                                if (ia>0)
                                {
                                    imedia=iAtomMed[ia];
                                }
                                else
                                {
                                    cout <<"assigning arbitrary epsilon in cavity" << endl;
                                    imedia=iAtomMed[1];
                                }// if
                            }
                            else
                            {
                                imedia=iAtomMed[iac];
                            }// if

                    //cout << "i,j,imap[4][j]: " << i << " " << j << " " << imap[4][j] << endl;
                            switch (imap[4][j])
                            {
                            case 1:
                                iepsmp[iz+imap[3][j]][iy+imap[2][j]][ix+imap[1][j]].nX=eps[j]+imedia*epsdim;
                                //cout << "i,j,ix+imap[1][j],iy+imap[2][j],iz+imap[3][j],iepsmp: " << i << " " << j << " " << ix+imap[1][j] << " " <<iy+imap[2][j] << " " << iz+imap[3][j]<< " "  << iepsmp[iz+imap[3][j]][iy+imap[2][j]][ix+imap[1][j]] << endl;
                                break;
                            case 2:
                                iepsmp[iz+imap[3][j]][iy+imap[2][j]][ix+imap[1][j]].nY=eps[j]+imedia*epsdim;
                                break;
                            case 3:
                                iepsmp[iz+imap[3][j]][iy+imap[2][j]][ix+imap[1][j]].nZ=eps[j]+imedia*epsdim;
                                break;
                            default:
                                cout <<"????? flag1" << endl;
                            }// select;

                            eps2[j]=imedia;

/**
 * not assigning the owner but only the medium, the former will be assigned in the fScale routine
 * check to see if the nearest neighbour status has been changed..
*/
                            ix2=ix;
                            iy2=iy;
                            iz2=iz;

/**
 * if the nearest neighbour is a box boundary point { skip this since box boundary
 * points can not also be dielctric boundary points
 * 2011-05-18 Multiple IFs replaced by SELECT CASE
*/
                            switch (j)
                            {
                            case 1:
                                ix2=ix+1;
                                if(ix2==iGrid) continue;
                                break;
                            case 2:
                                iy2=iy+1;
                                if(iy2==iGrid) continue;
                                break;
                            case 3:
                                iz2=iz+1;
                                if(iz2==iGrid) continue;
                                break;
                            case 4:
                                ix2=ix-1;
                                if(ix2==1) continue;
                                break;
                            case 5:
                                iy2=iy-1;
                                if(iy2==1) continue;
                                break;
                            case 6:
                                iz2=iz-1;
                                if(iz2==1) continue;
                            }// select;

/**
 * once again one distinguishes between internal,external,internal bgp and external bgp
*/
                            iext=0;
                            ibgp=0;

//2011-05-18 Changed to i nt_coord derived type
                            itmp[1]=abs(iepsmp[iz2][iy2][ix2].nX)/epsdim;
                            itmp[2]=abs(iepsmp[iz2][iy2][ix2].nY)/epsdim;
                            itmp[3]=abs(iepsmp[iz2][iy2][ix2].nZ)/epsdim;
                            itmp[4]=abs(iepsmp[iz2][iy2][ix2-1].nX)/epsdim;
                            itmp[5]=abs(iepsmp[iz2][iy2-1][ix2].nY)/epsdim;
                            itmp[6]=abs(iepsmp[iz2-1][iy2][ix2].nZ)/epsdim;

                            if(itmp[1]==0) iext=1;
                            if(itmp[1]!=itmp[6]) ibgp=1;

                            for(cont=2; cont<=6; cont++)
                            {
                                if(itmp[cont]==0) iext=1;
                                if(itmp[cont]!=itmp[cont-1]) ibgp=1;
                            }// do

#ifdef DEBUG
                            nt=0;
                            if((iepsmp[iz2][iy2][ix2].nX/epsdim)>0) nt=nt+1;
                            if((iepsmp[iz2][iy2][ix2].nY/epsdim)>0) nt=nt+1;
                            if((iepsmp[iz2][iy2][ix2].nZ/epsdim)>0) nt=nt+1;
                            if((iepsmp[iz2][iy2][ix2-1].nX/epsdim)>0) nt=nt+1;
                            if((iepsmp[iz2][iy2-1][ix2].nY/epsdim)>0) nt=nt+1;
                            if((iepsmp[iz2-1][iy2][ix2].nZ/epsdim)>0) nt=nt+1;
                            if(nbe[nt]!=(ibgp==1&&iext==1))
                            {
                                cout <<"PROBLEMS3" << ix2 << iy2 << iz2 << endl;
                            }// if
#endif //DEBUG
                            if ((ibgp==0)&&(bndeps[ix2][iy2][iz2][1]!=0))
                            {
//reset bndeps for that point (i.e. remove
//bgp flag).
//a bgp become internal
                                iBoundNumsurf=iBoundNumsurf-bndeps[ix2][iy2][iz2][2];
                                bndeps[ix2][iy2][iz2][1]=0;
                                bndeps[ix2][iy2][iz2][2]=0;
                                mr=mr+1;
                            }
                            else
                            {
                                if (ibgp==1&&iext==0&&bndeps[ix2][iy2][iz2][2]==1)
                                {
//an ext bgp is turned into an internal bgp
                                    iBoundNumsurf=iBoundNumsurf-1;
                                    bndeps[ix2][iy2][iz2][2]=0;
                                }// if
                            }// if

                            if (ibgp==1&&bndeps[ix2][iy2][iz2][1]==0)
                            {
//create a new boundary point..
                                m=m+1;
                                bndeps[ix2][iy2][iz2][1]=n2+m;
				if(n2+m > ibmx)
				{
				   cout << "this case is too big, ibmx need to be encreased." << endl; // Lin Li
				   exit(0);
				}
                                ibnd[n2+m]=int_coord(ix2,iy2,iz2);
                                bndeps[ix2][iy2][iz2][2]=iext;
                                iBoundNumsurf=iBoundNumsurf+bndeps[ix2][iy2][iz2][2];
                            }// if
                        }// if

//now jump to the next midpoint of the same grid point

                    }// do D200;

//remap iepsmp in case there have been changes.. (that is some 0's became -1's)
//in other words: midpoint must remain external to objects
                    for(jj=1; jj<=6; jj++)

                    {
//in this way I can deal with eps[jj]<0
                        isign=1;

//iord=owner of the midpoint jj before change or after eps=1
                        iiord=optComp(iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]],imap[4][jj]);
                        iord=(iiord%epsdim);
                        //if(i==562){
                        //    cout << "ix+imap[1][jj],iy+imap[2][jj],iz+imap[3][jj]" << ix+imap[1][jj] << " " << iy+imap[2][jj] << " " << iz+imap[3][jj] << endl;
                        //    cout << "jj,iiord,iepsmp,imap[4][jj]: " << jj << " " << iiord << " " << iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]] << " " << imap[4][jj] << endl;
                        //}
                        //if (iord<0){
                        //    cout << "i,jj,iord,iiord,epsdim: " << i << " " << jj << " " << iord << " " << iiord << " " << epsdim << endl;
                        //}


                        //cout << "1: jj,imap[4][jj]:" << jj << " " << imap[4][jj] << endl;
//the last changed for sure has not iord<0 there can be iord<0 due to nearest neighbors already changed
                        if (iord<0) continue;
                        //cout << "2: jj,imap[4][jj]:" << jj << " " << imap[4][jj] << endl;
//if it has changed at previous step, dont change anymore
                        if (eps[jj]<0)
                        {
                            isign=-1;
                            if (iord==0) iord=1;
                        }// if

                        jjj=abs(iiord)/epsdim;
                        //cout << "i,j,jj,iord, jjj: " << i << " " << j << " " << jj << " " << iord << " " << jjj << endl;

                        switch (imap[4][jj])
                        {
                        case 1:
                            //cout << "case 1:" << endl;
                            iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nX=isign*(iord+jjj*epsdim);
                            //cout << "i,jj,ix+imap[1][jj],iy+imap[2][jj],iz+imap[3][jj],iepsmp: " << i << " " << jj << " " << ix+imap[1][jj] << " " <<iy+imap[2][jj] << " " << iz+imap[3][jj]<< " "  << iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]] << endl;
                                break;
                        case 2:
                            //cout << "case 2:" << endl;
                            iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nY=isign*(iord+jjj*epsdim);
                            break;
                        case 3:
                            //cout << "case 3:" << endl;
                            iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nZ=isign*(iord+jjj*epsdim);
                            break;
                        default:
                            if(debug_space) cout <<"????? flag2" << endl;
                        }// select;

//left iord definition with mod since if it is <> 0, it keeps its identity
                    }// do

//at this point one still can trace what changes have been made check to see if this is still a boundary point
//once again one distinguishes between internal,external,internal bgp and external bgp
                    iext=0;
                    ibgp=0;

                    itmp[1]=abs(iepsmp[iz][iy][ix].nX)/epsdim;
                    itmp[2]=abs(iepsmp[iz][iy][ix].nY)/epsdim;
                    itmp[3]=abs(iepsmp[iz][iy][ix].nZ)/epsdim;
                    itmp[4]=abs(iepsmp[iz][iy][ix-1].nX)/epsdim;
                    itmp[5]=abs(iepsmp[iz][iy-1][ix].nY)/epsdim;
                    itmp[6]=abs(iepsmp[iz-1][iy][ix].nZ)/epsdim;

                    if(itmp[1]==0) iext=1;
                    if(itmp[1]!=itmp[6]) ibgp=1;

                    for(cont=2; cont<=6; cont++)
                    {
                        if(itmp[cont]==0) iext=1;
                        if(itmp[cont]!=itmp[cont-1]) ibgp=1;
                    }// do

#ifdef DEBUG

                    nt=0;
                    if((abs(iepsmp[iz][iy][ix].nX)/epsdim)>0)nt=nt+1;
                    if((abs(iepsmp[iz][iy][ix].nY)/epsdim)>0)nt=nt+1;
                    if((abs(iepsmp[iz][iy][ix].nZ)/epsdim)>0)nt=nt+1;
                    if((abs(iepsmp[iz][iy][ix-1].nX)/epsdim)>0) nt=nt+1;
                    if((abs(iepsmp[iz][iy-1][ix].nY)/epsdim)>0) nt=nt+1;
                    if((abs(iepsmp[iz-1][iy][ix].nZ)/epsdim)>0) nt=nt+1;

                    if (nbe[nt] != (ibgp==1&&iext==1))
                    {
                        cout <<"PROBLEMS4" << ix << iy << iz << endl;
                        cout <<"epsdim=" << epsdim << "ibgp=" << ibgp << "iext=" << iext << endl;
                        cout <<"itmp" << itmp << endl;
                        /*
                        cout <<iepsmp[ iy ][x ][ix << iy << iz<< iz].nX << endl;
                        cout <<iepsmp[ iy ][x ][ix << iy << iz<< iz].nY << endl;
                        cout <<iepsmp[ iy ][x ][ix << iy << iz<< iz].nZ << endl;
                        cout <<iepsmp[<< i][x-][ix-1 << iy << izy << iz].nX << endl;
                        cout <<iepsmp[ iy-][x ][ix << iy-1 << iz1 << iz].nY << endl;
                        cout <<iepsmp[ iy ][x ][ix << iy << iz-1<< iz-1].nZ << endl;
                        */
                    }// if

#endif //DEBUG
//if not now a boundary element change bndeps

                    if ((iext==0)||(ibgp==0))
                    {
                        iBoundNumsurf=iBoundNumsurf-bndeps[ix][iy][iz][2];
                        if(ibgp==1) bndeps[ix][iy][iz][2]=iext;

                        if (ibgp==0)
                        {
                            bndeps[ix][iy][iz][1]=0;
                            bndeps[ix][iy][iz][2]=0;
                            mr=mr+1;
                            if(iext==1)cout <<"//!!born a new external point!!!" << endl;
                        }// if
                    }// if

                }// if //if end for whether bndeps is nonzero;

            }// do //next boundary point FINISH;

            n1=n2+1;
            n2=n2+m;
#ifdef VORBOSE
            cout <<"bgp added m=" << m << " bgp removed mr =" << mr << endl;
#endif // VERBOSE

            cout <<"bgp added m=" << m << " bgp removed mr =" << mr << endl;
            if (m>mpr)
            {
                ndv=ndv+1;
                if (ndv>20)   //Lin Li: the value used to be 2,;
                {
// sometimes not enough
                    cout <<"surface iteration did not converge" << endl;
                    exit (0);
                }// if
            }
            else
            {
                ndv=0;
            }// if

//2011-05-18 Replaced goto 100 statement
            //cout << "m,mpr: " << m << " " << mpr << endl;
            if(m<=0)
            {

                break;
            } //exit D100;
        }// do D100;


//end big loop!!!!!!!!!





    if(false)
    {
        // testing iepsmp and idebmap:
        ofstream epsfile;
        ofstream debfile;
        epsfile.open ("epsmap_win_vw.txt");
        debfile.open ("debmap_win_vw.txt");
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









        if (n2>ibmx)
        {
            cout <<"ibnd upper bound " << n2 << " exceeds ibmx" << endl;
            exit (0);
        }//if

//2011-05-18 Memory cleanup is done with DEALLOCATE F95
//statement

        /*
        if(allocated(cbn1_v)) deallocate(cbn1_v);
        if(allocated(cbn2_v)) deallocate(cbn2_v);
        if(allocated(cbal)) deallocate(cbal);
        */

//        delete cbn1_v;
//        delete cbn2_v;
//        delete cbal;


        if(cbn1_v.size()>0) vector <delphi_integer>().swap(cbn1_v);
        if(cbn2_v.size()>0) vector <delphi_integer>().swap(cbn2_v);
        if(cbal.size()>0) vector <delphi_integer>().swap(cbal);

        if(cbn1 != NULL) free_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
        if(debug_space) cout << "### freed cbn1 ###" << endl;
        if(debug_space) cout << "### cbn1: " << cbn1 << endl;
        if(cbn2 != NULL) free_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);


#ifdef VERBOSE
        cout <<"no. cavity mid-points inaccessible to solvent = " << ncav << endl;
#endif // VERBOSE


//consolidate the list, removing dead boundary points, adding new ones..
        j=0;
        //ncms=0;

//2011-05-19 Array is re-sized keeping old values in the
//memory.


	if(debug_space)cout << "Lin Li: ibgrd_v.size(): " << ibgrd_v.size() << endl;

        if(ibgrd_v.size() > 0)
        {

            Nar=sizeof(ibgrd_v);
            if (Nar<ibmx)
            {
                ibgrd_v.resize(ibmx);
                ibgrd=&ibgrd_v[-1];

                /* don't need this because of using vectors
                allocate(ibgrd_temp(Nar));
                ibgrd_temp=ibgrd;
                deallocate(ibgrd);
                allocate(ibgrd(ibmx));
                ibgrd(1:Nar)=ibgrd_temp;
                deallocate(ibgrd_temp);
                */

                //Lin Li: disscuss
            }// if
        }
        else
        {
            //ibgrd_v.assign(ibmx+1,{0,0,0});
            ibgrd_v.assign(ibmx+1,sgrid_temp_int);
            ibgrd=&ibgrd_v[-1];

        }// if


        for(i=1; i<=n2; i++)
        {

            ixyz=ibnd[i];
            ix=ixyz.nX;
            iy=ixyz.nY;
            iz=ixyz.nZ;

            if (bndeps[ix][iy][iz][1]!=0)
            {
                j=j+1;
                bndeps[ix][iy][iz][1]=j;


//2011-05-19 Precaution not to exceed array size (see above comment)
                if (j<=ibmx)
                {
                    ibgrd[j]=ixyz;
                }
                else
                {
                    cout << "j=" << j << " is larger than ibmx= " << ibmx << " << thus stopped..." << endl;
                    exit (0);
                }// if
            }// if


            for(jj=1; jj<=6; jj++)
            {

                ntt=optComp(iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]],imap[4][jj]);
                nt=(ntt%epsdim);

                if (nt<0)
                {
                    ntt=-ntt;

                    switch (imap[4][jj])
                    {
                    case 1:
                        iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nX=ntt;
                        break;
                    case 2:
                        iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nY=ntt;
                        break;
                    case 3:
                        iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nZ=ntt;
                        break;
                    default:
                        if(debug_space) cout <<"????? flag3" << endl;
                    }// select;

                    if (nt==-1)
                    {
                        ntt=ntt-1;

                        switch (imap[4][jj])
                        {
                        case 1:
                            iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nX=ntt;
                            break;
                        case 2:
                            iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nY=ntt;
                            break;
                        case 3:
                            iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]].nZ=ntt;
                            break;
                        default:
                            if(debug_space) cout <<"????? flag4" << endl;
                        }// select;
                    }// if

                }// if

#ifdef DEBUG
                ntt=optComp(iepsmp[iz+imap[3][jj]][iy+imap[2][jj]][ix+imap[1][jj]],imap[4][jj]);
                nt=(ntt%epsdim);
                jjj=ntt/epsdim;
                if (nt==0&&jjj>0)
                {
                    cout <<"PROBLEMS 5" << ix << iy << iz << jj << endl;
                }// if

#endif // DEBUG

            }// do

        }// do



        if (j>ibmx)
        {
            cout << "no. ms points exceeds ibmx" << endl;
            exit (0);
        }// if

        iBoundNum=j;

#ifdef VERBOSE
        cout <<"after surface elaboration iBoundNum= " << iBoundNum << endl;
        cout <<" and iBoundNumsurf= " << iBoundNumsurf << endl;
#endif // VERBOSE

    }// if

    /*
    if (bDebug)
    {
        open(52,file='bgpDurante.log',form='formatted');
        write(52,*)iGrid,iBoundNum;
        for(iiii=1; iiii<=iBoundNum; iiii++)
        {
            write(52,*)iiii,ibnd(iiii);
            write(52,*)iiii,ibnd(3*iiii-2),ibnd(3*iiii-1),ibnd(3*iiii);
        }// do
        close (52);
    }//if
    */



    /*
    if(allocated(bndeps)) deallocate(bndeps);
    if(allocated(ibnd)) deallocate(ibnd);
    */
    //delete bndeps;
    delete [] ibnd;
    if(bndeps != NULL) free_pt4d<delphi_integer>(bndeps,iGrid+1,iGrid+1,iGrid+1,3);

    if(debug_space) cout << "free memery done" << endl;

    //fScale bondary grid point positions relative to acc data
        if (isolv&&(irea||logs||lognl||isen||isch))
        {
#ifdef VERBOSE
            cout <<"scaling boundary grid points ..." << endl;
#endif // VERBOSE



            //allocate(scspos(iBoundNum));
            //scspos_v.assign(iBoundNum,{0.,0.,0.});
            scspos_v.assign(iBoundNum,sgrid_temp_real);
            scspos=&scspos_v[-1];

            //cout << "Lin Li: 0" << endl;
            //float ***a=get_pt3d<float>(1000,1000,1);

            for(j=1; j<=iBoundNum; j++)
            {
                scspos[j]=optCast<delphi_real,delphi_integer>(ibgrd[j]);
                //cout << "j,scspos[j]: " << j << " " << scspos[j] << endl;
            }// do



            //allocate(scsnor(iBoundNum),atsurf(iBoundNum),atndx(iBoundNum));

            //scsnor_v.assign(iBoundNum,{0.,0.,0.});
            scsnor_v.assign(iBoundNum,sgrid_temp_real);
            scsnor=&scsnor_v[-1];
            //scsnor = new  SGrid <delphi_real>  [iBoundNum+1];

            atsurf_v.assign(iBoundNum,0);
            atsurf=&atsurf_v[-1];
            //atsurf = new delphi_integer [iBoundNum+1];

            atndx_v.assign(iBoundNum,0);
            atndx=&atndx_v[-1];

            //atndx = new delphi_integer [iBoundNum+1];


            sclbp();

#ifdef VERBOSE
                cout << iall << " points had to be assigned by global comparison" << endl;
#endif // VERBOSE

            if (!isite && scsnor_v.size() >0 ) vector<SGrid <delphi_real> >().swap(scsnor_v); //no need to deallocate for vectors
        }// if

    //if (isrf&&!ivertnorm) //ivertnorm: only used in vdtms once, without initialized, removed
    if (isrf)
    {
        if (bOnlyMol)
        {
//2011-05-19 Parameters transfered via module architecture
            //allocate(egrid(iGrid,iGrid,iGrid));

            //egrid = get_pt3d <SGrid <delphi_integer> > (iGrid+1,iGrid+1,iGrid+1);
            get_pt3d <SGrid <delphi_integer> > (egrid,iGrid+1,iGrid+1,iGrid+1);

           // msrf();


            //if(allocated(egrid) )deallocate(egrid);
            if(egrid != NULL) free_pt3d(egrid,iGrid+1,iGrid+1,iGrid+1);
        }
        else
        {
            cout <<"msrf routine cannot be run" << endl;
            cout <<"because there are also geometric objects" << endl;
        }// if
    }// if
        //cout << "Lin Li: free atsurf??" << endl;
        //cout << isitsf << " " << isite << " " << isch << " " << scrgfrm << endl;
        if (!isitsf&&!isite&&!(isch&&scrgfrm!=0))
        {
            //if(allocated(atndx)) deallocate(atndx);
            //if(allocated(atsurf)) deallocate(atsurf);
            //cout << "Lin Li: free atsurf..." << endl;

            if(atndx_v.size()>0) vector <delphi_integer>().swap(atndx_v);
            if(atsurf_v.size()>0) vector <delphi_integer>().swap(atsurf_v);

        }// if
    /*
        if(allocated(iab1)) deallocate(iab1);
        if(allocated(iab2)) deallocate(iab2);
        if(allocated(icume)) deallocate(icume);
        if(allocated(r0)) deallocate(r0);
        if(allocated(r02)) deallocate(r02);
        if(allocated(rs2)) deallocate(rs2);
        if(allocated(ast)) deallocate(ast);
    */
        if(iab1 != NULL) free_pt3d<delphi_integer>(iab1,lcb1+1,mcb1+1,ncb1+1);
        if(iab2 != NULL) free_pt3d<delphi_integer>(iab2,lcb1+1,mcb1+1,ncb1+1);


        if(cbn1 != NULL) free_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
        if(debug_space) cout << "### freed cbn1 ###" << endl;
        if(cbn2 != NULL) free_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);

        if(r0.size()>0) vector <delphi_real>().swap(r0);
        if(r02.size()>0) vector <delphi_real>().swap(r02);
        if(rs2.size()>0) vector <delphi_real>().swap(rs2);
        if(ast.size()>0) vector <delphi_integer>().swap(ast);



        ibgrd_v.resize(iBoundNum);
        ibgrd=&ibgrd_v[-1];

}// void vwtms;
