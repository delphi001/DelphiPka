#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"

#include "space.h"

using namespace std;
//#####################################################################
//ATTENTION! This file is part of epsmakmod module.
// Do not compile separately!
//
//2011-05-20 Parameters to void are transfered via qlog and
// pointers modules
//#####################################################################
void CDelphiSpace::msrf()
{


    delphi_integer itot, vtot;
    string fixyn;
    string vdat1, line[6], fname;
    //SGrid <delphi_real> xo,xo2;
    bool fix;
    delphi_integer epsdim;
    SGrid <delphi_real> v1,v2,vxyz;
    //bool out,nbe[7],exists;
    SGrid <delphi_integer> iv123;
//2011-05-27 Declarations added due to IMPLICIT NONE
    delphi_integer mxvtx,i,ia1,ia2,ia3,iate,it,imxvtx;
    delphi_integer ib,iv1,iv2,iv3,j,k,mxtri,ntot2;
    delphi_real aa,area,cc,bb,rad,rad2,ss;
    delphi_real tar,tne4,vmg;

    iv1=0;
    iv2=0;
    iv3=0;
    vtot=0;
    itot=0;
    ntot2=0;

    if(debug_space) cout << "### in msrf: ###" <<endl;

    mxvtx=iBoundNum*2.2;
    mxtri=2*mxvtx;
    epsdim=iNatom+iNObject+2;

//2011-05-20 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
    for(k=1; k<=iGrid; k++)
    {
        for(j=1; j<=iGrid; j++)
        {
            for(i=1; i<=iGrid; i++)
            {
//cambiato da mod a div, mi dovrebbe servire solo il mezzo
                //egrid[i][j][k]=iEpsMap[i][j][k]/epsdim;
                egrid[i][j][k]=iepsmp[i][j][k]/epsdim;
            }// do
        }// do
    }// do

    for(i=1; i<=iGrid; i++)
    {
        for(j=1; j<=iGrid; j++)
        {
            egrid[i][1][j].nX = egrid[i][1][j].nY;
            egrid[i][j][1].nX = egrid[i][j][1].nZ;
        }// do
    }// do

    for(k=2; k<=iGrid; k++)
    {
        for(j=2; j<=iGrid; j++)
        {
            for(i=2; i<=iGrid; i++)
            {
                iate = 0;
                if (egrid[i][j][k].nX > 0) iate = iate + 1;
                if (egrid[i][j][k].nY > 0) iate = iate + 1;
                if (egrid[i][j][k].nZ > 0) iate = iate + 1;
                if (egrid[i-1][j][k].nX > 0) iate = iate + 1;
                if (egrid[i][j-1][k].nY > 0) iate = iate + 1;
                if (egrid[i][j][k-1].nZ > 0) iate = iate + 1;

                if (iate<=3)
                {
                    egrid[i][j][k].nX = 0;
                }
                else
                {
                    egrid[i][j][k].nX = 1;
                }// if
            }// do
        }// do
    }// do

//2011-05-20 Arrays are allocated with ALLOCATE F95 statement
// Initial values are added (just in case)
    /*
        allocate(vindx(mxtri),vert(mxvtx));
        vindx=int_coord(0,0,0);
        vert=coord(0,0,0);
      */
    //vindx.assign(mxtri+1, {0,0,0});
    //vert.assign(mxvtx+1, {0.,0.,0.});
    vindx.assign(mxtri+1, sgrid_temp_int);
    vert.assign(mxvtx+1, sgrid_temp_real);

    //print*,'mxtri,mxvtx=',mxtri,mxvtx
    cout << "mxtri,mxvtx= " << " " << mxtri << " " << mxvtx << endl;
//bin //bioengr //boot //cgroup //collective //common //common1 //cugi //dev //etc //fast //home //lib //lib64 //local_scratch //lost+found //media //misc //mnt //newscratch //opt //proc //root //sbin //selinux //smlc //software //srv //sys //tmp //usr //vaar //var //xcatpost
//2011-05-20 In void ex egrid array is three-dimensional while in this void this array was four-dimensional.
// Also array vindx was 2D here, but 1D in the EX void.
// Also it should be clarified why file in /usr/local/bin is opened and read in ex void.

    vdat1="./";

    //call ex(vtot, itot, vdat1, 2 );

    if (vtot>mxvtx)
    {
        cout <<"vtot = " << vtot << " > mxvtx = " << mxvtx << endl;
        cout <<"increase mxvtx in msrf.f" << endl;
        exit(0);
    }// if

    for(ib=1; ib<=vtot; ib++)
    {
        vert[ib]=vert[ib]/2.;
    }// do

    itot = itot/3;

//fScale boundary grid pointeger positions relative to acc data
    cout <<"scaling vertices" << endl;
    //allocate(vnorm(vtot),vnorm2(vtot));

    //vnorm.assign(vtot+1, {0.,0.,0.});
    //vnorm2.assign(vtot+1, {0.,0.,0.});
    vnorm.assign(vtot+1, sgrid_temp_real);
    vnorm2.assign(vtot+1, sgrid_temp_real);

    //call sclbp(vert,vnorm,vtot,iab1,iab2);

//fix holes and make vertex to triangle arrays allocate hole-fixing arrays next hole-fixing variables
    if (vtot < mxvtx/2)
    {
        imxvtx = vtot*2;
    }
    else
    {
        imxvtx = mxvtx;
    }// if


//    if (itot < mxtri/2)
//    {
//        imxtri = itot*2;
//    }
//    else
//    {
//        imxtri = mxtri;
//    }// if


    //allocate(vtlen(imxvtx),vtlst(6*imxvtx));
    //allocate(vtpnt(imxvtx),tmlst(9,imxvtx));
    vtlen.assign(imxvtx+1,0);
    vtpnt.assign(imxvtx+1,0);
    vtlst.assign(6*imxvtx+1,0);

    //tmlst=get_pt2d <delphi_integer> (9+1,imxvtx+1);
    get_pt2d <delphi_integer> (tmlst,9+1,imxvtx+1);

//2011-05-25 Arrays to void are transfered via pointers
// module
    //call mkvtl(1,vtot, 1,itot, imxtri, imxvtx);

    //call fxhl(1,vtot,1,itot,ntot, itot, imxvtx, imxtri);

    //call fxhl(1,vtot,1,itot,ntot2, itot, imxvtx, imxtri);

//write(*,*) "number of triangles added to fix holes= ",ntot+ntot2

    if (ntot2 > 0)
    {
        fix = true;
    }
    else
    {
        fix = false;
    }// if

    while(fix)
    {
        //call fxhl(1,vtot,1,itot,ntot2, itot, imxvtx, imxtri);

//write(*,*) "number of triangles added to fix holes= ",ntot2

        if (ntot2 > 0)
        {
            fix = true;
        }
        else
        {
            fix = false;
        }// if
    }// do

    if (itot>mxtri)
    {
        cout <<"itot = " << itot << " > mxtri = " << mxtri << endl;
        cout <<"increase mxtri in msrf.f" << endl;
        exit(0);
    }// if


    //if(allocated(vtlen)) deallocate(vtlen);
    //if(allocated(vtlst)) deallocate(vtlst);
    //if(allocated(tmlst)) deallocate(tmlst);
    //if(allocated(vtpnt)) deallocate(vtpnt);


    if(vtlen.size()>0) vector <delphi_integer>().swap(vtlen);
    if(vtlst.size()>0) vector <delphi_integer>().swap(vtlst);
    if(vtpnt.size()>0) vector <delphi_integer>().swap(vtpnt);
    if(tmlst != NULL) free_pt2d(tmlst,9+1,imxvtx+1);



    cout << "number of vertices = " << vtot << endl;
    cout << "number of triangles = " << itot << endl;

    //vnorm2={0.,0.,0}; //already initialized

//calculate area
    area=0.0;
    //areas=0.0;
    //areac=0.0;
    //arear=0.0;


    for(it=1; it<=itot; it++)
    {
        iv123=vindx[it];
        v1=vert[iv2]-vert[iv1];
        v2=vert[iv3]-vert[iv1];

//2011-05-20 Vector product defined in operators_on_coordinates module
        vxyz=optCross(v1,v2);
        vmg=sqrt(optDot(vxyz,vxyz));
        tar=vmg/2.;
        vxyz=vnorm[iv1]+vnorm[iv2]+vnorm[iv3];
        vmg=sqrt(optDot(vxyz,vxyz));

        vnorm2[iv1]=vnorm2[iv1]+(vxyz/vmg);
        vnorm2[iv2]=vnorm2[iv2]+(vxyz/vmg);
        vnorm2[iv3]=vnorm2[iv3]+(vxyz/vmg);

//calculate spherical triangle area if appropriate
        ia1=atndx[iv1];
        ia2=atndx[iv2];
        ia3=atndx[iv3];

        if (ia1>0)
        {
            if (ia1==ia2&&ia1==ia3)
            {
                rad=sDelPhiPDB[ia1].radius;
                rad2=rad*rad;
                aa=optSum((vert[iv2]-vert[iv1])*(vert[iv2]-vert[iv1]));
                bb=optSum((vert[iv3]-vert[iv2])*(vert[iv3]-vert[iv2]));
                cc=optSum((vert[iv1]-vert[iv3])*(vert[iv1]-vert[iv3]));

                aa=acos(1.-aa/(2.*rad2));
                bb=acos(1.-bb/(2.*rad2));
                cc=acos(1.-cc/(2.*rad2));
                ss=(aa+bb+cc)*.5;
                tne4=sqrt(tan(ss*.5)*tan((ss-aa)*.5)*tan((ss-bb)*.5)*tan((ss-cc)*.5));
                tar=4.*atan(tne4)*rad2;
            }// if
        }// if

        area=area+tar;

    }// do

    for(i=1; i<=vtot; i++)
    {
        vmg=sqrt( optDot( vnorm2[i],vnorm2[i] ) );
        vnorm2[i]=vnorm2[i]/vmg;
    }// do

    cout <<"MS area = " << area << endl;

//2011-05-26 Subroutine wrtsurf is short and called only once,
// therefore put the body of the void here
//-------- Begin of wrtsurf body -----------------------------


    if (!ibem)
    {

        ofstream surfile;
        surfile.open ("grasp.srf");

        fname="grasp.srf";
        cout << "writing GRASP file to " << fname << endl;
        //line=" ";
        line[1]="format=2";
        line[2]="vertices,normals,triangles";
        //write(line[4],'(3i6,f12.6)') vtot,itot,iGrid,fScale;
        //write(line[5],'(3f12.6)') cOldMid;
        surfile << "format=2" << endl;
        surfile << "vertices,normals,triangles" << endl;
        surfile << endl;
        surfile <<  setw(6) << vtot << setw(6) << itot << setw(6) << iGrid << setw(12) << setprecision(6) << fScale;
        surfile << setw(6) << setw(12) << setprecision(6) << cOldMid.nX
                << setw(6) << setw(12) << setprecision(6) << cOldMid.nY
                << setw(6) << setw(12) << setprecision(6) << cOldMid.nZ;

        cout << "writing data for" << vtot << " vertices and" << itot << " triangles" << endl;
        for(i=1;i<=mxvtx;i++){
            surfile << vert[i] << endl;
        }
        for(i=1;i<=vtot;i++){
            surfile << vnorm[i] << endl;
        }
        for(i=1;i<=mxtri;i++){
            surfile << vindx[i] << endl;
        }

        surfile.close();

        cout << "finished writing " << fname << endl;
    }
    else
    {
        ofstream surfile;
        surfile.open ("bem.srf");

        //open(7,file="bem.srf");
        //write(7,*)vtot,itot;
        surfile << vtot << " " << itot << endl;

        for(i=1; i<=vtot; i++)
        {
            surfile << vert[i] << endl;
        }// do

        for(i=1; i<=itot; i++)
        {
            surfile << vindx[i] << endl;
        }// do

        for(i=1; i<=vtot; i++)
        {
            surfile << vnorm[i] << endl;
        }// do

        surfile.close();
    }// if
//-------- End of wrtsurf body --------------------------------

}// void msrf;
