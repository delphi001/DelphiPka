#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "../space/space.h"

using namespace std;
//#####################################################################
//ATTENTION! This file is part of epsmakmod module.
// Do not compile separately!
//
//2011-05-26 Other parameters transfered via module architecture
//#####################################################################

//coi contains circle of intersection data for pairs when objects are involved
    struct ccoi
    {
        SGrid <delphi_real> xyz; //vector applied in the center of the sphere and pointing to the center of the coi
        delphi_real rad; //coi radius;
        delphi_integer is; //number of the intersections between an atom and an object
    };// type ccoi;

void CDelphiSpace::sas()
{

    delphi_integer nver=520,nedge=1040;
//2011-05-26 Arrays declared in pointers module and allocated in
//vwtms sub

    SGrid <delphi_real> ver[nver+1];
    delphi_integer edgv[2+1][nedge+1],edg[nedge+1],oti[nver+1],st[nedge+1];
//a third field in pls takes into account the correspondence nprobj-npr

//2011-05-26 Array declared in pointers module and allocated here
// delphi_integer*4 pls(3,1)
//2011-05-26 Removed variables declared in qlog module and
//redeclared some variables as coord type
    //delphi_integer objecttype,itmp,jprec,kk,nprtobj,nprobj,inter;
    delphi_integer jprec,nprtobj,nprobj;
    delphi_integer ii,jj,dim,dim1;
    string strtmp,strtmp1;
    //delphi_real side,dist;
    //SGrid <delphi_real> omin,omax,vectx,xmin,xmax,xa,xb, xc, xd, xq;

    SGrid <delphi_real> xyzm,x123, dx123, tij123;
    SGrid <delphi_real> rmv[3+1],cf123, dy123;
    SGrid <delphi_integer> ix123, ic123;
    vector < SGrid <delphi_integer> > pls; // pls is changed to be local variable
    vector < SGrid <delphi_integer> > plstemp;
    /*
    type(ccoi), allocatable :: coi(:),coitemp(:);
    type(int_coord), allocatable :: plstemp(:);
    type(coord), allocatable :: expostemp(:);
    */
    //ccoi * coi,coitemp;
    ccoi coi_1;
    vector <ccoi> coi;
    vector <ccoi> coitemp;
    vector < SGrid <delphi_real> > expostemp;


    //delphi_real tmp,tmp1;
    //delphi_real radius,modul,mod2;
    //delphi_real alpha;
    //delphi_real dx,dy,dz,seno,cose,modx,mody,modul2;
    delphi_real rad2;
//2011-05-26 Declaration due IMPLICIT NONE
    delphi_integer nacc,nacct,i,j,k,ie,ia2,ie2,ie1,ilvl,ia1,iv,iv1,iv2;
    delphi_integer ip,liml,limu,ne,nlvl,npr,nprp,nprx,nst,nprt,nvo;
    delphi_integer nxa,nvi, nv;
    delphi_real ctf,ctf2,cst,dctf,d2,del,dx1,dx2,dx3,ds2,pre,rad, radj;
    delphi_real rij,rdn,rv1,rv2,rvmg,dctf2,sm1,sm2,snt,tta,vmg;
    delphi_real cbln,csp,dmg,tm;

    if(debug_space) cout << "######## in sas: ############" << endl;
    coi_1.xyz= sgrid_temp_real;
    coi_1.rad=0.;
    coi_1.is=0;
    //sgrid_temp_real= {0.,0.,0.};

    fill_n(ver,nver+1,sgrid_temp_real);

#ifdef PARALLEL_OMP
        #pragma omp for schedule(auto)
#endif


    for (i=0; i<3; i++)
    {
        for (j=0; j<nedge+1; j++)
        {
            edgv[i][j]=0;
        }
    }

//    ver=ver=coord(0.,0.,0.);
//    edgv=0.;
    nacc=extot;
    nacct=0;
    nprt=0;
    nprtobj=0;
    nlvl=5;
    nvi=12;

    radpmax=max(fRadPrb[1],fRadPrb[2]);
    if(debug_space) cout <<"radpmax, fRadPrb " << radpmax << fRadPrb[1] << fRadPrb[2] << endl;

    cbln=2.*(fRMax+radpmax); //cube length
    if(debug_space)  cout << "cbln: " << cbln << endl;

    cubedata(1.0,cbln);

    dim=(lcb+1)*(mcb+1)*(ncb+1); // dim: how many cubes

//2011-05-26 Arrays are allocated with ordinary ALLOCATE F95 statement
// WARNING! Allocation may be incorrect since memalloc
// keeps old values of the array in memory, if exist,
// while ALLOCATE statement kills all old values. Arrays
// are allocated again later in sas voids and in
// vwtms void

    //allocate(cbn1_v[dim]][cbn2_v[dim]);
    cbn1_v.assign(dim+1,0);
    cbn2_v.assign(dim+1,0);

    dim1=27;

    if ((iNObject-numbmol)>0) dim1=max(dim,27);

    //allocate(cbal(dim1*(natom+nobject-numbmol)))
    cbal.assign(dim1*(iNatom+iNObject-numbmol)+1,0);

    if(debug_space) cout << "dim, dim1,(natom+nobject-numbmol): " << dim << " " << dim1 << " " << (iNatom+iNObject-numbmol) << endl;

    cube();


//generate a template of vertices....
    tta=2.*fPi/nvi;
#ifdef PARALLEL_OMP
        #pragma omp for schedule(auto)
#endif



    for(i=1; i<=nvi; i++)
    {
        rdn=(i-1)*tta;
        //ver[i]= {cos(rdn),sin(rdn),0.};
        ver[i].nX=cos(rdn);
        ver[i].nY=sin(rdn);
        ver[i].nZ=0.;

        j=i+1;
        if(i==nvi)j=1;
        edgv[1][i]=i;
        edgv[2][i]=j;
    }// do

    nv=nvi;
    ne=nvi;
    ie1=1;
    ie2=0;


    for(ilvl=1; ilvl<=nlvl; ilvl++)
    {
        ie1=ie2+1;
        ie2=ne;
        for(ie=ie1; ie<=ie2; ie++)
        {
            iv1=edgv[1][ie];
            iv2=edgv[2][ie];

//2011-05-26 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
            xyzm=ver[iv1]+ver[iv2];
            vmg=sqrt(optDot(xyzm,xyzm));

            nv=nv+1;
            ver[nv]=xyzm/vmg;
            //cout << "ie,iv1,iv2,nv,ver[nv]: " << ie << " " << iv1 << " " << iv2 << " " << nv << " " << ver[nv] << endl;
            ne=ne+1;
            edg[ie]=ne;
            edgv[1][ne]=iv1;
            edgv[2][ne]=nv;
            ne=ne+1;
            edgv[1][ne]=nv;
            edgv[2][ne]=iv2;
        }// do
    }// do

    ne=ie2;
    //edg[ie1:ne]=-1;
#ifdef PARALLEL_OMP
        #pragma omp for schedule(auto)
#endif


    for (i=ie1; i<=ne; i++)
    {
        edg[i]=-1;
    }
    cout <<"# of vertices = " << nv << " # of edges = " << ne << endl;
    //ast=1;
    //for (i=0; i<=iNatom+1; i++)
    for (i=0; i<iNatom+1; i++)
    {
        ast[i]=1;
    }

    nacc=0;
    npr=0;
    nprobj=0;
    nprp=0;


//it finds pairs.......
    for(i=1; i<=iNatom; i++)
    {
        rad=r0[i];
        rad2=r02[i];
        if(sDelPhiPDB[i].radius==0.) continue;
        x123=xn1[i];
        ix123=optCast <delphi_integer,delphi_real> ((x123-xyzo)*cbai);
        //Lin Li: 1 dimensional operation is abandoned
        //liml=cbn1_v[ix123.nX+1+(lcb+1)*ix123.nY+(lcb+1)*(mcb+1)*ix123.nZ];
        //limu=cbn2_v[ix123.nX+1+(lcb+1)*ix123.nY+(lcb+1)*(mcb+1)*ix123.nZ];

        liml=cbn1[ix123.nX][ix123.nY][ix123.nZ];
        limu=cbn2[ix123.nX][ix123.nY][ix123.nZ];

        //cout << "i,ix123.nX+1+(lcb+1)*ix123.nY+(lcb+1)*(mcb+1)*ix123.nZ,cbn1[]: " << cbn1[ix123.nX][ix123.nY][ix123.nZ] << endl;
//2011-06-17 Resizing of arrays keeping old values intact
        if ((npr+limu-liml+1)>nprt)
        {
            nprt=nprt+5000;
            //if (allocated(pls))
            if(pls.size() > 0)
            {
                /*
                //allocate(plstemp(nprt-5000));
                plstemp.assign(nprt-5000+1, sgrid_temp_int);
                plstemp=pls;
                //deallocate(pls);
                vector < SGrid <delphi_integer> >().swap(pls);
                //allocate(pls(nprt));
                pls.assign(nprt+1, {0,0,0});
                //pls(1:nprt-5000)=plstemp;
                for(j=1; j<=nprt-5000; j++)
                {
                    pls[j]=plstemp[j];
                }
                //deallocate(plstemp);
                vector < SGrid <delphi_integer> >().swap(plstemp);
                */
                pls.resize(nprt+1);
            }
            else
            {
                //allocate(pls(nprt));
                pls.assign(nprt+1, sgrid_temp_int);
            }// if
        }// if

        if ((nprobj+limu-liml+1)>nprtobj)
        {
            nprtobj=nprtobj+1000;
            //if (allocated(coi))
            if(coi.size() >0)
            {
                //allocate(coitemp(nprtobj-1000));
                //coi contains circle of intersection data for pairs when objects are involved

                /*
                coitemp.assign(nprtobj-1000+1,coi_1);
                coitemp=coi;
                //deallocate(coi);
                vector <ccoi> ().swap(coi);
                //allocate(coi(nprtobj));
                coi.assign(nprobj+1,coi_1);
                //coi(1:nprtobj-1000)=coitemp;
                for(j=1; j<=nprtobj-1000; j++)
                {
                    coi[j]=coitemp[j];
                }
                //deallocate(coitemp);
                vector <ccoi> ().swap(coitemp);
                */
                coi.resize(nprobj+1);
            }
            else
            {
                //allocate(coi(nprtobj));
                coi.assign(nprobj+1,coi_1);
            }// if
        }// if

        jprec=0;
        //cout << "liml, limu: " << liml << " " << limu << endl;
        for(jj=liml; jj<=limu; jj++)
        {
            j=cbal[jj];
            //cout << "jj,cbal[jj]: " << jj << " " << cbal[jj] << endl;
            if (j<=iNatom)
            {
                radj=r0[j];
                //cout << "sDelPhiPDB[j].radius,j,i: " << sDelPhiPDB[j].radius << " " << j << " " << i << endl;
                if (sDelPhiPDB[j].radius>0.&&j>i)
                {
                    //cout << "sDelPhiPDB[j].radius,j,i: " << sDelPhiPDB[j].radius << " " << j << "" << i << endl;
                    ctf=rad+radj;
                    ctf2=ctf*ctf;
                    dctf=abs(rad-radj);
                    dctf2=dctf*dctf;
                    dx123=xn1[j]-x123;
                    d2=optDot(dx123,dx123);
                    del=ctf2-d2;

                    //cout << "npr,del,d2,dctf2: " << npr << " " << del << " " << d2 << " " << dctf2 << endl;
                    if (del>0.01&&d2>dctf2)
                    {
                        npr=npr+1;
                        //pls[npr]= {i,j,0};
                        pls[npr].nX=i;
                        pls[npr].nY=j;
                        pls[npr].nZ=0;
                        //cout << "npr,del,d2,dctf2: " << npr << " " << del << " " << d2 << " " << dctf2 << endl;
                    }// if
                }// if
            }
            else //objects are abandoned:
            {
                /*
                   if (j!=jprec)
                   {
                //it finds out if there is intersection between i and
                //kk and it generates suitable parameters
                //kk= objectnumber
                       kk=j-iNatom;
                       strtmp=dataobject[kk][1];
                       strtmp1=dataobject[kk][2];
                       inter=0;
                       xq=xn1[i];
                       read(strtmp(8:10),'(I3)')objecttype;

                       switch (objecttype)
                       {
                       case 1 //if (objecttype==1) {
                //dealing with a sphere
                               call distobj(xq][dist][dxyz][kk][fRadPrb[1]][false);

                           if (abs(dist)<r0[i])
                           {
                               inter=inter+1;
                               npr=npr+1;
                               nprobj=nprobj+1;
                               pls(npr)=int_coord(i,j,nprobj);
                               coi(nprobj)=ccoi((-dxyz)*dist][sqrt(r02[i]-dist**2)][inter);
                           }// if
                        break;

                       case 2 //if (objecttype==2) {
                //dealing with a cylinder
                               call distobj(xq][dist][dxyz][kk][fRadPrb[1]][true);

                           if (abs(dist)<r0[i])
                           {
                               read(strtmp(20:80),*)xa,xb,radius;
                               read(strtmp1,'(5f8.3)')vectz,modul,modul2;
                               tmp=zeta+fRadPrb[1];

                               if (abs(tmp)<r0[i])  //side B;
                               {
                                   inter=inter+1;
                                   npr=npr+1;
                                   nprobj=nprobj+1;
                                   pls(npr)=int_coord(i,j,nprobj);
                                   coi(nprobj)=ccoi(vectz*(tmp/modul)][sqrt(r02[i]-tmp**2)][inter);
                               }// if

                               tmp=zeta-fRadPrb[1];

                               if (abs(tmp)<r0[i])   //side A;
                               {
                                   inter=inter+1;
                                   npr=npr+1;
                                   nprobj=nprobj+1;
                                   pls(npr)=int_coord(i,j,nprobj);
                                   coi(nprobj)=ccoi(vectz*(tmp/modul)][sqrt(r02[i]-tmp**2)][inter);
                               }// if

                               tmp=axdist-radius-fRadPrb[1];

                               if (abs(tmp)<r0[i])   //lateral][closest][;
                               {
                //approximating
                //with planes
                                   if (axdist==0.)
                                   {
                                       cout <<"cannot use planar approximation" << i << j << endl;
                                   }
                                   else
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi((tmp/axdist)*&(xq-xb-vectz*(zeta/modul))][sqrt(r02[i]-tmp**2)][inter);
                                   }// if
                               }// if

                               tmp=axdist+radius+fRadPrb[1];

                               if (abs(tmp)<r0[i])   //lateral][farthest][;
                               {
                //approximating
                //with planes
                                   cout <<"planar approx very poor << sphere too big" << i << j << endl;
                                   if (axdist==0.)
                                   {
                                       cout <<"cannot use planar approximation" << i << j << endl;
                                   }
                                   else
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi((tmp/axdist)*&(xq-xb-vectz*(zeta/modul))][sqrt(r02[i]-tmp**2)][inter);
                                   }// if
                               }// if
                           }// if
                           break;
                       case 3 //if (objecttype==3) {
                //dealing with a cone
                               call distobj(xq][dist][dxyz][kk][fRadPrb[1]][true);

                           if (abs(dist)<r0[i])
                           {
                               read(strtmp(20:80),*)xa,xb,alpha;
                               alpha=alpha*pi/180;
                               seno=sin(alpha);
                               cose=cos(alpha);

                               read(strtmp1,'(5f8.3)')vectz,modul,modul2;
                               tmp=zeta+fRadPrb[1];

                               if (abs(tmp)<r0[i])   //side B;
                               {
                                   inter=inter+1;
                                   npr=npr+1;
                                   nprobj=nprobj+1;
                                   pls(npr)=int_coord(i,j,nprobj);
                                   coi(nprobj)=ccoi(vectz*(tmp/modul)][sqrt(r02[i]-tmp**2)][inter);
                               }// if

                               tmp=axdist*cose-fRadPrb[1]-seno*(modul-zeta);

                               if (abs(tmp)<r0[i])   //lateral][closest][;
                               {
                //approximating
                //with planes
                                   if (axdist==0.)
                                   {
                                       cout <<"no planar approx since sphere too large" << i << j << endl;
                                   }
                                   else
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       tmp1=tmp*cose/axdist;
                                       tv=tmp1*(xb-xq+vectz*((zeta-axdist*seno/cose)/modul));
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi(tv][sqrt(r02[i]-tmp**2)][inter);
                                   }// if
                               }// if

                               tmp=seno*(modul-zeta)+cose*axdist+fRadPrb[1];

                               if (tmp<r0[i])
                               {
                                   if (zeta>modul+fRadPrb[1]||axdist==0.)
                                   {
                                       cout <<"cannot use planar approx in this pos" << i << j << endl;
                                       else  //lateral,farthest, approximating with;
                                       {
                //planes
                                           inter=inter+1;
                                           npr=npr+1;
                                           nprobj=nprobj+1;
                                           pls(npr)=int_coord(i,j,nprobj);
                                           tv=(tmp/axdist)*(-cose*(xq-xb)+((cose*zeta+axdist*seno)/modul)*vectz);
                                           coi(nprobj)=ccoi(tv][sqrt(r02[i]-tmp**2)][inter);
                                       }// if
                                   }// if

                //beyond the tip, close to the axis
                                   if (zeta>modul+fRadPrb[1]&&axdist<1.0e-6)
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       tmp=(dist-r0[i])*(1.-cose)/modul;
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi(tmp*vectz][r0[i]*seno][inter);
                                   }// if
                               }// if
                               break;
                           case 4 //if (objecttype==4) {
                //dealing with a parallelepiped
                                   read(strtmp(20:80),*)xa,xb,xc,xd;
                               read(strtmp1,'(12f6.2)')vectz,modul,vecty,mody,vectx,modx;

                //now using the newest notation for box vertices
                //new notation:vectx=B-A,vecty=C-A,vectz=D-A

                // xq=P-A

                //chiamate da togliere!!!!!!!!!!!!!!!!!!!!!!!!!!
                               xq=xq-xa;

                //working on z interval
                               dot=vectz.dot.xq;
                               tmp=dot/modul;
                               tmp1=fRadPrb[1]+modul;

                               if (tmp>-fRadPrb[1]-rad&&tmp<rad+tmp1)
                               {
                                   tmp2=tmp1-tmp;
                                   rdx2=rad2-tmp2**2;
                                   tmp3=fRadPrb[1]+tmp;
                                   rsx2=rad2-tmp3**2;

                                   if (rdx2>0.0)
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       tmp2=tmp2/modul;
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi(vectz*tmp2,sqrt(rdx2),inter);
                                   }// if

                                   if (rsx2>0.0)
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       tmp3=-tmp3/modul;
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi(vectz*tmp3,sqrt(rsx2),inter);
                                   }// if

                //working on y interval
                                   dot=vecty.dot.xq;
                                   tmp=dot/mody;
                                   tmp1=fRadPrb[1]+mody;
                                   tmp2=tmp1-tmp;
                                   rdx2=rad2-tmp2**2;
                                   tmp3=fRadPrb[1]+tmp;
                                   rsx2=rad2-tmp3**2;

                                   if (rdx2>0.0)
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi(vecty*(tmp2/mody),sqrt(rdx2),inter);
                                   }// if

                                   if (rsx2>0.0)
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi(vecty*(-tmp3/mody),sqrt(rsx2),inter);
                                   }// if

                //working on x interval
                                   dot=vectx.dot.xq;
                                   tmp=dot/modx;
                                   tmp1=fRadPrb[1]+modx;
                                   tmp2=tmp1-tmp;
                                   rdx2=rad2-tmp2**2;
                                   tmp3=fRadPrb[1]+tmp;
                                   rsx2=rad2-tmp3**2;

                                   if (rdx2>0.0)
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi(vectx*(tmp2/modx),sqrt(rdx2),inter);
                                   }// if

                                   if (rsx2>0.0)
                                   {
                                       inter=inter+1;
                                       npr=npr+1;
                                       nprobj=nprobj+1;
                                       pls(npr)=int_coord(i,j,nprobj);
                                       coi(nprobj)=ccoi(vectx*(-tmp3/modx),sqrt(rsx2),inter);
                                   }// if
                               }// if
                           }// select;
                       }// if
                       */
            }// if //end of j<=iNatom;

            jprec=j;

        }// do

        if (npr==nprp)
        {
            ast[i]=0;
        }// if

        nprp=npr;

    }// do

//2011-05-26 Added deallocation in order to avoid confusion
//with later allocation of the same arrays
    //if(allocated(cbn1_v)) deallocate(cbn1_v);
    //if(allocated(cbn2_v)) deallocate(cbn2_v);
    //if(allocated(cbal)) deallocate(cbal);

    if(cbn1_v.size()>0) vector <delphi_integer> ().swap(cbn1_v);
    if(cbn2_v.size()>0) vector <delphi_integer> ().swap(cbn2_v);
    if(cbal.size()>0) vector <delphi_integer> ().swap(cbal);

    if(cbn1 != NULL) free_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
    if(debug_space) cout << "### freed cbn1 ###" << endl;
    if(cbn2 != NULL) free_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);


    cout <<"# of pairs = " << npr << endl;


//2011-05-26 Temporarly removed time calculations
    cbln=fRMax+radpmax;

    cubedata(2.0,cbln);



    //dim=(lcb+1)*(mcb+1)*(ncb+1);
    //allocate(cbn1_v[dim]);
    //allocate(cbn2_v[dim]);


    dim1=27;
    if ((iNObject-numbmol)>0) dim1=max(dim,27);

    //allocate(cbal[dim1*(iNatom+iNObject-numbmol)]);
    cbal.assign(dim1*(iNatom+iNObject-numbmol)+1,0);
    cube();



    nprx=0;
    for(ip=1; ip<=npr; ip++)
    {
        i=pls[ip].nX;
        j=pls[ip].nY;

        if (j<=iNatom)
        {
            dx123=xn1[j]-xn1[i];
            d2=optDot(dx123,dx123);
            dmg=sqrt(d2);
            pre=1.+(r02[i]-r02[j])/d2;
            tij123=xn1[i]+((0.5*pre)*dx123);
            rij=0.5 * sqrt((r0[i]+r0[j])*(r0[i]+r0[j])-d2) * sqrt(d2-(r0[i]-r0[j])*(r0[i]-r0[j]))/dmg;
            //cout <<"i,j,rij: " << i << " " << j << " " << rij << endl;
        }
        else // never go to this else statement
        {
            cout << "### Warning: if j < iNatom else:" << endl;
            nprobj=pls[ip].nZ;
            rij=coi[nprobj].rad;

//pay attention, here dx has a different meaning from previous one
            dx123=coi[nprobj].xyz;
            d2=optDot(dx123,dx123);
            tij123=xn1[i]+dx123;
            dmg=sqrt(d2);
        }// if

        dx1=dx123.nX;
        dx2=dx123.nY;
        dx3=dx123.nZ;
        rvmg=sqrt(dx1*dx1+dx2*dx2);
        //cout << "rvmg: " << rvmg << endl;

        if (rvmg>1.0e-8)
        {
            rv1=-dx2/rvmg;
            rv2=dx1/rvmg;
            cst=dx3/dmg;

            snt=sqrt(1.-cst*cst);
//snt=rvmg/dmg !doesn't lead to any improved performance

            csp=1.0-cst;
            tm=csp*rv1;
            sm1=snt*rv1;
            sm2=snt*rv2;
            //rmv[1]= {tm*rv1+cst,tm*rv2,sm2};
            rmv[1].nX=tm*rv1+cst;
            rmv[1].nY=tm*rv2;
            rmv[1].nZ=sm2;

            //rmv[2]= {tm*rv2,csp*rv2*rv2+cst,-sm1};
            rmv[2].nX=tm*rv2;
            rmv[2].nY=csp*rv2*rv2+cst;
            rmv[2].nZ=-sm1;

            //rmv[3]= {-sm2,sm1,cst};
            rmv[3].nX=-sm2;
            rmv[3].nY=sm1;
            rmv[3].nZ=cst;
        }
        else
        {
            //rmv[1]= {1.,0.,0.};
            rmv[1].nX=1.0;
            rmv[1].nY=0.0;
            rmv[1].nZ=0.0;

            //rmv[2]= {0.,1.,0.};
            rmv[2].nX=0.0;
            rmv[2].nY=1.0;
            rmv[2].nZ=0.0;


            //rmv[3]= {0.,0.,1.};
            rmv[3].nX=0.0;
            rmv[3].nY=0.0;
            rmv[3].nZ=1.0;


        }// if

        nvo=0;
        //nbv=0;

/**
 * assign memory to expos if needed
*/
//2011-06-17 Re-sizing array keeping old value
        if ((nacc+nv)>nacct)
        {
            nacct=nacct+1000;
            //if (allocated(expos))
            if (expos.size()>0)
            {
                /*
                //allocate(expostemp(nacct-1000));
                expostemp.assign(nacct-1000+1,sgrid_temp_real);
                expostemp=expos;
                //deallocate(expos);
                //allocate(expos[nacct]);
                expos.assign(nacct+1,sgrid_temp_real);
                //expos[1:nacct-1000]=expostemp;
                for(itemp=1; itemp<=nacct-1000; itemp++)
                {
                    expos[itemp]=expostemp[itemp];
                }
                //deallocate(expostemp);
                */
                expos.resize(nacct+1);
            }
            else
            {
                //allocate(expos[nacct]);
                expos.assign(nacct+1,sgrid_temp_real);
            }// if
        }// if

        //do iv=1,nvi; //D10
        iv=1;
        //for (iv=1; iv<=nvi; iv++)
        //cout << "nvi: " << nvi << endl;
D10:

        while(iv<=nvi)
        {
/**
 * +rm(7)*ver[3][iv] has been removed because it is always zero
*/
            cf123.nX=rmv[1].nX*ver[iv].nX+rmv[1].nY*ver[iv].nY;
            cf123.nY=rmv[2].nX*ver[iv].nX+rmv[2].nY*ver[iv].nY;
            cf123.nZ=rmv[3].nX*ver[iv].nX+rmv[3].nY*ver[iv].nY;
            cf123=tij123+(cf123*rij);
            //cout << "iv,cf123: " << iv << " " << cf123 << endl;

            /*
                    if (j>iNatom) //for objects, abandoned
                    {
                        inter=coi[nprobj].is;

                        if (inter>1)
                        {
            //if inter>1 we are close to tips in object, thus, false vertices might have been previously
            //generated, so now, if a vertex is outside the object it is fictiously checked for occlusion by
            //oti but discarded afterwards kk= objectnumber
                            kk=j-iNatom;
                            xq=cf123;

                            //call distobj(xq][dist][dxyz][kk][fRadPrb[1]][true);

                            if (dist>5.0e-4)
                            {
                                oti[iv]=j;
                                cycle D10;
                            }// if
                        }// if
                    }// if
            */
            ic123=optCast <delphi_integer,delphi_real> ((cf123-xyzo)*cbai);

            //cout << "ic123: " << ic123 << endl;
            //liml=cbn1_v[ic123.nX+1+(lcb+1)*ic123.nY+(lcb+1)*(mcb+1)*ic123.nZ];
            //limu=cbn2_v[ic123.nX+1+(lcb+1)*ic123.nY+(lcb+1)*(mcb+1)*ic123.nZ];
            liml=cbn1[ic123.nX][ic123.nY][ic123.nZ];
            limu=cbn2[ic123.nX][ic123.nY][ic123.nZ];


            //cout << "iv,liml: " << iv << " " << liml << endl;
            ii=liml;
D05:
            //do ii=liml,limu;
            //for (ii=liml; ii<=limu; ii++)
            while(ii<=limu)
            {

                //cout << "flag1: ii: " << ii << endl;
                k=cbal[ii];
                if (k>iNatom)
                {
                    oti[iv]=k;
                    ii++;
                    goto D05;
                }// if
                //cout << "flag2: ii: " << ii << endl;
                dy123=xn1[k]-cf123;
                ds2=optDot(dy123,dy123);

                if (ds2<rs2[k])
                {
                    oti[iv]=k;
                    iv++;
                    goto D10;
                }// if
                //cout << "flag3: ii: limu: " << ii << " " << limu << endl;
                ii++;
            }// do D05;

            //cout << "iv,ii: " << iv << " " << ii << endl;
            nvo=nvo+1;
            nacc=nacc+1;
            expos[nacc]=cf123;
            //cout << "Lin Li1 : expos: " << nacc << " " << expos[nacc] << endl;
            oti[iv]=0;


            iv++;
        }// do D10;
                //cout << "### flag3: "   << endl;
        nst=0;
        if (nlvl>0)
        {

            for(ie=nvi; ie>=1; ie--)
            {
                //cout << "### flag4: ie: " << ie  << endl;

                ia1=oti[edgv[1][ie]];
                ia2=oti[edgv[2][ie]];

                if(ia1>0&&ia1==ia2) continue;
                nst=nst+1;
                st[nst]=ie;
            }// do
        }// if
        //cout << "### flag5: " << endl;

        if (nst>0)
        {
D030:
            while(true)
            {
                ie=st[nst];
                nst=nst-1;
                ia1=oti[edgv[1][ie]];
                ia2=oti[edgv[2][ie]];

                if ((ia1>iNatom)||(ia2>iNatom))
                {
                    if(nst>0) goto D030;
                    //exit D030;
                    //cout << "break D030: " << endl;
                    break;
                }// if

                iv=ie+nvi;

                //rm(7)*ver[3][iv] has been removed because it is always
                //zero
                //cf123= {optDot(rmv[1],ver[iv]), optDot(rmv[2],ver[iv]), optDot(rmv[3],ver[iv])};
                cf123.nX=optDot(rmv[1],ver[iv]);
                cf123.nY=optDot(rmv[2],ver[iv]);
                cf123.nZ=optDot(rmv[3],ver[iv]);

                cf123=tij123+(cf123*rij);

                if (ia1!=0)
                {
                    dy123=xn1[ia1]-cf123;
                    ds2=optDot(dy123,dy123);

                    if (ds2<rs2[ia1])
                    {
                        oti[iv]=ia1;

                        if (edg[ie]>0)
                        {
                            nst=nst+1;
                            st[nst]=edg[ie]+1;
                        }// if

                        if(nst>0) goto D030;
                        //exit D030;
                        //cout << "break D030: " << endl;
                        break;
                    }// if
                }// if

                if (ia2!=0)
                {
                    dy123=xn1[ia2]-cf123;
                    ds2=optDot(dy123,dy123);

                    if (ds2<rs2[ia2])
                    {
                        oti[iv]=ia2;

                        if (edg[ie]>0)
                        {
                            nst=nst+1;
                            st[nst]=edg[ie];
                        }// if

                        if(nst>0) goto D030;
                        //exit D030;
                        //cout << "break D030: " << endl;
                        break;
                    }// if
                }// if

                ic123=optCast <delphi_integer,delphi_real> ((cf123-xyzo)*cbai);
                //cout << "ic123: " << ic123 << endl;
                //liml=cbn1_v[ic123.nX+1+(lcb+1)*ic123.nY+(lcb+1)*(mcb+1)*ic123.nZ];
                //limu=cbn2_v[ic123.nX+1+(lcb+1)*ic123.nY+(lcb+1)*(mcb+1)*ic123.nZ];
                liml=cbn1[ic123.nX][ic123.nY][ic123.nZ];
                limu=cbn2[ic123.nX][ic123.nY][ic123.nZ];


                ii=liml;
D055:
                //do ii=liml,limu;
                //for(ii=liml; ii<=limu; ii++)
                while(ii<=limu)
                {
                    k=cbal[ii];

                    if (k>iNatom)
                    {
                        oti[iv]=k;
                        ii++;
                        goto D055;
                    }// if

                    dy123=xn1[k]-cf123;
                    ds2=optDot(dy123,dy123);

                    if (ds2<rs2[k])
                    {
                        oti[iv]=k;

                        if (edg[ie]>0)
                        {
                            nst=nst+1;
                            st[nst]=edg[ie]+1;
                            nst=nst+1;
                            st[nst]=edg[ie];
                        }// if

                        if(nst>0) goto D030;
                        goto END030; //exit D030;
                    }// if
                    ii++;
                }// do D055;

                nvo=nvo+1;
                nacc=nacc+1;
                expos[nacc]=cf123;


                oti[iv]=0;
                if (edg[ie]>0)
                {
                    if ( edg[edg[ie]+1] > 0 || ia2>0 )
                    {
                        nst=nst+1;
                        st[nst]=edg[ie]+1;
                    }// if

                    if (edg[edg[ie]]>0||ia1>0)
                    {
                        nst=nst+1;
                        st[nst]=edg[ie];
                    }// if
                }// if

                if(nst<=0) goto END030; //exit D030;
            }// do D030;
END030:

        cout << "";
        //cout << "END030" << endl;
        }// if



        if (nvo>0)
        {
            //considering pairs also where one 'partner' is an object
            nprx=nprx+1;
            ast[i]=0;
            if (j<=iNatom) ast[j]=0;
        }// if

    }// do

    //if(allocated(pls)) deallocate(pls);
    //if(allocated(coi)) deallocate(coi);
    //if(allocated(cbn1_v)) deallocate(cbn1_v);
    //if(allocated(cbn2_v)) deallocate(cbn2_v);
    //if(allocated(cbal)) deallocate(cbal);
    if(pls.size()>0) vector < SGrid <delphi_integer> >().swap(pls);
    if(coi.size()>0)  vector <ccoi> ().swap(coi);
    if(cbn1_v.size()>0)  vector <delphi_integer> ().swap(cbn1_v);
    if(cbn2_v.size()>0)  vector <delphi_integer> ().swap(cbn2_v);
    if(cbal.size()>0)  vector <delphi_integer> ().swap(cbal);

    if(cbn1 != NULL) free_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
    if(debug_space) cout << "### freed cbn1 ###" << endl;
    if(cbn2 != NULL) free_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);


    nxa=0;
    for(i=1; i<=iNatom; i++)
    {
        if(ast[i]==0)nxa=nxa+1;
    }// do

/*
    if (iNObject-numbmol>1) // for object, abandoned
    {
        cout <<"now calculating object-object exposed vertices" << endl;
        jj=0;
        for(ii=1; ii<=iNObject; ii++)
        {
            //strtmp=dataobject[ii][1];
            strtmp=dataobject_v[(ii-1)*2];

            //read(strtmp(16:18),*)kind;
            kind = atoi(strtmp.substr(15,3).c_str());
            if (strtmp.substr(0,4)!="is a" && kind!=2)

            //if (strtmp(1:4)!='is a'&&kind!=2)
            {
                if (jj==0)
                {
                    omin=sLimObject[ii].nMin;
                    omax=sLimObject[ii].nMax;
                }// if

                jj=1;
                omin=min(omin,sLimObject[ii].nMin);
                omax=min(omax,sLimObject[ii].nMax);
            }// if
        }// do

        //make cubedata
        //nside = number of initial subdivisions
        nside=4;
        //for objects, only water probes involved

        //OBS: it is advisable to improve the sideinter estimation as a
        //function of smallest object dimensions and of radprobe as
        //well, here probably is too small!!
        tmp=max(omax.nX-omin.nX,omax.nY-omin.nY);
        tmp=max(tmp,omax.nZ-omin.nZ)+2*fRadPrb[1];

        side=tmp/Int2Float(nside);
        xmin=omin-(fRadPrb[1]-0.5*side);
        xmax=omax+(fRadPrb[1]-0.5*side);
        xmax=omin+side*(0.5+optCast <delphi_real,delphi_integer> (optCast <delphi_integer,delphi_real> (0.999+((xmax-xmin)/side))));

        h=1./fScale;
        sideinter=max(fRadPrb[1]/4.,0.05);
        sidemin=sideinter/4;
        cout <<"Generating vertices between objects:" << endl;
        cout <<"Finite difference grid spacing:" << h << endl;
        cout <<"Initial side:" << side << endl;
        cout <<"Intermediate side:" << sideinter << endl;
        cout <<"Minimum side:" << sidemin << xmin << xmax << endl;
        nacct=nacc+(side/sidemin)*(side/sidemin)*(side/sidemin);
        cout <<"threshold for number of exposed vertices:" << nacct << endl;

        //2011-05-26 Array is already allocated
        //2011-05-26 Other parameters are transfered via qlog module
        //call objvertices(nacc,xmin,xmax,side,h, (numbmol>0));
    }// if
*/
    cout <<"# pairs analyzed (atom-atom and atom-object)= " << npr << endl;
    cout <<"# exposed pairs (atom-atom and atom-object)= " << nprx << endl;
    cout <<"no. arc points = " << nacc << endl;
    cout <<"no. surface atoms = " << nxa << " nbur = " << iNatom-nxa << endl;
//if(nacc>exmax)stop 'nacc limit exceeded'

    extot=nacc;

}// void sas;
