#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"
#include "space.h"

using namespace std;

//#####################################################################

void CDelphiSpace::crgarr()
{


    delphi_integer epsdim;
    SGrid <delphi_real> rxyz;
    SGrid <delphi_integer> jxyz;

//2011-05-30 Declarations added due to IMPLICIT NONE
    delphi_integer ic1,ic2,i,ix,ii,imed,jx,jy,jz;
    delphi_real chrg,rgrid;
    bool ipassed;
    SGridValue<delphi_real> temp;
    //SGrid <delphi_real> sgrid_temp_real;
    delphi_real radpolext=1.0;


/**
 * 2011-05-30 First determined ic1=nqass - number of assigned
 *charges in order to avoid re-sizing of arrays atmcrg,chgpos and
 *crgatn.
*/

//Part1 (easy) : from molecules
    //temp.nGrid={0.,0.,0.};
    temp.nGrid.nX=0.;
    temp.nGrid.nY=0.;
    temp.nGrid.nZ=0.;
    temp.nValue=0;



    ic1=0;
    for(i=1; i<=iNatom; i++)
    {
        if(abs(sDelPhiPDB[i].charge)>1.e-6) ic1=ic1+1;
    }// do

/**
 *2011-05-30 Part 2 (more difficult) : from charge distributions
 *running distrTOpointeger without assigning values to  arrays not allocated yet
*/
    if (ndistr>0) {
        cout << "Warning: call distrTOpoFloat2Int(ic1,false)..." << endl;
            //call distrTOpoFloat2Int(ic1,false);
    }
    nqass=ic1;

/**
 *atmcrg contains grid positions of all charges AND the charge in
 *the 4th field atmeps6 contains 6*epsilon/epkt as a function of
 *ic2-th charge internal atmeps6 to the grid NO LONGER USED
 *nqgrdtonqass maps ic2 to ic1 atmeps contains epsilon/epkt
 *as a funcion of ic1-th general charge
*/
    //allocate(atmcrg[nqass],chgpos[nqass]);
    //allocate(crgatn[nqass],atmeps[nqass]);

    atmcrg_v.assign(nqass,temp);
    //chgpos_v.assign(nqass,{0.,0.,0.});
    chgpos_v.assign(nqass,sgrid_temp_real);

    crgatn_v.assign(nqass,0);
    atmeps_v.assign(nqass,0.);
    atmcrg=&atmcrg_v[-1];
    chgpos=&chgpos_v[-1];
    crgatn=&crgatn_v[-1];
    atmeps=&atmeps_v[-1];


    epsdim=iNObject+iNatom+2;

/**
 * find charge moments for dipole approximation
*/
    qnet=0.0;
    qplus=0.0;
    qmin=0.0;
    //cqplus={0.,0.,0.};
    cqplus.nX=0.;
    cqplus.nY=0.;
    cqplus.nZ=0.;

    //cqmin={0.,0.,0.};
    cqmin.nX=0.;
    cqmin.nY=0.;
    cqmin.nZ=0.;


    ic1=0;
    for(ix=1; ix<=iNatom; ix++)
    {
        if (abs(sDelPhiPDB[ix].charge)>1.e-6)
        {
            ic1=ic1+1;
            chrg=sDelPhiPDB[ix].charge;
            atmcrg[ic1].nGrid=xn2[ix];
            atmcrg[ic1].nValue=chrg;
            chgpos[ic1]=xn1[ix];
            crgatn[ic1]=ix;
            qnet=qnet + chrg;

            if (chrg>0.)
            {
                qplus=qplus + chrg;

//2011-05-30 Using operations on coord type variables
// defined in module operators_on_coordinates
                cqplus=cqplus+(chrg*atmcrg[ic1].nGrid);
            }
            else
            {
                qmin=qmin + chrg;
                cqmin=cqmin+(chrg*atmcrg[ic1].nGrid);
            }// if
        }// if
    }// do

#ifdef VERBOSE
        if(!iGaussian==1&&inhomo==1&&logs)cout <<"number of charges coming from molecules " << ic1 << endl;
#endif // VERBOSE


//insert charges from charge distributions
    if (ndistr>0) {
        cout << "Warning: call distrTOpoFloat2Int(ic1,false)..." << endl;
            //call distrTOpoFloat2Int(ic1,false);
    }
//bDebug++++++++++++++++++++++++WWW
/*
    if (false)
    {
        open(52,file='Charges.txt',form='formatted');
        for(iiii=1; iiii<=ic1; iiii++)
        {
            write(52][*) iiii][chgpos[iiii];
        }// do
        close (52);
    }// if
*/
//}// bDebug+++++++++++++++++++

//assign charges for boundary conditions
//ic1 = number of charges

//divide by charge totals
    if (qplus>1.e-6) cqplus=cqplus/qplus;

    if (abs(qmin)>1.e-6) cqmin=cqmin/qmin;

//select those charges which will be charging the grid
//Arrays of correct size are already allocated
    rgrid=iGrid;
    ic2=0;

//2011-06-02 First determine correct size of arrays chgrv2 and
// nqgrdtongass
    for(ix=1; ix<=nqass; ix++)
    {
        if ( optANDGT(atmcrg[ix].nGrid,1.0) && optANDLT(atmcrg[ix].nGrid,rgrid) ) ic2=ic2+1;
    }// do

    nqgrd=ic2;
    ic2=0;

    //allocate(chrgv2(nqgrd),nqgrdtonqass(nqgrd));
    chrgv2_v.assign(nqgrd,temp);
    nqgrdtonqass_v.assign(nqgrd,0);
    chrgv2=&chrgv2_v[-1];
    nqgrdtonqass=&nqgrdtonqass_v[-1];

    for(ix=1; ix<=nqass; ix++)
    {
//crgatn[crg number]=atom number or iNatom+objectnumber or -
//distr.number
        ii=crgatn[ix];

        if (ii<0)
        {
/**
 * now we have to consider charge distributions
 *in this case the distribution is not linked to any object
 *[jx,jy,jz]=coordinates of closest grid pointeger to charge
 *(rx,ry,rz)=coordinates of the charge relatives to the
 *current grid point
*/

//2011-06-02 Using operations on coord and int_coord type
// variables defined in module
// operators_on_coordinates
            jxyz=optCast<delphi_integer,delphi_real>(atmcrg[ix].nGrid+0.5);
            rxyz=atmcrg[ix].nGrid -optCast<delphi_real,delphi_integer>(jxyz);
            jx=jxyz.nX;
            jy=jxyz.nY;
            jz=jxyz.nZ;

            if (rxyz.nZ>rxyz.nX)
            {
                if (rxyz.nZ>-rxyz.nX)
                {
                    if (rxyz.nZ>rxyz.nY)
                    {
                        if (rxyz.nZ>-rxyz.nY)
                        {
                            imed=iepsmp[jz][jy][jx].nZ;
                        }
                        else
                        {
                            imed=iepsmp[jz][jy-1][jx].nY;
                        }// if
                    }
                    else
                    {
                        imed=iepsmp[jz][jy][jx].nY;
                    }// if
                }
                else
                {
                    if (rxyz.nY>rxyz.nX)
                    {
                        if (rxyz.nY>-rxyz.nX)
                        {
                            imed=iepsmp[jz][jy][jx].nY;
                        }
                        else
                        {
                            imed=iepsmp[jz][jy][jx-1].nX;
                        }// if
                    }
                    else
                    {
                        imed=iepsmp[jz][jy-1][jx].nY;
                    }// if
                }// if
            }
            else
            {
                if (rxyz.nZ>-rxyz.nX)
                {
                    if (rxyz.nY>rxyz.nX)
                    {
                        imed=iepsmp[jz][jy][jx].nY;
                    }
                    else
                    {
                        if (rxyz.nY>-rxyz.nX)
                        {
                            imed=iepsmp[jz][jy][jx].nX;
                        }
                        else
                        {
                            imed=iepsmp[jz][jy-1][jx].nY;
                        }// if
                    }// if
                }
                else
                {
                    if (rxyz.nZ>rxyz.nY)
                    {
                        imed=iepsmp[jz][jy-1][jx].nY;
                    }
                    else
                    {
                        if (rxyz.nZ>-rxyz.nY)
                        {
                            imed=iepsmp[jz][jy][jx].nY;
                        }
                        else
                        {
                            imed=iepsmp[jz-1][jy][jx].nZ;
                        }// if
                    }// if
                }// if
            }// if

            imed=imed/epsdim;
        }
        else
        {
            imed=iAtomMed[ii];
        }// if

        atmeps[ix]=medeps[imed];

        if ( optANDGT(atmcrg[ix].nGrid,1.0) && optANDLT(atmcrg[ix].nGrid,rgrid) )
        {
            ic2=ic2+1;
            chrgv2[ic2].nGrid=atmcrg[ix].nGrid;
            chrgv2[ic2].nValue=atmcrg[ix].nValue;
            nqgrdtonqass[ic2]=ix;
        }// if

    }// do

    ipassed=false;

    for(i=1; i<=nqass; i++)
    {
        ii=crgatn[i];

        if (ii>0&&ii<=iNatom)
        {
            if (sDelPhiPDB[ii].radius<=0.)
            {
                ipassed=true;
                cout << ii << " " << delphipdb[ii-1].getAtInf() << " is charged! Radius moved from zero to " << radpolext << endl;
                sDelPhiPDB[ii].radius=radpolext;
            }// if
        }// if
    }// do

    if(ipassed) cout <<"BE CAREFUL!! A WRONG ASSIGNMENT FOR THE RADIUS MIGHT LEAD TO INACCURATE REACTION FIELD ENERGY !!!" << endl;

}// void crgarr;
