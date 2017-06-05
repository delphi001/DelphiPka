/*
 * Space_Run.cpp
 *
 */

#include "space.h"

void CDelphiSpace::run()
{

    debug_space=false; //debug_space=true will print a lot of debug information
    if(debug_space) cout << "############### run in Space module.... ##################" << endl;






    //cout << "###: cutoff: " << cutoff << endl;
    //cout << "###: sigma: " << sigma << endl;
    //cout << "###: inhomo: " << inhomo << endl;
    //cout << "###: srfcut: " << srfcut << endl;
    //cout << "###: iGaussian: " << iGaussian << endl;
    //cout << "###: radipz: " << radipz << endl;


//    SGrid <delphi_integer> *** testEpsmap=pdc->getKey_Ptr("iepsmap",iNatom,iNatom,iNatom);
    /*
        cout << "iNatom: " << iNatom << endl;
        cout << "xn1[0].nX: " << xn1[0].nX << endl;
        cout << "delphipdb[0].getRadius: " << delphipdb[0].getRadius() << endl;
        cout << "delphipdb[0].getPose: " << delphipdb[0].getPose().nX << endl;

        //bDebMap.assign(iNatom*iNatom*iNatom,false);
        //iEpsMap.assign(iNatom*iNatom*iNatom, {0,0,0});

        iEpsMap[1][2][3].nX=54;

        egrid_v[0].nX=66;

        cout << "iEpsMap_v[iGrid*iGrid-1].nX " << iEpsMap_v[iGrid*iGrid+2*iGrid+3].nX << endl;
        cout << "iEpsMap[1][0][0].nX " << iEpsMap[1][2][3].nX << endl;
        cout << "bDebMap_v[0]" << bDebMap_v[0] << endl;
        cout << "iEpsMap_v[0]" << iEpsMap_v[0] << endl;
        cout << "egrid_v[0].nX: " << egrid_v[0].nX << endl;
    */

    delphi_integer i,j,k,ix,iy,iz,ic,ico,natom2;
    delphi_real fRMid;
    SGrid <delphi_integer> epstmp;
    SGrid <delphi_real> xl,xr;
    fRMid=(iGrid+1)/2.0;

    sgrid_temp_real.nX=0.;
    sgrid_temp_real.nY=0.;
    sgrid_temp_real.nZ=0.;

    sgrid_temp_int.nX=0;
    sgrid_temp_int.nY=0;
    sgrid_temp_int.nZ=0;

    //iEpsMap=pdc->getKey_Ptr < SGrid <delphi_integer> > ( "iepsmp",iGrid,iGrid,iGrid);

    //##### navigate iepsmp and idebmap, xn1,xn2,fRadProb: ##########
    //idebmap=Move_index_3d <bool> (bDebMap,iGrid,iGrid,iGrid);
    //iepsmp=Move_index_3d <SGrid <delphi_integer> > (iEpsMap,iGrid,iGrid,iGrid);

    //Move_index_3d <bool> (idebmap,bDebMap,iGrid,iGrid,iGrid);


    get_pt3d <bool> (idebmap,iGrid+1,iGrid+1,iGrid+1);
    //######### initialize bDebMap: #########

    for (ix=1; ix<=iGrid; ix++)
    {
        for (iy=1; iy<=iGrid; iy++)
        {
            for (iz=1; iz<=iGrid; iz++)
            {

                idebmap[ix][iy][iz]=true;
            }
        }
    }

    //Move_index_3d <SGrid <delphi_integer> > (iepsmp,iEpsMap,iGrid,iGrid,iGrid);
    if(iGaussian==0)
    {
        get_pt3d <SGrid <delphi_integer> > (iepsmp,iGrid+1,iGrid+1,iGrid+1);
    }
    else
    {
        get_pt3d <SGrid <delphi_real> > (gepsmp,iGrid+1,iGrid+1,iGrid+1);
        get_pt3d <SGrid <delphi_real> > (gepsmp2,iGrid+1,iGrid+1,iGrid+1);
    }

    if(iGaussian==1&&inhomo==0&&logs) //for 2nd Gaussian run
    {

        for(i=1; i<=iGrid; i++)
        {
            //for(i=1;i<=iGrid;i++){
            for(j=1; j<=iGrid; j++)
            {
                for(k=1; k<=iGrid; k++)
                {
                    //gepsmp2[i][j][k]=fGepsMap2_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1];
                    //cout << "LinLi: " << i << " " << j << " " << k << " " << gepsmp2[i][j][k] << endl;
                    gepsmp[i][j][k]=fGepsMap_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1];
                   //if(gepsmp[i][j][k].nX>0.) cout << "LinLi: " << i << " " << j << " " << k << " " << gepsmp[i][j][k] << endl;
                }
            }
        }
    }

    xn1=&xn1_v[-1];
    xn2=&xn2_v[-1];
    fRadPrb=&fRadPrb_v[-1];
    //ibgrd=&ibgrd_v[-1];
    //atmcrg=&atmcrg_v[-1];
    //chgpos=&chgpos_v[-1];
    //crgatn=&crgatn_v[-1];
    //atmeps=&atmeps_v[-1];
    //chrgv2=&chrgv2_v[-1];
    //nqgrdtonqass=&nqgrdtonqass_v[-1];

    //for (i=0; i<=iNatom-1; i++)
    //{
    //    cout << "xn2,xn2_v: " << xn2[i] << " " << xn2_v[i-1] << endl;
    //}

    //#########for focusing #######
    if(debug_space)cout << "######ibctyp #### " << ibctyp << endl;
    if(ibctyp==3)
    {

        delphi_real halfl;
        SGrid <delphi_real> edge_low,edge_high,xyz_temp;
        vector <CAtomPdb> delphipdb2;              //delphipdb2 for focusing


        if(debug_space)cout << "######focusing method: detecting atoms outside box ####" << endl;

        if(debug_space)cout << fScale << " " << iGrid << " " << acenter << endl;
        halfl=(iGrid-1)/fScale/2.0;
        //edge_low=acenter-halfl;
        //edge_high=acenter+halfl;

        edge_low=acenter-halfl-3.5;
        edge_high=acenter+halfl+3.5;
        if(debug_space)cout << "halfl,edge_low,edge_high: " << halfl<< " " <<edge_low << " " << edge_high << endl;

        natom2=0;

        for (i=0; i<iNatom; i++)
        {
            xyz_temp=delphipdb[i].getPose();
            if( optORGT(xyz_temp,edge_high) || optORLT(xyz_temp,edge_low) )
                //if( optORGT(xyz_temp,xr) || optORLT(xyz_temp,xl) )
            {
                //---------- outside box ----------
                //cout << "xyz_temp: " << xyz_temp << endl;
            }
            else
            {
                delphipdb2.push_back(delphipdb[i]);
                natom2++;
                //xyz_temp=delphipdb2[natom2-1].getPose();
                //cout << xyz_temp << " " << delphipdb2.size() << endl;
            }

        }



        if(debug_space)cout << "halfl,edge_low,edge_high: " << halfl<< " " <<edge_low << " " << edge_high << endl;

        if(debug_space)cout << iNatom << " " << natom2 << endl ;

        delphipdb.resize(natom2);
        for (i=0; i<natom2; i++)
        {
            delphipdb[i]=delphipdb2[i];
        }

        iNatom=delphipdb.size();

        if(debug_space)cout << iNatom << " "  << delphipdb2.size() << endl;

        //exit(0);
        vector <CAtomPdb> ().swap(delphipdb2);
    }

//######### initialize sDelPhiPDB: #########

    for (i=0; i<=iNatom-1; i++)
    {
        sDelPhiPDB[i+1].radius=delphipdb[i].getRadius();
        sDelPhiPDB[i+1].xyz=delphipdb[i].getPose();
        sDelPhiPDB[i+1].charge=delphipdb[i].getCharge();

    }
//######## from qdiff4v:  ########
    for (i=1; i<=iNatom; i++)
    {
        xn1[i]=sDelPhiPDB[i].xyz;
        xn2[i]=(sDelPhiPDB[i].xyz-cOldMid)*fScale+fRMid;
        //cout << "cOldMid, fScale, fRMid: " << cOldMid.nX  << " " << cOldMid.nY  << " " << cOldMid.nZ << " " << fScale << " " << fRMid << endl;
        //cout << "i,xn1[i],xn2[i]: " <<i<< xn1[i] << xn2[i] << endl;
    }

    epsmak(); //Lin Li reset


// Now start crgarr


    if (isolv)
    {
//increased the starting dimension of crgatn and..+++++
        extracrg=0;
        if (ndistr>0) extracrg=iGrid*iGrid*iGrid;
        if(debug_space) cout << "ndistr: " << ndistr << endl;
//2011-05-30 Allocation of the arrays below is moved to the body of crgarr void, arrays are accessible
// via pointers module. Sizes of arrays are determined before allocation inside the crgarr void

        crgarr(); //Lin Li reset

        xl=cOldMid-(1.0/fScale)*(iGrid+1)*0.5;
        xr=cOldMid+(1.0/fScale)*(iGrid+1)*0.5;

        if (logs||lognl)
        {
            ico=0;
            for(ic=1; ic<=nqass; ic++)
            {
                if ( optORLT(chgpos[ic],xl) || optORGT(chgpos[ic],xr) )
                {
                    if (crgatn[ic]<0)
                    {
                        cout << "!WARNING: distribution " << -crgatn[ic] << "outside the box" << endl;
                    }
                    else
                    {
                        if (crgatn[ic]>iNatom)
                        {
                            cout << "WARNING:crg " << ic << "object " << crgatn[ic]-iNatom << "outside the box " << chgpos[ic] << endl;
                        }
                        else
                        {
                            cout << "//!! WARNING : charge " << delphipdb[crgatn[ic]-1].getAtInf() << "outside the box" << endl;
                            //if(debug_space) cout << "ico: " << ico  << endl;
                        }// if
                    }// if
                    ico=1;
                }// if
            }// do
            //if(debug_space) cout << "ico, ibctype: " << ico << " " << ibctyp << endl;
            if (ico>0&&ibctyp!=3)
            {
                cout <<"CHARGES OUTSIDE THE BOX AND NOT DOING FOCUSSING << THEREFORE STOP" << endl;
                exit(0);
            }// if
        }// if
    }// if

    if(debug_space) cout << "#### Lin Li: coverting the matrix...." << endl;
    /*
        for(i=1;i<=iGrid-1;i++){
        //for(i=1;i<=iGrid;i++){
           for(j=1;j<=iGrid;j++){
              for(k=i+1;k<=iGrid;k++){
              //for(k=1;k<=iGrid;k++){
    	     //cout <<"i,j,k,iepsmp: "<< i <<" " << j << " " << " " << k << " "  << iepsmp[k][j][i] << endl;
    	     //cout <<"i,j,k,idebmap: "<< i <<" " << j << " " << " " << k << " "  << idebmap[i][j][k] << " "  << idebmap[k][j][i] << endl;

                 epstmp=iepsmp[k][j][i];
                 iepsmp[k][j][i]=iepsmp[k][j][i];
                 iepsmp[i][j][k]=epstmp;

                 //debtmp=idebmap[i][j][k];
                 //idebmap[i][j][k]=idebmap[k][j][i];
                 //idebmap[k][j][i]=debtmp;

    	     //if(i==5&&j==6&&k==1) cout << "i,j,k,iepsmp: " << i << " " << j << " " << k << " " << iepsmp[k][j][i] << endl;
              }
           }
        }
    */
    /*
        // testing..........
        for(i=1; i<=iGrid; i++)
        {        //for(i=1;i<=iGrid;i++){
            for(j=1; j<=iGrid; j++)
            {
                for(k=1; k<=iGrid; k++)
                {
                    //if(iepsmp[i][j][k] != iEpsMap_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1])cout << "iepsmp: " << i << " "<< j<< " "<< k << " " << iepsmp[i][j][k] << " " << iEpsMap_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1] << endl;
                    if(idebmap[i][j][k] != bDebMap_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1])cout << "bdebmap: " << i << " "<< j<< " "<< k << " " << idebmap[i][j][k] << " " << bDebMap_v[(k-1)*iGrid*iGrid+(j-1)*iGrid+i-1] << endl;


                }
            }
        }
    */

    //######### initialize iepsmp: #########
    if(iGaussian==0)
    {
        iEpsMap_v.assign(iGrid*iGrid*iGrid, sgrid_temp_int);
    }
    else
    {
        fGepsMap_v.assign(iGrid*iGrid*iGrid, sgrid_temp_real);
        fGepsMap2_v.assign(iGrid*iGrid*iGrid, sgrid_temp_real);
    }


    if(iGaussian==0)
    {
        for(i=1; i<=iGrid; i++)
        {
            //for(i=1;i<=iGrid;i++){
            for(j=1; j<=iGrid; j++)
            {
                for(k=1; k<=iGrid; k++)
                {
                    iEpsMap_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1]=iepsmp[i][j][k];
                }
            }
        }
    }
    else
    {
        for(i=1; i<=iGrid; i++)
        {
            //for(i=1;i<=iGrid;i++){
            for(j=1; j<=iGrid; j++)
            {
                for(k=1; k<=iGrid; k++)
                {
                    fGepsMap2_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1]=gepsmp2[i][j][k];
                    fGepsMap_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1]=gepsmp[i][j][k];
                }
            }
        }
    }

    for(i=1; i<=iGrid; i++)
    {
        //for(i=1;i<=iGrid;i++){
        for(j=1; j<=iGrid; j++)
        {
            for(k=1; k<=iGrid; k++)
            {
#ifdef IKJ
                bDebMap_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+(k-1)]=idebmap[k][j][i];
#endif // IKJ
                bDebMap_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+(k-1)]=idebmap[i][j][k];

            }
        }
    }


    //free_pt3d <bool> (bDebMap,iGrid,iGrid,iGrid);

    if(iepsmp != NULL) free_pt3d(iepsmp,iGrid+1,iGrid+1,iGrid+1);
    if(idebmap != NULL) free_pt3d(idebmap,iGrid+1,iGrid+1,iGrid+1);
    if(gepsmp != NULL) free_pt3d(gepsmp,iGrid+1,iGrid+1,iGrid+1);
    if(gepsmp2 != NULL) free_pt3d(gepsmp2,iGrid+1,iGrid+1,iGrid+1);

    // free_pt3d_p <bool> (idebmap,iGrid+1,iGrid+1);

    // free_pt3d_p <SGrid <delphi_integer> > (iEpsMap,iGrid+1,iGrid+1);
    //fGepsMap2_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1]=gepsmp2[i][j][k];
    //cout << "gepsmp2[6][7][8]:" << gepsmp2[5][7][8] << endl;
    if(debug_space) cout << "#### going out of space_run... ####" <<endl;

}

