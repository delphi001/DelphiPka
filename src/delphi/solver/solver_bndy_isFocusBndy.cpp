/*
 * solver_bndy_isFocusBndy.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: chuan
 */

#include "../misc/misc_interpl.h"
#include "solver_fastSOR.h"

bool CDelphiFastSOR::isFocusBndy(delphi_real *** phimap)
{
    char toblbl[20],label[10],title[60],botlbl[16];
    delphi_real fScale1,goff,gmid,goff1;
    SGrid<delphi_real> fgBoxCenter1,fgXYZ,fcXYZ,fgTemp;
    delphi_integer iGrid1,isgrid;
    int i,ix,iy,iz,iout,flag_s,flag_f;//flag_s: start; flag_f: finish.
    delphi_real *** phimap_pre;
    string strLine,strSubLine;
    SGrid<delphi_real> oldmid_1;

#if !defined(MCCE) && !defined(PRIME)

    if(debug_solver)cout << "#### in isFocusBndy: phiintype: " << phiintype << endl;

    if(phiintype==0)
    {
        ifstream ifPhiiFileStream;
        ifPhiiFileStream.open(strPhiiFile.c_str());
        if (!ifPhiiFileStream.is_open())
        {
            throw CUnknownPhiiFile(strPhiiFile);
            return false;
        }

        cout << "\n";
        cout << "focussing boundary condition\n";
        cout << "read from file\n";
        cout << strPhiiFile << endl << endl;

        ifPhiiFileStream.read(toblbl,21);
        ifPhiiFileStream.read(label ,10);
        ifPhiiFileStream.read(title ,54);

        //cout << "###title: " << title << endl;

        for (iz = 0; iz < iGrid; iz++)
        {
            for (iy = 0; iy < iGrid; iy++)
            {
                for (ix = 0; ix < iGrid; ix++)
                {
                    ifPhiiFileStream.read(reinterpret_cast<char*>(&phimap[iz][iy][ix]),sizeof(delphi_real));
                    //cout << "iz,iy,iz,phi: " << iz << iy << ix << phimap[iz][iy][ix] << endl;
                }
            }
        }

        ifPhiiFileStream.read(botlbl,16);

        ifPhiiFileStream.read(reinterpret_cast<char*>(&fScale1),sizeof(fScale1));
        ifPhiiFileStream.read(reinterpret_cast<char*>(&fgBoxCenter1),sizeof(fgBoxCenter1));
        ifPhiiFileStream.read(reinterpret_cast<char*>(&iGrid1),sizeof(iGrid1));

        ifPhiiFileStream.close();

    }
    else if(phiintype==1) // cube format
    {
        ifstream ifPhiiFileStream;
        ifPhiiFileStream.open(strPhiiFile.c_str());
        if (!ifPhiiFileStream.is_open())
        {
            throw CUnknownPhiiFile(strPhiiFile);
            return false;
        }

        cout << "\n";
        cout << "focussing boundary condition\n";
        cout << "read from file\n";
        cout << strPhiiFile << endl << endl;

        getline(ifPhiiFileStream,strLine);
        //cout << strLine << endl;


        strSubLine = strLine.substr(0,9);
        //cout << strSubLine << endl;
        //exit(0);
        fScale1 = atof( strSubLine.c_str() );

        strSubLine = strLine.substr(10,15);
        iGrid1=atoi( strSubLine.c_str() );

        strSubLine = strLine.substr(16,25);
        oldmid_1.nX = atof( strSubLine.c_str() );

        strSubLine = strLine.substr(26,35);
        oldmid_1.nY = atof( strSubLine.c_str() );

        strSubLine = strLine.substr(36,45);
        oldmid_1.nZ = atof( strSubLine.c_str() );

        fgBoxCenter1=oldmid_1;

        cout << fScale1 << " " << iGrid1 << " " << oldmid_1.nX << " " << oldmid_1.nY << " " << oldmid_1.nZ<< endl;
        cout << fScale << " " << iGrid << " " << fgBoxCenter.nX << " " << fgBoxCenter.nY << " " << fgBoxCenter.nZ<< endl;

        //--------- allocate memory for previous phimap -------
        phimap_pre_v.assign(iGrid1*iGrid1*iGrid1,0.0);
        phimap_pre = pdc->getKey_Ptr<delphi_real>("phimap_pre",iGrid1,iGrid1,iGrid1); // pointer to 3D phimap

        for(i=1; i<=6; i++)
        {
            getline(ifPhiiFileStream,strLine);
            cout << strLine << endl;
        }


        flag_s=0;
        for (iz = 0; iz < iGrid1; iz++)
        {
            for (iy = 0; iy < iGrid1; iy++)
            {
                flag_s=0;
                while (flag_s<iGrid1)
                {
                    getline(ifPhiiFileStream,strLine);
                    if(ifPhiiFileStream.eof()) break;
                    if(flag_s+5 > iGrid1-1)
                    {
                        flag_f=iGrid1-1;
                    }
                    else
                    {
                        flag_f=flag_s+5;
                    }
                    i=0;//for ith number in a line
                    for(ix=flag_s; ix<=flag_f; ix++)
                    {
                        strSubLine = strLine.substr(13*i,13);
                        i++; //for ith number in a line
                        //cout << "strSubLine: " << strSubLine << endl;
                        //phimap_pre_v[iz*iGrid1*iGrid1+iy*iGrid1+ix] = atof( strSubLine.c_str() );
                        phimap_pre[iz][iy][ix] = atof( strSubLine.c_str() );
                        //cout << phimap_pre_v[iz*iGrid1*iGrid1+iy*iGrid1+ix] << " ";
                    }
                    flag_s=flag_s+6;
                    //cout << endl;

                }
            }
        }

        ifPhiiFileStream.close();

    }
    else
    {
        cout << "Focusing is not CUBE or PHI" << endl;
        exit(0);
    }
#endif

#ifdef MCCE

    phimap_pre_v = pmcce->phimap;
    fScale1      = pmcce->scale1;
    fgBoxCenter1 = pmcce->oldmid1;
    iGrid1       = pmcce->igrid1;

    phimap_pre   = pdc->getKey_Ptr<delphi_real>("phimap_pre",iGrid1,iGrid1,iGrid1);

#endif

#ifdef PRIME

    phimap_pre_v   = pPrime->phimap;
    fScale1        = pPrime->scale1;
    fgBoxCenter1   = pPrime->oldmid1;
    iGrid1         = pPrime->igrid1;

    phimap_pre = pdc->getKey_Ptr<delphi_real>("phimap_pre",iGrid1,iGrid1,iGrid1);

#endif

    if (fZero > abs(fScale1-fScale))
    {
        cout << "scales are the same.\n";
        cout << "therefore assuming this to be a continuence\n";
    }
    else
    {
        cout << "\n";
        cout << " focussing potential map:\n";
        if (debug_solver)   cout << title << endl;
        cout << "original scale (grids/A)      : " << fScale1 << endl;
        cout << "object centre at (A) : " << fgBoxCenter1.nX << " "
             << fgBoxCenter1.nY << " " << fgBoxCenter1.nZ << endl;

        //----- check to see that new grid lies within old one that is going to provide bc's
        iout = 0;
        goff = (iGrid+1.0)/2.0;
        goff1 = (iGrid1+1.0)/2.0;

        for (iz = 0; iz < iGrid; iz += iGrid-1)
        {
            for (iy = 0; iy < iGrid; iy += iGrid-1)
            {
                for (ix = 0; ix < iGrid; ix += iGrid-1)
                {
                    fgXYZ.nX = (delphi_real)(ix+1);
                    fgXYZ.nY = (delphi_real)(iy+1);
                    fgXYZ.nZ = (delphi_real)(iz+1);

                    //for each new grid corner, calculate old grid coords
                    fcXYZ = (fgXYZ-goff)/fScale+fgBoxCenter;
                    fgTemp = (fcXYZ-fgBoxCenter1)*fScale1+goff1;

                    if (optORLE<delphi_real>(fgTemp,1.0) || optORGE<delphi_real>(fgTemp,(delphi_real)iGrid1)) iout = 1;
                }
            }
        }

        if (0 != iout) throw COutsideFocus(fScale1,fgBoxCenter1,fScale,fgBoxCenter);

        /*
         * for each boundary point
         *    convert to delphi_real coordinates
         *    convert to old grid coordinates
         *    interpolate potential
         * note that can use same potential array for boundaries since old potentials at boundary are not used for new ones
         */

        //----- save new grid size, and set temporarily to 65
        isgrid = iGrid;
        iGrid  = iGrid1;
        gmid = (isgrid+1.0)/2.0;
        cout << "pulling boundary values out of old potential map...\n";
//        cout << "phimap_pre[0]: " << phimap_pre[0][0][0] << endl;
//        cout << "phimap_pre[1][1][2]: " << phimap_pre[0][0][1] << endl;

        for (iz = 0; iz < isgrid; iz++)
        {
            for (iy = 0; iy < isgrid; iy++)
            {
                for (ix = 0; ix < isgrid; ix += isgrid-1)
                {
                    fgXYZ.nX = (delphi_real)(ix+1);
                    fgXYZ.nY = (delphi_real)(iy+1);
                    fgXYZ.nZ = (delphi_real)(iz+1);
                    fcXYZ = (fgXYZ-gmid)/fScale+fgBoxCenter;
                    fgTemp = (fcXYZ-fgBoxCenter1)*fScale1+goff1;

                    //for each new grid side, calculate old grid coords find potential
                    //phimap[iz][iy][ix] = interpl(iGrid,phimap_pre,fgTemp);
                    phimap[ix][iy][iz] = interpl(iGrid,phimap_pre,fgTemp);
                }
            }
        }

        for (iz = 0; iz < isgrid; iz++)
        {
            for (iy = 0; iy < isgrid; iy += isgrid-1)
            {
                for (ix = 0; ix < isgrid; ix++)
                {
                    fgXYZ.nX = (delphi_real)(ix+1);
                    fgXYZ.nY = (delphi_real)(iy+1);
                    fgXYZ.nZ = (delphi_real)(iz+1);
                    fcXYZ = (fgXYZ-gmid)/fScale+fgBoxCenter;
                    fgTemp = (fcXYZ-fgBoxCenter1)*fScale1+goff1;

                    //for each new grid side, calculate old grid coords find potential
                    //phimap[iz][iy][ix] = interpl(iGrid,phimap_pre,fgTemp);
                    phimap[ix][iy][iz] = interpl(iGrid,phimap_pre,fgTemp);
                }
            }
        }

        for (iz = 0; iz < isgrid; iz += isgrid-1)
        {
            for (iy = 0; iy < isgrid; iy++)
            {
                for (ix = 0; ix < isgrid; ix++)
                {
                    fgXYZ.nX = (delphi_real)(ix+1);
                    fgXYZ.nY = (delphi_real)(iy+1);
                    fgXYZ.nZ = (delphi_real)(iz+1);
                    fcXYZ = (fgXYZ-gmid)/fScale+fgBoxCenter;
                    fgTemp = (fcXYZ-fgBoxCenter1)*fScale1+goff1;

                    //for each new grid side, calculate old grid coords find potential
                    //phimap[iz][iy][ix] = interpl(iGrid,phimap_pre,fgTemp);
                    phimap[ix][iy][iz] = interpl(iGrid,phimap_pre,fgTemp);
                }
            }
        }

        iGrid = isgrid;
    }

    if (debug_solver)
    {
        cout << "phimap[0][0][0]: " << phimap[0][0][0]<< endl;
        cout << "phimap[1][1][1]: " << phimap[1][1][1]<< endl;
        cout << "phimap[1][1][0]: " << phimap[1][1][0]<< endl;
        cout << "phimap[1][0][0]: " << phimap[1][0][0]<< endl;

        cout << "iGrid: " << iGrid << endl;
        cout << "phimap[iGrid-1][iGrid-1][iGrid-1]: " << phimap[iGrid-1][iGrid-1][iGrid-1]<< endl;
        cout << "phimap[1][1][1]: " << phimap[1][1][1]<< endl;
        cout << "phimap[1][1][iGrid-1]: " << phimap[1][1][iGrid-1]<< endl;
        cout << "phimap[1][iGrid-1][iGrid-1]: " << phimap[1][iGrid-1][iGrid-1]<< endl;
    }

    /** remove phimap_pre pointers*/
    for(int i = 0; i != iGrid1; ++i)
    {
        delete[] phimap_pre[i];
    }
    delete[] phimap_pre;

    //exit(0);

    return true;
}
