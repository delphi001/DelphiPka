/**
 * \return:
 *    delphi_integer                  iCrgedGridEven    icount1a         number of charged grid points in odd/even arries, to be used in making qval and iqpos
 *    vector<delphi_integer>          prgiCrgPose       iqpos(icount1b)
 *    vector<delphi_real>             prgfCrgValA       qval(icount1b)   charge value in electron units assigned to each grid point
 *    vector<delphi_real>             prgfCrgValG       gval(icount1b)   charge value of each such charge on the grid
 *
 *    delphi_integer                  iCrgedGridSum     icount1b         number of charged grid points in odd/even arries, to be used in making qval and iqpos
 *    vector<delphi_real>             prgfGridCrg       gchrg(icount1b)  gchrg is the fractional charge in electron units assigned to each grid point
 *    vector< SGrid<delphi_integer> > prgigGridCrgPose  gchrgp(icount1b) gchrgp is the position of each such charge on the grid
 *    delphi_integer                  iCrgBndyGridNum   ibc              number of charged boundary grid points
 *    vector<SDoubleGridValue> prgdgvCrgBndyGrid cgbp(ibc)        information on the charged boundary grid points
 *
 */

#include "solver_fastSOR.h"

//-----------------------------------------------------------------------//
void CDelphiFastSOR::setCrg()
{
    //++++++++++ INPUT:
    //const delphi_integer& nqgrd = iCrg2GridNum;
    const vector< SGridValue<delphi_real> >& chrgv2 = prggvCrg2Grid;
    const vector< SGrid<delphi_integer> >& iepsmp = prgigEpsMap;
    const vector<bool>& idebmap = prgbDielecMap;
    const int& idirectalg = iDirectEpsMap;
    const vector<delphi_real>& medeps = prgfMediaEps;
    vector<delphi_integer>::const_iterator nqgrdtonqass = prgiCrg2GridMap.cbegin();
    const vector<delphi_real>& atmeps = prgfAtomEps;

    //++++++++++ LOCAL:
    int kb1,kb2,kb3,i,j,k,m,ix,iy,iz;
    delphi_integer iw;
    delphi_real cg1,cg2,cg3,re,deb;
    SGrid<delphi_integer> igCoord;
    vector<delphi_real> phimap(iGrid*iGrid*iGrid,0.0); // phimap(1:igrid,1:igrid,1:igrid) = 0.0

    delphi_integer inearest[6],iext,ibgp,cont,itmp,ico;
    delphi_real temp;
    vector<delphi_real> gchrgd,gchrgtmp;
    SDoubleGridValue fdbGridVal;
    vector<delphi_integer> gchrg2;

    //++++++++++ OUTPUT:
    delphi_integer& icount1a = iCrgedGridEven;
    delphi_integer& icount1b = iCrgedGridSum;
    delphi_integer& ibc      = iCrgBndyGridNum;
    vector<delphi_integer>& iqpos = prgiCrgPose;
    vector<delphi_real>& qval = prgfCrgValA, &gval = prgfCrgValG, &gchrg = prgfGridCrg;
    vector< SGrid<delphi_integer> >& gchrgp = prgigGridCrgPose;
    vector<SDoubleGridValue>& cgbp = prgdgvCrgBndyGrid;


    iext=0;
    ibgp=0;

    for (vector< SGridValue<delphi_real> >::const_iterator it = chrgv2.begin(); it != chrgv2.end(); ++it)
    {
        kb1 = it->nGrid.nX;
        kb2 = it->nGrid.nY;
        kb3 = it->nGrid.nZ;

        for (ix = 0; ix <= 1; ix++)
        {
            i = kb1 + ix;
            cg1 = kb1 - it->nGrid.nX + 1 - ix;
            for (iy = 0; iy <= 1; iy++)
            {
                j = kb2 + iy;
                cg2 = kb2 - it->nGrid.nY + 1 - iy;
                for (iz = 0; iz <= 1; iz++)
                {
                    k = kb3 + iz;
                    cg3 = kb3 - it->nGrid.nZ + 1 - iz;
                    re = abs(cg1*cg2*cg3)*(it->nValue);
                    phimap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)] += re;
                }
            }
        }
    }

    for (k = 2; k < iGrid; k++)
    {
        for (j = 2; j < iGrid; j++)
        {
            for (i = 2; i < iGrid; i++)
            {
                m = (k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1);
                if(abs(phimap[m]) > fZero) // phimap(i,j,k).ne.0.
                {
                    igCoord.nX = i;
                    igCoord.nY = j;
                    igCoord.nZ = k;
                    gchrgp.push_back(igCoord);
                    gchrg.push_back(phimap[m]);
                    phimap[m] = 0.0;

                    //---------- determine how many charged grid points are even (odd in Fortran)
                    if (0 != optSum<delphi_integer>(igCoord)%2) icount1a += 1;
                }
            }
        }
    }

    for (vector< SGridValue<delphi_real> >::const_iterator it = chrgv2.begin(); it != chrgv2.end(); ++it)
    {
        kb1 = it->nGrid.nX;
        kb2 = it->nGrid.nY;
        kb3 = it->nGrid.nZ;

        for (ix = 0; ix <= 1; ix++)
        {
            i = kb1 + ix;
            cg1 = kb1 - it->nGrid.nX + 1 - ix;
            for (iy = 0; iy <= 1; iy++)
            {
                j = kb2 + iy;
                cg2 = kb2 - it->nGrid.nY + 1 - iy;
                for (iz = 0; iz <= 1; iz++)
                {
                    k = kb3 + iz;
                    cg3 = kb3 - it->nGrid.nZ + 1 - iz;
                    re = abs(cg1*cg2*cg3)*(it->nValue);
                    deb = 0.0;
                    if (idebmap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)]) deb = 1.0;
                    phimap[(k-1)*iGrid*iGrid+(j-1)*iGrid+(i-1)] += re/(6.0*atmeps[*nqgrdtonqass-1]+fDebFct*deb); // atmeps[*nqgrdtonqass-1]
                }
            }
        }

        nqgrdtonqass++; // ic1=nqgrdtonqass(ig)
    }

    icount1b = gchrg.size();

    for (vector< SGrid<delphi_integer> >::iterator it = gchrgp.begin(); it != gchrgp.end(); ++it)
    {
        ix = it->nX;
        iy = it->nY;
        iz = it->nZ;
        iw = (iz-1)*iGrid*iGrid+(iy-1)*iGrid+(ix-1);
        gchrgtmp.push_back(phimap[iw]);
    }

    vector<delphi_real>().swap(phimap); // remove the vector phimap. No need in below

    //---------------------------- Begin of wrtcrg.f -----------------------------//
    if (bGridCrgOut)
    {
        string strCrgFile = "crg.dat"; // charge file name
        ofstream ofCrgFileStream;
        ofCrgFileStream.open(strCrgFile.c_str());

        ofCrgFileStream << fixed << setprecision(3);

        ofCrgFileStream << "DELPHI OUTPUT FILE: GRID CHARGE" << endl;
        ofCrgFileStream << "FORMAT NUMBER= 1" << endl;
        ofCrgFileStream << "NUMBER OF CHARGES=" << setw(8) << right << icount1b << endl;

        for (int i = 0; i <= iMediaNum; i++)
            ofCrgFileStream << "DIELECTRIC IN MEDIUM NUMBER " << setw(3) << right << i << " : "
                            << setw(8) << right << prgfMediaEps[i]*fEPKT << endl;

        ofCrgFileStream << "GRID SCALE=" << fScale << endl;

        for (int i = 0; i < icount1b; i++)
        {
            ofCrgFileStream << scientific << setprecision(17) << setw(25) << right << gchrg[i];
            ofCrgFileStream << fixed << setw(8) << right << gchrgp[i].nX
                            << setw(8) << right << gchrgp[i].nY << setw(8) << right << gchrgp[i].nZ << endl;
        }

        ofCrgFileStream.close();
    }
    //---------------------------- End of wrtcrg.f -----------------------------//

    //---------- set up odd/even pointer array, to be used in making qval and iqpos
    i = 1;
    j = icount1a+1;
    for (vector< SGrid<delphi_integer> >::iterator it = gchrgp.begin(); it != gchrgp.end(); ++it)
    {
        if (0 != optSum<delphi_integer>(*it)%2)
        {
            gchrg2.push_back(i);
            i++;
        }
        else
        {
            gchrg2.push_back(j);
            j++;
        }
    }

    //---------- determine denominator at all charged grid points
    if (0 != ibc) ibc = 0;
    ico = 0;
    i = 0;



    for (vector< SGrid<delphi_integer> >::iterator it = gchrgp.begin(); it != gchrgp.end(); ++it)
    {
        ix = it->nX;
        iy = it->nY;
        iz = it->nZ;

        if(iGaussian==0)
        {
            iext=0;
            ibgp=0;

            iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);
            inearest[0] = iepsmp[iw].nX/iEpsDim;
            inearest[1] = iepsmp[iw].nY/iEpsDim;
            inearest[2] = iepsmp[iw].nZ/iEpsDim;

            iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1-1);
            inearest[3] = iepsmp[iw].nX/iEpsDim;

            iw = (iz-1)*iGrid*iGrid + (iy-1-1)*iGrid + (ix-1);
            inearest[4] = iepsmp[iw].nY/iEpsDim;

            iw = (iz-1-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);
            inearest[5] = iepsmp[iw].nZ/iEpsDim;

            if (0 == inearest[0]) iext = 1;

            if (inearest[0] != inearest[5]) ibgp = 1;

            for (cont = 1; cont < 6; cont++)
            {
                if (0 == inearest[cont]) iext = 1;
                if (inearest[cont-1] != inearest[cont]) ibgp = 1;
            }
        } //iGaussian
        deb = 0.0;

        if (idebmap[(iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1)]) deb = 1.0;

        if (0 == idirectalg)
        {
            throw CDirectEpsilonMap(idirectalg);
        }
        else
        {
            if(iGaussian==0)
            {
                temp=0.0;
                iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);
                itmp = iepsmp[iw].nX/iEpsDim;
                temp += medeps[itmp];
                itmp = iepsmp[iw].nY/iEpsDim;
                temp += medeps[itmp];
                itmp = iepsmp[iw].nZ/iEpsDim;
                temp += medeps[itmp];

                iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1-1);
                itmp = iepsmp[iw].nX/iEpsDim;
                temp += medeps[itmp];

                iw = (iz-1)*iGrid*iGrid + (iy-1-1)*iGrid + (ix-1);
                itmp = iepsmp[iw].nY/iEpsDim;
                temp += medeps[itmp];

                iw = (iz-1-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);
                itmp = iepsmp[iw].nZ/iEpsDim;
                temp += medeps[itmp];

                temp += fDebFct*deb;

                gchrgd.push_back(temp);
            }
            else if(iGaussian==1)
            {
                //cout << setprecision(6) << "fEpsDiff,fEpsIn,fEpsOut: " << fEpsDiff << " " << fEpsIn << " " << fEpsOut << endl;

                if(uniformdiel)
                {
                    //temp=epsin*6
                    temp=fEpsIn*6;
                    //gchrgd(i)=temp+debfct*deb
                    gchrgd.push_back(temp+fDebFct*deb);
                }

                else
                {

                    temp=0;
                    //temp=gepsmp2(ix,iy,iz)%x+gepsmp2(ix,iy,iz)%y+gepsmp2(ix,iy,iz)%z &
                    //& +gepsmp2(ix-1,iy,iz)%x+gepsmp2(ix,iy-1,iz)%y+gepsmp2(ix,iy,iz-1)%z
                    iw = (ix-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1);
                    temp+=gepsmp2[iw].nX;
                    temp+=gepsmp2[iw].nY;
                    temp+=gepsmp2[iw].nZ;
                    iw = (ix-1-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1);
                    temp+=gepsmp2[iw].nX;
                    iw = (ix-1)*iGrid*iGrid + (iy-1-1)*iGrid + (iz-1);
                    temp+=gepsmp2[iw].nY;
                    iw = (ix-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1-1);
                    temp+=gepsmp2[iw].nZ;

                    //temp=temp/epkt
                    temp=temp/fEPKT;
                    //gchrgd(i)=temp+ debfct*deb
                    gchrgd.push_back(temp+fDebFct*deb);
                    //cout << "gchrgd: " << temp+fDebFct*deb << endl;
                } //endif

            }

        }


        if (1 == ibgp || 1 == iext)
        {
            ibc += 1;

            if (0 == ibgp) ico += 1;

            fdbGridVal.fgCoord.nX = (delphi_real)ix;
            fdbGridVal.fgCoord.nY = (delphi_real)iy;
            fdbGridVal.fgCoord.nZ = (delphi_real)iz;
            fdbGridVal.fVal1 = gchrgtmp[i]*f4Pi*fScale;
            fdbGridVal.fVal2 = gchrg2[i];

            cgbp.push_back(fdbGridVal);
        }

        i++;
    } //----- end of for (vector< SGrid<delphi_integer> >::iterator it = gchrgp.begin(); it != gchrgp.end(); ++it)

#ifdef VERBOSE
    cout << "no. charged boundary grid points = " << ibc << endl;
    if (0 != ibc) CCrgedPtsInSolution waring(ico);
#endif

    //---------- make qval, fpoh term so potentials will be in kt/e
    qval.assign(icount1b,0.0);
    gval.assign(icount1b,0.0);
    for (i = 0; i < icount1b; i++)
    {
        j = gchrg2[i];
        qval[j-1] = gchrg[i]*(f4Pi*fScale/gchrgd[i]);
        gval[j-1] = gchrg[i];
    }

    vector<delphi_real>().swap(gchrgd); // remove gchrgd, need in below

    iqpos.assign(icount1b,0);
    for (i = 0; i < icount1b; i++)
    {
        j  = gchrg2[i];
        ix = gchrgp[i].nX;
        iy = gchrgp[i].nY;
        iz = gchrgp[i].nZ;
        iw = (iz-1)*iGrid*iGrid+(iy-1)*iGrid+(ix+1);
        iqpos[j-1] = iw/2;
    }

#ifdef DEBUG_DELPHI_SOLVER_SETCRG
    {
        string strTestFile = "test_setcrg.dat";
        ofstream ofTestStream(strTestFile.c_str());
        ofTestStream << boolalpha;
        ofTestStream << fixed << setprecision(7);

        int index;

        ofTestStream << "icount1a = " << setw(6) << right << icount1a << " icount1b = " << setw(6) << right << icount1b
                     << " ibc = " << setw(6) << right << ibc << endl;

        index = 0;
        for (vector< SGrid<delphi_integer> >::iterator it = gchrgp.begin(); it != gchrgp.end(); ++it)
        {
            index += 1;
            ofTestStream << "gchrgp(" << setw(6) << right << index << ") = "
                         << setw(6) << right << it->nX << setw(6) << right << it->nY << setw(6) << right << it->nZ << endl;
        }

        index = 0;
        for (vector<delphi_real>::iterator it = gchrg.begin(); it != gchrg.end(); ++it)
        {
            index += 1;
            ofTestStream << "gchrg(" << setw(6) << right << index << ") = " << setw(12) << right << *it << endl;
        }

        index = 0;
        for (vector<SDoubleGridValue>::iterator it = cgbp.begin(); it != cgbp.end(); ++it)
        {
            index += 1;
            ofTestStream << "cgbp(" << setw(6) << right << index << ") = "
                         << setw(12) << right << it->fgCoord.nX << setw(12) << right << it->fgCoord.nY << setw(12) << right << it->fgCoord.nZ
                         << setw(12) << right << it->fVal1 << setw(12) << right << it->fVal2 << endl;
        }

        index = 0;
        for (vector<delphi_integer>::iterator it = iqpos.begin(); it != iqpos.end(); ++it)
        {
            index += 1;
            ofTestStream << "iqpos(" << setw(6) << right << index << ") = " << setw(6) << right << *it << endl;
        }

        index = 0;
        for (vector<delphi_real>::iterator it = qval.begin(); it != qval.end(); ++it)
        {
            index += 1;
            ofTestStream << "qval(" << setw(6) << right << index << ") = " << setw(12) << right << *it << endl;
        }

        index = 0;
        for (vector<delphi_real>::iterator it = gval.begin(); it != gval.end(); ++it)
        {
            index += 1;
            ofTestStream << "gval(" << setw(6) << right << index << ") = " << setw(12) << right << *it << endl;
        }

        ofTestStream.close();
    }

#endif // DEBUG_DELPHI_SOLVER_SETCRG
}
