/*
 * solver_fastSOR_itit.cpp
 *
 *  Created on: Feb 10, 2014
 *      Author: chuan
 */

#include "solver_fastSOR.h"

void CDelphiFastSOR::itit()
{
    delphi_real rmsch,rmsch2,rmxch,rmxch2;
    int itr,ires;
    delphi_real grden,grdn[5] = {0.0,0.0,0.0,0.0,0.0};
    delphi_integer ix,iy,iz;
    delphi_real maxres = (fRmsc>fMaxc)?fRmsc:fMaxc;
    maxres = (maxres>fGridConverge)?maxres:fGridConverge;
    vector<delphi_real> rmsl(nxran,0.0),rmaxl(nxran,0.0);

    //int flag_itr; //Lin: for debugging

#ifdef VERBOSE
    if (0.0 < fGridConverge) cout << "  rms-change     max change    grid energy    #iterations" << endl;
    else                     cout << "  rms-change     max change       #iterations" << endl;
#endif


    if (0 == iConvergeFract)
    {
        iIterateInterval = 10;
        iConvergeFract = 1;
    }

    if (iIterateInterval > iLinIterateNum) iIterateInterval = iLinIterateNum;

    initOddEvenItr(1); // forWhom = 1

#ifdef DEBUG_DELPHI_SOLVER_ITIT
    {
        string strTestFile = "test_itit.dat";
        ofstream ofTestStream(strTestFile.c_str());
        ofTestStream << boolalpha;
        ofTestStream << fixed << setprecision(7);

        ix = 1;
        for (vector<delphi_integer>::iterator it = ibndx.begin(); it != ibndx.end(); ++it)
        {
            ofTestStream << "ibndx[" << setw(6) << right << ix << "] = " << setw(8) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_integer>::iterator it = ibndy.begin(); it != ibndy.end(); ++it)
        {
            ofTestStream << "ibndy[" << setw(6) << right << ix << "] = " << setw(8) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_integer>::iterator it = ibndz.begin(); it != ibndz.end(); ++it)
        {
            ofTestStream << "ibndz[" << setw(6) << right << ix << "] = " << setw(8) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = bndx1.begin(); it != bndx1.end(); ++it)
        {
            ofTestStream << "bndx1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = bndx2.begin(); it != bndx2.end(); ++it)
        {
            ofTestStream << "bndx2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = bndx3.begin(); it != bndx3.end(); ++it)
        {
            ofTestStream << "bndx3[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = bndx4.begin(); it != bndx4.end(); ++it)
        {
            ofTestStream << "bndx4[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
        {
            ofTestStream << "phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ix = 1;
        for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
        {
            ofTestStream << "phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
            ix++;
        }

        ofTestStream.close();
    }
#endif // DEBUG_DELPHI_SOLVER_ITIT

    //solver_pdc->showMap("test_delphicpp_itit0.dat");
        //for (int ii=1;ii<=700;ii++)
        //{
        //    cout << "Bphimap2: " << setw(6) << right << ii << " " << setw(11) << right << phimap2[ii-1] << endl;
        //}
    itr = 1;
    ires = 0;
    do
    {
            //cout    << " LinLi: 7,7,7; igrid-1,7,7 : " << scientific << prgfPhiMap[(7-1)*iGrid*iGrid+(7-1)*iGrid+7-1] << "   "
            //<< prgfPhiMap[(7-1-1)*iGrid*iGrid+(7-1)*iGrid+(iGrid-1)] << endl;

        rmsch = 0.0;
        rmxch = 0.0;

        /*
         * iterate over odd points
         */
        itrOddPoints(1,itr); // forWhom = 1

        if (bFixedRelaxParam)
        {
            int itr2 = 2*itr-1;
            om3 = 1.0/(1.0-om2*fSpec*0.25);
            if (fZero > om1) om3 = 1.0/(1.0-om2*fSpec*0.5);
            om4 = om3/om2;
            om2 = om3;
            om1 = 1.0-om2;

            if (0.0 < fIonStrength)
            {
                if (1 == itr2%2)
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
                        *it = (*it)*om4;
                }
                else
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
                        *it = (*it)*om4;
                }
            }

            for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
                *it = (*it)*om4;

            for (delphi_integer iy = 0; iy < iDielecBndyOdd; iy++)
                for (delphi_integer ix = 0; ix < 6; ix++)
                    prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix]*om4;

            sixth = sixth*om4;
        }

#ifdef DEBUG_DELPHI_SOLVER_ITIT
        if (1 == itr)
        {
            string strTestFile = "test_itit.dat";
            ofstream ofTestStream(strTestFile.c_str());
            ofTestStream << boolalpha;
            ofTestStream << fixed << setprecision(7);

            ix = 1;
            for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
            {
                ofTestStream << "phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
                ix++;
            }

            ix = 1;
            for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
            {
                ofTestStream << "phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
                ix++;
            }

            ofTestStream.close();
        }
#endif // DEBUG_DELPHI_SOLVER_ITIT

        /*
         * Next update phimap2 using the new phimap1
         */
        itrEvenPoints(1,itr); // forWhom = 1

        if (bFixedRelaxParam)
        {
            int itr2 = 2*itr;
            om3 = 1.0/(1.0-om2*fSpec*0.25);
            if (fZero > om1) om3 = 1.0/(1.0-om2*fSpec*0.5);
            om4 = om3/om2;
            om2 = om3;
            om1 = 1.0-om2;

            if (0.0 < fIonStrength)
            {
                if (1 == itr2%2)
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap1.begin(); it != prgfSaltMap1.end(); ++it)
                        *it = (*it)*om4;
                }
                else
                {
                    for (vector<delphi_real>::iterator it = prgfSaltMap2.begin(); it != prgfSaltMap2.end(); ++it)
                        *it = (*it)*om4;
                }
            }

            for (vector<delphi_real>::iterator it = prgfCrgValA.begin(); it != prgfCrgValA.end(); ++it)
                *it = (*it)*om4;

            for (delphi_integer iy = 0; iy < iDielecBndyOdd; iy++)
                for (delphi_integer ix = 0; ix < 6; ix++)
                    prgfBndyDielec[iy][ix] = prgfBndyDielec[iy][ix]*om4;

            sixth = sixth*om4;
        }

#ifdef DEBUG_DELPHI_SOLVER_ITIT
        if (1 == itr)
        {
            string strTestFile = "test_itit.dat";
            ofstream ofTestStream(strTestFile.c_str());
            ofTestStream << boolalpha;
            ofTestStream << fixed << setprecision(7);

            ix = 1;
            for (vector<delphi_real>::iterator it = phimap1.begin(); it != phimap1.end(); ++it)
            {
                ofTestStream << "phimap1[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
                ix++;
            }

            ix = 1;
            for (vector<delphi_real>::iterator it = phimap2.begin(); it != phimap2.end(); ++it)
            {
                ofTestStream << "phimap2[" << setw(6) << right << ix << "] = " << setw(11) << right << *it << endl;
                ix++;
            }

            ofTestStream.close();
        }
#endif // DEBUG_DELPHI_SOLVER_ITIT

        /*
         * we also save time by only checking convergence every 10 iterations, rather than every single iteration.
         * store phi2 in phi3 to compare against next iteration
         */
        if (iIterateInterval-1 == itr%iIterateInterval) // itr = 9,19,29,...
        {
            for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract)
            {
                prgfPhiMap[ix] = phimap2[ix];
            }
        }

        if (0.0 < fGridConverge)
        {
            grden = 0.0;

            for(ix = 0; ix < iCrgedGridEven; ix++)
            {
                iy = prgiCrgPose[ix];
                grden += phimap1[iy-1]*prgfCrgValG[ix];
            }

            for(ix = iCrgedGridEven; ix < iCrgedGridSum; ix++)
            {
                iy = prgiCrgPose[ix];
                grden += phimap2[iy-1]*prgfCrgValG[ix];
            }

            grdn[itr%5] = grden/2.0;
            if (10 < itr)
            {
                bool igt = true;
                for (int i = 0; i < 5; i++)
                    for (int j = 0; j < 5; j++)
                        if (abs(grdn[j]-grdn[i]) > fGridConverge) igt =false;
                if (igt)
                {
                    cout << grdn[0] << " " << grdn[1] << " " << grdn[2] << " " << grdn[3] << " " << grdn[4] << endl;
                    ires = 1;
                }
            }
        }

        if (0 == itr%iIterateInterval || 1 == ires) //----- check to see if accuracy is sufficient
        {
            delphi_real rnorm2=0.0, temp2;

            for (ix = 1; ix < iHalfGridNum; ix += iConvergeFract)
            {
                temp2   = prgfPhiMap[ix] - phimap2[ix];
                rnorm2 += temp2*temp2;

                //if(rmxch<abs(temp2)) flag_itr=ix; //Lin Li: to see which grid is the key grid
                rmxch   = max(rmxch,abs(temp2));

            }

            rmsch   = sqrt((delphi_real)iConvergeFract*rnorm2/((iGrid-2)*(iGrid-2)*(iGrid-2)));
            //rnormch = sqrt(rnorm2);
            rmsch2  = rmsch;
            rmxch2  = rmxch;

#ifdef VERBOSE
            if (0.0 < fGridConverge)
                cout << scientific << rmsch2 << "  "  << rmxch2 << "  " << grden << "  at  " << setw(5) << left << itr << " iterations\n";
            else
                cout << scientific << rmsch2 << "  "  << rmxch2 << "  at  " << setw(5) << left << itr << " iterations\n";
                //cout << scientific << "grid: " << flag_itr << " " << rmsch2 << "  "  << rmxch2 << "  at  " << setw(5) << left << itr << " iterations\n";
#endif

            if (fRmsc > rmsch || fMaxc > rmxch) ires = 1;

            if (bLogGraph)
            {
                int ibin;
                for (int j = itr-9; j <= itr; j++)
                {
                    ibin = (j-1)*(60-1)/(iLinIterateNum-1)+1;
                    rmsl[ibin-1]  = rmsch;
                    rmaxl[ibin-1] = rmxch;
                }
            }
        }

        //cout << "itr: " << itr << endl;
        /*
        if(itr==2) solver_pdc->showMap("test_delphicpp_itit2.dat");
        if(itr==3) solver_pdc->showMap("test_delphicpp_itit3.dat");
        if(itr==4) solver_pdc->showMap("test_delphicpp_itit4.dat");
        if(itr==5) solver_pdc->showMap("test_delphicpp_itit5.dat");
        if(itr==6) solver_pdc->showMap("test_delphicpp_itit6.dat");
        if(itr==7) solver_pdc->showMap("test_delphicpp_itit7.dat");
        if(itr==8) solver_pdc->showMap("test_delphicpp_itit8.dat");
        if(itr==9) solver_pdc->showMap("test_delphicpp_itit9.dat");
*/


        itr++;

        /*
         * check to see if accuracy is sufficient
         */
        if (1.0e-7 > maxres)
        {
            if (iLinIterateNum >= itr && 0 == ires)
                continue;
            else
                break;
        }
        else
        {
            if ((iLinIterateNum >= itr || bAutoConverge) && (0 == ires))
                continue;
            else
                break;
        }

    }
    while(true);

    postItr(rmaxl,rmsl);

    /*
     * code phimap corner, for use in transference from irises to convex and via versa
     */
    {
        delphi_real ap1,ap2,ap3,ap4;
        ap1 = prgfPhiMap[0];
        ap2 = ap1*10000.0;
        ap3 = (int)ap2;
        if (0 < ap3) ap4 = (ap3+0.8)/10000.0;
        else         ap4 = (ap3-0.8)/10000.0;
        prgfPhiMap[0] = ap4;
    }

#ifdef DEBUG_DELPHI_SOLVER_ITIT
    {
        string strTestFile = "test_itit.dat";
        ofstream ofTestStream(strTestFile.c_str());
        ofTestStream << boolalpha;
        ofTestStream << fixed << setprecision(7);

        const delphi_real *** phimap = pdc->getKey_constPtr<delphi_real>("phimap",iGrid,iGrid,iGrid); // const pointer to 3D phimap
        for (iz = 0; iz < iGrid; iz++)
        {
            for (iy = 0; iy < iGrid; iy++)
            {
                for (ix = 0; ix < iGrid; ix++)
                {
                    ofTestStream << "phimap[" << setw(6) << right << iz+1 << "," << setw(6) << right << iy+1 << ","
                                 << setw(6) << right << ix+1 << "] = " << setw(12) << right << phimap[iz][iy][ix] << endl;
                }
            }
        }

        ofTestStream.close();
    }
#endif // DEBUG_DELPHI_SOLVER_ITIT
}
