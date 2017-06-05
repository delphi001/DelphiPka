/**
 * In:
 *    delphi_integer ibnum = iBndyGridNum;
 *    vector< SGrid<delphi_integer> > iepsmp = prgigEpsMap;
 *    vector<bool> idebmap = prgbDielecMap;
 *    vector< SGrid<delphi_integer> > ibgrd = prgigBndyGrid;
 *
 * return:
 *    delphi_integer                iDielecBndyEven     icount2a   used for realigning idpos and db,compressing to contingous space
 *    delphi_integer                iDielecBndyOdd      icount2b   used for realigning idpos and db,compressing to contingous space
 *    vector<delphi_integer>        prgiBndyDielecIndex idpos(nsp) indices of the dielectric boundary points used to recover the boundary values
 *    vector< vector<delphi_real> > prgfBndyDielec      db(6,nsp)  dielectric values on the boundary
 *    vector<delphi_real>           prgfSaltMap1        sf1(nhgp)  coefficients used in GS/SOR iterations
 *    vector<delphi_real>           prgfSaltMap2        sf2(nhgp)  coefficients used in GS/SOR iterations
 */

/*
 * sample code for creating multi-dimensional vector:
 * int num_of_col = 5;
 * int num_of_row = 9;
 * double init_value = 3.14;
 *
 * vector< vector<double> > matrix;
 * //now we have an empty 2D-matrix of size (0,0). Resizing it with one single command:
 * matrix.resize( num_of col , vector<double>( num_of_row , init_value ) );
 * // and we are good to go ...
 */

#include "solver_fastSOR.h"

//-----------------------------------------------------------------------//
void CDelphiFastSOR::setDielecBndySaltMap()
{
    //++++++++++ INPUT:
    vector< SGrid<delphi_integer> >::const_iterator iepsmp = prgigEpsMap.cbegin();
    vector<bool>::const_iterator idebmap = prgbDielecMap.cbegin();
    //vector< SGrid<delphi_integer> >::const_iterator ibgrd = prgigBndyGrid.cbegin();

    //++++++++++ LOCAL:
    vector<int> it(6,0);          // it(6) rgiDBIndex(6,0)
    vector<delphi_real> vecttemp(6,0.0); // vecttemp rgfDBTemp(6,0.0)
    delphi_integer ieps; // iEps
    delphi_real temp; // fTempVal
    int  deb; // iDielec
    delphi_integer idbs = 0; // iDielecBndyGridNum
    vector<delphi_integer> iepsv; // iepsv(nsp) prgiEpsv

    //++++++++++ OUTPUT:
    delphi_integer& icount2a = iDielecBndyEven;
    delphi_integer& icount2b = iDielecBndyOdd;
    vector<delphi_integer>& idpos = prgiBndyDielecIndex;
    vector< vector<delphi_real> >& db = prgfBndyDielec;
    vector<delphi_real>& sf1 = prgfSaltMap1;
    vector<delphi_real>& sf2 = prgfSaltMap2;

    //-------------------------------- dbsfd -----------------------------//
    /*
     *  create rgfDielecBndyValue(dbval) and rgfSaltFuncDiff(sfd) to store the values
     *  later to be assigned later in db(db) and prgSaltMap(sf1-2)
     */
    delphi_real rgfDielecBndyValue[2][7][2];       // dbval(0:1,0:6,0:1)
    delphi_real rgfSaltFuncDiff[6][2];             // sfd(5,0:1)
    delphi_real fDenom;                            // denom


    //fGepsMap2_v[(i-1)*iGrid*iGrid+(j-1)*iGrid+k-1]=gepsmp2[i][j][k];

    //cout <<   "gepsmp2[5][7][8]:"
    //<< " " << gepsmp2[(5-1)*iGrid*iGrid+(7-1)*iGrid+8-1].nX
    //<< " " << gepsmp2[(5-1)*iGrid*iGrid+(7-1)*iGrid+8-1].nY
    //<< " " << gepsmp2[(5-1)*iGrid*iGrid+(7-1)*iGrid+8-1].nZ
    //<< endl;

    if(debug_solver) cout << "##### in dbsf:  ###" << endl;
    if(debug_solver) cout << "iGaussian " << iGaussian << endl;
    if(iGaussian == 1) bDbOut = false; //for Gaussian

    if (0.0 < fIonStrength)
    {
        for (int iz = 0; iz <= 1; iz += 1)
        {
            for (int iy = 1; iy <= 3; iy += 1)
            {
                fDenom = fSixEps + iy*fEpsDiff + iz*fDebFct;
                rgfDielecBndyValue[0][iy][iz] = 0.0;
                rgfDielecBndyValue[1][iy][iz] = fEpsDiff/fDenom;
                rgfSaltFuncDiff[iy][iz] = fEpsOut/fDenom;
            }
        }

        for (int iz = 0; iz <= 1; iz += 1)
        {
            for (int iy = 4; iy <= 5; iy += 1)
            {
                fDenom = fSixEps + iy*fEpsDiff + iz*fDebFct;
                rgfDielecBndyValue[0][iy][iz] = -fEpsDiff/fDenom;
                rgfDielecBndyValue[1][iy][iz] = 0.0;
                rgfSaltFuncDiff[iy][iz] = fEpsOut/fDenom;
            }
        }
    }
    else
    {
        for (int iz = 0; iz <= 1; iz += 1)
        {
            for (int iy = 1; iy <= 5; iy += 1)
            {
                fDenom = fSixEps + iy*fEpsDiff;
                rgfDielecBndyValue[0][iy][iz] = fEpsOut/fDenom - fSixth;
                rgfDielecBndyValue[1][iy][iz] = fEpsIn/fDenom  - fSixth;
            }
        }
    }
    //---------------------------- End of dbsfd.f -----------------------//

    //-------------------------------- mkdbsf ----------------------------//
    int ix,iy,iz;
    delphi_integer iw;

    //------ Due to the nature of vector in C++, realignment of idpos, db, sf1 and sf2 is unnecessary.
    vector<delphi_integer>::iterator idposEven = idpos.end();
    vector< vector<delphi_real> >::iterator dbEven = db.end();


    vector<delphi_integer>::iterator iepsvEven = iepsv.end();

    string strDbFile = "db.dat"; // db file name
    ofstream ofDbFileStream;

    if (bDbOut) //in fortran: idbwrt
    {
        ofDbFileStream.open(strDbFile.c_str());
        ofDbFileStream << fixed << setprecision(3);
        ofDbFileStream << "DELPHI DB FILE"  << endl;
        ofDbFileStream << "FORMAT NUMBER=1" << endl;
        ofDbFileStream << "NUMBER OF BOUNDARY POINTS= " << iBndyGridNum << endl;
    }


    for (vector< SGrid<delphi_integer> >::const_iterator const_itr = prgigBndyGrid.cbegin(); const_itr != prgigBndyGrid.cend(); ++const_itr)
    {

        ix=const_itr->nX;
        iy=const_itr->nY;
        iz=const_itr->nZ;



        if (0 == iDirectEpsMap)
        {
            it.assign(6,0);

            iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);
            if (0 != (iepsmp+iw)->nX/iEpsDim) it[0] = 1;
            if (0 != (iepsmp+iw)->nY/iEpsDim) it[1] = 1;
            if (0 != (iepsmp+iw)->nZ/iEpsDim) it[2] = 1;

            iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1-1);
            if (0 != (iepsmp+iw)->nX/iEpsDim) it[3] = 1;

            iw = (iz-1)*iGrid*iGrid + (iy-1-1)*iGrid + (ix-1);
            if (0 != (iepsmp+iw)->nY/iEpsDim) it[4] = 1;

            iw = (iz-1-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);
            if (0 != (iepsmp+iw)->nZ/iEpsDim) it[5] = 1;

            ieps = it[0]+it[1]+it[2]+it[3]+it[4]+it[5];
        }
        else //Here is the key point for Gaussian:
        {
            ieps = 0;
            temp = 0.0;

            delphi_integer iIndex;

            if(iGaussian == 0)
            {
                iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);
                iIndex = (iepsmp+iw)->nX/iEpsDim;
                temp += prgfMediaEps[iIndex];
                vecttemp[0] = prgfMediaEps[iIndex];

                iIndex = (iepsmp+iw)->nY/iEpsDim;
                temp += prgfMediaEps[iIndex];
                vecttemp[1] = prgfMediaEps[iIndex];

                iIndex = (iepsmp+iw)->nZ/iEpsDim;
                temp += prgfMediaEps[iIndex];
                vecttemp[2] = prgfMediaEps[iIndex];

                iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1-1);
                iIndex = (iepsmp+iw)->nX/iEpsDim;
                temp += prgfMediaEps[iIndex];
                vecttemp[3] = prgfMediaEps[iIndex];

                iw = (iz-1)*iGrid*iGrid + (iy-1-1)*iGrid + (ix-1);
                iIndex = (iepsmp+iw)->nY/iEpsDim;
                temp += prgfMediaEps[iIndex];
                vecttemp[4] = prgfMediaEps[iIndex];

                iw = (iz-1-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);
                iIndex = (iepsmp+iw)->nZ/iEpsDim;
                temp += prgfMediaEps[iIndex];
                vecttemp[5] = prgfMediaEps[iIndex];
            }

            else //Gaussian:
            {
                //cout << ix << " " << iy << " " << iz << endl;
                iw = (ix-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1);
                temp+=gepsmp2[iw].nX;
                vecttemp[0]=gepsmp2[iw].nX;

                temp+=gepsmp2[iw].nY;
                vecttemp[1]=gepsmp2[iw].nY;

                temp+=gepsmp2[iw].nZ;
                vecttemp[2]=gepsmp2[iw].nZ;

                iw = (ix-1-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1);
                temp+=gepsmp2[iw].nX;
                vecttemp[3]=gepsmp2[iw].nX;

                iw = (ix-1)*iGrid*iGrid + (iy-1-1)*iGrid + (iz-1);
                temp+=gepsmp2[iw].nY;
                vecttemp[4]=gepsmp2[iw].nY;

                iw = (ix-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1-1);
                temp+=gepsmp2[iw].nZ;
                vecttemp[5]=gepsmp2[iw].nZ;
                //cout << "ix,iy,iz,temp,vecttemp: " << ix << " " << iy << " " << iz << " " << temp << " "
                //     << vecttemp[0] << " " << vecttemp[1] << " "<< vecttemp[2] << " " << vecttemp[3] << " "
                //     << vecttemp[4] << " " << vecttemp[5] << endl;
            }
        }

        deb = 0;

        //iw = (iz-1-1)*iGrid*iGrid + (iy-1-1)*iGrid + (ix-1);
        //iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);
        iw = (ix-1)*iGrid*iGrid + (iy-1)*iGrid + (iz-1);
        if ( *(idebmap+iw) )
        {
            deb = 1;
            idbs += 1;
        }

        vector<delphi_real> dbrow;
        if (0 == iDirectEpsMap)
        {
            dbrow.push_back( rgfDielecBndyValue[it[3]][ieps][deb] );
            dbrow.push_back( rgfDielecBndyValue[it[0]][ieps][deb] );
            dbrow.push_back( rgfDielecBndyValue[it[4]][ieps][deb] );
            dbrow.push_back( rgfDielecBndyValue[it[1]][ieps][deb] );
            dbrow.push_back( rgfDielecBndyValue[it[5]][ieps][deb] );
            dbrow.push_back( rgfDielecBndyValue[it[2]][ieps][deb] );
        }
        else
        {
            fDenom = temp+deb*fDebFct;

            if (0.0 == fIonStrength)
            {
                dbrow.push_back( vecttemp[3]/fDenom-fSixth );
                dbrow.push_back( vecttemp[0]/fDenom-fSixth );
                dbrow.push_back( vecttemp[4]/fDenom-fSixth );
                dbrow.push_back( vecttemp[1]/fDenom-fSixth );
                dbrow.push_back( vecttemp[5]/fDenom-fSixth );
                dbrow.push_back( vecttemp[2]/fDenom-fSixth );
            }
            else
            {
                dbrow.push_back( vecttemp[3]/fDenom );
                dbrow.push_back( vecttemp[0]/fDenom );
                dbrow.push_back( vecttemp[4]/fDenom );
                dbrow.push_back( vecttemp[1]/fDenom );
                dbrow.push_back( vecttemp[5]/fDenom );
                dbrow.push_back( vecttemp[2]/fDenom );
            }
        }

        iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1); // iw=isgrid*(k-1) + igrid*(j-1) + i

        if (0 == iw%2) // iw is even
        {
            idposEven = idpos.begin() + icount2a;
            idposEven = idpos.insert(idposEven,iw/2+1);

            iepsvEven = iepsv.begin() + icount2a;
            iepsvEven = iepsv.insert(iepsvEven,ieps);
            iepsvEven++;

            dbEven    = db.begin() + icount2a;
            dbEven    = db.insert(dbEven,dbrow);
            dbEven++;

            icount2a += 1;
        }
        else // iw is odd
        {
            idpos.push_back( (iw+1)/2 );
            iepsv.push_back(ieps);
            db.push_back(dbrow);
        }

        // cout << "Lin Li dbsf: ix: " << iw << " db[iw][0] " << db[iw][0] << " " << db[iw][1] << " "
        //    << db[iw][2] << " " << db[iw][3] << " "
        //    << db[iw][4] << " " << db[iw][5] << endl;

        icount2b += 1;
        //cout << "###icount2a,icount2b: " << icount2a << " " << icount2b << endl;

        if (bDbOut)
        {
            ofDbFileStream << setw(3) << left << ix << " " << setw(3) << left << iy << " " << setw(3) << left << iz << " "
                           << setw(8) << left << dbrow[0] << setw(8) << left << dbrow[1] << setw(8) << left << dbrow[2]
                           << setw(8) << left << dbrow[3] << setw(8) << left << dbrow[4] << setw(8) << left << dbrow[5] << endl;
        }

        // dbs=dbs+db(1,ibnum3) is NOT used

    } //---------- end of for (delphi_integer ix = 0; ix < iBndyGridNum; ix += 1)

    if (bDbOut) ofDbFileStream.close();

#ifdef VERBOSE
    cout << "no. dielectric boundary points in salt = " << idbs << endl;
#endif

    //---------- realign idpos and db,compressing to contingous space

    //---------- set saltmaps 1 and 2, i.e., sf1 and sf2.
    if (0.0 < fIonStrength)
    {
        const delphi_real fSixSalt = fSixth*(1.0/(1.0+fDebFct/fSixEps)-1.0);// sixsalt

        iw = 0;

        for (vector<bool>::const_iterator const_itr = prgbDielecMap.cbegin(); const_itr != prgbDielecMap.cend(); ++const_itr)
        {
            if (0 == iw%2) // even pts
            {
                if (*const_itr) sf1.push_back(fSixth+fSixSalt);
                else            sf1.push_back(fSixth);
            }
            else // odd pts
            {
                if (*const_itr) sf2.push_back(fSixth+fSixSalt);
                else            sf2.push_back(fSixth);
            }

            iw += 1;
        }

        sf2.push_back(0.0); //sf[iGrid3]

#ifdef DEBUG_DELPHI_SOLVER_MKDBSF1
        {
            string strTestFile = "test_mkdbsf1.dat";
            ofstream ofTestStream(strTestFile.c_str());
            ofTestStream << boolalpha;
            ofTestStream << fixed << setprecision(7);

            int count = 1;
            for (vector<delphi_real>::const_iterator const_itr = sf1.cbegin(); const_itr != sf1.cend(); ++const_itr)
            {
                ofTestStream << "sf1(" << setw(6) << right << count << ") = " << setw(11) << right << *const_itr << endl;
                count += 1;
            }

            count = 1;
            for (vector<delphi_real>::const_iterator const_itr = sf2.cbegin(); const_itr != sf2.cend(); ++const_itr)
            {
                ofTestStream << "sf2(" << setw(6) << right << count << ") = " << setw(11) << right << *const_itr << endl;
                count += 1;
            }

            ofTestStream.close();
        }
#endif // DEBUG_DELPHI_SOLVER_MKDBSF1

        for (int i = 0; i < icount2a; i++)
        {
            if (0 != iDirectEpsMap)
                sf1[idpos[i]-1] = 0.0;
            else
            {
                if (abs(sf1[idpos[i]-1]-fSixth) < fZero)
                    sf1[idpos[i]-1] = rgfSaltFuncDiff[ iepsv[idpos[i]-1] ][0];
                else
                    sf1[idpos[i]-1] = rgfSaltFuncDiff[ iepsv[idpos[i]-1] ][1];
            }
        }

        for (int i = icount2a; i < icount2b; i++)
        {
            if (0 != iDirectEpsMap)
                sf2[idpos[i]-1] = 0.0;
            else
            {
                if (abs(sf2[idpos[i]-1]-fSixth) < fZero)
                    sf2[idpos[i]-1] = rgfSaltFuncDiff[ iepsv[ idpos[i]-1] ][0];
                else
                    sf2[idpos[i]-1] = rgfSaltFuncDiff [iepsv[ idpos[i]-1] ][1];
            }
        }
    }
    //---------------------------- end of mkdbsf.f -----------------------//

#ifdef DEBUG_DELPHI_SOLVER_MKDBSF
    {
        string strTestFile = "test_mkdbsf.dat";
        ofstream ofTestStream(strTestFile.c_str());
        ofTestStream << boolalpha;
        ofTestStream << fixed << setprecision(7);

        ofTestStream << "icount2a = " << setw(6) << right << icount2a << " icount2b = " << setw(6) << right << icount2b << endl;

        int count = 1;
        for (vector<delphi_integer>::const_iterator const_itr = idpos.cbegin(); const_itr != idpos.cend(); ++const_itr)
        {
            ofTestStream << "idpos(" << setw(6) << right << count << ") = "  << setw(6) << right << *const_itr << endl;
            count += 1;
        }

        count = 1;
        for (vector< vector<delphi_real> >::const_iterator const_itr = db.cbegin(); const_itr != db.cend(); ++const_itr)
        {
            ofTestStream << "db("  << setw(6) << right << count << ") = "
                         << setw(11) << right << const_itr->at(0) << setw(11) << right << const_itr->at(1)
                         << setw(11) << right << const_itr->at(2) << setw(11) << right << const_itr->at(3)
                         << setw(11) << right << const_itr->at(4) << setw(11) << right << const_itr->at(5) << endl;
            count += 1;
        }

        count = 1;
        for (vector<delphi_real>::const_iterator const_itr = sf1.cbegin(); const_itr != sf1.cend(); ++const_itr)
        {
            ofTestStream << "sf1(" << setw(6) << right << count << ") = " << setw(11) << right << *const_itr << endl;
            count += 1;
        }

        count = 1;
        for (vector<delphi_real>::const_iterator const_itr = sf2.cbegin(); const_itr != sf2.cend(); ++const_itr)
        {
            ofTestStream << "sf2(" << setw(6) << right << count << ") = " << setw(11) << right << *const_itr << endl;
            count += 1;
        }

        ofTestStream.close();
    }
#endif // DEBUG_DELPHI_SOLVER_MKDBSF

    if(debug_solver) cout << "### out of dbsf... ### " << endl;

}

