#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"

#include "space.h"


using namespace std;

void CDelphiSpace::epsmak()
{
    if(debug_space) cout << "######## This is in epsmak: ##########" << endl;

    SGrid <delphi_real> amin,amax;
    //SGrid <delphi_real> sgrid_temp_real;
    SExtrema<delphi_real> sextrema_temp;

    delphi_real fRMid,fRMaxTemp;
    delphi_integer i,ix,iy,iz;




    sextrema_temp.nMax=sgrid_temp_real;
    sextrema_temp.nMin=sgrid_temp_real;

    i=iGrid;

    //sLimGridUnit.assign(iNObject, { {0.,0.,0.} , {0.,0.,0.} } );
    sLimGridUnit.assign(iNObject, sextrema_temp );
/**
 * here limobject is expressed in grid units
*/
    fRMid=(iGrid+1)/2.0;

    //cout << "iNObject: " << iNObject << endl;

    //cout << "sLimGridUnit: " << sLimGridUnit[0].nMin << endl;
    //cout << "sLimObject: " << sLimObject[0] << endl;

    for(i=0; i<=iNObject-1; i++)
    {
        sLimGridUnit[i].nMin=(sLimObject[i].nMin-cOldMid)*fScale+fRMid;
        sLimGridUnit[i].nMax=(sLimObject[i].nMax-cOldMid)*fScale+fRMid;
    }

    if(bUniformDiel)
    {
        cout << "not going to calculate boundary elements since" << endl;
        cout << "uniform dielectric" << endl;
        iBoundNum=0;
        return;
    }

/**
 * lepsx.y.z and uepsx.y.z should be the upper and lower limits of
 * the expanded box. if the molecule is smaller than this then
 * reduce leps and upeps accordingly note leps/ueps not yet defined..
*/
    //2011-05-10 Converted to SGrid <float> derived type
    amin=sLimGridUnit[0].nMin;
    amax=sLimGridUnit[0].nMax;
//    cout << "amaxf,RMaxTemp: "<< amax << " " << fRMaxTemp << endl;


/**
 * find global limits IN GRID UNITS, both, molecule and objects, are considered
*/
    if(iNObject > 1)
    {
        for(i=1; i<=iNObject-1; i++)
        {
            cout << "to be finished using SGrid <float> operations:" << endl;
        }
    }
    fRMaxTemp=rdmx;


    if(fIonStrenth !=0 )
    {
        fRMaxTemp=max(fRMaxTemp,fExternRadius);
    }

    fRMaxTemp=fRMaxTemp*fScale;

    //Using operations on SGrid <float> type variables defined
    //in module operators_on_coordinates

    amin=amin-fRMaxTemp;
    amax=amax+fRMaxTemp;
//    cout << "amaxf,RMaxTemp: "<< amax << " " << fRMaxTemp << endl;

    //Using operations on SGrid <float> and int_coord
    //type variables defined in module operators_on_coordinates
    //SGrid<int> test = optCast<delphi_real>(amin);
//    cout << " LimEps.nMax,amax: "<< LimEps.nMax << " " << amax;
//    cout << "amax,amin: " << amax << amin << endl;


//    cout << "optCast <delphi_integer,delphi_real> (amin)-2: " << iamin << endl;

    //determine the limit of limEps, no larger than igrdd and no less than 1.
    LimEps.nMin=optMax(optCast <delphi_integer,delphi_real> (amin)-2,1);
    LimEps.nMax=optMin(optCast <delphi_integer,delphi_real> (amax)+2,iGrid);
//    cout << "LimEps" << LimEps << endl;


/**
 *point is out of any kind of object (probably in solution)
 **if radprb is less than half of grid spacing, then use old
 *algorithm (sri 29March 93)
 *The new algorithm should be able to handle all scales; still
 *remains to be tested (Sri Apr 95)
*/
    if(iGaussian==0)
    {
        if(debug_space) cout << "go to setout..." << endl;
        setout();
    }
    else
    {
        if(debug_space) cout << "go to setgaussian..." << endl;
        setGaussian();
    }


    if(iGaussian==0)
    {
        if(debug_space) cout <<"going to VdwToMs..." << endl;

         VdwToMs();

    }

/*
    cout << "print time 3:" << endl;
    cout << "amin.nX:" << amin.nX << endl;
    cout << "amin.nY:" << amin.nY << endl;
    cout << "iGrid:" << iGrid << endl;
    cout << "fRMid:" << fRMid << endl;
    cout << "sLimGridUnit[0].nMax.nX:" << sLimGridUnit[0].nMax.nX << endl;
    cout << "sLimGridUnit[0].nMin.nY:" << sLimGridUnit[0].nMin.nY << endl;
    cout << "cOldMid.nX:" << cOldMid.nX << endl;
    cout << "bUniformDiel:" << bUniformDiel << endl;
*/


    if(debug_space && iGaussian==0)
    {
        // testing iepsmp and idebmap:
        ofstream epsfile;
        ofstream debfile;


        epsfile.open ("epsmap_win.txt");
        debfile.open ("debmap_win.txt");
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

#ifdef LIN
        epsfile.open ("epsmap_win_kji.txt");
        debfile.open ("debmap_win_kji.txt");
        for (iz=1; iz<=iGrid; iz++)
        {
            for (iy=1; iy<=iGrid; iy++)
            {
                for (ix=1; ix<=iGrid; ix++)
                {
                    epsfile << setw(5) << iepsmp[ix][iy][iz].nX
                            << setw(5) << iepsmp[ix][iy][iz].nY
                            << setw(5) << iepsmp[ix][iy][iz].nZ
                            << endl;
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
#endif // LIN

        epsfile.close();
        debfile.close();

    }


/*
#ifdef DEBUG
        {
            //open(52,file='iepsmapnewfirst',form='formatted');
            ofstream file1;
            file1.open("iepsmap.txt");

            //write (52,*)iGrid,iGrid*iGrid;
            file1 << "TEST" << endl;
            for(ix=1; ix<=iGrid; ix++)
            {
                for(iy=1; iy<=iGrid; iy++)
                {
                    iz=(iGrid+1)/2;
                    file1 << ix << " " << iy << " " << iz << " " << iEpsMap[ix][iy][iz].nZ << endl;
                }// do
            }// do
            //close (52);
            file1.close();
#endif
*/


}





