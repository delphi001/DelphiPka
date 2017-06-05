#ifndef SOLVER_FASTSOR_H
#define SOLVER_FASTSOR_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <memory>
#include <cmath>      // std::abs
//#include <deque>

#include "../interface/interface_abstractmodule.h"
#include "../misc/misc_timer.h"
#include "../delphi/delphi_constants.h"
#include "../io/io.h"
#include "solver_exceptions.h"

using namespace std;

class CDelphiFastSOR:virtual public IAbstractModule
{
private:                                            // In DATA CONTAINER
    shared_ptr<CTimer> pTimer;

    /*********************************************************************************************
     *                                                                                           *
     *              references to the variables obtained from the data container                 *
     *                                                                                           *
     ********************************************************************************************/

    //++++++++++++++ const references to read-only variables from data container +++++++++++++++//
    //----- uniform parameters
    const string&  strBioModel;                      // biomodel
    const string&  strNumSolver;                     // solver
    //----- set by Statements
    const delphi_integer& iAtomNum;                  // natom
    const bool&    bCrgInterplateType;               // isph true-spherical false-cubic
    const bool&    bSpectralRadius;                  // iuspec
    const delphi_real&    fSpectralRadius;           // uspec
    const vector<bool>& rgbPeriodicBndy;             // iper
    const int&     iNonIterateNum;                   // nnit
    const delphi_real&    fIonRadius;                       // exrad
    const bool&    bFixedRelaxParam;                 // icheb
    const bool&    bAutoConverge;                    // iautocon
    const delphi_real&    fGridConverge;             // gten
    const delphi_real&    fRmsc;                     // res1
    const delphi_real&    fMaxc;                     // res2
    const bool&    bLogGraph;                        // igraph
    const bool&    bLogPotential;                    // ipoten
    const bool&    bManualRelaxParam;                // imanual
    //----- io file names
    const string&  strEpsFile;                       // epsnam
    const string&  strPhiiFile;                      // phiinam
    //----- set by functions

    const bool&    bGridCrgOut;                      // iwgcrg
    const bool&    bEpsOut;                          // epswrt
    const bool&    bIonsEng;                         // logions
    const int&     iGaussian;                        // Gaussian
    const int&     inhomo;                           // inhomo
    //----- set by DelPhi
    const delphi_real&    fEpsOut;                   // epsout
    const delphi_real&    fDebyeLength;              // deblen
    const delphi_real&    fScale;                    // scale
    const delphi_real&    fEpsIn;                    // epsin
    const delphi_real&    fIonStrength;              // rionst
    const int&     iDirectEpsMap;                    // idirectalg
    const delphi_real&    fEPKT;                     // epkt
    const SGrid<delphi_real>& fgBoxCenter;           // oldmid
    const SGrid<delphi_real>& fgPotentialDrop;       // vdrop
    const delphi_real&    fTaylorCoeff2;             // chi2
    const delphi_real&    fTaylorCoeff3;             // chi3
    const delphi_real&    fTaylorCoeff4;             // chi4
    const delphi_real&    fTaylorCoeff5;             // chi5
    const bool&           uniformdiel;               //uniformdiel for Gaussian
    //----- set by IO class
    const delphi_integer& iMediaNum;                        // nmedia
    const delphi_integer& iObjectNum;                       // nobject
    const vector<delphi_real>& prgfMediaEps;                // medeps(0:nmediamax)
    //----- set by Surface class
    const delphi_integer& iCrgGridNum;                      // nqass
    const delphi_real&    fMinusCrg;                        // qmin
    const delphi_real&    fPlusCrg;                         // qplus
    const SGrid<delphi_real>& fgPlusCrgCenter;              // cqplus
    const SGrid<delphi_real>& fgMinusCrgCenter;             // cqmin
    const delphi_integer& iBndyGridNum;                     // ibnum
    const vector< SGrid<delphi_integer> >& prgigBndyGrid;   // ibgrd(ibnum)
    const vector< SGrid<delphi_integer> >& prgigEpsMap;     // iepsmp(igrid,igrid,igrid)
    const vector< SGrid<delphi_real> >& gepsmp2;        // gepsmap(igrid,igrid,igrid) LinLi: Gaussian
    const vector<bool>& prgbDielecMap;                      // idebmap(igrid,igrid,igrid)
    const delphi_integer& iCrg2GridNum;                     // nqgrd
    const vector< SGridValue<delphi_real> >& prggvCrg2Grid; // chrgv2(nqgrd)
    const vector<delphi_integer>& prgiCrg2GridMap;          // nqgrdtonqass(nqgrd)
    const vector<delphi_real>& prgfAtomEps;                 // atmeps(nqass)
    const vector< SGridValue<delphi_real> >& prggvCrgedAtom;// atmcrg(nqass)
    //const vector< SGrid<delphi_real> >& prgfgSurfCrgA;    // scspos(ibnum)

    //++++++++++++++++ reference to read-and-write variables from data container +++++++++++++++//
    delphi_integer& iGrid;                                  // igrid (modified in setFocusBndy)
    int&     iBndyType;                              // ibctyp (modified in setBndy)
    delphi_real&    fNetCrg;                                // qnet (modified in setCoulombBndy and setDipolarBndy)
    int&     iLinIterateNum;                         // nlit
    int&     iIterateInterval;                       // icon1
    int&     iConvergeFract;                         // icon2
    bool&    bDbOut;                           // idbwrt
    delphi_real&    fRelaxParam;                            // relpar
    //----- out into data container
    delphi_integer& iDielecBndyOdd;                         // icount2b
    delphi_integer& iCrgedGridSum;                          // icount1b
    vector<delphi_real>& prgfGridCrg;                       // gchrg(icount1b)
    vector< SGrid<delphi_integer> >& prgigGridCrgPose;      // gchrgp(icount1b)
    delphi_integer& iCrgBndyGridNum;                        // ibc
    vector< SDoubleGridValue >& prgdgvCrgBndyGrid;   // cgbp(ibc)
    vector<delphi_real>& prgfPhiMap;                        // phimap(igrid,igrid,igrid)
    vector<delphi_real>& phimap_pre_v;                    //phimap_pre
    /*********************************************************************************************
     *                                                                                           *
     *                            variables defined in this class                                *
     *                                                                                           *
     ********************************************************************************************/

    //+++++++++++++++++++++++++++++++++ const class variables ++++++++++++++++++++++++++++++++++//
    const delphi_integer iTotalGridNum;                     // ngp=igrid*igrid*igrid+1
    const delphi_integer iHalfGridNum;                      // nhgp=ngp/2
    const delphi_integer iBndyDielecNum;                    // nsp=2*(ibnum+1)
    const delphi_real    fDebFct;                           // debfct
    const delphi_real    fEpsDiff;                          // difeps
    const delphi_real    fSixEps;                           // sixeps
    const delphi_integer iEpsDim;                           // epsdim
    const int     nxran;                             // nxran
    const int     nyran;                             // nyran
    bool debug_solver;
    int phiintype;                                   // for focusing
    //+++++++++++++++++++++++++++++ local variables in this class ++++++++++++++++++++++++++++++//
    //------ Using either std:vector or std:deque we can construct idpos, db etc. w/o realignment
    delphi_integer iDielecBndyEven;                         // icount2a
    vector<delphi_integer> prgiBndyDielecIndex;             // idpos(nsp)
    vector< vector<delphi_real> > prgfBndyDielec;           // db(6,nsp) <--> prgfBndyDielec[iBndyDielecNum][6]
    vector<delphi_real> prgfSaltMap1;                       // sf1(nhgp)
    vector<delphi_real> prgfSaltMap2;                       // sf2(nhgp)

    delphi_integer iCrgedGridEven;                          // icount1a
    vector<delphi_integer> prgiCrgPose;                     // iqpos(icount1b)
    vector<delphi_real> prgfCrgValA;                        // qval(icount1b)
    vector<delphi_real> prgfCrgValG;                        // gval(icount1b)

    delphi_real fSpec;                                      // spec
    vector<delphi_integer> ibndx,ibndy,ibndz;
    delphi_integer idif1x,idif2x,inc1xa,inc1xb,inc2xa,inc2xb;
    delphi_integer idif1y,idif2y,inc1ya,inc1yb,inc2ya,inc2yb;
    delphi_integer idif1z,idif2z,inc1za,inc1zb,inc2za,inc2zb;
    vector<delphi_integer> sta1,sta2,fi1,fi2;
    delphi_integer lat1,lat2,long1,long2;
    vector<delphi_real> phimap1,phimap2;
    vector<delphi_real> bndx1,bndx2,bndx3,bndx4;
    vector<delphi_real> qmap1,qmap2;
    vector<delphi_real> debmap1,debmap2;
    delphi_real om1,om2,om3,om4,sixth;

#ifdef MCCE
    SMCCE* pmcce;
#endif
    
#ifdef PRIME
    SPrime* pPrime;
#endif
    
    void setDielecBndySaltMap(); // subroutine mkdbsf()

    void setCrg(); // subroutine setcrg()

    void setBndy(); // subroutine setbc()

    //+++++ choices of boundary conditions
    bool isDipolarBndy(delphi_real *** phimap); // ibctyp = 2

    bool isFocusBndy(delphi_real *** phimap);   // ibctyp = 3

    bool isCoulombBndy(delphi_real *** phimap); // ibctyp = 4

    void initOddEvenItr(const int& forWhom);

    void itrOddPoints(const int& forWhom, const int&);

    void itrEvenPoints(const int& forWhom, const int&);

    void conplt(const vector<delphi_real>& array,const string& title,const int& iclr,const int& iscl,const int& imk,
                const int& iplt,const char symb,const int& ixmin,const int& ixmax,vector<string>& iplot,delphi_real& ymin,delphi_real& ymax);

    void postItr(const vector<delphi_real>& rmaxl,const vector<delphi_real>& rmsl);

    delphi_real calculateRelaxFactor(); // subroutine relfac()

    void itit(); // subroutine itit()

    void nitit(const delphi_real& qfact); // subroutine nitit(qfact)

    shared_ptr<IDataContainer> solver_pdc;

public:
    CDelphiFastSOR(shared_ptr<IDataContainer> pdc,shared_ptr<CTimer> pt):
/*********************************************************************************************
 *                                                                                           *
 *              references to the variables obtained from the data container                 *
 *                                                                                           *
 ********************************************************************************************/

        //++++++++++++++ const references to read-only variables from data container +++++++++++++++//
        IAbstractModule(pdc),
        pTimer(pt),
        //----- uniform parameters
        strBioModel(pdc->getKey_constRef<string>("biomodel")),
        strNumSolver(pdc->getKey_constRef<string>("solver")),
        //----- set by Statments
        iAtomNum(pdc->getKey_constRef<delphi_integer>("natom")),
        bCrgInterplateType(pdc->getKey_constRef<bool>("isph")),
        bSpectralRadius(pdc->getKey_constRef<bool>("iuspec")),
        fSpectralRadius(pdc->getKey_constRef<delphi_real>("uspec")),
        rgbPeriodicBndy(pdc->getKey_constRef< vector<bool> >("iper")),
        iNonIterateNum(pdc->getKey_constRef<int>("nnit")),
        fIonRadius(pdc->getKey_constRef<delphi_real>("exrad")),
        bFixedRelaxParam(pdc->getKey_constRef<bool>("icheb")),
        bAutoConverge(pdc->getKey_constRef<bool>("iautocon")),
        fGridConverge(pdc->getKey_constRef<delphi_real>("gten")),
        fRmsc(pdc->getKey_constRef<delphi_real>("res1")),
        fMaxc(pdc->getKey_constRef<delphi_real>("res2")),
        bLogGraph(pdc->getKey_constRef<bool>("igraph")),
        bLogPotential(pdc->getKey_constRef<bool>("ipoten")),
        bManualRelaxParam(pdc->getKey_constRef<bool>("imanual")),
        phiintype(pdc->getKey_constRef<int>("phiintype")),
        //----- io file names
        strEpsFile(pdc->getKey_constRef<string>("epsnam")),
        strPhiiFile(pdc->getKey_constRef<string>("phiinam")),
        //----- set by functions
        bDbOut(pdc->getKey_Ref<bool>("idbwrt")),
        bGridCrgOut(pdc->getKey_constRef<bool>("iwgcrg")),
        bEpsOut(pdc->getKey_constRef<bool>("epswrt")),
        bIonsEng(pdc->getKey_constRef<bool>("logions")),
        iGaussian(pdc->getKey_constRef<int>("gaussian")),
        inhomo(pdc->getKey_constRef<int>("inhomo")),
        //----- set by DelPhi
        fEpsOut(pdc->getKey_constRef<delphi_real>("epsout")),
        fDebyeLength(pdc->getKey_constRef<delphi_real>("deblen")),
        fScale(pdc->getKey_constRef<delphi_real>("scale")),
        fEpsIn(pdc->getKey_constRef<delphi_real>("epsin")),
        fIonStrength(pdc->getKey_constRef<delphi_real>("rionst")),
        iDirectEpsMap(pdc->getKey_constRef<int>("idirectalg")),
        fEPKT(pdc->getKey_constRef<delphi_real>("epkt")),
        fgBoxCenter(pdc->getKey_constRef< SGrid<delphi_real> >("oldmid")),
        fgPotentialDrop(pdc->getKey_constRef< SGrid<delphi_real> >("vdrop")),
        fTaylorCoeff2(pdc->getKey_constRef<delphi_real>("chi2")),
        fTaylorCoeff3(pdc->getKey_constRef<delphi_real>("chi3")),
        fTaylorCoeff4(pdc->getKey_constRef<delphi_real>("chi4")),
        fTaylorCoeff5(pdc->getKey_constRef<delphi_real>("chi5")),
        uniformdiel(pdc->getKey_constRef<bool>("uniformdiel")),
        //----- set by IO class
        iMediaNum(pdc->getKey_constRef<delphi_integer>("nmedia")),
        iObjectNum(pdc->getKey_constRef<delphi_integer>("nobject")),
        prgfMediaEps(pdc->getKey_constRef< vector<delphi_real> >("medeps")),
        //----- set by Surface class
        iCrgGridNum(pdc->getKey_constRef<delphi_integer>("nqass")),
        fMinusCrg(pdc->getKey_constRef<delphi_real>("qmin")),
        fPlusCrg(pdc->getKey_constRef<delphi_real>("qplus")),
        fgPlusCrgCenter(pdc->getKey_constRef< SGrid<delphi_real> >("cqplus")),
        fgMinusCrgCenter(pdc->getKey_constRef< SGrid<delphi_real> >("cqmin")),
        iBndyGridNum(pdc->getKey_constRef<delphi_integer>("ibnum")),
        prgigBndyGrid(pdc->getKey_constRef< vector< SGrid<delphi_integer> > >("ibgrd")),
        prgigEpsMap(pdc->getKey_constRef< vector< SGrid<delphi_integer> > >("iepsmp")),
        gepsmp2(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("gepsmp2")),
        prgbDielecMap(pdc->getKey_constRef< vector<bool> >("idebmap")),
        iCrg2GridNum(pdc->getKey_constRef<delphi_integer>("nqgrd")),
        prggvCrg2Grid(pdc->getKey_constRef< vector< SGridValue<delphi_real> > >("chrgv2")),
        prgiCrg2GridMap(pdc->getKey_constRef< vector<delphi_integer> >("nqgrdtonqass")),
        prgfAtomEps(pdc->getKey_constRef< vector<delphi_real> >("atmeps")),
        prggvCrgedAtom(pdc->getKey_constRef< vector<SGridValue<delphi_real> > >("atmcrg")),
        //prgfgSurfCrgA(pdc->getKey_constRef< vector< SGrid<delphi_real> > >("scspos")),

        //++++++++++++++++ reference to read-and-write variables from data container +++++++++++++++//
        iGrid(pdc->getKey_Ref<delphi_integer>("igrid")), // modified in setFocusBndy
        iBndyType(pdc->getKey_Ref<int>("ibctyp")), // modified in setBndy
        fNetCrg(pdc->getKey_Ref<delphi_real>("qnet")), // modified in setCoulombBndy and setDipolarBndy
        iLinIterateNum(pdc->getKey_Ref<int>("nlit")),
        iIterateInterval(pdc->getKey_Ref<int>("icon1")),
        iConvergeFract(pdc->getKey_Ref<int>("icon2")),
        fRelaxParam(pdc->getKey_Ref<delphi_real>("relpar")),
        //----- out into data container
        iDielecBndyOdd(pdc->getKey_Ref<delphi_integer>("icount2b")),
        iCrgedGridSum(pdc->getKey_Ref<delphi_integer>("icount1b")), // modified in setcrg
        prgfGridCrg(pdc->getKey_Ref< vector<delphi_real> >("gchrg")), // modified in setcrg
        prgigGridCrgPose(pdc->getKey_Ref< vector< SGrid<delphi_integer> > >("gchrgp")), // modified in setcrg
        iCrgBndyGridNum(pdc->getKey_Ref<delphi_integer>("ibc")), // modified in setcrg
        prgdgvCrgBndyGrid(pdc->getKey_Ref< vector<SDoubleGridValue> >("cgbp")), // modified in setcrg
        prgfPhiMap(pdc->getKey_Ref< vector<delphi_real> >("phimap")),
        phimap_pre_v(pdc->getKey_Ref< vector<delphi_real> >("phimap_pre")),

/*********************************************************************************************
 *                                                                                           *
 *                            variables defined in this class                                *
 *                                                                                           *
 ********************************************************************************************/

        //+++++++++++++++++++++++++++++++++ const class variables ++++++++++++++++++++++++++++++++++//
        iTotalGridNum(iGrid*iGrid*iGrid+1),
        iHalfGridNum(iTotalGridNum/2),
        iBndyDielecNum(2*(iBndyGridNum+1)),
        fDebFct(fEpsOut/(fDebyeLength*fScale*fDebyeLength*fScale)),
        fEpsDiff(fEpsIn-fEpsOut),
        fSixEps(fEpsOut*6.0),
        iEpsDim(iAtomNum+iObjectNum+2),
        nxran(60),
        nyran(20)
    {
#ifdef DEBUG_OBJECT
        cout << endl;
        cout << "****************************************************************\n";
        cout << "*               CDelphiFastSOR is constructed                  *\n";
        cout << "****************************************************************\n";
#endif

        if (0 != strBioModel.compare("PBE") && 0 != strNumSolver.compare("DELPHI"))
            throw CUnknownBioModelSolver(strBioModel,strNumSolver);

        //++++++++++++++++++++++++++++ local variables in this class ++++++++++++++++++++++++++++//
        iDielecBndyEven = 0;
        iCrgedGridEven  = 0;
        fSpec           = 0.0;
        prgfPhiMap.assign(iGrid*iGrid*iGrid,0.0);
        qmap1.assign(iHalfGridNum,0.0);
        qmap2.assign(iHalfGridNum,0.0);

        solver_pdc=pdc;
    };

    ~CDelphiFastSOR()
    {
#ifdef DEBUG_OBJECT
        cout << endl;
        cout << "****************************************************************\n";
        cout << "*                CDelphiFastSOR is destroyed                   *\n";
        cout << "****************************************************************\n";
#endif
    };

#ifdef MCCE
    void getMCCE(SMCCE* mcce_data)
    {
        pmcce = mcce_data;
    }
#endif
    
#ifdef PRIME
    void getPRIME(shared_ptr<SPrime> param)
    {
        pPrime = &*param;
    }
#endif

    virtual void validateInput();

    virtual void run();
};

#endif // SOLVER_FASTSOR_H
