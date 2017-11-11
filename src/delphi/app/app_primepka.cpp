#include "../delphi/delphi_data.h"
#include "../space/space.h"
#include "../solver/solver_fastSOR.h"
#include "../energy/energy.h"
#include "../site/site.h"
#include "../../delphiPKa/prime_environment.h"

class StreamRedirector
{
public:
    explicit StreamRedirector(std::ios& stream, std::streambuf* newBuf) : savedBuf_(stream.rdbuf()), stream_(stream)
    {
        stream_.rdbuf(newBuf);
    }
    
    ~StreamRedirector()
    {
        stream_.rdbuf(savedBuf_);
    }
    
private:
    std::streambuf* savedBuf_;
    std::ios& stream_;
};



void runDelphi(shared_ptr<SPrime> param)
{
    
#ifdef DELPHI_OUTPUT
    
    string delphi_ntime = "delphi_run_" + to_string(param->ntimes) + ".log";
    
#endif
    
#ifndef DELPHI_OUTPUT
    
    string delphi_ntime = "delphi_run_" + to_string(param->ntimes) + ".log";
    
#endif
    
    ofstream logFile(delphi_ntime);
    StreamRedirector redirect_cout(cout,logFile.rdbuf());
    StreamRedirector redirect_cerr(cerr,logFile.rdbuf());

    cout << boolalpha;
    
    cerr << boolalpha;
    
    
    try
    {
        CTimer * pTester =  new CTimer; // record execution time
        
        pTester->start();
        
        shared_ptr<CTimer> pTimer( new CTimer); // record execution time
        
        
        shared_ptr<IDataContainer> pDataContainer( new CDelphiData(param,pTimer) );

//        pDataContainer->showMap("showmap_atbegin.dat");
        
        int& inhomo(pDataContainer->getKey_Ref<int>("inhomo"));
        const int& iGaussian(pDataContainer->getKey_constRef<int>("gaussian"));
        bool& logs (pDataContainer->getKey_Ref<bool>("logs"));
        bool& logg (pDataContainer->getKey_Ref<bool>("logg"));
        
        inhomo=0;
        
        
        if( iGaussian==1 && logs )
        {
            logg=true; //for gaussian
            inhomo=1;
        }
        
        //********************************************************************************//
        //                                                                                //
        //    realize an object of CDelphiSpace class to construct molecular surfaces     //
        //                                                                                //
        //********************************************************************************//
        
        
        
        unique_ptr<IAbstractModule> pSpace( new CDelphiSpace(pDataContainer,pTimer) );
        
        pSpace->run();
        
        pSpace.reset();
        
//        pDataContainer->showMap("showmap_aftersuf.dat");
        
        if( !(iGaussian==1&&inhomo==0&&logs) )
        {
            cout << " number of atom coordinates read  :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_integer>("natom") << endl;
        }
        
        if (pDataContainer->getKey_constRef<bool>("isolv"))
        {
            if( !(iGaussian==1&&inhomo==0&&logs) )
            {
                cout << " total number of assigned charges :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_integer>("nqass") << endl;
                cout << " net assigned charge              :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_real>("qnet") << endl;
                cout << " assigned positive charge         :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_real>("qplus") << endl;
                cout << " centred at (gu)                  :" << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqplus").nX << " "
                    << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqplus").nY << " "
                    << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqplus").nZ << endl;
                cout << " assigned negative charge         :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_real>("qmin") << endl;
                cout << " centred at (gu)                  :" << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqmin").nX << " "
                    << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqmin").nY << " "
                    << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqmin").nZ << endl;
                cout << "\nnumber of dielectric boundary points" << right << setw(10) << pDataContainer->getKey_constRef<delphi_integer>("ibnum") << endl;
            }
            if (pDataContainer->getKey_constRef<bool>("iexun") && 0 == pDataContainer->getKey_constRef<delphi_integer>("ibnum"))
                throw CNoBndyAndDielec(pTimer);
            
            //********************************************************************************//
            //                                                                                //
            //   realize an object of CDelphiFastSOR class to calculate potentials on grids   //
            //                                                                                //
            //********************************************************************************//

            
            unique_ptr<CDelphiFastSOR> pSolver( new CDelphiFastSOR(pDataContainer,pTimer) );
            
            if (param->bndcon == 3) {
                
               pSolver->getPRIME(param);

            }
            
            pSolver->run();
            
            pSolver.reset();
 
//            pDataContainer->showMap("showmap_afteritr.dat");
            
            //********************************************************************************//
            //                                                                                //
            //          realize an object of CDelphiEnergy class to calculate energies        //
            //                                                                                //
            //********************************************************************************//

            
            unique_ptr<IAbstractModule> pEnergy( new CDelphiEnergy(pDataContainer,pTimer) );
            
            pEnergy->run();
            
            pEnergy.reset();
            
//            pDataContainer->showMap("showmap_aftereng.dat");
            
            if(iGaussian==1&&inhomo==1&&logs) //second run for Gaussian
            {
                inhomo=0;
                
                unique_ptr<IAbstractModule> pSpace( new CDelphiSpace(pDataContainer,pTimer) );
                pSpace->run();
                pSpace.reset();
                
                unique_ptr<IAbstractModule> pSolver( new CDelphiFastSOR(pDataContainer,pTimer) );
                pSolver->run();
                pSolver.reset();
                

                unique_ptr<IAbstractModule> pEnergy( new CDelphiEnergy(pDataContainer,pTimer) );
                pEnergy->run();
                pEnergy.reset();
                
            }

            
            //********************************************************************************//
            //                                                                                //
            //               realize an object of CSite class to write site info              //
            //                                                                                //
            //********************************************************************************//
            
            
            unique_ptr<CSite> pSite( new CSite(pDataContainer,pTimer) );
            
            if (pDataContainer->getKey_constRef<bool>("isite"))
            {
                int iisitsf = 0;
                if (pDataContainer->getKey_Ref<bool>("isitsf")) iisitsf = 1;
                pSite->writeSite(iisitsf);
            }
            
            if (pDataContainer->getKey_constRef<bool>("phiwrt")) pSite->writePhi();
            
//            pDataContainer->showMap("showmap_aftersite.dat");
            
            /*
             * equivalent to out(frc,file="filename") in the parameter file
             
            if (0 == pSite->prime_grdphiv.size())
            {
                pSite.reset();
                pTimer->exit(); pTimer.reset();
            }
            */
            
            param->strAtomDes = pSite->prime_atomdes;
            param->vecGridPot = pSite->prime_grdphiv;
            param->vecSiteCrg = pSite->prime_crhgv;
            
            
            
            
#ifdef PRIME
            pSite->clearIO();
#endif
            
            pSite.reset();

            
        }
        
        //********************************************************************************//
        //                                                                                //
        //       retrieve the solvation energy and grid energy from data container        //
        //                                                                                //
        //********************************************************************************//
        
        param->ergs = pDataContainer->getKey_Val<delphi_real>("ergs");
        param->ergg = pDataContainer->getKey_Val<delphi_real>("ergg");
        
        if(pDataContainer->getKey_constRef<int>("ibctyp") == 2)
        {
            param->igrid1 = pDataContainer->getKey_constRef<delphi_integer>("igrid");
            param->scale1 = pDataContainer->getKey_constRef<delphi_real>("scale");
            param->oldmid1 = pDataContainer->getKey_constRef<SGrid<delphi_real> >("oldmid");
            param->phimap = pDataContainer->getKey_constRef<vector<delphi_real> >("phimap");
            
//            pDataContainer->showMap("test_phimap");
        }
        
        pDataContainer.reset();
        
        pTimer->exit(); pTimer.reset();
        
        delete pTester;
        
#ifdef DELPHI_OUTPUT
        param->ntimes++;
#endif
        remove(&delphi_ntime[0]);
    } // ---------- end of try block
    catch (CException&)
    {
        cerr << "\n\n ......... PROGRAM ABORTS WITH AN EXCEPTION AND " << CWarning::iWarningNum << " WARNING(S) ........\n\n";
    }
    
    cout << "\n\n .........  PROGRAM EXITS SUCCESSFULLY : WITH TOTAL " << CWarning::iWarningNum << " WARNING(S) ........\n\n";
    
    cout.unsetf(ios_base::floatfield); // return to cout default notation 
    
}




