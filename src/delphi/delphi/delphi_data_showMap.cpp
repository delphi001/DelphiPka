/*
 * delphi_data_showMap.cpp
 *
 *  Created on: Apr 20, 2014
 *      Author: chuan
 */

#include "delphi_data.h"

//-----------------------------------------------------------------------//
void CDelphiData::showMap(const string& strMapFile)
{
   //string strMapFile = "delphicpp.dat";

   ofstream ofMapStream(strMapFile.c_str());
   ofMapStream << boolalpha;
   ofMapStream << fixed << setprecision(7);

//   ofMapStream << endl;
//   ofMapStream << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";
//   ofMapStream << "+                                                        + \n";
//   ofMapStream << "+  DELPHI DATA CONTAINER size = " <<myData.size() << ")  + \n";
//   ofMapStream << "+                                                        + \n";
//   ofMapStream << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";

   ofMapStream << "------------------- uniform parameters -----------------"    << endl;
   ofMapStream << "          biomodel : " << getKey_constRef<string>("biomodel")     << endl;
   ofMapStream << "            solver : " << getKey_constRef<string>("solver")       << endl;
   ofMapStream << "------------------- set by statements ------------------"    << endl;
   ofMapStream << "   iautocon(AUTOC) : " << (getKey_constRef<bool>("iautocon")?'T':'F') << endl;
   ofMapStream << "    ibctyp(BNDCON) : " << getKey_constRef<int>("ibctyp")          << endl;
   ofMapStream << "    perfil(PERFIL) : " << getKey_constRef<delphi_real>("perfil")         << endl;
   ofMapStream << "     icheb(CHEBIT) : " << (getKey_constRef<bool>("icheb")?'T':'F')<< endl;
   ofMapStream << "      isrf(CLCSRF) : " << (getKey_constRef<bool>("isrf")?'T':'F') << endl;
   ofMapStream << "     icon2(CONFRA) : " << getKey_constRef<int>("icon2")           << endl;
   ofMapStream << "     icon1(CONINT) : " << getKey_constRef<int>("icon1")           << endl;
   ofMapStream << "     iexun(EXITUN) : " << (getKey_constRef<bool>("iexun")?'T':'F')<< endl;
   ofMapStream << "     repsout(EXDI) : " << getKey_constRef<delphi_real>("repsout")        << endl;
   ofMapStream << "        isph(FCRG) : " << (getKey_constRef<bool>("isph")?'T':'F') << endl;
   ofMapStream << "      gten(GRDCON) : " << getKey_constRef<delphi_real>("gten")           << endl;
   ofMapStream << "      igrid(GSIZE) : " << getKey_constRef<delphi_integer>("igrid")       << endl;
   ofMapStream << "      repsin(INDI) : " << getKey_constRef<delphi_real>("repsin")         << endl;

   {
      const vector<delphi_real> vctfSalt = getKey_constRef< vector<delphi_real> >("conc");
      ofMapStream << "  conc(SALT,SALT2) : " << vctfSalt[0] << " " << vctfSalt[1] << endl;
   }

   ofMapStream << "     exrad(IONRAD) : " << getKey_constRef<delphi_real>("exrad")          << endl;
   ofMapStream << "       nlit(LINIT) : " << getKey_constRef<int>("nlit")            << endl;
   ofMapStream << "    igraph(LOGGRP) : " << (getKey_constRef<bool>("igraph")?'T':'F')<< endl;
   ofMapStream << "    ipoten(LOGPOT) : " << (getKey_constRef<bool>("ipoten")?'T':'F')<< endl;
   ofMapStream << "        res2(MAXC) : " << getKey_constRef<delphi_real>("res2")           << endl;
   ofMapStream << "       nnit(NONIT) : " << getKey_constRef<int>("nnit")            << endl;

   {
      const vector<bool> vctbPeriodicBndy = getKey_constRef< vector<bool> >("iper");
      ofMapStream << "    iper1-3(PBX-Z) : " << (vctbPeriodicBndy[0]?'T':'F') << " "
                                             << (vctbPeriodicBndy[1]?'T':'F') << " "
                                             << (vctbPeriodicBndy[2]?'T':'F') << endl;
      ofMapStream << " iper4-6(VDROPX-Z) : " << (vctbPeriodicBndy[3]?'T':'F') << " "
                                             << (vctbPeriodicBndy[4]?'T':'F') << " "
                                             << (vctbPeriodicBndy[5]?'T':'F') << endl;
   }

   ofMapStream << "     iconc(PHICON) : " << (getKey_constRef<bool>("iconc")?'T':'F') << endl;

   {
      const vector<delphi_real> vctfProbeRadius = getKey_constRef< vector<delphi_real> >("radprb");
      ofMapStream << "    radprb(PRBRAD) : " << vctfProbeRadius[0] << " " << vctfProbeRadius[1] << endl;
   }

   ofMapStream << "     uspec(RELFAC) : " << getKey_constRef<delphi_real>("uspec")          << endl;
   ofMapStream << "    relpar(RELPAR) : " << getKey_constRef<delphi_real>("relpar")         << endl;
   ofMapStream << "        res1(RMSC) : " << getKey_constRef<delphi_real>("res1")           << endl;
   ofMapStream << "      scale(SCALE) : " << getKey_constRef<delphi_real>("scale")          << endl;
   ofMapStream << "     isolv(SOLVPB) : " << (getKey_constRef<bool>("isolv")?'T':'F')<< endl;

   {
      const vector<int> vctiValence1 = getKey_constRef< vector<int> >("ival");
      ofMapStream << " ival(VAL+1,VAL-1) : " << vctiValence1[0] << " " << vctiValence1[1] << endl;
   }

   {
      const vector<int> vctiValence2 = getKey_constRef< vector<int> >("ival2");
      ofMapStream << "ival2(VAL+2,VAL-2) : " << vctiValence2[0] << " " << vctiValence2[1] << endl;
   }

   ofMapStream << "atompotdist(ATPODS): " << getKey_constRef<delphi_real>("atompotdist")    << endl;
   ofMapStream << "temperature(TEMPER): " << getKey_constRef<delphi_real>("temperature")    << endl;
   ofMapStream << "      vdrop(VDROP) : " << getKey_constRef< SGrid<delphi_real> >("vdrop").nX << " "
                                           << getKey_constRef< SGrid<delphi_real> >("vdrop").nY << " "
                                           << getKey_constRef< SGrid<delphi_real> >("vdrop").nZ << endl;
   ofMapStream << "            iuspec : " << (getKey_constRef<bool>("iuspec")?'T':'F') << endl;
   ofMapStream << "           imanual : " << (getKey_constRef<bool>("imanual")?'T':'F') << endl;
   ofMapStream << "------------------- io file names ------------------"        << endl;
   //ofMapStream << "           prmnam : " << getKey_constRef<string>("prmnam")     << endl;
   //ofMapStream << "          centnam : " << getKey_constRef<string>("centnam")    << endl;
   ofMapStream << "            siznam : " << getKey_constRef<string>("siznam")       << endl;
   ofMapStream << "            crgnam : " << getKey_constRef<string>("crgnam")       << endl;
   ofMapStream << "            pdbnam : " << getKey_constRef<string>("pdbnam")       << endl;
   ofMapStream << "            phinam : " << getKey_constRef<string>("phinam")       << endl;
   ofMapStream << "           frcinam : " << getKey_constRef<string>("frcinam")      << endl;
   ofMapStream << "            frcnam : " << getKey_constRef<string>("frcnam")       << endl;
   ofMapStream << "            epsnam : " << getKey_constRef<string>("epsnam")       << endl;
   ofMapStream << "           phiinam : " << getKey_constRef<string>("phiinam")      << endl;
   ofMapStream << "           mpdbnam : " << getKey_constRef<string>("mpdbnam")      << endl;
   ofMapStream << "           updbnam : " << getKey_constRef<string>("updbnam")      << endl;
   ofMapStream << "           ufrcnam : " << getKey_constRef<string>("ufrcnam")      << endl;
   ofMapStream << "            srfnam : " << getKey_constRef<string>("srfnam")       << endl;
   ofMapStream << "            nrgnam : " << getKey_constRef<string>("nrgnam")       << endl;
   ofMapStream << "           scrgnam : " << getKey_constRef<string>("scrgnam")      << endl;
   ofMapStream << "------------------ set by functions ----------------"             << endl;
   // set by CENTER or CENT function:
   ofMapStream << "            offset : " << getKey_constRef< SGrid<delphi_real> >("offset").nX << " "
                                          << getKey_constRef< SGrid<delphi_real> >("offset").nY << " "
                                          << getKey_constRef< SGrid<delphi_real> >("offset").nZ << endl;
   // set by ACENTER or ACENT function
   ofMapStream << "             acent : " << getKey_constRef< SGrid<delphi_real> >("acent").nX << " "
                                          << getKey_constRef< SGrid<delphi_real> >("acent").nY<< " "
                                          << getKey_constRef< SGrid<delphi_real> >("acent").nZ<< endl;
   ofMapStream << "            iacent : " << (getKey_constRef<bool>("iacent")?'T':'F') << endl;
   // set by READ or IN function
   ofMapStream << "            pdbfrm : " << getKey_constRef<int>("pdbfrm")          << endl;
   ofMapStream << "            ipdbrd : " << (getKey_constRef<bool>("ipdbrd")?'T':'F') << endl;
   // set by WRITE or OUT function
   ofMapStream << "            phiwrt : " << (getKey_constRef<bool>("phiwrt")?'T':'F') << endl;
   ofMapStream << "            phifrm : " << getKey_constRef<int>("phifrm")          << endl;
   ofMapStream << "             ibios : " << (getKey_constRef<bool>("ibios")?'T':'F')<< endl;
   ofMapStream << "              ibem : " << (getKey_constRef<bool>("ibem")?'T':'F') << endl;
   ofMapStream << "             isite : " << (getKey_constRef<bool>("isite")?'T':'F')<< endl;
   ofMapStream << "            frcfrm : " << getKey_constRef<int>("frcfrm")          << endl;
   ofMapStream << "            epswrt : " << (getKey_constRef<bool>("epswrt")?'T':'F') << endl;
   ofMapStream << "            iatout : " << (getKey_constRef<bool>("iatout")?'T':'F')<< endl;
   ofMapStream << "           mpdbfrm : " << getKey_constRef<int>("mpdbfrm")         << endl;
   ofMapStream << "           ipdbwrt : " << (getKey_constRef<bool>("ipdbwrt")?'T':'F')<< endl;
   ofMapStream << "           ifrcwrt : " << (getKey_constRef<bool>("ifrcwrt")?'T':'F')<< endl;
   ofMapStream << "           inrgwrt : " << (getKey_constRef<bool>("inrgwrt")?'T':'F')<< endl;
   ofMapStream << "            iwgcrg : " << (getKey_constRef<bool>("iwgcrg")?'T':'F')<< endl;
   ofMapStream << "              iacs : " << (getKey_constRef<bool>("iacs")?'T':'F') << endl;
   ofMapStream << "            idbwrt : " << (getKey_constRef<bool>("idbwrt")?'T':'F') << endl;
   ofMapStream << "              isen : " << (getKey_constRef<bool>("isen")?'T':'F') << endl;
   ofMapStream << "              isch : " << (getKey_constRef<bool>("isch")?'T':'F') << endl;
   ofMapStream << "           scrgfrm : " << getKey_constRef<int>("scrgfrm")         << endl;
   // set by ENERGY function
   ofMapStream << "              logg : " << (getKey_constRef<bool>("logg")?'T':'F') << endl;
   ofMapStream << "              logs : " << (getKey_constRef<bool>("logs")?'T':'F') << endl;
   ofMapStream << "             logas : " << (getKey_constRef<bool>("logas")?'T':'F')<< endl;
   ofMapStream << "              loga : " << (getKey_constRef<bool>("loga")?'T':'F') << endl;
   ofMapStream << "           logions : " << (getKey_constRef<bool>("logions")?'T':'F') << endl;
   ofMapStream << "              logc : " << (getKey_constRef<bool>("logc")?'T':'F') << endl;
   // set by SITE function: all MUST be initialized to to false
   ofMapStream << "             isita : " << (getKey_constRef<bool>("isita")?'T':'F')<< endl;
   ofMapStream << "             isitq : " << (getKey_constRef<bool>("isitq")?'T':'F')<< endl;
   ofMapStream << "             isitp : " << (getKey_constRef<bool>("isitp")?'T':'F')<< endl;
   ofMapStream << "            isitap : " << (getKey_constRef<bool>("isitap")?'T':'F') << endl;
   ofMapStream << "           isitdeb : " << (getKey_constRef<bool>("isitdeb")?'T':'F') << endl;
   ofMapStream << "             isitf : " << (getKey_constRef<bool>("isitf")?'T':'F')<< endl;
   ofMapStream << "             isitr : " << (getKey_constRef<bool>("isitr")?'T':'F')<< endl;
   ofMapStream << "             isitc : " << (getKey_constRef<bool>("isitc")?'T':'F')<< endl;
   ofMapStream << "             isitx : " << (getKey_constRef<bool>("isitx")?'T':'F')<< endl;
   ofMapStream << "             isiti : " << (getKey_constRef<bool>("isiti")?'T':'F')<< endl;
   ofMapStream << "             isitt : " << (getKey_constRef<bool>("isitt")?'T':'F')<< endl;
   ofMapStream << "            isitrf : " << (getKey_constRef<bool>("isitrf")?'T':'F')<< endl;
   ofMapStream << "            isitcf : " << (getKey_constRef<bool>("isitcf")?'T':'F')<< endl;
   ofMapStream << "            isitmd : " << (getKey_constRef<bool>("isitmd")?'T':'F')<< endl;
   ofMapStream << "            isitsf : " << (getKey_constRef<bool>("isitsf")?'T':'F')<< endl;
   ofMapStream << "            isittf : " << (getKey_constRef<bool>("isittf")?'T':'F')<< endl;
   ofMapStream << "           isitpot : " << (getKey_constRef<bool>("isitpot")?'T':'F') << endl;
   ofMapStream << "              irea : " << (getKey_constRef<bool>("irea")?'T':'F') << endl;
   ofMapStream << "             iself : " << (getKey_constRef<bool>("iself")?'T':'F')<< endl;
   // set by BUFFZ function
   ofMapStream << "              bufz : " << getKey_constRef< SExtrema<delphi_integer> >("buffz").nMin.nX << " "
                                          << getKey_constRef< SExtrema<delphi_integer> >("buffz").nMin.nY << " "
                                          << getKey_constRef< SExtrema<delphi_integer> >("buffz").nMin.nZ   << " "
                                          << getKey_constRef< SExtrema<delphi_integer> >("buffz").nMax.nX   << " "
                                          << getKey_constRef< SExtrema<delphi_integer> >("buffz").nMax.nY   << " "
                                          << getKey_constRef< SExtrema<delphi_integer> >("buffz").nMax.nZ   << endl;
   ofMapStream << "             ibufz : " << (getKey_constRef<bool>("ibufz")?'T':'F') << endl;
   // set by SURFACE function
   ofMapStream << "         isurftype : " << getKey_constRef<int>("isurftype")       << endl;
   ofMapStream << "----------------------- DelPhi ---------------------"        << endl;
   ofMapStream << "            deblen : " << getKey_constRef<delphi_real>("deblen")         << endl;
   ofMapStream << "            epsout : " << getKey_constRef<delphi_real>("epsout")         << endl;
   ofMapStream << "              cran : " << getKey_constRef< SGrid<delphi_real> >("cran").nX << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cran").nY << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cran").nZ << endl;
   ofMapStream << "              pmid : " << getKey_constRef< SGrid<delphi_real> >("pmid").nX << " "
                                          << getKey_constRef< SGrid<delphi_real> >("pmid").nY << " "
                                          << getKey_constRef< SGrid<delphi_real> >("pmid").nZ << endl;
   ofMapStream << "            oldmid : " << getKey_constRef< SGrid<delphi_real> >("oldmid").nX << " "
                                          << getKey_constRef< SGrid<delphi_real> >("oldmid").nY << " "
                                          << getKey_constRef< SGrid<delphi_real> >("oldmid").nZ << endl;
   ofMapStream << "            rionst : " << getKey_constRef<delphi_real>("rionst")         << endl;
   ofMapStream << "            chi1-5 : " << getKey_constRef<delphi_real>("chi1") << " "
                                          << getKey_constRef<delphi_real>("chi2") << " "
                                          << getKey_constRef<delphi_real>("chi3") << " "
                                          << getKey_constRef<delphi_real>("chi4") << " "
                                          << getKey_constRef<delphi_real>("chi5") << endl;
   ofMapStream << "             lognl : " << (getKey_constRef<bool>("lognl")?'T':'F')<< endl;
   ofMapStream << "              epkt : " << getKey_constRef<delphi_real>("epkt")           << endl;
   ofMapStream << "             epsin : " << getKey_constRef<delphi_real>("epsin")          << endl;
   ofMapStream << "            ifrcrd : " << (getKey_constRef<bool>("ifrcrd")?'T':'F') << endl;
   ofMapStream << "        idirectalg : " << getKey_constRef<int>("idirectalg")      << endl;
   ofMapStream << "           numbmol : " << getKey_constRef<delphi_integer>("numbmol")     << endl;
   ofMapStream << "              rdmx : " << getKey_constRef<delphi_real>("rdmx")           << endl;
   ofMapStream << "       uniformdiel : " << (getKey_constRef<bool>("uniformdiel")?'T':'F') << endl;

   {
      const vector< SExtrema<delphi_real> >* vctefExtrema = getKey_constPtr< vector< SExtrema<delphi_real> > >("limobject");

      ofMapStream << "         limobject : " << endl;
      if (0 != vctefExtrema->size())
      {
         for (unsigned int i = 0; i < vctefExtrema->size(); i++)
         {
            ofMapStream<<setw(6)<<right<<i<<" : ";

            //ofMapStream<<fixed<<setprecision(3);
            ofMapStream << setw(8) << right << vctefExtrema->at(i).nMin.nX << "  "
                        << setw(8) << right << vctefExtrema->at(i).nMin.nY << "  "
                        << setw(8) << right << vctefExtrema->at(i).nMin.nZ << " || "
                        << setw(8) << right << vctefExtrema->at(i).nMax.nX << "  "
                        << setw(8) << right << vctefExtrema->at(i).nMax.nY << "  "
                        << setw(8) << right << vctefExtrema->at(i).nMax.nZ << endl;
         }
      }
   }

   {
      const vector< SGrid<delphi_real> >* vctgfAtomCoordA = getKey_constPtr< vector< SGrid<delphi_real> > >("xn1");

      ofMapStream << "               xn1 : " << endl;
      if (0 != vctgfAtomCoordA->size())
      {
         for (unsigned int i = 0; i < vctgfAtomCoordA->size(); i++)
         {
            ofMapStream << setw(6) << right << i << " : " << setw(8) << right << vctgfAtomCoordA->at(i).nX << "  "
                        << setw(8) << right << vctgfAtomCoordA->at(i).nY << "  "
                        << setw(8) << right << vctgfAtomCoordA->at(i).nZ << endl;
         }
      }
   }

   {
      const vector< SGrid<delphi_real> >* vctgfAtomCoordG = getKey_constPtr< vector< SGrid<delphi_real> > >("xn2");

      ofMapStream << "               xn2 : " << endl;
      if (0 != vctgfAtomCoordG->size())
      {
         for (unsigned int i = 0; i < vctgfAtomCoordG->size(); i++)
         {
            ofMapStream << setw(6) << right << i << " : " << setw(8) << right << vctgfAtomCoordG->at(i).nX << "  "
                        << setw(8) << right << vctgfAtomCoordG->at(i).nY << "  "
                        << setw(8) << right << vctgfAtomCoordG->at(i).nZ << endl;
         }
      }
   }

   ofMapStream << "------------------ set by IO class -----------------"        << endl;
   ofMapStream << "         resnummax : " << getKey_constRef<delphi_integer>("resnummax")   << endl;
   ofMapStream << "            nmedia : " << getKey_constRef<delphi_integer>("nmedia")      << endl;

   {
      const vector<delphi_real>* vctfMediaEps = getKey_constPtr< vector<delphi_real> >("medeps");

      ofMapStream << "            medeps : " << endl;
      if (0 != vctfMediaEps->size())
      {
         for (unsigned int i = 0; i < vctfMediaEps->size(); i++)
         {
            ofMapStream << setw(6) << right << i << " : " << setw(8) << right << vctfMediaEps->at(i) << endl;
         }
      }
   }

   ofMapStream << "           nobject : " << getKey_constRef<delphi_integer>("nobject")     << endl;

   {
      const vector<string>* vctstrObject = getKey_constPtr< vector<string> >("dataobject");

      ofMapStream << "        dataobject : " << endl;
      if (0 != vctstrObject->size())
      {
         for (unsigned int i = 0; i < vctstrObject->size(); i=i+2)
         {
            ofMapStream << setw(6)  << right << i << " : "
                        << setw(15) << right << vctstrObject->at(i)   << " || "
                        << setw(15) << right << vctstrObject->at(i+1) << endl;
         }
      }
   }

   ofMapStream << "             natom : " << getKey_constRef<delphi_integer>("natom")       << endl;

   {
      const vector<CAtomPdb>* vctapAtomPdb = getKey_constPtr< vector<CAtomPdb> >("delphipdb");

      ofMapStream << "         delphipdb : " << endl;
      if (0 != vctapAtomPdb->size())
      {
         for (unsigned int i = 0; i < vctapAtomPdb->size(); i++)
         {
            ofMapStream<<setw(6)<< right << i << " : " << setw(8) << right << vctapAtomPdb->at(i).getPose().nX << "  "
                                << setw(8) << right << vctapAtomPdb->at(i).getPose().nY << "  "
                                << setw(8) << right << vctapAtomPdb->at(i).getPose().nZ <<" || "
                                << setw(8) << right << vctapAtomPdb->at(i).getRadius()  <<" || "
                                << setw(8) << right << vctapAtomPdb->at(i).getCharge()  <<" || "
                                << setw(15)<< right << vctapAtomPdb->at(i).getAtInf()   << endl;
         }
      }
   }

   {
      const vector<delphi_integer>* vctiAtomMediaNum = getKey_constPtr< vector<delphi_integer> >("iatmmed");

      ofMapStream << "           iatmmed : " << endl;
      if (0 != vctiAtomMediaNum->size())
      {
         for (unsigned int i=0; i < vctiAtomMediaNum->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : " << setw(8) << right << vctiAtomMediaNum->at(i) << endl;
         }
      }
   }

   ofMapStream << "          ionlymol : " << (getKey_constRef<bool>("ionlymol")?'T':'F') << endl;

   ofMapStream << "---------------------- Surface ---------------------"             << endl;
   ofMapStream << "             nqass : " << getKey_constRef<delphi_integer>("nqass")       << endl;
   ofMapStream << "              qnet : " << getKey_constRef<delphi_real>("qnet")           << endl;
   ofMapStream << "              qmin : " << getKey_constRef<delphi_real>("qmin")           << endl;
   ofMapStream << "             qplus : " << getKey_constRef<delphi_real>("qplus")          << endl;
   ofMapStream << "            cqplus : " << getKey_constRef< SGrid<delphi_real> >("cqplus").nX << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cqplus").nY << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cqplus").nZ << endl;
   ofMapStream << "             cqmin : " << getKey_constRef< SGrid<delphi_real> >("cqmin").nX  << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cqmin").nY  << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cqmin").nZ  << endl;
   ofMapStream << "              cmin : " << getKey_constRef< SGrid<delphi_real> >("cmin").nX   << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cmin").nY   << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cmin").nZ   << endl;
   ofMapStream << "              cmax : " << getKey_constRef< SGrid<delphi_real> >("cmax").nX   << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cmax").nY   << " "
                                          << getKey_constRef< SGrid<delphi_real> >("cmax").nZ   << endl;
   ofMapStream << "             ibnum : " << getKey_constRef<delphi_integer>("ibnum")           << endl;

   {
      const delphi_integer iGrid = getKey_constRef<delphi_integer>("igrid");

      ofMapStream << "            iepsmp : " << endl;
      if (0 != getKey_constRef< vector< SGrid<delphi_integer> > >("iepsmp").size())
      {
         vector< SGrid<delphi_integer> >::iterator it = getKey_Ref< vector< SGrid<delphi_integer> > >("iepsmp").begin();

         for (unsigned int k = 0; k < (size_t)iGrid; k++)
         {   for (unsigned int j = 0; j < (size_t)iGrid; j++)
             {  for (unsigned int i = 0; i < (size_t)iGrid; i++)
                {
                   ofMapStream << setw(3) << right << i << " " << setw(3) << right << j << " " << setw(3) << right << k << " : "
                               << setw(6) << right << it->nX << setw(6) << right << it->nY << setw(6) << right << it->nZ << endl;
                   it++;
                }
             }
          }
       }
   }

   /*------------------------------ test of getKey_constPtr overloading --------------------------//
   {
      const delphi_integer iGrid = getKey_constRef<delphi_integer>("igrid");
      const SGrid<delphi_integer> *** vctgiEpsMap = getKey_constPtr< SGrid<delphi_integer> >("iepsmp",iGrid,iGrid,iGrid);

      ofMapStream << "      iepsmp2 : " << endl;
      if (0 != iGrid)
      {
         for (unsigned int i = 0; i < (size_t)iGrid; i++)
         {   for (unsigned int j = 0; j < (size_t)iGrid; j++)
             {  for (unsigned int k = 0; k < (size_t)iGrid; k++)
                {
                   ofMapStream << setw(3) << right << i << " " << setw(3) << right << j << " " << setw(3) << right << k << " : "
                               << setw(6) << right << vctgiEpsMap[i][j][k].nX
                               << setw(6) << right << vctgiEpsMap[i][j][k].nY
                               << setw(6) << right << vctgiEpsMap[i][j][k].nZ << endl;
                }
             }
          }
       }
   }
   //------------------------------ test of getKey_constPtr overloading --------------------------*/

   {
      const delphi_integer iGrid = getKey_constRef<delphi_integer>("igrid");

      ofMapStream << "           idebmap : " << endl;

      vector<bool>::iterator it = getKey_Ref< vector<bool> >("idebmap").begin();
      if (0 != getKey_constRef< vector<bool> >("idebmap").size())
      {
         for (unsigned int k = 0; k < (size_t)iGrid; k++)
         {
            for (unsigned int j = 0; j < (size_t)iGrid; j++)
            {
               for (unsigned int i = 0; i < (size_t)iGrid; i++)
               {
                  ofMapStream << setw(3) << right << i << "," << setw(3) << right << j << "," << setw(3) << right << k << " : "
                              << setw(8) << right << (*it?'T':'F') << endl;
                  it++;
               }
            }
         }
      }
   }

   {
      const vector< SGrid<delphi_integer> >* vctgiBndyGrid = getKey_constPtr< vector< SGrid<delphi_integer> > >("ibgrd");

      ofMapStream << "             ibgrd : " << endl;
      if (0 != vctgiBndyGrid->size())
      {
         for (unsigned int i=0; i < vctgiBndyGrid->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : "
                       << setw(6) << right << vctgiBndyGrid->at(i).nX << " "
                       << setw(6) << right << vctgiBndyGrid->at(i).nY << " "
                       << setw(6) << right << vctgiBndyGrid->at(i).nZ << endl;
         }
      }
   }

   ofMapStream << "             nqgrd : " << getKey_constRef<delphi_integer>("nqgrd") << endl;

   {
      const vector< SGridValue<delphi_real> >* vctgvfCrg2Grid = getKey_constPtr< vector< SGridValue<delphi_real> > >("chrgv2");

      ofMapStream << "            chrgv2 : " << endl;
      if (0 != vctgvfCrg2Grid->size())
      {
         for (unsigned int i=0; i < vctgvfCrg2Grid->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : "
                   << setw(8) << right << vctgvfCrg2Grid->at(i).nGrid.nX << " "
                   << setw(8) << right << vctgvfCrg2Grid->at(i).nGrid.nY << " "
                   << setw(8) << right << vctgvfCrg2Grid->at(i).nGrid.nZ << " "
                   << setw(8) << right << vctgvfCrg2Grid->at(i).nValue   << endl;
         }
      }
   }

   {
      const vector<delphi_integer>* vctiCrg2GridMap = getKey_constPtr< vector<delphi_integer> >("nqgrdtonqass");

      ofMapStream << "      nqgrdtonqass : " << endl;
      if (0 != vctiCrg2GridMap->size())
      {
         for (unsigned int i=0; i < vctiCrg2GridMap->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : " << setw(6) << right << vctiCrg2GridMap->at(i) << endl;
         }
      }
   }

   {
      const vector< SGridValue<delphi_real> >* vctgvfAtomCrg = getKey_constPtr< vector< SGridValue<delphi_real> > >("atmcrg");

      ofMapStream << "            atmcrg : " << endl;
      if (0 != vctgvfAtomCrg->size())
      {
         for (unsigned int i=0; i < vctgvfAtomCrg->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : "
                       << setw(8) << right << vctgvfAtomCrg->at(i).nGrid.nX << " "
                       << setw(8) << right << vctgvfAtomCrg->at(i).nGrid.nY << " "
                       << setw(8) << right << vctgvfAtomCrg->at(i).nGrid.nZ << " "
                       << setw(8) << right << vctgvfAtomCrg->at(i).nValue << endl;
         }
      }
   }

   {
      const vector< SGrid<delphi_real> >* vctgfCrgPoseA = getKey_constPtr< vector< SGrid<delphi_real> > >("chgpos");

      ofMapStream << "            chgpos : " << endl;
      if (0 != vctgfCrgPoseA->size())
      {
         for (unsigned int i = 0; i < vctgfCrgPoseA->size(); i++)
         {
            ofMapStream << setw(6) << right << i << " : "
                        << setw(8) << right << vctgfCrgPoseA->at(i).nX << "  "
                        << setw(8) << right << vctgfCrgPoseA->at(i).nY << "  "
                        << setw(8) << right << vctgfCrgPoseA->at(i).nZ << endl;
         }
      }
   }

   {
      const vector< SGrid<delphi_real> >* vctgfSurfCrgA = getKey_constPtr< vector< SGrid<delphi_real> > >("scspos");

      ofMapStream << "            scspos : " << endl;
      if (0 != vctgfSurfCrgA->size())
      {
         for (unsigned int i = 0; i < vctgfSurfCrgA->size(); i++)
         {
            ofMapStream << setw(6) << right << i <<" : "
                      << setw(8) << right << vctgfSurfCrgA->at(i).nX << "  "
                      << setw(8) << right << vctgfSurfCrgA->at(i).nY << "  "
                      << setw(8) << right << vctgfSurfCrgA->at(i).nZ << endl;
         }
      }
   }

   {
      const vector<delphi_integer>* vctiCrgAt = getKey_constPtr< vector<delphi_integer> >("crgatn");

      ofMapStream << "            crgatn : " << endl;
      if (0 != vctiCrgAt->size())
      {
         for (unsigned int i=0; i < vctiCrgAt->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : " << setw(6) << right << vctiCrgAt->at(i) << endl;
         }
      }
   }

   {
      const vector<delphi_integer>* vctiAtSurf = getKey_constPtr< vector<delphi_integer> >("atsurf");

      ofMapStream << "            atsurf : " << endl;
      if (0 != vctiAtSurf->size())
      {
         for (unsigned int i=0; i < vctiAtSurf->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : " << setw(6) << right << vctiAtSurf->at(i) << endl;
         }
      }
   }

   {
      const vector<delphi_integer>* vctiAtNdx = getKey_constPtr< vector<delphi_integer> >("atndx");

      ofMapStream << "             atndx : " << endl;
      if (0 != vctiAtNdx->size())
      {
         for (unsigned int i=0; i < vctiAtNdx->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : " << setw(6) << right << vctiAtNdx->at(i) << endl;
         }
      }
   }

   {
      const vector< SGrid<delphi_real> >* vctgfSurfCrgE = getKey_constPtr< vector< SGrid<delphi_real> > >("scsnor");

      ofMapStream << "            scsnor : " << endl;
      if (0 != vctgfSurfCrgE->size())
      {
         for (unsigned int i = 0; i < vctgfSurfCrgE->size(); i++)
         {
            ofMapStream << setw(6) << right << i <<" : "
                        << setw(8) << right << vctgfSurfCrgE->at(i).nX << "  "
                        << setw(8) << right << vctgfSurfCrgE->at(i).nY << "  "
                        << setw(8) << right << vctgfSurfCrgE->at(i).nZ << endl;
         }
      }
   }

   {
      const vector<delphi_real>* vctfAtomEps = getKey_constPtr< vector<delphi_real> >("atmeps");

      ofMapStream << "            atmeps : " << endl;
      if (0 != vctfAtomEps->size())
      {
         for (unsigned int i=0; i < vctfAtomEps->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : " << setw(8) << right << vctfAtomEps->at(i) << endl;
         }
      }
   }

   ofMapStream << "---------------------- Solver ----------------------"        << endl;
   ofMapStream << "          icount2b : " << getKey_constRef<delphi_integer>("icount2b")    << endl;
   ofMapStream << "          icount1b : " << getKey_constRef<delphi_integer>("icount1b")    << endl;

   {
      const vector<delphi_real>* vctfGridCrg = getKey_constPtr< vector<delphi_real> >("gchrg");

      ofMapStream << "             gchrg : " << endl;
      if (0 != vctfGridCrg->size())
      {
         for (unsigned int i=0; i < vctfGridCrg->size(); i++)
         {
           ofMapStream << setw(6) << right << i << " : " << setw(8) << right << vctfGridCrg->at(i) << endl;
         }
      }
   }

   {
      const vector< SGrid<delphi_integer> >* vctgiGridCrgPose = getKey_constPtr< vector< SGrid<delphi_integer> > >("gchrgp");

      ofMapStream << "            gchrgp : " << endl;
      if (0 != vctgiGridCrgPose->size())
      {
         for (unsigned int i = 0; i < vctgiGridCrgPose->size(); i++)
         {
            ofMapStream << setw(6) << right << i << " : "
                        << setw(8) << right << vctgiGridCrgPose->at(i).nX << "  "
                        << setw(8) << right << vctgiGridCrgPose->at(i).nY << "  "
                        << setw(8) << right << vctgiGridCrgPose->at(i).nZ << endl;
         }
      }
   }

   ofMapStream << "               ibc : " << getKey_constRef<delphi_integer>("ibc")    << endl;

   {
      const vector<SDoubleGridValue>* vctdgvCrgBndyGrid = getKey_constPtr< vector<SDoubleGridValue> >("cgbp");

      ofMapStream << "              cgbp : " << endl;
      if (0 != vctdgvCrgBndyGrid->size())
      {
         for (unsigned int i = 0; i < vctdgvCrgBndyGrid->size(); i++)
         {
            ofMapStream << setw(6) << right << i << " : "
                        << setw(8) << right << vctdgvCrgBndyGrid->at(i).fgCoord.nX << "  "
                        << setw(8) << right << vctdgvCrgBndyGrid->at(i).fgCoord.nY << "  "
                        << setw(8) << right << vctdgvCrgBndyGrid->at(i).fgCoord.nZ << "  "
                        << setw(8) << right << vctdgvCrgBndyGrid->at(i).fVal1      << "  "
                        << setw(8) << right << vctdgvCrgBndyGrid->at(i).fVal2      << endl;
         }
      }
   }

   {
      const delphi_integer iGrid = getKey_constRef<delphi_integer>("igrid");
      const vector<delphi_real>* vctfPhiMap = getKey_constPtr< vector<delphi_real> >("phimap");

      ofMapStream << "            phimap : " << endl;
      if (0 != vctfPhiMap->size())
      {
         for (unsigned int k = 0; k < (size_t)iGrid; k++)

         {
            for (unsigned int j = 0; j < (size_t)iGrid; j++)
            {
               for (unsigned int i = 0; i < (size_t)iGrid; i++)
               {
                  ofMapStream << setw(3) << right << i << "," << setw(3) << right << j << "," << setw(3) << right << k << " : "
                              << setw(11) << right << vctfPhiMap->at(k*iGrid*iGrid+j*iGrid+i) << endl;
               }
            }
         }
      }
   }

   ofMapStream << "---------------------- Energy ----------------------"        << endl;
   {
      const vector<delphi_real>* vctfSurfCrgE = getKey_constPtr< vector<delphi_real> >("schrg");

      ofMapStream << "             schrg : " << endl;

      if(0 != vctfSurfCrgE->size())
      {
         for (unsigned int i = 0; i < vctfSurfCrgE->size(); i++)
         {
            ofMapStream << setw(6) << right << i << " : " << setw(8) << right << vctfSurfCrgE->at(i) << endl;
         }
      }
   }

   ofMapStream << "                  total grid energy (ergg) : " << setw(21) << right
               << getKey_constRef<delphi_real>("ergg") << endl;
   ofMapStream << "                   coulombic energy (ergc) : " << setw(21) << right
               << getKey_constRef<delphi_real>("ergc") << endl;
   ofMapStream << "    corrected reaction field energy (ergs) : " << setw(21) << right
               << getKey_constRef<delphi_real>("ergs") << endl;
   ofMapStream << "        total reaction field energy (ergr) : " << setw(21) << right
               << getKey_constRef<delphi_real>("ergr") << endl;
   ofMapStream << " total ionic direct contribution (ergions) : " << setw(21) << right
               << getKey_constRef<delphi_real>("ergions") << endl;

//   ofMapStream << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";
//   ofMapStream << "+                                                        + \n";
//   ofMapStream << "+                    DELPHI DATA CONTAINER               + \n";
//   ofMapStream << "+                                                        + \n";
//   ofMapStream << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";

   ofMapStream.close();

}
