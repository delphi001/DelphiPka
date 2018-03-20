/*
 * delphi_datamarshal_showParameters.cpp
 *
 *  Created on: Mar 07, 2014
 *      Author: chuan
 */

#include "delphi_datamarshal.h"
#include <boost/lexical_cast.hpp>

//-----------------------------------------------------------------------//
void CDelphiDataMarshal::showParameters() const
{

   cout << fixed; // set the floatfield format flag for the str stream to fixed
   cout << "   " << endl;

   size_t MAXWIDTH=45;
   std::string time_string = " Time>";


#ifdef VERBOSE
   if (1 < iMediaNum)
   {
      cout << "Attention, many dielectrics! Not all the surface charges are facing the solution!!" << endl;
   }
#endif

   cout << time_string << left << setw(MAXWIDTH) << " Time to read in and/or assign rad/chrg" << ":"; pTimer->showElapse();
   cout << endl;

   cout << left << setw(MAXWIDTH) << " Direct mapping of epsilon" << " : " << (iDirectEpsMap?"1/Y":"0/N") << endl;

   if (1000000 < fDebyeLength)
   {
      delphi_real fDebyeNum = (100.0/fPercentageFill-1.0)*fMaxDimension/fDebyeLength;
#ifdef VERBOSE
      cout << left << setw(MAXWIDTH) << " Debye Lengths contained in the finite diff. box" << " : " << fDebyeNum << endl;
#endif
   }


//#ifdef VERBOSE
   cout << left << setw(MAXWIDTH) << " Grid size" << " : " << iGrid << endl;
   cout << left << setw(MAXWIDTH) << " Percent of box occupied" << " : " << fPercentageFill << endl;
   cout << left << setw(MAXWIDTH) << " Scale,in grids (1/A)" << " : " << setw(10) << fScale << endl;
   cout << left << setw(MAXWIDTH) << " xmin,xmax (A)"   << " : " << gfMinCoordinate.nX << " " << right << setw(10) << gfMaxCoordinate.nX << endl;
   cout << left << setw(MAXWIDTH) << " ymin,ymax (A)"   << " : " << gfMinCoordinate.nY << " " << right << setw(10) << gfMaxCoordinate.nY << endl;
   cout << left << setw(MAXWIDTH) << " zmin,zmax (A)"   << " : " << gfMinCoordinate.nZ << " " << right << setw(10) << gfMaxCoordinate.nZ << endl;
   cout << left << setw(MAXWIDTH) << " x,y,z range (A)" << " : " << gfCoordinateRange.nX << " "
                                << right << setw(10) << gfCoordinateRange.nY << " "
                                << right << setw(10) << gfCoordinateRange.nZ << endl;
   cout << left << setw(MAXWIDTH) << " System geometric center (A)" << " : " << gfGeometricCenter.nX << " "
                                             << right << setw(10) << gfGeometricCenter.nY << " "
                                             << right << setw(10) << gfGeometricCenter.nZ << endl;
   cout << left << setw(MAXWIDTH) << " Grid box is centered (A)" << " : " << gfBoxCenter.nX << " "
                                             << right << setw(10) << gfBoxCenter.nY << " "
                                             << right << setw(10) << gfBoxCenter.nZ << endl;
   cout << left << setw(MAXWIDTH) << " Object centre offset (gu)" << " : " << gfOffCenter.nX << " "
                                             << right << setw(10) << gfOffCenter.nY << " "
                                             << right << setw(10) << gfOffCenter.nZ  << endl;
//#endif


   if (bSolvePB)
   {
//#ifdef VERBOSE
      cout << left << setw(MAXWIDTH) << " Outer dielectric" << " : " << fExDielec << endl;
//#endif
      for (delphi_integer i = 1; i <= iMediaNum; i++)
      {
         delphi_real fEspInMedium = vctfMediaEps[i]*fEPKT;

         if (1.0 > fEspInMedium) throw CInvalidEpsInMedium(i);

         std::string dielec_info = " Dielectric in Medium " + boost::lexical_cast<std::string>(i);
         cout << left << setw(MAXWIDTH) << dielec_info << " : " << fEspInMedium << endl;

      }


// #ifdef VERBOSE

      cout << left << setw(MAXWIDTH) << " First kind salt [C] (M) " << " : " << vctfSalt[0] << endl;
      cout << left << setw(MAXWIDTH) << " Valences salt 1 are"      << " : " << vctiValence1[0] <<  right << setw(10) << vctiValence1[1] << endl;

		if(vctfSalt[1]!=0.0f) {
			cout << left << setw(MAXWIDTH) << " Second kind salt [C] (M) " << " : " << vctfSalt[1] << endl;
			cout << left << setw(MAXWIDTH) << " Valences salt 2 are"       << " : " << vctiValence2[0] <<  right << setw(10) << vctiValence2[1] << endl;
		}
      cout << left << setw(MAXWIDTH) << " Ionic strength (M)"        << " : " << fIonStrength << endl;
      cout << left << setw(MAXWIDTH) << " Debye length (A)"          << " : " << fDebyeLength << endl;
      cout << left << setw(MAXWIDTH) << " Absolute temperature (K)"  << " : " << fTemper      << endl;
      cout << left << setw(MAXWIDTH) << " Ion exclusion [r] (A)"    << " : " << fIonRadius   << endl;
      cout << left << setw(MAXWIDTH) << " Probe[r] facing water (A)" << " : " << vctfProbeRadius[0] << endl;
      cout << left << setw(MAXWIDTH) << " Probe[r] internal (A)" << " : " << vctfProbeRadius[1] << endl;
// #endif


      switch (iBndyType)
      {
         case 1: cout << left << setw(MAXWIDTH) << " Boundary conditions" << " : " << "ZERO\n";         break;
         case 2: cout << left << setw(MAXWIDTH) << " Boundary conditions" << " : " << "DIPOLAR\n";      break;
         case 3: cout << left << setw(MAXWIDTH) << " Boundary conditions" << " : " << "FOCUSSING\n";    break;
         case 4: cout << left << setw(MAXWIDTH) << " Boundary conditions" << " : " << "COULOMBIC\n";    break;
         case 5: cout << left << setw(MAXWIDTH) << " Boundary conditions" << " : " << "VOLTAGE DROP\n"; break;
      }

      //cout << "Gaussian space module     :" << right << setw(10) << (iGaussian?"ON":"OFF") << endl;

#ifdef VERBOSE
      cout << left << setw(MAXWIDTH) << " x,y,z PBC. and volt drop" << " : " << vctbPeriodicBndy[0] << right << setw(6) << vctbPeriodicBndy[1]
                                                          << right << setw(6) << vctbPeriodicBndy[2] << right << setw(6) << vctbPeriodicBndy[3]
														                              << right << setw(6) << vctbPeriodicBndy[4] << right << setw(6) << vctbPeriodicBndy[5] << endl;



      if (vctbPeriodicBndy[3] || vctbPeriodicBndy[4] || vctbPeriodicBndy[5])
         cout << left << setw(MAXWIDTH) << " Voltage drops along x,y,z" << " : " << gfPotentialDrop.nX
                                                                                 << right << setw(10) << gfPotentialDrop.nY
                                                                                 << right << setw(10) << gfPotentialDrop.nZ << endl;
#endif
      if (bAutoConverge)
      {
         if (0.0 < fGridConverge)
            cout << left << setw(MAXWIDTH) << " Convergence by grid energy" << " : " << fGridConverge << " kt\n";
         else
            cout << left << setw(MAXWIDTH) << " # of linear iterations" << " : " << "Automatic\n";
      }
      else
         cout << left << setw(MAXWIDTH) << " # of linear iterations" << " : " << iLinIterateNum << endl;

      if (0.0 < fRmsc || 0.0 < fMaxc)
      {
         cout << left << setw(MAXWIDTH) << " Convergence by rms change" << " : " << scientific << fRmsc << " kT\n";
         cout << left << setw(MAXWIDTH) << " Convergence by max change" << " : " << scientific << fMaxc << " kT\n";
      }

      if (fZero > fIonStrength && 0 < iNonIterateNum)
         cout << left << setw(MAXWIDTH) << " Ionic strength=0 ==> only linear iterations \n";
      else
      {
         cout << left << setw(MAXWIDTH) << " # of non-linear iterations" << " : " << iNonIterateNum << endl;
         cout << left << setw(MAXWIDTH) << " Non-linear energy calc."    << " : " << bNonlinearEng << endl;
         cout << left << setw(MAXWIDTH) << " Manual relaxation para."    << " : " << bManualRelaxParam << endl;
      }

      cout << left << setw(MAXWIDTH) << " GAUSSIAN module" <<  " : " << (iGaussian?"ON":"OFF") << endl;
      if ( iGaussian ) cout << left << setw(MAXWIDTH) << " GAUSSIAN sigma" <<  " : " << fSigma << endl;
      if ( iGaussian ) cout << left << setw(MAXWIDTH) << " GAUSSIAN surface cut off" <<  " : " << fSrfcut << endl;

      cout << left << setw(MAXWIDTH) << " CONVOLUTE module" << " : " << (iConvolute?"ON":"OFF") << endl;
      if ( iConvolute ) cout << left << setw(MAXWIDTH) << " CONVOLUTE kernel sigma"      <<  " : " << fksigma << endl;
      if ( iConvolute ) cout << left << setw(MAXWIDTH) << " CONVOLUTE heavyside epsilon" <<  " : " << fhvsd_eps << endl;

      cout << left << setw(MAXWIDTH) << " Surface potential calculations" << " : " << (zetaOn?"ON":"OFF") << endl;
      if ( zetaOn ) cout << left << setw(MAXWIDTH) << " Surface requested at distance (A)"      <<  " : " << zetaDistance << endl;


      cout << endl;

#ifdef VERBOSE
      cout << left << setw(MAXWIDTH) << " Ionic direct energy"      << " : " << bIonsEng << endl;
      cout << left << setw(MAXWIDTH) << " Concentration map output" << " : " << bOutCrgDensity << endl;
      cout << left << setw(MAXWIDTH) << " Spherical charge distbn." << " : " << bCrgInterplateType << endl;
      cout << left << setw(MAXWIDTH) << " INSIGHT format output"    << " : " << bBiosystemOut << endl;
      cout << left << setw(MAXWIDTH) << " Ionic direct energy"      << " : " << bSiteOut << endl;
#endif
   } // ---------- end of if (bSolvePB)
#ifdef VERBOSE
   cout << left << setw(MAXWIDTH) << " Modified atom file output" << " : " << bModPdbOut << endl;
   cout << left << setw(MAXWIDTH) << " Map file label"            << " : " << rgcFileMap << endl;

   if (bPdbUnformatIn)  cout << " Set to read  unformatted pdb file\n";
   if (bUnformatPdbOut) cout << " Set to write unformatted pdb file\n";
   if (bFrcUnformatIn)  cout << " Set to read  unformatted frc.pdb file\n";
   if (bUnformatFrcOut) cout << " Set to write unformatted frc.pdb file\n";
   if (!bLogGraph)      cout << " Convergence graph turned off\n";
   if (!bLogPotential)  cout << " Potential listings turned off\n";

   if (10 != iIterateInterval || 1 != iConvergeFract)
   {
      cout << " Convergence test interval is every" << right << setw(6) << iIterateInterval << " loops \n";
      cout << " Testing" << right << setw(6) << 100/iConvergeFract << "% \n";
   }

   cout << endl;
#endif
}
