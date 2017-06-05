/*
 * delphi_datamarshal_getStatement.cpp
 *
 *  Created on: May 04, 2014
 *      Author: chuan
 *
 * Reading statements is recode using istringstream as well.
 *
 * The order included headers in this file is crucial in order to avoid ambiguous reference to "real" when
 * compiling the code in Mac system.
 */

#include "delphi_datamarshal.h"
#include <boost/lexical_cast.hpp>

//-----------------------------------------------------------------------//
bool CDelphiDataMarshal::getStatement(string strLineNoSpace)
{
   using boost::lexical_cast;
   using boost::bad_lexical_cast;

   /*
    * It is known that newlines in DOS and Windows end with the combination of two characters, namely '\r\n',
    * while they end with a single '\n' to indicate a new line in Unix and a single '\r' in Mac.
    * Need to check if "\r" exists (for Windows machine). If so, remove \r in the string
    *
    * (chuan 2014May17) this operation is moved to IDataMarshal::read(strFileName)
    */
    //if('\r' == strLineNoSpace[strLineNoSpace.size()-1])
    //{
    //    strLineNoSpace.erase( strLineNoSpace.size()-1);
    //}

   string strLineUpperCase, strStatement, strArgument;

   locale loc;

   /*
    * transform the statement to upper case
    */
   for (size_t ii = 0; ii < strLineNoSpace.length(); ++ii)
		strLineUpperCase += toupper(strLineNoSpace[ii],loc);

   /*
    * find the equal sign which is used to separate the statement and argument
    */
   size_t found = strLineUpperCase.find_first_of("=");

   /*
    * get the statement and the argument of the statement
    */
   strStatement = strLineUpperCase.substr(0,found);
   strArgument  = strLineUpperCase.substr(found+1,strLineUpperCase.size()-found-1);
   //cout << "?????#### strLineNoSpace: " << strLineNoSpace << endl; // LinLi
   //cout << "?????#### strStatement: " << strStatement << endl; // LinLi

   /*
    * initialize to be 0 as an error indicator
    */
   size_t typearg = 0;

   /*
    * strStatment is a 2-letter short statement. One in strLineUpperCase
    */
   if (2 == found)
   {
      for (int ii = 1; ii < iStatementNum; ii++)
      {
         if (rgstrStatement_2lAbbre[ii] == strStatement)
         { typearg = ii; break; }
      }
   }
   /*
    * strStatment is a full statement shown in rgstrStatement_ShortForm. The string length is btw 2 and 7.
    */
   else if ( (2 < found) && (7 > found) )
   {
      for (int ii = 1; ii < iStatementNum; ii++)
      {
         if (rgstrStatement_ShortForm[ii] == strStatement)
         { typearg = ii; break; }
      }
   }
   /*
    * use statements in long form which are described in section 4.4. Index of statements and their
    * shorthand, delphi manual version 5.1, pp 19 - 21.
    */
   else
   {
      char c1stLetter[1];
      strStatement.copy(c1stLetter,1,0);


      switch (c1stLetter[0])
      {
         case 'A':
            if (       "AUTOMATICCONVERGENCE" == strStatement) typearg = 21;
            if (            "AUTOCONVERGENCE" == strStatement) typearg = 21;
            if (                    "AUTOCON" == strStatement) typearg = 21;
            if (                "ATOMPOTDIST" == strStatement) typearg = 43;
            break;
         case 'B':
            if (                    "BOXFILL" == strStatement) typearg =  3;
            if (          "BOUNDARYCONDITION" == strStatement) typearg =  9;
            if (         "BOUNDARYCONDITIONS" == strStatement) typearg =  9;
            break;
         case 'C':
            if (        "CONVERGENCEINTERVAL" == strStatement) typearg = 16;
            if (        "CONVERGENCEFRACTION" == strStatement) typearg = 17;
            if (                     "CUTOFF" == strStatement) typearg = 45;
            break;
         case 'E':
            if (         "EXTERIORDIELECTRIC" == strStatement) typearg =  5;
            if (         "EXTERNALDIELECTRIC" == strStatement) typearg =  5;
            if (      "EXITUNIFORMDIELECTRIC" == strStatement) typearg = 22;
            break;
         case 'F':
            if (                "FANCYCHARGE" == strStatement) typearg = 13;
            break;
         case 'G':
            if (                   "GRIDSIZE" == strStatement) typearg =  1;
            if (            "GRIDCONVERGENCE" == strStatement) typearg = 23;
            if (                   "GAUSSIAN" == strStatement) typearg = 49; //Lin Li: gaussian
            break;
         case 'I':
            if (                  "IONRADIUS" == strStatement) typearg =  7;
            if (         "INTERIORDIELECTRIC" == strStatement) typearg =  4;
            if (                 "ITERATIONS" == strStatement) typearg = 10;
            if (                  "ITERATION" == strStatement) typearg = 10;
            if (                     "INHOMO" == strStatement) typearg = 47;
            break;
         case 'L':
            if (            "LINEARITERATION" == strStatement) typearg = 10;
            if (           "LINEARITERATIONS" == strStatement) typearg = 10;
            if (          "LOGFILEPOTENTIALS" == strStatement) typearg = 14;
            if (         "LOGFILECONVERGENCE" == strStatement) typearg = 15;
            break;
         case 'M':
            //if (               "MEMBRANEDATA" == strStatement) typearg = 12; // OBSOLELE. REMOVED FROM THE LIST.
            if (             "MAXCONVERGENCE" == strStatement) typearg = 38;
            break;
         case 'N':
            if (         "NONLINEARITERATION" == strStatement) typearg = 11;
            if (        "NONLINEARITERATIONS" == strStatement) typearg = 11;
            //if (            "NORMCONVERGENCE" == strStatement) typearg = 39; // UNUSED. REMOVED FROM THE LIST.
            break;
         case 'P':
            if (          "PERIODICBOUNDARYX" == strStatement) typearg = 18;
            if (          "PERIODICBOUNDARYY" == strStatement) typearg = 19;
            if (          "PERIODICBOUNDARYZ" == strStatement) typearg = 20;
            if (                "PROBERADIUS" == strStatement) typearg =  6;
            if (             "PERCENTBOXFILL" == strStatement) typearg =  3;
            if (                "PERCENTFILL" == strStatement) typearg =  3;
            break;
         case 'R':
            if (           "RELAXATIONFACTOR" == strStatement) typearg = 24;
            if (                  "RADPOLEXT" == strStatement) typearg = 29;
            if (                     "RELPAR" == strStatement) typearg = 30;
            if (             "RMSCONVERGENCE" == strStatement) typearg = 37;
            if (                     "RADIPZ" == strStatement) typearg = 50;
            break;
         case 'S':
            if (                   "SALTCONC" == strStatement) typearg =  8;
            if (          "SALTCONCENTRATION" == strStatement) typearg =  8;
            if ("SPHERICALCHARGEDISTRIBUTION" == strStatement) typearg = 13;
            if (                      "SIGMA" == strStatement) typearg = 46;
            if (                     "SRFCUT" == strStatement) typearg = 48;
            break;
         case 'T':
            if (                "TEMPERATURE" == strStatement) typearg = 44;
            break;
      } // end of switch (c1stLetter[1])

   }

   /*
    * return if invalid statement
    */
   if (0 == typearg) return false;

   /*
    * The try-catch block tries to catch any exception thrown by lexical_cast during converting the string to numbers.
    * A better solution when comparing to atoi and atof which return zeroes silencely.
    */
   try
   {
      int iboolint;

      switch (typearg)
      {
         case 1:
            iGrid = boost::lexical_cast<delphi_integer>(strArgument);
            if (5 > iGrid || 2000 < iGrid) COutOfRange_GSIZE warning(iGrid);
            break;
         case 2:
            fScale = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 >= fScale || 40.0 <= fScale) COutofRange_SCALE warning(fScale);
            break;
         case 3:
            fPercentageFill = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (20.0 > fPercentageFill || 100.0 < fPercentageFill) COutOfRange_PERFIL warning(fPercentageFill);
            break;
         case 4:
            fInDielec = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (1000.0 < fInDielec) COutOfRange_INDI warning(fInDielec);
            break;
         case 5:
            fExDielec = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 >= fExDielec || 1000.0 < fExDielec) COutOfRange_EXDI warning(fExDielec);
            break;
         case 6:
            vctfProbeRadius[0] = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > vctfProbeRadius[0] || 10.0 <= vctfProbeRadius[0]) COutOfRange_PRBRAD warning(vctfProbeRadius[0]);
            break;
         case 7:
            fIonRadius = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > fIonRadius || 10.0 < fIonRadius) COutOfRange_IONRAD warning(fIonRadius);
            break;
         case 8:
            vctfSalt[0] = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > vctfSalt[0] || 10.0 < vctfSalt[0]) COutOfRange_SALT warning(vctfSalt[0]);
            break;
         case 9:
            iBndyType = boost::lexical_cast<int>(strArgument);
            break;
         case 10:
            iLinIterateNum = boost::lexical_cast<int>(strArgument);
            if (3 >= iLinIterateNum || 10000 < iLinIterateNum) COutOfRange_LINIT warning(iLinIterateNum);
            bAutoConverge  = false;
            break;
         case 11:
            iNonIterateNum = boost::lexical_cast<int>(strArgument);
            if (30 > iNonIterateNum || 10000 < iNonIterateNum) COutOfRange_NONIT warning(iNonIterateNum);
            break;
         case 12:
         {
            /*
             * MEMBRANEDATA has been decided to be removed after 1st of tests
             */
            //if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
            //   imem = iboolint;

            CObsoleteStatement warning(strLineNoSpace);
            break;
         }
         case 13:
            if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
               bCrgInterplateType = iboolint;
            if (bCrgInterplateType) CSphericalCrgIntelp warning(bCrgInterplateType);
            break;
         case 14:
            if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
               bLogPotential = iboolint;
            break;
         case 15:
         {
            /*
             * LOGGRP has been decided to be removed after 1st of tests
             */
            //if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
            //   bLogGraph = iboolint;

            CObsoleteStatement warning(strLineNoSpace);
            break;
         }
         case 16:
            iIterateInterval = boost::lexical_cast<int>(strArgument);
            if (1 > iIterateInterval || 100 < iIterateInterval) COutOfRange_CONINT warning(iIterateInterval);
            break;
         case 17:
            iConvergeFract = boost::lexical_cast<int>(strArgument);
            if (1 > iConvergeFract || 100 < iConvergeFract) COutOfRange_CONFRA warning(iConvergeFract);
            break;
         case 18:
            if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
               vctbPeriodicBndy[0] = (bool)iboolint;
            break;
         case 19:
            if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
               vctbPeriodicBndy[1] = (bool)iboolint;
            break;
         case 20:
            if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
               vctbPeriodicBndy[2] = (bool)iboolint;
            break;
         case 21:
            if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
               bAutoConverge = (bool)iboolint;
            break;
         case 22:
            if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
               bExitUniformDielect = (bool)iboolint;
            break;
         case 23:
         {
            /*
             * GRDCON has been decided to be removed after 1st of tests
             */
            //fGridConverge = boost::lexical_cast<delphi::delphi_real>(strArgument);
            //if (0.0 > fGridConverge || 100.0 < fGridConverge) COutOfRange_GRDCON warning(fGridConverge);
            //
            //bAutoConverge = true;

            CObsoleteStatement warning(strLineNoSpace);
            break;
         }
         case 24:
            fSpectralRadius = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 >= fSpectralRadius || 1.0 <= fSpectralRadius) COutOfRange_RELFAC warning(fSpectralRadius);
            bSpectralRadius = true;
            break;
         case 25:
            if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
               bFixedRelaxParam = (bool)iboolint;
            break;
         case 26:
            if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
               bSolvePB = (bool)iboolint;
            break;
         case 27:
         {
            /*
             * CLCSRF has been decided to be removed after 1st of tests
             */
            //if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
            //   bOutGraspSurf = (bool)iboolint;

            CObsoleteStatement warning(strLineNoSpace);
            break;
         }
         case 28:
         {
            /*
             * PHICON has been decided to be removed after 1st of tests
             */
            //if ( -1 < (iboolint = yesno(strArgument, strStatement)) )
            //   bOutCrgDensity = (bool)iboolint;

            CObsoleteStatement warning(strLineNoSpace);
            break;
         }
         case 29:
         {
            /*
             * RADPOL has been decided to be removed after 1st of tests
             */
            //fRadPolExt = boost::lexical_cast<delphi::delphi_real>(strArgument);

            CObsoleteStatement warning(strLineNoSpace);
            break;
         }
         case 30:
            fRelaxParam = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > fRelaxParam || 2.0 < fRelaxParam) COutOfRange_RELPAR warning(fRelaxParam);
            bManualRelaxParam = true;
            break;
         case 31:
            vctfSalt[1] = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > vctfSalt[1] || 10.0 < vctfSalt[1]) COutOfRange_SALT2 warning(vctfSalt[1]);
            break;
         case 32:
            vctfProbeRadius[1] = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > vctfProbeRadius[1] || 10.0 <= vctfProbeRadius[1]) COutOfRange_PRBRAD2 warning(vctfProbeRadius[1]);
            break;
         case 33:
            vctiValence1[0] = boost::lexical_cast<int>(strArgument);
            if (1 > vctiValence1[0] || 10 < vctiValence1[0]) COutOfRange_VALPLUS1 warning(vctiValence1[0]);
            break;
         case 34:
            vctiValence1[1] = boost::lexical_cast<int>(strArgument);
            if (1 > vctiValence1[1] || 10 < vctiValence1[1]) COutOfRange_VALMINUS1 warning(vctiValence1[1]);
            break;
         case 35:
            vctiValence2[0] = boost::lexical_cast<int>(strArgument);
            if (1 > vctiValence2[0] || 10 < vctiValence2[0]) COutOfRange_VALPLUS2 warning(vctiValence2[0]);
            break;
         case 36:
            vctiValence2[1] = boost::lexical_cast<int>(strArgument);
            if (1 > vctiValence2[1] || 10 < vctiValence2[1]) COutOfRange_VALMINUS2 warning(vctiValence2[1]);
            break;
         case 37:
            fRmsc = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 >= fRmsc || 1.0 <= fRmsc) COutOfRange_RMSC warning(fRmsc);
            break;
         case 38:
            fMaxc = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 >= fMaxc || 1.0 <= fMaxc) COutOfRange_MAXC warning(fMaxc);
            break;
         case 39:
         {
            /*
             * NORMC has been decided to be removed after 1st of tests
             */
            //fNormc = boost::lexical_cast<delphi::delphi_real>(strArgument) ;

            CObsoleteStatement warning(strLineNoSpace);
            break;
         }
         case 40:
            gfPotentialDrop.nX = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > gfPotentialDrop.nX || 1000.0 < gfPotentialDrop.nX)
               COutOfRange_VDROP warning(gfPotentialDrop.nX);
            else
               vctbPeriodicBndy[3] = true;
            break;
         case 41:
            gfPotentialDrop.nY = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > gfPotentialDrop.nY || 1000.0 < gfPotentialDrop.nY)
               COutOfRange_VDROP warning(gfPotentialDrop.nY);
            else
               vctbPeriodicBndy[4] = true;
            break;
         case 42:
            gfPotentialDrop.nZ = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > gfPotentialDrop.nZ || 1000.0 < gfPotentialDrop.nZ)
               COutOfRange_VDROP warning(gfPotentialDrop.nZ);
            else
               vctbPeriodicBndy[5] = true;
            break;
         case 43:
            fPotentialUpperBond = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 >= fPotentialUpperBond || 10.0 <= fPotentialUpperBond) COutOfRange_ATPODS warning(fPotentialUpperBond);
            break;
         case 44:
            fTemper = boost::lexical_cast<delphi::delphi_real>(strArgument);
            if (0.0 > fTemper || 1000.0 < fTemper) COutOfRange_TEMPER warning(fTemper);
            fTemper -= dAbsoluteZero;
            break;
         case 45:
            fCutoff = boost::lexical_cast<float>(strArgument);
            break;
        case 46:
            fSigma = boost::lexical_cast<float>(strArgument);
            break;
         case 47:
            iInhomo = boost::lexical_cast<int>(strArgument);
            break;
         case 48:
            fSrfcut = boost::lexical_cast<float>(strArgument);
            break;
         case 49:
            iGaussian = boost::lexical_cast<int>(strArgument);
            break;
        case 50:
            fRadipz = boost::lexical_cast<float>(strArgument);
            break;

      } // end of switch (typearg)

      return true;
   }
   catch(bad_lexical_cast&)
   {
      return false;
   }
}

//-----------------------------------------------------------------------//
int CDelphiDataMarshal::yesno(string strArgument, string strStatement)
{

   if (0 == strArgument.compare("TRUE"))  return 1;

   if (0 == strArgument.compare("YES"))   return 1;

   if (0 == strArgument.compare("ON"))    return 1;

   if (0 == strArgument.compare("T"))     return 1;

   if (0 == strArgument.compare("FALSE")) return 0;

   if (0 == strArgument.compare("OFF"))   return 0;

   if (0 == strArgument.compare("NO"))    return 0;

   if (0 == strArgument.compare("F"))     return 0;

   CBadStatementAssignment warning(strArgument, strStatement);

   return -1;
}

