/*
 * delphi_datamarshal_getFunction.cpp
 *
 *  Created on: May 03, 2014
 *      Author: chuan
 *
 * The way how these functions are coded is confusing and unsafe. I re-coded them using istringstream
 * (see below). Hope it's better in C++ language...
 *
 * The order included headers in this file is crucial in order to avoid ambiguous reference to "real" when
 * compiling the code in Mac system
 */

#include "delphi_datamarshal.h"

#include <sstream>
//#include <iterator>
//#include <algorithm>
#include <boost/lexical_cast.hpp>

//-----------------------------------------------------------------------//
bool CDelphiDataMarshal::getFunction(string strLineNoSpace)
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;

    locale loc;

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

    /*
     * Find the position of the bracket. Report error if there is no '(' or ')'
     */
    if (string::npos == strLineNoSpace.find_first_of('(') || string::npos == strLineNoSpace.find_first_of(')') )
        return false;

    /*
     * get the function name and transform the function name to upper case
     */
    size_t pos = strLineNoSpace.find_first_of('(');

    string strFunc = strLineNoSpace.substr(0,pos);

    for (size_t ii = 0; ii < strFunc.length(); ++ii)
        strFunc[ii] = toupper(strFunc[ii], loc);

    /*
     * determine function type based on its name
     */
    size_t typearg = 0; // initialize to be 0 as an error indicator

    for (int ii = 1; ii < iFunctionNum_FullName; ii++)
    {
        if (rgstrFunction_FullForm[ii] == strFunc)
        {
            typearg = ii;
            break;
        }
    }

    for (int ii = 1; ii < iFunctionNum_ShortName; ii++)
    {
        if (rgstrFunction_ShortForm[ii] == strFunc)
        {
            typearg = ii;
            break;
        }
    }

    if (0 == typearg) return false; // invalid function name

    /*
     * get the arguments list of the function w/0 '(' and ')', convert the characters to upper case
     */
    if (2 == strLineNoSpace.substr(pos, strLineNoSpace.size()-pos).length()) // function has no argument
    {
        CEnmptyParameter_FUNCTION warning(strFunc);
        return false;
    }

    string strArgs_fromInput = strLineNoSpace.substr(pos+1, strLineNoSpace.find_first_of(')')-pos-1); // take off '(' and ')'

    string strArgs_UpperCase; // inputs in upper case

    for (size_t ii = 0; ii < strArgs_fromInput.length(); ++ii)
        strArgs_UpperCase += toupper(strArgs_fromInput[ii], loc);

    /*
     * find the unit or form argument
     */
    //bool bFile = false;
    //if (string::npos != strArgs_UpperCase.find("UNIT=") || string::npos != strArgs_UpperCase.find("FILE=") )
    //   bFile = true;

    /*
     * determine if there is a frm, form or format argument
     */
    //bool bFormat = false;
    //if (string::npos != strArgs_UpperCase.find("FRM=") || string::npos != strArgs_UpperCase.find("FORM=") || string::npos != strArgs_UpperCase.find("FORMAT="))
    //   bFormat = true;

    /*
     * pharse strArgs_UpperCase into takens for easy handling
     */
    std::vector<std::string> prgstrArgTokens_UpperCase = getArguments(strArgs_UpperCase);
    std::vector<std::string> prgstrArgTokens_fromInput = getArguments(strArgs_fromInput);

    /*
     * The try-catch block tries to catch any exception thrown by lexical_cast during converting the string to numbers.
     * A better solution when comparing to atoi and atof which return zeroes silencely.
     */
    try
    {
        switch (typearg)
        {
        /*
         * CENTER or CENT function
         * =======================
         * This function is used to specify the offset (expressed in grid units) with respect to the lattice center
         * at which the center of the molecule [pmid(3)] is placed. This will influence what point in the real space
         * (expressed in Angstroms) is placed at the center of the grid [oldmid(3)]. The relationship between real
         * space r(i) and grid g(i) coordinates for a grid size of igrid, with a scale of gpa grids/angstrom is as follows.
         *
         * The centre of the grid is:
         *    midg = (igrid+1)/2
         *    oldmid(i) = pmid(i) - OFFSET(i)/gpa
         *    g(i) = (r(i) - pmid(i))*gpa + midg + OFFSET(i)
         *    r(i) = (g(i) - midg)/gpa + oldmid(3)
         *
         * The scale, the system center and the shift are printed in the logfile.
         *
         * Note that a certain error inevitably results from the mapping of the molecule onto the grid. By moving the
         * molecule slightly (changing CENTER offset between 0,0,0 and 1,1,1) and repeating the calculations, it is possible
         * to see whether the results are sensitive to the particular position on the grid, and if so, to improve the accuracy
         * by averaging (this is related to rotational averaging, discussed in the J. Comp Chem paper of Gilson et al.).
         * However using a larger scale is a more effective way of improving accuracy than averaging. 
         */
        case 1:
        {
            vector<string>::iterator itt = prgstrArgTokens_fromInput.begin();

            /*
             * read center info from a file
             */
            if (string::npos != prgstrArgTokens_UpperCase[0].find("UNIT=") || string::npos != prgstrArgTokens_UpperCase[0].find("FILE="))
            {
                strCentFile = getFile_Name_or_Format(prgstrArgTokens_UpperCase[0],prgstrArgTokens_fromInput[0]);

                gfOffCenter.nX = 999.0;
                gfOffCenter.nY = 0.0;
                gfOffCenter.nZ = 0.0; // center(999,0,0)

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing input file name
                        if (string::npos != it->find("AN=1"))
                        {
                            gfOffCenter.nX = 999.0;
                            gfOffCenter.nY = 999.0;
                            gfOffCenter.nZ = 0.0; // center(999,999,0)
                        }
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * read center coordinates
             */
            else if (0 < prgstrArgTokens_UpperCase.size() && 4 > prgstrArgTokens_UpperCase.size())
            {
                gfOffCenter.nX = boost::lexical_cast<delphi::delphi_real>(prgstrArgTokens_UpperCase[0]);

                if (1 < prgstrArgTokens_UpperCase.size())
                    gfOffCenter.nY = boost::lexical_cast<delphi::delphi_real>(prgstrArgTokens_UpperCase[1]);


                if (2 < prgstrArgTokens_UpperCase.size())
                    gfOffCenter.nZ = boost::lexical_cast<delphi::delphi_real>(prgstrArgTokens_UpperCase[2]);
            }
            else
                return false;

            break;
        }

        /*
         * ACENTER or ACENT function
         * =========================
         * Acenter takes three absolute coordinates, i.e. in Å and uses those as the center
         */
        case 2:
        {
            if (0 < prgstrArgTokens_UpperCase.size() && 4 > prgstrArgTokens_UpperCase.size())
            {
                gfAcent.nX = boost::lexical_cast<delphi::delphi_real>(prgstrArgTokens_UpperCase[0]);

                if (1 < prgstrArgTokens_UpperCase.size())
                    gfAcent.nY = boost::lexical_cast<delphi::delphi_real>(prgstrArgTokens_UpperCase[1]);

                if (2 < prgstrArgTokens_UpperCase.size())
                    gfAcent.nZ = boost::lexical_cast<delphi::delphi_real>(prgstrArgTokens_UpperCase[2]);
            }
            else
                return false;

            bIsAcent = true;

            break;
        }

        /*
         * READ or IN function
         * ===================
         * This function allows files to be read as input. It comes with several specifiers, namely:
         *    SIZ    : for the radius file  CRG: for the charge file
         *    PDB    : for the pdb structure file (possible alternative formats: frm=UN and frm=MOD)
         *    MODPDB4: for the modified pdb structure file (possible alternative formats: frm=PQR and frm=MOD) which
         *             contain charges and radii values with 4 digits precession after decimal points.
         *    FRC    : for positions of site potentials.
         *    PHI    : for the phimap used in focusing
         *
         * The main use, at present will be to give the user flexibility to specify the file name or unit number of
         * any of these files. Note that the default files for all read (and write) operations are the standard DelPhi
         * files.
         */
        case 3:
        {
            vector<string>::iterator itt = prgstrArgTokens_fromInput.begin();

            /*
             * size file
             */
            if (0 == prgstrArgTokens_UpperCase[0].compare("SIZ"))
            {
                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing input file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strSizeFile = getFile_Name_or_Format(*it,*itt);
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * charge file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("CRG"))
            {
                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing input file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strCrgFile = getFile_Name_or_Format(*it,*itt);
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * standard PDB file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("PDB"))
            {
                iPdbFormatIn = STDPDB;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing input file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strPdbFile = getFile_Name_or_Format(*it,*itt);
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * modified pdb file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("MODPDB"))
            {
                iPdbFormatIn = MODPDB;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing input file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strPdbFile = getFile_Name_or_Format(*it,*itt);
                        // a token describing input file format
                        else if (string::npos != it->find("FRM=") || string::npos != it->find("FORM=") || string::npos != it->find("FORMAT="))
                        {
                            string strFormat = getFile_Name_or_Format(*it,*itt);

                            if (0 == strFormat.compare("MOD"))
                                iPdbFormatIn = MODPDB;
                            else if (0 == strFormat.compare("PQR"))
                                iPdbFormatIn = PQRPDB;
                            else
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                        }
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * modified pdb (4-digit) file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("MODPDB4"))
            {
                iPdbFormatIn = MOD4PDB;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing input file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strPdbFile = getFile_Name_or_Format(*it,*itt);
                        // a token describing input file format
                        else if (string::npos != it->find("FRM=") || string::npos != it->find("FORM=") || string::npos != it->find("FORMAT="))
                        {
                            string strFormat = getFile_Name_or_Format(*it,*itt);

                            if (0 == strFormat.compare("MOD"))
                                iPdbFormatIn = MOD4PDB;
                            else if (0 == strFormat.compare("PQR"))
                                iPdbFormatIn = PQR4PDB;
                            else
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                        }
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * a shortcut to read PQR file (removed after 1st round of tests)
             */
            //else if (0 == prgstrArgTokens_UpperCase[0].compare("PQR"))
            //{
            //   iPdbFormatIn = PQRPDB;
            //
            //   if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
            //   {
            //      for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
            //      {
            //         itt++;
            //
            //         // a token describing input file name
            //         if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
            //            strPdbFile = getFile_Name_or_Format(*it,*itt);
            //         // invalid token
            //         else
            //            CUnknownParameter_FUNCTION warning(strFunc,*itt);
            //      }
            //   }
            //}
            /*
             * a shortcut to read PQR4 file (removed after 1st round of tests)
             */
            //else if (0 == prgstrArgTokens_UpperCase[0].compare("PQR4"))
            //{
            //   iPdbFormatIn = PQR4PDB;
            //
            //   if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
            //   {
            //      for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
            //      {
            //         itt++;
            //
            //         // a token describing input file name
            //         if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
            //            strPdbFile = getFile_Name_or_Format(*it,*itt);
            //         // invalid token
            //         else
            //            CUnknownParameter_FUNCTION warning(strFunc,*itt);
            //      }
            //   }
            //}
            /*
             * FRC file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("FRC"))
            {
                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing input file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strFrciFile = getFile_Name_or_Format(*it,*itt);
                        // a token describing input file format
                        else if (string::npos != it->find("FRM=") || string::npos != it->find("FORM=") || string::npos != it->find("FORMAT="))
                        {
                            string strFormat = getFile_Name_or_Format(*it,*itt);

                            if (0 == strFormat.compare("SELF"))
                                bPDB2FRCInSite = true;
                            else
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                        }
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * PHI file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("PHI"))
            {

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing input file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                        {
                            strPhiiFile = getFile_Name_or_Format(*it,*itt);
                            //cout << "get: " <<strPhiiFile<< endl;
                            //exit(0);
                        }
                        else if (string::npos != it->find("FORMAT=") || string::npos != it->find("FOR=") ) //for focusing cube
                        {
                            if( getFile_Name_or_Format(*it,*itt) == "CUBE")
                                iPhiInType=1;
                        }


                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * invalid token prgstrArgTokens_UpperCase[0]
             */
            else
                CUnknownParameter_FUNCTION warning(strFunc,prgstrArgTokens_fromInput[0]);

            break;
        }

        /*
         * WRITE or OUT function
         * =====================
         * Equally obviously this deals with output. The specifiers are:
         *    PHI    : for phimaps (possible other formats: frm=BIOSYM, frm=GRASP, frm=CUBE; see Unit14 in Files)
         *    FRC    : for site potentials (possible other formats: frm=RC, frm=R, frm=UN;)
         *    EPS    : for epsmaps
         *    MODPDB : for modified pdb files
         *    MODPDB4: modified pdb files that contain charges and radii with 4 digits precision after decimal points
         *    UNPDB  : for unformatted pdb file
         *    UNFRC  : for unformatted frc files
         *    ENERGY : writes the file "energy.dat" containing energy data. (Example: out(energy))
         *             Note that this is different from the Energy function!
         */
        case 4:
        {
            vector<string>::iterator itt = prgstrArgTokens_fromInput.begin();

            /*
             * PHI file
             */
            if (0 == prgstrArgTokens_UpperCase[0].compare("PHI"))
            {
                bPhimapOut = true;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing output file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strPhiFile = getFile_Name_or_Format(*it,*itt);
                        // a token describing output file format
                        else if (string::npos != it->find("FRM=") || string::npos != it->find("FORM=") || string::npos != it->find("FORMAT="))
                        {
                            string strFormat = getFile_Name_or_Format(*it,*itt);

                            if (0 == strFormat.compare("BIOSYM"))
                            {
                                iPhiFormatOut = 1;
                                bBiosystemOut = true;
                            }
                            else if (0 == strFormat.compare("GRASP"))
                                iPhiFormatOut = 2;
                            else if (0 == strFormat.compare("CCP4"))
                                iPhiFormatOut = 3;
                            else if (0 == strFormat.compare("DIFF"))
                                iPhiFormatOut = 4;
                            else if (0 == strFormat.compare("CUBE"))
                                iPhiFormatOut = 5;
                            else
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                        }
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * SRF file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("SRF"))
            {
                /*
                 * out(srf,...) has been decided to be removed after 1st of tests
                 */

                /*
                bOutGraspSurf = true;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                   for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                   {
                      itt++;

                      // a token describing output file name
                      if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                         strGraspFile = getFile_Name_or_Format(*it,*itt);
                      // a token describing output file format
                      else if (string::npos != it->find("FRM=") || string::npos != it->find("FORM=") || string::npos != it->find("FORMAT="))
                      {
                         string strFormat = getFile_Name_or_Format(*it,*itt);

                         if (0 == strFormat.compare("BEM"))
                            bBemSrfOut = true;
                         else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                      }
                      // invalid token
                      else
                         CUnknownParameter_FUNCTION warning(strFunc,*itt);
                   }
                }
                */

                CObsoleteFunction warning(strLineNoSpace);

            }
            /*
             * FRC file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("FRC"))
            {
                bSiteOut = true;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing output file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strFrcFile = getFile_Name_or_Format(*it,*itt);
                        // a token describing output file format
                        else if (string::npos != it->find("FRM=") || string::npos != it->find("FORM=") || string::npos != it->find("FORMAT="))
                        {
                            string strFormat = getFile_Name_or_Format(*it,*itt);

                            if (0 == strFormat.compare("RC"))
                            {
                                /*
                                 * format = RC has been decided to be removed after 1st round of tests
                                 */
                                //iFrcFormatOut = 1; bReactFieldInFRC = true;
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                            }
                            else if (0 == strFormat.compare("R"))
                            {
                                /*
                                 * format = RC has been decided to be removed after 1st round of tests
                                 */
                                //iFrcFormatOut = 2; bReactFieldInFRC = true;
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                            }
                            else if (0 == strFormat.compare("UN"))
                            {
                                /*
                                 * format = RC has been decided to be removed after 1st round of tests
                                 */
                                //iFrcFormatOut = 3;
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                            }
                            else
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                        }
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * EPS file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("EPS"))
            {
                bEpsOut = true;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing output file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strEpsFile = getFile_Name_or_Format(*it,*itt);
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * MODPDB file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("MODPDB"))
            {
                bModPdbOut = true;

                iModPdbFormatOut = MODPDB;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing output file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strModifiedPdbFile = getFile_Name_or_Format(*it,*itt);
                        // a token describing output file format
                        else if (string::npos != it->find("FRM=") || string::npos != it->find("FORM=") || string::npos != it->find("FORMAT="))
                        {
                            string strFormat = getFile_Name_or_Format(*it,*itt);

                            if (0 == strFormat.compare("PQR"))
                                iModPdbFormatOut = PQRPDB;
                            else
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                        }
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * MODPDB4 file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("MODPDB4"))
            {
                bModPdbOut = true;

                iModPdbFormatOut = MOD4PDB;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing output file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strModifiedPdbFile = getFile_Name_or_Format(*it,*itt);
                        // a token describing output file format
                        else if (string::npos != it->find("FRM=") || string::npos != it->find("FORM=") || string::npos != it->find("FORMAT="))
                        {
                            string strFormat = getFile_Name_or_Format(*it,*itt);

                            if (0 == strFormat.compare("PQR"))
                                iModPdbFormatOut = PQR4PDB;
                            else
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                        }
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * UNPDB file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("UNPDB"))
            {
                /*
                 * out(unpdb,...) has been decided to be removed after 1st round of tests
                 */
                //bUnformatPdbOut = true;
                //
                //if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                //{
                //   for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                //   {
                //      itt++;
                //
                //      // a token describing output file name
                //      if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                //         strUnformatPdbFile = getFile_Name_or_Format(*it,*itt);
                //      // invalid token
                //      else
                //         CUnknownParameter_FUNCTION warning(strFunc,*itt);
                //   }
                //}
                CObsoleteFunction warning(strLineNoSpace);
            }
            /*
             * UNFRC file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("UNFRC"))
            {
                /*
                 * out(unfrc,...) has been decided to be removed after 1st round of tests
                 */
                //bUnformatFrcOut = true;
                //
                //if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                //{
                //   for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                //   {
                //      itt++;
                //
                //      // a token describing output file name
                //      if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                //         strUnformatFrcFile = getFile_Name_or_Format(*it,*itt);
                //      // invalid token
                //      else
                //         CUnknownParameter_FUNCTION warning(strFunc,*itt);
                //   }
                //}
                CObsoleteFunction warning(strLineNoSpace);
            }
            /*
             * ENERGY file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("ENERGY"))
            {
                bEngOut = true;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing output file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strEnergyFile = getFile_Name_or_Format(*it,*itt);
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * GCRG file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("GCRG"))
            {
                bGridCrgOut = true;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;
                        CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }

                }
            }
            /*
             * HSURF file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("HSURF"))
            {
                bHsurf2DatOut = true;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;
                        CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }

                }
            }
            /*
             * DB file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("DB"))
            {
                bDbOut = true;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;
                        CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }

                }
            }
            /*
             * SURFEN file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("SURFEN"))
            {
                bSurfEngOut = true;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;
                        CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }

                }
            }
            /*
             * SCRG file
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("SCRG"))
            {
                bSurfCrgOut = true;

                iModPdbFormatOut = 0;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;

                        // a token describing output file name
                        if (string::npos != it->find("UNIT=") || string::npos != it->find("FILE=") )
                            strScrgFile = getFile_Name_or_Format(*it,*itt);
                        // a token describing output file format
                        else if (string::npos != it->find("FRM=") || string::npos != it->find("FORM=") || string::npos != it->find("FORMAT="))
                        {
                            string strFormat = getFile_Name_or_Format(*it,*itt);

                            if (0 == strFormat.compare("PDB"))
                                iSurfCrgFormatOut = 1;
                            else if (0 == strFormat.compare("PDBA"))
                                COutOBSOLETEPDBAFormat warning(strFormat);
                            else
                                CUnknownParameter_FUNCTION warning(strFunc,*itt);
                        }
                        // invalid token
                        else
                            CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }
                }
            }
            /*
             * unconditionally turn off write function for all
             */
            else if (0 == prgstrArgTokens_UpperCase[0].compare("OFF"))
            {
                bPhimapOut      = false;
                bOutGraspSurf   = false;
                bSiteOut        = false;
                bEpsOut         = false;
                bModPdbOut      = false;
                bUnformatPdbOut = false;
                bUnformatFrcOut = false;
                bEngOut         = false;
                bGridCrgOut     = false;
                bHsurf2DatOut   = false;
                bDbOut          = false;
                bSurfEngOut     = false;
                bSurfCrgOut     = false;

                if (1 < prgstrArgTokens_UpperCase.size()) // more than one argument
                {
                    for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin()+1; it != prgstrArgTokens_UpperCase.end(); ++it)
                    {
                        itt++;
                        CUnknownParameter_FUNCTION warning(strFunc,*itt);
                    }

                }

                CUncondionalOff_WRITE warning;
            }
            /*
             * invalid token prgstrArgTokens_UpperCase[0]
             */
            else
                CUnknownParameter_FUNCTION warning(strFunc,*itt);

            break;
        }

        /*
         * ENERGY function
         * ===============
         * At present it takes as its argument any of the following:
         *    G   or GRID                      : for the grid energy,
         *    S   or SOL     or SOLVATION      : for the corrected reaction field energy
         *    AS  or ANASURF or ANALYTICSURFACE: for the analytical surface energy
         *    AG  or ANAGRID or ANALYTICGRID   : for the analytical grid energy
         *    C   or COU     or COULOMBIC      : for the coulombic energy
         *    ION or IONIC   or IONIC_C for the direct ionic contribution (see Ionic direct Contribution)
         * separated by commas. (As always there is no case sensitivity here.)
         *
         * Note that the calculation of the non linear contributions are automatically turned on whenever non-linear PBE
         * solver is invoked.
         *
         * For the energy definition we recommend the Rocchia et al. J. Phys. Chem, however a brief explanation is given below:
         *
         * The grid energy is obtained from the product of the potential at each point on the grid and the charge at that point,
         * summed over all points on the grid. However, the potential computed for each charge on the grid includes not only the
         * potentials induced by all other charges, but also the "self" potential. The effect is caused by the partitioning of
         * the real charges into the grid points. Thus, two neighboring grid points might have partial charges that originate
         * from the same real charge. Since the product of a charge with its own potential is not a true physical quantity, the
         * grid energy should not be taken as a physically meaningful number by itself. Instead, the grid energy is only meaningful
         * when comparing two DelPhi runs with exactly the same grid conditions (e.g constant structure and constant scale). The
         * difference can then be used to extract solvation energies, salt effects, and others.
         *
         * The coulombic energy is calculated using Coulomb's law. It is defined as the energy required to bring charges from
         * infinite distance to their resting positions within the dielectric specified for the molecule. This term has been revised
         * in the new DelPhi to be consistent with the new multiple dielectric model. For the most recent definition, we again refer
         * the reader to the previously mentioned paper.
         *
         * The reaction field energy (also called the solvation energy) is obtained from the product of the potential due to induced
         * surface charges with all fixed charges of the solute molecule. This includes any fixed charge in the molecule that happens
         * to be outside of the grid box. The induced surface charges are calculated at each point on the boundary between two
         * dielectrics, e.g. the surface of the molecule. If the entire molecule lies within the box and salt is absent, this energy is
         * the energy of transferring the molecule from a medium equal to the interior dielectric of the molecule into a medium of
         * external dielectric of the solution. Depending on the physical process being described, this may be the actual solvation
         * energy, but in general the solvation energy is obtained by taking the difference in reaction field energies between suitable
         * reference states - hence we make the distinction between this physical process and our calculated energy term.
         *
         * The osmotic and electrostatic stress terms, defined in the same paper, are calculated whenever the program is asked to run
         * nonlinear iterations; their contribution arises from the solvent inside the box.
         *
         * The calculation of the external ion contribution is toggled by the flag ion in the energy statement, it calculates the direct
         * interaction of ions to real charges.
         *
         * This direct calculation needs some comments: it calculates the coulombic interaction between the ions in solution (screened
         * by the surrounding polarized water). First of all, it can be computationally expensive, especially if the grid size is very
         * high. Secondly, this energy neglects all of the ionic contribution due to the ions located outside the box. This can result in
         * a underestimation of the real contribution, especially if the percentage filling is high. In the linear case this problem has
         * been reduced almost completely by a mixed numerical/analytical technique. The analytical contribution of the first term in a
         * spherical expansion of the molecular system outside the box is calculated. So the direct ionic contribution is better estimated
         * by the sum of two terms: the one internal to the box and the one external to the box.
         *
         * For the non-linear case this is not so easy, because an analytical solution of the non-linear PBE for the spherical case is not
         * available. The use of the linear approximation routine would result in an overestimation of the ionic contribution.
         *
         * Another way to get an estimation of the same energy term is to subtract the grid and reaction field energy in two cases: with
         * and without salt. This means running the program twice on the system.
         *
         * In order to help the user decide if these two runs are necessary, the program estimates and reports the number of Debye lengths
         * bounded within the grid. This number, together with the ionic strength, aids the user in estimating the error due to neglecting
         * the ionic contribution outside the box.
         *
         * Note: If periodic boundary conditions are used, the estimated external ion contribution will be taken only in the directions
         * where periodicity is NOT invoked. This is done in order to give an estimate of energy per periodic cell.
         *
         * The self-reaction field energy is calculated whenever the usual reaction field energy is calculated. It is an attempt to
         * estimate the energy that the system gains when a point charge is moved from vacuum to a dielectric medium. In fact, when a
         * charge moves from vacuum to a dielectric a polarization charge of opposite sign builds up around it. This results in an
         * energetically favorable process. The higher the dielectric constant of the medium, the higher is the amount of polarization
         * charge and the more energetically favorable is the process. In our model, it is assumed that polarization charge builds up
         * spherically around any real point charge, at a distance the user can decide via the statement radpolext=[radius], which by
         * default is set to be 1Å. This information is no attempt to give a quantitative estimation of this energy but just to take into
         * account in a more detailed way the fact that the most favorable place for a charge is within a high dielectric medium. This
         * necessity arises when a system with many different dielectric regions is studied.
         */
        case 5:
        {
            vector<string>::iterator itt = prgstrArgTokens_fromInput.begin();

            for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin(); it != prgstrArgTokens_UpperCase.end(); ++it)
            {
                if      (0 == it->compare("G") || 0 == it->compare("GRID"))
                    bGridEng  = true;
                else if (0 == it->compare("S") || 0 == it->compare("SOL") || 0 == it->compare("SOLVATION"))
                    bSolvEng  = true;
                else if (0 == it->compare("AS") || 0 == it->compare("ANASURF") || 0 == it->compare("ANALYTICSURFACE"))
                    bAnalySurfEng = true;
                else if (0 == it->compare("AG") || 0 == it->compare("ANAGRID") || 0 == it->compare("ANALYTICGRID"))
                    bAnalyEng = true;
                else if (0 == it->compare("ION") || 0 == it->compare("IONIC") || 0 == it->compare("IONIC_C"))
                    bIonsEng = true;
                else if (0 == it->compare("C") || 0 == it->compare("COU") || 0 == it->compare("COULOMBIC"))
                    bCoulombEng = true;
                else
                    CUnknownParameter_FUNCTION warning(strFunc,*itt);

                itt++;
            }

            break;
        }

        /* SITE function
         * =============
         * Reports the potentials and electrostatic field components at the positions of the subset of atoms specified in the frc file.
         * The atoms specified in frc file should not be charged in the delphi run.
         *
         * The argument is a list of identifiers that can be:
         *    ATOM              or A
         *    CHARGE            or Q
         *    POTENTIAL         or P
         *    ATOMICPOT         or ATPO
         *    DEBYEFRACTION     or DEB
         *    FIELD             or F
         *    REACTION          or R
         *    COULOMB           or C
         *    COORDINATES       or X
         *    SALT              or I
         *    TOTAL             or T
         *    REACTIONFORCE     or RF
         *    COULOMBFORCE      or CF
         *    MOLECULARDYNAMICS or MDF
         *    SURFACECHARGE     or SF
         *    TOTALFORCE        or TF
         */
        case 6:
        {
            vector<string>::iterator itt = prgstrArgTokens_fromInput.begin();

            for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin(); it != prgstrArgTokens_UpperCase.end(); ++it)
            {
                if      (0 == it->compare("ATOM") || 0 == it->compare("A"))
                    bAtomInSite = true;
                else if (0 == it->compare("CHARGE") || 0 == it->compare("Q"))
                    bCrgInSite = true;
                else if (0 == it->compare("POTENTIAL") || 0 == it->compare("P"))
                    bGridPotentialInSite = true;
                else if (0 == it->compare("ATOMICPOT") || 0 == it->compare("ATPO"))
                    bAtomPotentialInSite = true;
                else if (0 == it->compare("DEBYEFRACTION") || 0 == it->compare("DEB"))
                    bDebyeFractionInSite = true;
                else if (0 == it->compare("FIELD") || 0 == it->compare("F"))
                    bFieldInSite = true;
                else if (0 == it->compare("REACTION") || 0 == it->compare("R"))
                    bReactPotentialInSite = true;
                else if (0 == it->compare("COULOMB") || 0 == it->compare("C"))
                    bCoulombPotentialInSite = true;
                else if (0 == it->compare("COORDINATES") || 0 == it->compare("X"))
                    bAtomCoordInSite = true;
                else if (0 == it->compare("SALT") || 0 == it->compare("I"))
                    bSaltInSite = true;
                else if (0 == it->compare("TOTAL") || 0 == it->compare("T"))
                    bTotalPotentialInSite = true;
                else if (0 == it->compare("REACTIONFORCE") || 0 == it->compare("RF")) // new
                    bReactForceInSite = true;
                else if (0 == it->compare("COULOMBFORCE") || 0 == it->compare("CF")) // new
                    bCoulombForceInSite = true;
                else if (0 == it->compare("MOLECULARDYNAMICS") || 0 == it->compare("MDF")) // new
                    bMDInSite = true;
                else if (0 == it->compare("SURFACECHARGE") || 0 == it->compare("SF")) // new
                    bSurfCrgInSite = true;
                else if (0 == it->compare("TOTALFORCE") || 0 == it->compare("TF")) // new
                    bTotalForceInSite = true;
                else
                    CUnknownParameter_FUNCTION warning(strFunc,*itt);

                itt++;
            }

            //if (!(bReactForceInSite||bCoulombForceInSite||bMDInSite||bSurfCrgInSite||bTotalForceInSite))
            //   bFieldInSite = true;

            if (bReactPotentialInSite||bReactForceInSite||bMDInSite||bTotalForceInSite||bTotalPotentialInSite||bSurfCrgInSite)
                bReactFieldInFRC = true;

            if (bSurfCrgInSite) bPDB2FRCInSite = true;

            break;
        }

        /*
         * BUFFZ function
         * ==============
         * Defines a box with sides parallel to grid unit vectors that the reaction field energy will then be calculated using
         * ONLY the polarization charges contained in that box. The fixed format is BUFFZ(6i3).
         */
        case 7:
        {
            /*
             * function BUFFZ has been decided to be removed after 1st round of tests
             */

            /*
            if (1 != prgstrArgTokens_UpperCase.size() || 18 != prgstrArgTokens_UpperCase[0].length())
            {
               CUnknownParameter_FUNCTION warning(strFunc,prgstrArgTokens_UpperCase[0]);
               return false;
            }

            eiBuffz.nMin.nX = boost::lexical_cast<delphi_integer>(prgstrArgTokens_UpperCase[0].substr(0,3));
            eiBuffz.nMin.nY = boost::lexical_cast<delphi_integer>(prgstrArgTokens_UpperCase[0].substr(3,3));
            eiBuffz.nMin.nZ = boost::lexical_cast<delphi_integer>(prgstrArgTokens_UpperCase[0].substr(6,3));

            eiBuffz.nMax.nX = boost::lexical_cast<delphi_integer>(prgstrArgTokens_UpperCase[0].substr(9,3));
            eiBuffz.nMax.nY = boost::lexical_cast<delphi_integer>(prgstrArgTokens_UpperCase[0].substr(12,3));
            eiBuffz.nMax.nZ = boost::lexical_cast<delphi_integer>(prgstrArgTokens_UpperCase[0].substr(15,3));

            bIsBuffz = true;
            */

            CObsoleteFunction warning(strLineNoSpace);

            break;
        }

        /*
         * QPREF function
         * ==============
         * Obsolete. Take no action.
         */
        case 8:
        {
            CObsoleteFunction warning(strLineNoSpace);
            break;
        }

        /*
         * INSOBJ function
         * ===============
         * ONLY FOR OBJECTS. REMOVED FROM THE LIST AND TAKE NO ACTION.
         */
        case 9:
        {
            CObsoleteFunction warning(strLineNoSpace);
            break;
        }

        /*
         * SURFACE function
         * ================
         */
        case 10:
        {
            vector<string>::iterator itt = prgstrArgTokens_fromInput.begin();

            for (vector<string>::iterator it = prgstrArgTokens_UpperCase.begin(); it != prgstrArgTokens_UpperCase.end(); ++it)
            {
                if (0 == it->compare("CONNOLLY"))
                    iTypeSurf = 0;
                else if (0 == it->compare("SKIN"))
                    iTypeSurf = 1;
                else if (0 == it->compare("BLOBBY"))
                    iTypeSurf = 2;
                else if (0 == it->compare("MESH"))
                    iTypeSurf = 3;
                else if (0 == it->compare("MSMS"))
                    iTypeSurf = 4;
                else
                    CUnknownParameter_FUNCTION warning(strFunc,*itt);

                itt++;
            }

            break;
        }

        default:
        {
            CUnknownFunction warning(strFunc);
            break;
        }

        } //--------------------- end of switch (typearg)

        return true;
    }
    catch(bad_lexical_cast&)
    {
        return false;
    }
}

//-----------------------------------------------------------------------//
string CDelphiDataMarshal::getFile_Name_or_Format(string strArg_UpperCase,string strArg_fromInput)
{
    size_t pos;
    string strFile_Name_or_Format;
    string strNumber = "0123456789";
    string strLetter = "ABCDEFGHIJKLMNOPQRSTUVWXYZ.:_-+=!@#$^1234567890abcdefghijklmnopqrstuvwxyz|/?><";

    int kk2 = 0, kk3 = 0, kk4 = 0, kk5 = 0, kk6 = 0;

    if (string::npos != strArg_UpperCase.find("UNIT="))
        kk2 = 5;
    else if (string::npos != strArg_UpperCase.find("FILE="))
        kk3 = 5;
    else if (string::npos != strArg_UpperCase.find("FRM="))
        kk4 = 4;
    else if (string::npos != strArg_UpperCase.find("FORM="))
        kk5 = 5;
    else if (string::npos != strArg_UpperCase.find("FORMAT="))
        kk5 = 7;

    if (0 != kk2)
    {
        strFile_Name_or_Format = "fort.";

        for (size_t ii = kk2; ii < strArg_fromInput.length(); ++ii)
        {
            pos = strNumber.find_first_of(strArg_fromInput[ii]);
            if (string::npos == strNumber.find_first_of(strArg_fromInput[ii])) break;
            strFile_Name_or_Format += strNumber[pos];
        }
    }

    if (0 != kk3)
    {
        for (size_t ii = kk3; ii < strArg_fromInput.length(); ++ii)
        {
            if (string::npos != strLetter.find_first_of(strArg_fromInput[ii]))
                strFile_Name_or_Format += strArg_fromInput[ii];
        }
    }

    if (0 != kk4)
    {
        for (size_t ii = kk4; ii < strArg_UpperCase.length(); ++ii)
        {
            if (string::npos != strLetter.find_first_of(strArg_UpperCase[ii]))
                strFile_Name_or_Format += strArg_UpperCase[ii];
        }
    }

    if (0 != kk5)
    {
        for (size_t ii = kk5; ii < strArg_UpperCase.length(); ++ii)
        {
            if (string::npos != strLetter.find_first_of(strArg_UpperCase[ii]))
                strFile_Name_or_Format += strArg_UpperCase[ii];
        }
    }

    if (0 != kk6)
    {
        for (size_t ii = kk6; ii < strArg_UpperCase.length(); ++ii)
        {
            if (string::npos != strLetter.find_first_of(strArg_UpperCase[ii]))
                strFile_Name_or_Format += strArg_UpperCase[ii];
        }
    }

    return strFile_Name_or_Format;
}

//-----------------------------------------------------------------------//
inline vector<string> CDelphiDataMarshal::getArguments(string strArgs)
{
    /*
     * take off "(" and ")" from strArgs if they present
     */
    if ('(' == strArgs.at(0) && ')' == strArgs.at(strArgs.length()-1))
        strArgs = strArgs.substr(1,strArgs.length()-2);

    std::vector<std::string> prgstrTokens;
    istringstream issBuffer(strArgs);
    for(string token; getline(issBuffer,token,',');) prgstrTokens.push_back(token);

    /*
     *  to test the vector prgstrTokens, must include <iterator> and <algorithm> at the beginning of this file
     */
    //copy(prgstrTokens.begin(),prgstrTokens.end(),ostream_iterator<string>(cout,"."));
    //cout << '\n';

    return prgstrTokens;
}


