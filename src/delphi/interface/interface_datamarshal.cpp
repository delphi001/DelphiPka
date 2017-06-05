/*
 * interface_datamarshal.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: chuan
 */

#include "interface_datamarshal.h"

using namespace std;

//-----------------------------------------------------------------------//
//                                                                       //
//                       protected members                               //
//                                                                       //
//-----------------------------------------------------------------------//

//-----------------------------------------------------------------------//
bool IDataMarshal::getBiomodel(const string &strLineNoSpace)
{
   string strLineUpperCase, strStatement, strArgument;

   size_t found;

   locale loc;

   static bool bNoBiomodel = true, bNoSolver = true;

   // transform the statement to upper case
   for (size_t ii = 0; ii < strLineNoSpace.length(); ++ii)
		strLineUpperCase += toupper(strLineNoSpace[ii],loc);

   // find the equal sign which is used to separate the statement and argument
   found = strLineUpperCase.find_first_of("=");

   // get the statement and the argument of the statement
   strStatement = strLineUpperCase.substr(0,found);
   strArgument  = strLineUpperCase.substr(found+1,strLineUpperCase.size()-found-1);

   if ("BIOMODEL" == strStatement)
   {
      strBioModel = strArgument; bNoBiomodel = false;
      return true;
   }

   if ("SOLVER" == strStatement)
   {
      strNumSolver = strArgument; bNoSolver = false;
      return true;
   }

   if (bNoBiomodel || bNoSolver)
   {
      CNoBiomodel warning(strBioModel,strNumSolver);

      bNoBiomodel = false; bNoSolver = false;

      return false;
   }

   return false;
}

//-----------------------------------------------------------------------//
bool IDataMarshal::getQinclude(const string &strLineNoSpace)
{
    if ((string::npos == strLineNoSpace.find("qinclude")) && (string::npos == strLineNoSpace.find("QINCLUDE")))
       return false;

    size_t pos = strLineNoSpace.length();

    string strFileName = strLineNoSpace.substr(9,pos-10);

    cout << "QINCLUDE found - reading file " << strFileName.c_str() << endl;

    this->read(strFileName);

    return true;
}

//-----------------------------------------------------------------------//
bool IDataMarshal::getParameter(const string &strLineNoSpace)
{
   bool bSuccess = false, bIsStatement = false;

   size_t bb=0;

   // both functions and statements can have equal signs, but only the
   // statements have them outside of brackets.
   for (size_t ii=0; ii<strLineNoSpace.size(); ++ii)
   {
      if ('(' == strLineNoSpace[ii])
         bb++;
	  else if (')' == strLineNoSpace[ii])
         bb--;
	  else if ( ('=' == strLineNoSpace[ii]) && (0 == bb) )
         {bIsStatement = true;
         //cout << "##in interface_datamarshal.cp: " << strLineNoSpace << " " << bIsStatement << endl;//LinLi
         }
   }

   if (bIsStatement)
      {
          //cout << "##in interface_datamarshal.cp: " << strLineNoSpace << endl;//LinLi
          bSuccess = getStatement(strLineNoSpace);
      }
   else
      bSuccess = getFunction(strLineNoSpace);

   return bSuccess;
}

//-----------------------------------------------------------------------//
//                                                                       //
//                          public members                               //
//                                                                       //
//-----------------------------------------------------------------------//

void IDataMarshal::read(string strFileName)
{
   ifstream ifFileHandle;

   // open the file with name strFileName
   ifFileHandle.open(strFileName.c_str());

   // if the file doesnt exists, exit
   if (!ifFileHandle.is_open()) throw CUnknownFile(strFileName);

   string strLineFromFile, strLineNoSpace; // strings used to read

   do // execute the loop until a break occurs
   {
      getline(ifFileHandle, strLineFromFile); // get a line from the file and store it at the string @strLineFromFile

      if (ifFileHandle.eof()) break; // if we have reached the end of the file, exit the loop

          //cout << "## in Read1: " << strLineFromFile << endl; //LinLi

      for (size_t ii = 0; ii<strLineFromFile.size(); ++ii) // for each letter of the line
      {
         if ('!' == strLineFromFile[ii])
            break; // if we have reached a comment, skip
         else if (' ' != strLineFromFile[ii])
            strLineNoSpace += strLineFromFile[ii]; // if the character is not a white space, add the character
      }
          //cout << "## in Read2: " << strLineNoSpace << endl; //LinLi

      if (0 == strLineNoSpace.size()) continue; // a comment line or a line of spaces, skip


      /*
       * Newlines in DOS and Windows end with the combination of two characters '\r\n',
       * while they end with a single '\n' to indicate a new line in Unix and a single '\r' in Mac.
       * Remove '\r' in the string if it exists (for Windows machine).
       */
       if('\r' == strLineNoSpace[strLineNoSpace.size()-1])
       {
           strLineNoSpace.erase( strLineNoSpace.size()-1);
       }

      if (getBiomodel(strLineNoSpace)) // a line indicating either biomodel or numerical solver
      { strLineNoSpace.clear(); continue; }

      if (getQinclude(strLineNoSpace)) // a line of qinclude
	   { strLineNoSpace.clear(); continue; }

      // call the function getParameter(). If the result is false, there was a problem, and an error message is displayed
      if (!getParameter(strLineNoSpace))
      {

         //throw CUnknownLine(strFileName, strLineFromFile);
         CUnknownLine warning(strFileName, strLineFromFile);
      }

      strLineNoSpace.clear();

   } while (true);

   ifFileHandle.close();

}

