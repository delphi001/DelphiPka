/**
 * @file interface_exceptions.h
 * @brief warnings/exceptions occurring in interface classes
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * These warnings and exceptions should inherit from the classes CWarning and CException in order to
 * keep track the number of warnings and maintain consistent format.
 */

#ifndef INTERFACE_EXCEPTIONS_H_
#define INTERFACE_EXCEPTIONS_H_

#include "exceptions.h"
	   
/**
 * user-specified key is not found in the data container
 */
class CInexistentKey : public CException
{
   public:
	   CInexistentKey(const string & strKey)
	   {
	      cerr << "KEY " << strKey << " DOES NOT EXIST" << endl;   
	   }
};

/**
 * user-specified file is not found
 */
class CUnknownFile : public CException
{
   public:
	CUnknownFile(const string & strFileName)
	{
	   cerr << "ERROR OCCURS WHILE READING " << strFileName << " FILE : FILE DOES NOT EXIST" << endl;
	}
};

/**
 * mismatch size occurs when generating a pointer to the vector-type entry in data container
 */
class CMismatchSize : public CException
{
   public:
   CMismatchSize(const string & strKey)
   {
      cerr << "ERROR OCCURS WHILE CONVERTING THE VECTOR " << strKey << " TO MULTI-DIMENSIONAL ARRAY : THE SPECIFIED DIMESION DOES NOT MATCH" << endl;
   }
};

/**
 * skip unrecognized line(s) in the parameter file
 */
class CUnknownLine : public CWarning
{
   public:
	   CUnknownLine(const string & strFileName, const string & strLine)
	   {
	      cerr << "UNKNOWN LINE \"" << strLine << "\" IN FILE " << strFileName << "(SKIP THIS LINE...)" << endl;
	   }
};

/**
 * no bio-model or solver is described in the parameter file
 */
class CNoBiomodel : public CWarning
{
   public:
      CNoBiomodel(const string & strBiomodel, const string & strSolver)
      {
         cerr << "UNDEFINED BIOMODEL OR SOLVER (USING DEFAULT: BIOMODEL=" << strBiomodel << ", SOLVER=" << strSolver << ")" << endl;
      }
};


//class CResetArray : public CWarning
//{
//   public:
//      CResetArray(const string & strKey)
//      {
//         cerr << "THE ARRAY POINTED BY KEY \"" << strKey << "\" IS NOT NULL. \n"
//              << "\n\tTHE MEMORY ALLOCATED BEFORE WILL BE RELEASED..." << endl;
//      }
//};

#endif // INTERFACE_EXCEPTIONS_H_
