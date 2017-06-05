/**
 * @file io_exceptions.h
 * @brief warnings/exceptions occurring in class CIO
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * These warnings and exceptions should inherit from the classes CWarning and CException in order to
 * keep track the number of warnings and maintain consistent format.
 */

#ifndef IO_EXCEPTIONS_H_
#define IO_EXCEPTIONS_H_

#include "../interface/exceptions.h"

/**
 * user-specified IO file does not exist
 */
class CUnknownIOFile : public CException
{
   public:
      CUnknownIOFile(const string & strFileName)
      {
         cerr << "ERROR OCCURS WHILE READING IO FILE " << strFileName << " : FILE DOES NOT EXIST" << endl;
      }
};

/**
 * unknown header in .siz or .crg file
 */
class CUnknownForceFileHeader : public CException
{
   public:
      CUnknownForceFileHeader(const string & strFile, const string & strHeader)
      {
         cerr << "UNKNOWN HEADER IN FILE " << strFile << " : " << strHeader << endl;
         cerr << "\t VALID HEADER FOR SIZE FILE IS \'atom__res_radius\' OR \'atom__resnumbc_radius_\'" << endl;
         cerr << "\t VALID HEADER FOR CHARGE FILE IS \'atom__res_radius\' OR \'atom__resnumbc_radius_\'" << endl;
      }
};

/**
 * unrecognized delphi .pdb file
 */
class CUnknownDelphiPdb : public CException
{
   public:
      CUnknownDelphiPdb(const string & strFile)
      {
         cerr << "THE FILE " << strFile << " IS NOT A DELPHI PDB FILE! " << endl;
         cerr << "\t CHECK THE FILE TYPE! " << endl;           
      }
};

/**
 * unrecognized unformatted .pdb file
 */
class CUnknownUnformattedPdb : public CException
{
   public:
      CUnknownUnformattedPdb(const string & strFile)
      {
         cerr << "ERROR OCCURS WHEN READING THE UNFORMATTED PDB FILE " << strFile << endl;
         cerr << "\t CHECK THE FILE AND TRY AGAIN! " << endl;           
      }
};

/**
 * unrecognized modified .pdb file
 */
class CUnknownModPdb : public CException
{
   public:
      CUnknownModPdb(const string & strFile, const int & iFormatID)
      {
         cerr << "THE FORMAT ID " << iFormatID << " SPECIFIED IN THE MOD PDB FILE " << strFile << " IS NO LONG SUPPORTED." << endl;
         cerr << "\t CURRENT SUPPORTED ID = 1, 4 FOR MOD PDB FILE" << endl;           
      }
};

/**
 * obsolete line in .pdb file
 */
class CUnsupportedCRGDST : public CException
{
   public:
      CUnsupportedCRGDST(const string & strLine)
      {
         cerr << "THE INPUT PDB FILE CONTAINS A LINE \"" << strLine << "\" WITH KEYWORD \"crgdst\" WHICH IS NO LONGER SUPPORTED " << endl;
         cerr << "\t PLEASE UPDATE THE PDB FILE AND TRY AGAIN..." << endl;           
      }
};

/**
 * no object or atom after reading user-specified files
 */
class CNoAtomObjectData : public CException
{
   public:
      CNoAtomObjectData()
      {
         cerr << "# OF ATOMS = 0 AND # OF OBJECTS = 0" << endl;
         cerr << "\t EXITING DUE T NON-EXISTENCE OF ATOM FILE NOR OBJECT DATA" << endl;
      }
};

/**
 * no radius record for specified atom
 */
class CUnknownRadius : public CWarning
{
   public:
      CUnknownRadius(const string & strLine)
      {
         cerr << "NO RADIUS RECORD FOUND FOR \"" << strLine << "\" " << endl;
         //cerr << "\t PLEASE CHECK THE .siz FILE AND RETRY..." << endl;
      }
};

/**
 * radius of specified heave atom is set to be zero by default
 */
class CZeroHeavyAtomRadius : public CWarning
{
   public:
      CZeroHeavyAtomRadius(const string & strLine)
      {
         cerr << "RADIUS OF HEAVEY ATOM IN \"" << strLine << "\" IS SET TO ZERO"<< endl;
      }
};

/**
 * charge of specified atom is set to be zero by default
 */
class CUnknownCharge : public CWarning
{
   public:
      CUnknownCharge(const string & strLine)
      {
         cerr << "NO CHARGE RECORD FOUND FOR \"" << strLine << "\" (CHARGE IS SET = 0)" << endl;
      }
};

/**
 * non-zero net charge of a residue
 */
class CNonZeroNetCrg : public CWarning
{
   public:
      CNonZeroNetCrg(const string & strResidue, const string & strResidueNum, const delphi_real & fNetCharge)
      {
         cerr << "\"" << strResidue << strResidueNum << "\" HAS A NET CHARGE OF " << fNetCharge << endl;
      }
};

/**
 * empty atom in .pdb file
 */
class CEmptyAtomsInFile : public CWarning
{
   public:
      CEmptyAtomsInFile(const string & strFile)
      {
         cerr << "EMPTY ATOMS FOR MIDPOINT DETERMINATION IN FILE "  << strFile << " (ASSUMING ZERO OFFSETS)" << endl;
      }	   
};	   

/**
 * using modified .pdb file for sizes and charges as user specifies
 */
class CReadSizeCrgFromPDB : public CWarning
{
   public:
      CReadSizeCrgFromPDB()
      {
         cerr << "IN PARAMETER FILE, USER SPECIFIES TO READ PDB FILE IN MODIFIED FORMAT "
              << "(CHARGE AND SIZE ARE READ FROM THE PDB FILE AS WELL) \n";
      }
};

#endif // IO_EXCEPTIONS_H_
