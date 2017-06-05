/**
 * @file io.h
 * @brief IO class for reading/writing various standard files
 *
 * @author Chuan Li, chuanli@clemson.edu
 */

#ifndef CIO_H_
#define CIO_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>   // std::string
#include <vector>   // std::vector
#include <locale>
#include <stdlib.h> // atof
#include <cmath>

#include "../interface/environment.h"
#include "../misc/misc_grid.h"
#include "io_datatype.h"
#include "io_exceptions.h"

using namespace std;

const int SIZEFILE   = 1;
const int CHARGEFILE = 2;

const int STDPDB     = 10;
const int MODPDB     = 20;
const int PQRPDB     = 21;
const int MOD4PDB    = 30;
const int PQR4PDB    = 31;

/* CommPDB format is updated and used for communication with DelPhi and PrimePKA without read/write files
 * Updated: Dec 15, 2014 by Lin Wang 
 */
const int CommPDB    = 40;


class CIO
{
   protected:
      /**
       * F95 var.: repsin
       */
      const delphi_real fDielec;

      /**
       * F95 var.: epkt
       */
      const delphi_real fEPKT;

      /*
       * ******************** miscellany functions defined in io_misc.cpp *******************
       */

      /**
       * function to convert a string to upper case
       * @param strLine input string
       * @return converted string in upper case
       */
      string toUpperCase(const string& strLine);

      /**
       * function to remove spaces from a string
       * @param strLine input string
       * @return converted string without spaces
       */
      string removeSpace(const string& strLine);

      /**
       * function to check if the input file is formatted or not
       * @param strFile name of the input file
       * @return true if the file formatted, false otherwise.
       */
      bool checkFileFormat(const string& strFile);
                 
      /*
       * ******************* for reading atom force (radii/size and charge) file *******************
       * ********** these static attributes are declared here and defined in io_force.cpp **********
       */

      /**
       * F95 var.: nrmax
       * # of entries in radius file
       */
      static delphi_integer iRadiusNum;

      /**
       * F95 var.: radii(:)
       * a vector containing all radius read from .siz file
       */
      static vector<CForce> prgas;

      /**
       * F95 var.: nchrec
       * # of entries in charge file
       */
      static delphi_integer iCrgNum;

      /**
       * F95 var.: charge(:)
       * a vector containing all charges read from .crg file
       */
      static vector<CForce> prgac;

      /**
       * read a .siz or .crg file in non-PK format
       * @param ifFileStream input file stream
       * @param iFileType input file type
       */
      void readFileInNotPKFormat(ifstream& ifFileStream, const int& iFileType);

      /**
       * read a .siz or .crg file in PK format
       * @param ifFileStream input file stream
       * @param iFileType input file type
       */
      void readFileInPKFormat(ifstream& ifFileStream, const int& iFileType);

#ifdef DEBUG_IO_FORCE
      /**
       * debug function to print vector storing values read from .siz or .crg file
       * @param iFileType input file type
       * @prgaf pointer pointing to the vector of values read from .siz or .crg file
       */
      void printForce(const int& iFileType, const vector<CForce>& prgaf);
#endif

      /**
       * function to get index of a specified atom in the record
       * @param strAtom atom name
       * @param strResidue residue name
       * @param strResidueNum residue number
       * @param strChain chain name
       * @param iRecordNum total number of records
       * @param prgaf pointer to the vector of records
       * @return -1 if no record is found, index otherwise.
       */
      delphi_integer FindRecordIndex(const string& strAtom,const string& strResidue,const string& strResidueNum,
                                     const string& strChain,const delphi_integer& iRecordNum,vector<CForce>* prgaf);

      /**
       * function to find a particular record matching the specified atom
       * @param strAtom atom name
       * @param strResidue residue name
       * @param strResidueNum residue number
       * @param strChain chain name
       * @param iFileType file type
       * @param fValue value (charge or radius) associated with this atom
       * @return -1 if unfounded, index otherwise
       */
      delphi_integer FindRecord(const string& strAtom,const string& strResidue,const string& strResidueNum,
                                const string& strChain,const int& iFileType, delphi_real& fValue);

      //void checkCharge(const bool& bSurfCrgInSite);

      /*
       * ********** functions defined in io_pdb.cpp for reading PDB file in various formats ********
       */

      /**
       * F95 var.: imedianumb
       * number of objects read from .pdb file
       */
      delphi_integer iObjectMediaNum;

      /**
       * function to read standard .pdb file
       * @param ifPdbFileStream input .pdb file stream
       */
      void readStdPdbFile(ifstream& ifPdbFileStream);

      /**
       * function to read modified .pdb file (type = 1)
       * @param ifPdbFileStream input .pdb file stream
       */
      void readModFile1(ifstream& ifPdbFileStream);

      /**
       * function to read modified .pdb file (type = 4)
       * @param ifPdbFileStream input .pdb file stream
       */
      void readModFile4(ifstream& ifPdbFileStream);

      /**
       * function to read .pqr file
       * @param ifPdbFileStream input .pqr file stream
       */
      void readPqrFile(ifstream& ifPdbFileStream);

      /**
       * function to read mod4 pdb file
       * @param ifPdbFileStream input .pdb file stream
       */
      void readMod4File(ifstream& ifPdbFileStream);

      /**
       * function to read pqr4 file
       * @param ifPdbFileStream input .pdb file stream
       */
      void readPqr4File(ifstream& ifPdbFileStream);

      /**
       * function to read unformatted pdb file
       * @param strPdbFile input .pdb file name
       * @param ifPdbFileStream input .pdb file stream
       * @param bPostProcess flag for post processing
       */
      void readUnformattedPdb(const string& strPdbFile,ifstream& ifPdbFileStream,bool& bPostProcess);
    
      void commVecPDB(const vector<string>& strCommPDB);
    
#ifdef DEBUG_IO_PDB
      /**
       * debug function to print values read from .pdb file
       */
      void printPDB();                        
#endif
           
      /**
       * F95 var.: iatrad
       * whether there is radius info
       */
      bool    bExistRadiiInfo;

      /**
       * F95 var.: tmpiatmmed(nobjectmax)
       * vector containing internal media-number per object
       */
      vector<delphi_integer> prgiObjectMediaNum;

   public:     

      /**
       * F95 var.: nmedia
       * # of media
       */
      delphi_integer    iMediaNum;

      /**
       * F95 var.: nobject
       * # of objects
       */
      delphi_integer    iObjectNum;

      /**
       * F95 var.: natom
       * # of atoms
       */
      delphi_integer    iAtomNum;

      /**
       * F95 var.: resnummax
       * maximum residue number
       */
      delphi_integer    iResidueNum;

      /**
       * F95 var.: ionlymol
       * true if there are only molecules in the system (no objects)
       */
      bool       bOnlyMolecule;

      /**
       * F95 var.: delphipdb(Natom)
       * array of structure to store info read from pdb file
       */
      vector<CAtomPdb> vctapAtomPdb;

      /**
       * F95 var.: medeps(0:nmediamax)
       * vector containing correspondence media<->epsilon/epkt
       */
      vector<delphi_real>    vctfMediaEps;

      /**
       * F95 var.: dataobject(nobjectmax,2)
       * vector containing string with object data, and pre-elab data
       * changed it to vctstrObject(2*nobjectmax)
       */
      vector<string>  vctstrObject;

      /**
       * F95 var.: iatmmed(Natom+Nobjectmax)
       * vector containing internal media-number per atom and object
       */
      vector<delphi_integer> vctiAtomMediaNum;

      /**
       * constructor I using default Dielectric constant and epkt values
       */
      CIO():fDielec(4.0),fEPKT(167100.9162872952/297.3342119)
      {     
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                     CIO is constructed                       *\n";
         cout << "****************************************************************\n";
#endif

         //iRadiusNum      = 0;
         //iCrgNum         = 0;
         iMediaNum       = 1;
         iObjectNum      = 1;  
         iAtomNum        = 0;
         iObjectMediaNum = 1;
         bOnlyMolecule   = true;
         bExistRadiiInfo = false;
         iResidueNum     = 0;        
      };

      /**
       * constructor II using user-specified Dielectric constant and epkt values
       */
      CIO(delphi_real fDielecIn,delphi_real fEPKTIn):fDielec(fDielecIn),fEPKT(fEPKTIn)
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                     CIO is constructed                       *\n";
         cout << "****************************************************************\n";
#endif

         //iRadiusNum      = 0;
         //iCrgNum         = 0;
         iMediaNum       = 1;
         iObjectNum      = 1;  
         iAtomNum        = 0;
         iObjectMediaNum = 1;
         bOnlyMolecule   = true;
         bExistRadiiInfo = false;
         iResidueNum     = 0;        
      };
      
      /**
       * destructor
       */
      ~CIO()
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                      CIO is destroyed                        *\n";
         cout << "****************************************************************\n";
#endif
      };


      /**
       * function for reading atom force (radii/size and charge) file
       * @param strFile file name
       */
      void readForceFile(const string& strFile);

      /**
       * function for reading PDB file in various formats
       * @param strPdbFile pdb file name
       * @param iPdbFormat pdb file format
       * @param bPdbUnformat true if unformatted pdb file
       */
      void readPdbFile(const string& strPdbFile,const int& iPdbFormat,const bool& bPdbUnformat, const vector<string>& strCommPDB);

      /**
       * function for writing unformatted pdb file
       * @param strPdbFile pdb file name
       */
      void writeUnformatPdb(const string& strPdbFile);

      /**
       * function for writing modified pdb file
       * @param strPdbFile pdb file name
       * @param iModPdbFormatOut file format
       */
      void writeModifiedPdb(const string& strPdbFile,const int& iModPdbFormatOut);
                       
      /*
       * miscellany
       */

      /**
       * function to set DelPhi-style atom list vctapAtomPdb(iAtomNum), i.e, (delphipdb(natom))
       * @param bSolvePB true if solving PB equation
       * @param bSurfCrgInSite true if writing surface chagre info into site file
       * @param strSizeFile name of .siz file
       * @param strCrgFile name of .crg file
       * @param strPdbFile name of .pdb file
       * @param iPdbFormat .pdb file format
       * @param bPdbUnformat true if unformatted pdb file
       */
      void setDelphiAtom(const bool& bSolvePB,const bool& bSurfCrgInSite,const string& strSizeFile,
                         const string& strCrgFile,const string& strPdbFile,const int& iPdbFormat,
                         const bool& bPdbUnformat, const vector<string>& strCommPDB);
                                                 
      /**
       * function for reading FRC file
       * @param strFrcFile name of frc file
       * @param fgOffCenter offset from box center
       * @param fScale scale used for calculations
       * @return position of box center
       */
      SGrid<delphi_real> readFrcFile(const string& strFrcFile,const SGrid<delphi_real>& fgOffCenter,
                                     const delphi_real& fScale);

      /**
       * function for writing EPS file
       * @param iAtomNumIn number of atoms
       * @param iObjectNumIn number of objects
       * @param iGrid number of grid points on each direction
       * @param fScale scale used for calculations
       * @param fgBoxCenter position of box center
       * @param vctigEpsMap vector of epsilon map
       * @param vctbDielecMap vector of dielectric constants map
       * @param strEpsFile output eps file name
       */
      void writeEpsMap(const delphi_integer& iAtomNumIn,const delphi_integer& iObjectNumIn,
                       const delphi_integer& iGrid,const delphi_real& fScale,const SGrid<delphi_real>& fgBoxCenter,
                       const vector< SGrid<delphi_integer> >& vctigEpsMap,const vector<bool>& vctbDielecMap,
                       const string& strEpsFile);
      
      /*
       * for reading/writing PHI file
       */

      /*
       *for writing ENERGY file
       */

      /*
       * for writing HSURF2 file
       */

      /*
       *for writing SURFEN file
       */

      /*
       * for writing SCRG file
       */
    

};

#endif // CIO_H_
