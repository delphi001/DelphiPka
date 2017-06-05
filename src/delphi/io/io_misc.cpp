/*
 * io_misc.cpp
 *
 *  Created on: Feb 4, 2014
 *      Author: chuan
 */

#include "io.h"

//-----------------------------------------------------------------------//
string CIO::toUpperCase(const string& strLine)
{
   const string strAlphabet = "abcdefghijklmnopqrstuvwxyz";
   string strNewLine;
   locale loc;
   size_t found;
   
   for (size_t i = 0; i < strLine.size(); i++)
   {
      found = strAlphabet.find_first_of(strLine[i]);
      
      if (string::npos != found)
         strNewLine += toupper(strLine[i],loc);
      else
         strNewLine += strLine[i];   
   }
  
   return strNewLine;
}

//-----------------------------------------------------------------------//
string CIO::removeSpace(const string& strLine)
{
    string strNewLine;
    
    for (size_t i = 0; i < strLine.size(); i++)
       if (' ' != strLine[i]) strNewLine += strLine[i];
    
    return strNewLine;
}

//-----------------------------------------------------------------------//
bool CIO::checkFileFormat(const string& strFile)
{
   bool bIsFormatted = true; // initially assume formatted file

   ifstream ifFileHandle;

   ifFileHandle.open(strFile.c_str());
    
   if (!ifFileHandle.is_open()) throw CUnknownIOFile(strFile);   
   
   char cTestChar[80];  
   
   ifFileHandle.read(cTestChar,80);

   string strASCI = "1234567890 .-+#,$asdfghjklzxcvbnmqwertyuiopASDFGHJKLZXCVBNMQWERTYUIOP)(}{][/";
   int iCount = 0;
   
   for (int i = 0; i < 80; i++)
      if (string::npos == strASCI.find(cTestChar[i]) ) iCount += 1;

   if (10 < iCount) bIsFormatted = false; // unformatted file
   
   ifFileHandle.close();   
   
   return bIsFormatted;
}


//-----------------------------------------------------------------------//
void CIO::setDelphiAtom(const bool& bSolvePB, const bool& bSurfCrgInSite, const string& strSizeFile, const string& strCrgFile, const string& strPdbFile, const int& iPdbFormat, const bool& bPdbUnformat, const vector<string>& strCommPDB)
{
   int    iFound  = -1;
   delphi_real   fValue  = 0.0; // radius or charge
   string strAtInf,strSub0,strSub1,strAtom,strResidue,strResidueNum,strChain; 

   cout << "\nassigning charges and radii... \n\n";

   /*
    * only standard PDB requires to read charge and size from separate files. Otherwise, charges and sizes are already
    * included in the input PDB file.
    */
   if (STDPDB == iPdbFormat)
   {
      this->readForceFile(strSizeFile); // read delphi size file
      if (bSolvePB) this->readForceFile(strCrgFile);  // read delphi charge file
   }
   else
   {
      CReadSizeCrgFromPDB warning;
   }
    
   this->readPdbFile(strPdbFile,iPdbFormat,bPdbUnformat, strCommPDB); // read pdb file

   if (!bExistRadiiInfo) // pdb file, such as std pdb, does not have radii info
   {
      for (delphi_integer iThisAtom = 0; iThisAtom < iAtomNum; iThisAtom++)
      {
         strAtInf = vctapAtomPdb[iThisAtom].getAtInf(); strAtInf = toUpperCase(strAtInf);

         /*
          * for certain pdb files such as 1AB1.pdb, position 5 (maybe and 9) are for special purpose
          */
         //strAtom       = strAtInf.substr(0,6);  strAtom       = removeSpace(strAtom);
         //strResidue    = strAtInf.substr(6,4);  strResidue    = removeSpace(strResidue);
         //strChain      = strAtInf.substr(10,1); strChain      = removeSpace(strChain);
         //strResidueNum = strAtInf.substr(11,4); strResidueNum = removeSpace(strResidueNum);

         strAtom       = strAtInf.substr(0,5);  strAtom       = removeSpace(strAtom);
         strResidue    = strAtInf.substr(6,3);  strResidue    = removeSpace(strResidue);
         strChain      = strAtInf.substr(10,1); strChain      = removeSpace(strChain);
         strResidueNum = strAtInf.substr(11,4); strResidueNum = removeSpace(strResidueNum);

         // assign radius, searching for decreasingly specific specification ending with generic atom type
         // note all atoms must have an assignment         
         fValue = 0.0;

         iFound = FindRecord(strAtom,strResidue,strResidueNum,strChain,SIZEFILE,fValue);
      
         if (-1 == iFound)
         {
            //throw  CUnknownRadius(strAtInf); // need stop here
            CUnknownRadius warning(strAtInf);
         }

         strSub0 = strAtom.substr(0,1); strSub1 = strAtom.substr(1,1);
         if (1.0e-6 > fValue && 0 != strSub0.compare("H") &&  0 != strSub1.compare("H"))
            CZeroHeavyAtomRadius warning(strAtInf);
      
         vctapAtomPdb[iThisAtom].setRadius(fValue);
         
         fValue = 0.0;
         
         if (bSolvePB)
         {
            iFound = FindRecord(strAtom,strResidue,strResidueNum,strChain,CHARGEFILE,fValue);
            if (-1 == iFound) CUnknownCharge warning(strAtInf); // need stop here
            vctapAtomPdb[iThisAtom].setCharge(fValue);
         }
      } //---------- end of loop over iThisAtom =0:iAtomNum-1 
   } //---------- end of if (!bExistRadiiInfo)

   // check charge 
   // inc_setrcmod_chkcrg but combined the two cases of isitsf = .true. and .false.
   string strLastResidue    = "    "; // a4
   string strLastResidueNum = "    "; // a4
   delphi_real   fCrgSum,fError; // total net charge

   fCrgSum = 0.0;
   
   for (delphi_integer iThisAtom = 0; iThisAtom < iAtomNum; iThisAtom++)
   {
      fValue = vctapAtomPdb[iThisAtom].getCharge();
         
      strAtInf      = vctapAtomPdb[iThisAtom].getAtInf();
      strResidue    = strAtInf.substr(6,4); 
      strResidueNum = strAtInf.substr(11,4);
         
      if ( 0 != strResidue.compare(strLastResidue) || 0 != strResidueNum.compare(strLastResidueNum) )
      {
         if (1.0e-4 < abs(fCrgSum))
         {
            fError = abs(fCrgSum) - 1.0;
               
            if (1.0e-4 < abs(fError)) 
               CNonZeroNetCrg waring(strLastResidue,strLastResidueNum,fCrgSum);
         }
         
         strLastResidue    = strResidue;
         strLastResidueNum = strResidueNum;
         fCrgSum           = fValue; 
      }
      else
         fCrgSum += fValue;
          

      // specific section to find maximum residue number. Goes here when isitsf = .true.
      if (bSurfCrgInSite)
      {
         delphi_integer iNumber =  atoi( strResidueNum.c_str() );
         if (iNumber > iResidueNum) iResidueNum = iNumber;
      }
   }
  
   if (1.0e-4 < abs(fCrgSum))
   {
      fError = abs(fCrgSum) - 1.0;
      if (1.0e-4 < abs(fError))
         CNonZeroNetCrg waring(strLastResidue,strLastResidueNum,fCrgSum);
   }  

#ifdef DEBUG_IO_PDB
   printPDB();
#endif

}                        
















