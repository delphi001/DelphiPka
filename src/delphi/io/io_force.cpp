/*
 * io_force.cpp
 *
 *  Created on: Feb 4, 2014
 *      Author: chuan
 */

#include "io.h"

/*
 * Declaring a static member variable is not enough, you need to also define it somewhere.
 * This should be in a .cpp file that includes the header (.h) file that declared them.
 */

delphi_integer CIO::iRadiusNum = 0;
vector<CForce> CIO::prgas;
delphi_integer CIO::iCrgNum    = 0;
vector<CForce> CIO::prgac;


//-----------------------------------------------------------------------//
// atom(a6),residue(a3),radius(f8.4) or charge(f8.4)
// aaaaaarrrfff.ffff
void CIO::readFileInNotPKFormat(ifstream& ifFileStream, const int& iFileType)
{   
   //getForceFileRecordNum(ifFileStream,iFileType);

   vector<CForce> * prgaf; // pointer to an array fo CForce  
 
   if      (SIZEFILE   == iFileType) { prgaf = &prgas; }
   else if (CHARGEFILE == iFileType) { prgaf = &prgac; }

   string strLine, strSubLine;   

   ifFileStream.clear(); ifFileStream.seekg(0); // rewind

   while (!ifFileStream.eof()) // skip comments and header
   {
      getline(ifFileStream,strLine);      

      strSubLine = removeSpace(strLine); 

      if (0 == strSubLine.compare(0,1,"")) continue; // ignore possible empty lines      

      if ('!' != strLine[0]) break; // header
   }

   CForce tmpForceObj;

   while (!ifFileStream.eof())
   {  
      getline(ifFileStream, strLine);

      if ('!' == strLine[0]) continue; // skip commented line

      // to ignore possible empty lines
      strSubLine = removeSpace(strLine); 
      if (0 == strSubLine.compare(0,1,"")) continue;
      
      strSubLine = strLine.substr(0,6); // atom name, 6 characters
      
      if (0 != strSubLine.compare(0,6,"      "))
      {        
         strSubLine = removeSpace(strSubLine); 
         strSubLine = toUpperCase(strSubLine);
      }
      else
         strSubLine = " ";

      tmpForceObj.setAtom(strSubLine);
      
      strSubLine = strLine.substr(6,3); // residue name, 3 characters

      if (0 != strSubLine.compare(0,3,"   "))
      {
         strSubLine = removeSpace(strSubLine); 
         strSubLine = toUpperCase(strSubLine); 
      }        
      else
         strSubLine = " ";

      tmpForceObj.setResidue(strSubLine);
         
      tmpForceObj.setResidueNum(" "); // residue number

      tmpForceObj.setChain(" "); // subunit name

      strSubLine = strLine.substr(9,8); // radius or charge
         
      strSubLine = removeSpace(strSubLine);

      tmpForceObj.setValue( atof(strSubLine.c_str()) );
      
      prgaf->push_back(tmpForceObj);
   } // ---------- end of while (!ifFileStream.eof())  
}


//-----------------------------------------------------------------------//
// atom(a6),residue(a3),residue_number(a4),subunit(a1),radius(f8.4)
// aaaaaarrrnnnncfff.ffff
void CIO::readFileInPKFormat(ifstream & ifFileStream, const int & iFileType)
{
   vector<CForce> * prgaf; // pointer to an array fo CForce  
 
   if      (SIZEFILE   == iFileType) { prgaf = &prgas; }
   else if (CHARGEFILE == iFileType) { prgaf = &prgac; }
 
   string strLine, strSubLine;
   
   ifFileStream.clear(); ifFileStream.seekg(0); // rewind

   while (!ifFileStream.eof()) // skip comments and header
   {
      getline(ifFileStream,strLine);      

      strSubLine = removeSpace(strLine); 

      if (0 == strSubLine.compare(0,1,"")) continue; // ignore possible empty lines      

      if ('!' != strLine[0]) break; // header
   }
  
   CForce tmpForceObj;

   while (!ifFileStream.eof())
   {  
      getline(ifFileStream, strLine);

      if ('!' == strLine[0]) continue; // skip commented line

      // to ignore possible empty lines
      strSubLine = removeSpace(strLine); 
      if (0 == strSubLine.compare(0,1,"")) continue;

      strSubLine = strLine.substr(0,6); // atom name, 6 characters
      
      if (0 != strSubLine.compare(0,6,"      "))
      {        
         strSubLine = removeSpace(strSubLine); 
         strSubLine = toUpperCase(strSubLine);
      }
      else
         strSubLine = " ";

      tmpForceObj.setAtom(strSubLine);

      strSubLine = strLine.substr(6,3); // residue name, 3 characters

      if (0 != strSubLine.compare(0,3,"   "))
      {
         strSubLine = removeSpace(strSubLine); 
         strSubLine = toUpperCase(strSubLine); 
      }        
      else
         strSubLine = " ";
            
      tmpForceObj.setResidue(strSubLine);      

      strSubLine = strLine.substr(9,4); // residue number, 4 characters
         
      if (0 != strSubLine.compare(0,4,"    "))
      {
         strSubLine = removeSpace(strSubLine); 
         strSubLine = toUpperCase(strSubLine); 
      }        
      else
         strSubLine = " ";         
         
      tmpForceObj.setResidueNum(strSubLine);      

      strSubLine = strLine.substr(13,1); // subunit name, 1 characters
         
      if (0 != strSubLine.compare(0,1," "))
      {
         strSubLine = removeSpace(strSubLine); 
         strSubLine = toUpperCase(strSubLine); 
      }        
      else
         strSubLine = " ";         
         
      tmpForceObj.setChain(strSubLine);    

      strSubLine = strLine.substr(14,8); // radius or charge
        
      strSubLine = removeSpace(strSubLine);
         
      tmpForceObj.setValue( atof(strSubLine.c_str()) ); 
      
      prgaf->push_back(tmpForceObj);
   } // ---------- end of while (!ifFileStream.eof())
}


//-----------------------------------------------------------------------//
delphi_integer CIO::FindRecordIndex(const string& strAtom, const string& strResidue, const string& strResidueNum, const string& strChain, const delphi_integer& iRecordNum, vector<CForce>* prgaf)
{
   delphi_integer iThisRecord;

   string strThisAtom,strThisResidue,strThisChain,strThisResidueNum;   

   for (iThisRecord = 0; iThisRecord < iRecordNum; iThisRecord++ )
   {
      strThisAtom       = prgaf->at(iThisRecord).getAtom();
      strThisResidue    = prgaf->at(iThisRecord).getResidue();
      strThisResidueNum = prgaf->at(iThisRecord).getResidueNum();
      strThisChain      = prgaf->at(iThisRecord).getChain();
                        
      if ( ("*" == strAtom || 0 == strThisAtom.compare(strAtom)) && 0 == strThisResidue.compare(strResidue) && 0 == strThisResidueNum.compare(strResidueNum) && 0 == strThisChain.compare(strChain) )
    	  return iThisRecord;
   }
   
   return -1; // -1 indicates no record found 
}


//-----------------------------------------------------------------------//
delphi_integer CIO::FindRecord(const string& strAtom, const string& strResidue, const string& strResidueNum, const string& strChain, const int& iFileType, delphi_real& fValue)
{
   delphi_integer iRecordNum = 0;

   vector<CForce> * prgaf; // pointer to an array fo CForce 

   if      (SIZEFILE   == iFileType) { iRecordNum = iRadiusNum; prgaf = &prgas; }
   else if (CHARGEFILE == iFileType) { iRecordNum = iCrgNum;    prgaf = &prgac; }

   delphi_integer iFound = -1; fValue = 0.0;

   string strEmptyAtom       = " "; // a6
   string strWildAtom        = "*"; // a6
   string strEmptyResidue    = " "; // a3
   string strEmptyResidueNum = " "; // a4
   string strEmptyChain      = " "; // a1

   //---------- Search for record (atm ; res ; rnum ; ch) 
   iFound = FindRecordIndex(strAtom,strResidue,strResidueNum,strChain,iRecordNum,prgaf);

   if (-1 == iFound)
   {
      //---------- Search for record (atm ; res ; rnum ; ' ')
      iFound = FindRecordIndex(strAtom,strResidue,strResidueNum,strEmptyChain,iRecordNum,prgaf);
      if (-1 == iFound)
      {
         //---------- Search for record (atm ; res ; ' ' ; ' ') 
         iFound = FindRecordIndex(strAtom,strResidue,strEmptyResidueNum,strEmptyChain,iRecordNum,prgaf);
         if (-1 == iFound)
         {
            //---------- Search for record (atm ; ' ' ; ' ' ; ' ') 
            iFound = FindRecordIndex(strAtom,strEmptyResidue,strEmptyResidueNum,strEmptyChain,iRecordNum,prgaf);
            if (-1 == iFound)
            {
               strEmptyAtom = strAtom.substr(0,1); strEmptyAtom.append("     ");
               
               //---------- Search for record (generic_atom_name;res; ' ';' ')
               iFound = FindRecordIndex(strEmptyAtom,strEmptyResidue,strEmptyResidueNum,strEmptyChain,iRecordNum,prgaf);
               if (-1 == iFound)
               {
                  //---------- Search for record (* ; ' ' ; ' ' ; ch)
                  iFound = FindRecordIndex(strWildAtom,strEmptyResidue,strEmptyResidueNum,strChain,iRecordNum,prgaf);
                  if (-1 == iFound)
                  {
                     //---------- Search for record (atm ; res ; ' ' ; ch)
                     iFound = FindRecordIndex(strAtom,strResidue,strEmptyResidueNum,strChain,iRecordNum,prgaf);
                     if (-1 == iFound)
                     {
                        //---------- Search for record (atm ; ' ' ; rnum ; ch)
                        iFound = FindRecordIndex(strAtom,strEmptyResidue,strResidueNum,strChain,iRecordNum,prgaf);
                        if (-1 == iFound)
                        {
                           //---------- Search for record (atm;' ';' ';ch)
                           iFound = FindRecordIndex(strAtom,strEmptyResidue,strEmptyResidueNum,strChain,iRecordNum,prgaf);
                           if (-1 == iFound)
                           {
                              //---------- Search for record (atm;' ';rnum;' ')
                              iFound = FindRecordIndex(strAtom,strEmptyResidue,strResidueNum,strEmptyChain,iRecordNum,prgaf);
                           }                    
                        }
                     }                    
                  }               
               }                    
            }                     
         }
      }                    
   }

   if (-1 != iFound) // found one in the record with specified atm. residue, chain and 
                     // residue #
   { fValue = prgaf->at(iFound).getValue(); return iFound; }
   
   fValue = 0.0; return iFound; // no record found, set the redii to be 0.0 by default 
}                       


//-----------------------------------------------------------------------//
//void CIO::checkCharge(const bool & bSurfCrgInSite)
//{
//   string strLastResidue    = "    "; // a4
//   string strLastResidueNum = "    "; // a4
//   string strAtInf, strResidue, strResidueNum;
//   delphi_real   fCharge, fChargeSum, fError;
//
//   fChargeSum = 0.0;
//   
//   if (bSurfCrgInSite)
//   {
//      for (delphi_integer i = 0; i < iAtomNum; i++)
//      {
//         fCharge = vctapAtomPdb[i].getCharge();
//         
//         strAtInf      = vctapAtomPdb[i].getAtInf();
//         strResidue    = strAtInf.substr(6,4); 
//         strResidueNum = strAtInf.substr(11,4);
//         
//         if ( 0 != strResidue.compare(strLastResidue) || 
//              0 != strResidueNum.compare(strLastResidueNum) )
//         {
//            if (1.0e-4 < abs(fChargeSum))
//            {
//               fError = abs(fChargeSum) - 1.0;
//               
//               if (1.0e-4 < abs(fError))
//               {
//                  CNonZeroNetCrg waring(strLastResidue,strLastResidueNum,fChargeSum);
//               }
//               
//               strLastResidue    = strResidue;
//               strLastResidueNum = strResidueNum;
//               fChargeSum        = fCharge; 
//            }
//            else
//               fChargeSum += fCharge;
//         }    
//      }
//   }
//   
//   if (1.0e-4 < abs(fChargeSum))
//   {
//      fError = abs(fChargeSum) - 1.0;
//      if (1.0e-4 < abs(fError))
//      {
//         CNonZeroNetCrg waring(strLastResidueNum,strLastResidue,fChargeSum);
//      }
//   }   
//}


//-----------------------------------------------------------------------//
void CIO::readForceFile(const string & strFile)
{
   ifstream ifFileStream;
   
   // open the file with name strParamFile
   ifFileStream.open(strFile.c_str());
   
   // if the file doesnt exists, exit 
   if (!ifFileStream.is_open()) throw CUnknownIOFile(strFile);   

   string strLine, strSubLine;

   while (!ifFileStream.eof()) // skip comments and empty lines till the header
   {
      getline(ifFileStream,strLine);      

      strSubLine = removeSpace(strLine); 

      if (0 == strSubLine.compare(0,1,"")) continue; // ignore possible empty lines      

      if ('!' != strLine[0]) break; // header
   } 

   if ( 0 == strLine.find("atom__res_radius") || 0 == strLine.find("atom__resnumbc_radius_") )
      cout << "\n" << "atom radii read from file " << strFile << "\n"; 
   else if ( 0 == strLine.find("atom__res_charge") || 0 == strLine.find("atom__resnumbc_charge_") )
      cout << "\n" << "atomic charges read from file " << strFile << "\n";
   else
      throw CUnknownForceFileHeader(strFile,strLine);

   ifFileStream.clear(); ifFileStream.seekg(0); // rewind

   while (!ifFileStream.eof()) // skip and print comments (start with !) till the header
   {
      getline(ifFileStream,strLine);      

      strSubLine = removeSpace(strLine); 

      if (0 == strSubLine.compare(0,1,"")) continue; // ignore possible empty lines      

      if ('!' != strLine[0]) break; // header
      
      cout << strLine << endl; // print comments
   }

   // determine file format 

   // isPK = false:
   // atom(a6),residue(a3),radius(f8.4): aaaaaarrrfff.ffff
   // isPK = true:
   // atom(a6),residue(a3),residue_number(a4),subunit(a1),radius(f8.4):
   // aaaaaarrrnnnncfff.ffff 
   if (0 == strLine.find("atom__res_radius") || 0 == strLine.find("atom__resnumbc_radius_"))
   {
      if (0 == strLine.find("atom__res_radius")) 
      {
         readFileInNotPKFormat(ifFileStream,SIZEFILE);  
   
#ifdef DEBUG_IO_SIZE
         printForce(SIZEFILE);
#endif   
      }
      else if (0 == strLine.find("atom__resnumbc_radius_")) 
      {
         cout << "reading pK style radius file \n";
         readFileInPKFormat(ifFileStream,SIZEFILE);

#ifdef DEBUG_IO_SIZE
         printForce(SIZEFILE);
#endif
      }
    
      iRadiusNum = prgas.size(); // # of entries in radius file
      cout << "# of radius parameter records: \t\t" << iRadiusNum << "\n";      
   }     
   else // 0 == strLine.find("atom__res_charge") ||  
        // 0 == strLine.find("atom__resnumbc_charge_")
   {
      if (0 == strLine.find("atom__res_charge"))
      {  
         readFileInNotPKFormat(ifFileStream,CHARGEFILE);     

#ifdef DEBUG_IO_CHARGE
         printForce(CHARGEFILE);
#endif
      }
      else if (0 == strLine.find("atom__resnumbc_charge_"))
      {
         cout << "reading pK style charge file \n";
         readFileInPKFormat(ifFileStream,CHARGEFILE);

#ifdef DEBUG_IO_CHARGE
         printForce(CHARGEFILE);
#endif
      }
      
      iCrgNum = prgac.size(); // # of entries in charge file
      cout << "# of charge parameter records: \t\t" << iCrgNum << "\n";         
   }   

   ifFileStream.close();
}   

#ifdef DEBUG_IO_FORCE
//-----------------------------------------------------------------------//
void CIO::printForce(const int & iFileType)
{
   vector<CForce> * prgaf; // pointer to an array fo CForce

   if (SIZEFILE == iFileType)
   {
      prgaf = &prgas;

      cout << endl;
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";
      cout << "+                                                        + \n";
      cout << "+                       ATOM RADII                       + \n";
      cout << "+                                                        + \n";
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";
   }
   else if (CHARGEFILE == iFileType)
   {
      prgaf = &prgac;

      cout << endl;
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";
      cout << "+                                                        + \n";
      cout << "+                       ATOM CHARGE                      + \n";
      cout << "+                                                        + \n";
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";
   }

   for (delphi_integer i = 0; i < iRecordNum; i++)
   {
      cout << setw(6) << i+1 << " : " << setw(6) << prgaf->at(i).getAtom() << " || " << setw(3) << prgaf->at(i).getResidue() << " || " << setw(4) << prgaf->at(i).getResidueNum() << " || "
    	   << setw(1) << prgaf->at(i).getChain() << " || " << setw(8) << prgaf->at(i).getValue() << endl;
   }

   cout << "========================================================== \n";
}
#endif

