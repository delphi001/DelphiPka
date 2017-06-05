/*
 * io_pdb.cpp
 *
 *  Created on: Feb 4, 2014
 *      Author: chuan
 */

#include "io.h"

//-----------------------------------------------------------------------//
void CIO::readStdPdbFile(ifstream& ifPdbFileStream)
{
   int  iObjectType = 0;   // objecttype, 1-sphere, 2-cylinder
   delphi_real fInDielec = 0.0;   // repsintmp
   string strLine,strSubLine,strHeader,strX,strY,strZ;
   delphi_integer iObjectIndex = 0,iAtomIndex = 0;
    
   ifPdbFileStream.clear(); ifPdbFileStream.seekg(0); // rewind   

   CAtomPdb tmpAtomObj;
   
   while (!ifPdbFileStream.eof())
   {   
      getline(ifPdbFileStream,strLine); 
      
      // to ignore possible empty lines
      strSubLine = removeSpace(strLine); 
      if (0 == strSubLine.compare(0,1,"")) continue;
      
      strHeader = strLine.substr(0,6); strHeader = toUpperCase(strHeader);
         
      if (0 == strHeader.compare("MEDIA ")) 
         // # of media has been read during first file reading           
         continue;    
      else if (0 == strHeader.compare("OBJECT"))
      {
         strSubLine = strLine.substr(7,3);  strSubLine = removeSpace(strSubLine);
         iObjectIndex = atoi( strSubLine.c_str() );  
         
         strSubLine = strLine.substr(11,3); strSubLine = removeSpace(strSubLine);
         iObjectType = atoi( strSubLine.c_str() ); 
         
         strSubLine = strLine.substr(15,3); strSubLine = removeSpace(strSubLine);         
         iObjectMediaNum = atoi( strSubLine.c_str() );
          
         strSubLine = strLine.substr(19,8); strSubLine = removeSpace(strSubLine);         
         fInDielec = atof( strSubLine.c_str() );
         
         // entries of vctfMediaEps have been initialized to be 1/epkt
         vctfMediaEps.push_back(fInDielec/fEPKT);
         
         if (0 != iObjectType)
         {
            bOnlyMolecule = false;
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(strLine);
            vctstrObject.push_back(string(" "));
         }
         else
         {
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(string("is a molecule    0"));
            vctstrObject.push_back(string(" "));
            prgiObjectMediaNum.push_back(iObjectMediaNum);
         }
      }
      else if (0 == strHeader.compare("CRGDST"))
         throw CUnsupportedCRGDST(strLine);
      else if (0 == strHeader.compare("ATOM  ") || 0 == strHeader.compare("HETATM")) 
      {
         tmpAtomObj.setAtInf(strLine.substr(11,15));

         strX = strLine.substr(30,8); strX = removeSpace(strX);
         strY = strLine.substr(38,8); strY = removeSpace(strY);
         strZ = strLine.substr(46,8); strZ = removeSpace(strZ);
         tmpAtomObj.setPose(atof( strX.c_str()),atof(strY.c_str()),atof(strZ.c_str()) );
         
         tmpAtomObj.setRadius((delphi_real)0.0);
         tmpAtomObj.setCharge((delphi_real)0.0);
         
         vctapAtomPdb.push_back(tmpAtomObj);
         
         vctiAtomMediaNum.push_back(iObjectMediaNum);
         
         iAtomIndex += 1;       
      }
      else // ignore the rest lines          
         continue;
   } // ---------- end of while (!ifPdbFileStream.eof())
}

//-----------------------------------------------------------------------//
void CIO::readModFile1(ifstream& ifPdbFileStream)
{
   int  iObjectType = 0;   // objecttype, 1-sphere, 2-cylinder
   delphi_real fInDielec = 0.0;   // repsintmp
   string strLine,strSubLine,strHeader,strX,strY,strZ;
   delphi_integer iObjectIndex = 0,iAtomIndex = 0;
     
   ifPdbFileStream.clear(); ifPdbFileStream.seekg(0); // rewind   

   CAtomPdb tmpAtomObj;

   while (!ifPdbFileStream.eof())
   {   
      getline(ifPdbFileStream,strLine); 
      
      // to ignore possible empty lines
      strSubLine = removeSpace(strLine); 
      if (0 == strSubLine.compare(0,1,"")) continue;
      
      strHeader = strLine.substr(0,6); strHeader = toUpperCase(strHeader);
         
      if (0 == strHeader.compare("MEDIA ")) 
         // # of media has been read during first file reading           
         continue;    
      else if (0 == strHeader.compare("OBJECT"))
      {
         strSubLine = strLine.substr(7,3);  strSubLine = removeSpace(strSubLine);
         iObjectIndex = atoi( strSubLine.c_str() );  
         
         strSubLine = strLine.substr(11,3); strSubLine = removeSpace(strSubLine);
         iObjectType = atoi( strSubLine.c_str() ); 
         
         strSubLine = strLine.substr(15,3); strSubLine = removeSpace(strSubLine);         
         iObjectMediaNum = atoi( strSubLine.c_str() );
          
         strSubLine = strLine.substr(19,8); strSubLine = removeSpace(strSubLine);         
         fInDielec = atof( strSubLine.c_str() );
         
         // entries of vctfMediaEps have been initialized to be 1/epkt
         vctfMediaEps.push_back(fInDielec/fEPKT);
         
         if (0 != iObjectType)
         {
            bOnlyMolecule = false;
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(strLine);
            vctstrObject.push_back(string(" "));
         }
         else
         {
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(string("is a molecule    0"));
            vctstrObject.push_back(string(" "));
            prgiObjectMediaNum.push_back(iObjectMediaNum);
         }
      }
      else if (0 == strHeader.compare("CRGDST"))
         throw CUnsupportedCRGDST(strLine);
      else if (0 == strHeader.compare("ATOM  ") || 0 == strHeader.compare("HETATM")) 
      {
         tmpAtomObj.setAtInf(strLine.substr(11,15));

         strX = strLine.substr(30,8); strX = removeSpace(strX);
         strY = strLine.substr(38,8); strY = removeSpace(strY);
         strZ = strLine.substr(46,8); strZ = removeSpace(strZ);
         tmpAtomObj.setPose( atof(strX.c_str()),atof(strY.c_str()),atof(strZ.c_str()) );
         
         strSubLine = strLine.substr(54,6); strSubLine = removeSpace(strSubLine);
         tmpAtomObj.setRadius( (delphi_real)atof(strSubLine.c_str()) );
         
         strSubLine = strLine.substr(60,7); strSubLine = removeSpace(strSubLine);         
         tmpAtomObj.setCharge( (delphi_real)atof(strSubLine.c_str()) );
         
         vctapAtomPdb.push_back(tmpAtomObj);
         
         vctiAtomMediaNum.push_back(iObjectMediaNum);
         
         iAtomIndex += 1;       
      }
      else
         continue;         

   } // ---------- end of while (!ifPdbFileStream.eof())
}

//-----------------------------------------------------------------------//
void CIO::readModFile4(ifstream & ifPdbFileStream)
{
   int  iObjectType = 0;   // objecttype, 1-sphere, 2-cylinder
   delphi_real fInDielec = 0.0;   // repsintmp
   string strLine,strSubLine,strHeader,strX,strY,strZ;
   delphi_integer iObjectIndex = 0,iAtomIndex = 0;
     
   ifPdbFileStream.clear(); ifPdbFileStream.seekg(0); // rewind   

   CAtomPdb tmpAtomObj;

   while (!ifPdbFileStream.eof())
   {   
      getline(ifPdbFileStream,strLine); 
      
      // to ignore possible empty lines
      strSubLine = removeSpace(strLine); 
      if (0 == strSubLine.compare(0,1,"")) continue;
      
      strHeader = strLine.substr(0,6); strHeader = toUpperCase(strHeader);
         
      if (0 == strHeader.compare("MEDIA ")) 
         // # of media has been read during first file reading           
         continue;    
      else if (0 == strHeader.compare("OBJECT"))
      {
         strSubLine = strLine.substr(7,3);  strSubLine = removeSpace(strSubLine);
         iObjectIndex = atoi( strSubLine.c_str() );  
         
         strSubLine = strLine.substr(11,3); strSubLine = removeSpace(strSubLine);
         iObjectType = atoi( strSubLine.c_str() ); 
         
         strSubLine = strLine.substr(15,3); strSubLine = removeSpace(strSubLine);         
         iObjectMediaNum = atoi( strSubLine.c_str() );
          
         strSubLine = strLine.substr(19,8); strSubLine = removeSpace(strSubLine);         
         fInDielec = atof( strSubLine.c_str() );
         
         // entries of vctfMediaEps have been initialized to be 1/epkt
         vctfMediaEps.push_back(fInDielec/fEPKT);
         
         if (0 != iObjectType)
         {
            bOnlyMolecule = false;
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(strLine);
            vctstrObject.push_back(string(" "));
         }
         else
         {
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(string("is a molecule    0"));
            vctstrObject.push_back(string(" "));
            prgiObjectMediaNum.push_back(iObjectMediaNum);
         }
      }
      else if (0 == strHeader.compare("CRGDST"))
         throw CUnsupportedCRGDST(strLine);
      else if (0 == strHeader.compare("ATOM  ") || 0 == strHeader.compare("HETATM")) 
      {
         tmpAtomObj.setAtInf(strLine.substr(11,15));

         strX = strLine.substr(30,8); strX = removeSpace(strX);
         strY = strLine.substr(38,8); strY = removeSpace(strY);
         strZ = strLine.substr(46,8); strZ = removeSpace(strZ);
         tmpAtomObj.setPose( atof(strX.c_str()),atof(strY.c_str()),atof(strZ.c_str()) );
         
         strSubLine = strLine.substr(54,8); strSubLine = removeSpace(strSubLine);
         tmpAtomObj.setCharge( (delphi_real)atof(strSubLine.c_str()) );

         strSubLine = strLine.substr(62,8); strSubLine = removeSpace(strSubLine);
         tmpAtomObj.setRadius( (delphi_real)atof(strSubLine.c_str()) );
         
         vctapAtomPdb.push_back(tmpAtomObj);
         
         vctiAtomMediaNum.push_back(iObjectMediaNum);
         
         iAtomIndex += 1;       
      }
      else
         continue;         

   } // ---------- end of while (!ifPdbFileStream.eof())
}

//-----------------------------------------------------------------------//
void CIO::readPqrFile(ifstream & ifPdbFileStream)
{
   int  iObjectType = 0;   // objecttype, 1-sphere, 2-cylinder
   delphi_real fInDielec = 0.0;   // repsintmp
   string strLine,strSubLine,strHeader,strX,strY,strZ;
   delphi_integer iObjectIndex = 0,iAtomIndex = 0;
     
   ifPdbFileStream.clear(); ifPdbFileStream.seekg(0); // rewind   

   CAtomPdb tmpAtomObj;

   while (!ifPdbFileStream.eof())
   {   
      getline(ifPdbFileStream,strLine); 
      
      // to ignore possible empty lines
      strSubLine = removeSpace(strLine); 
      if (0 == strSubLine.compare(0,1,"")) continue;
      
      strHeader = strLine.substr(0,6); strHeader = toUpperCase(strHeader);
         
      if (0 == strHeader.compare("MEDIA ")) 
         // # of media has been read during first file reading           
         continue;    
      else if (0 == strHeader.compare("OBJECT"))
      {
         strSubLine = strLine.substr(7,3);  strSubLine = removeSpace(strSubLine);
         iObjectIndex = atoi( strSubLine.c_str() );  
         
         strSubLine = strLine.substr(11,3); strSubLine = removeSpace(strSubLine);
         iObjectType = atoi( strSubLine.c_str() ); 
         
         strSubLine = strLine.substr(15,3); strSubLine = removeSpace(strSubLine);         
         iObjectMediaNum = atoi( strSubLine.c_str() );
          
         strSubLine = strLine.substr(19,8); strSubLine = removeSpace(strSubLine);         
         fInDielec = atof( strSubLine.c_str() );
         
         // entries of vctfMediaEps have been initialized to be 1/epkt
         vctfMediaEps.push_back(fInDielec/fEPKT);
         
         if (0 != iObjectType)
         {
            bOnlyMolecule = false;
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(strLine);
            vctstrObject.push_back(string(" "));
         }
         else
         {
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(string("is a molecule    0"));
            vctstrObject.push_back(string(" "));
            prgiObjectMediaNum.push_back(iObjectMediaNum);
         }
      }
      else if (0 == strHeader.compare("CRGDST"))
         throw CUnsupportedCRGDST(strLine);
      else if (0 == strHeader.compare("ATOM  ") || 0 == strHeader.compare("HETATM")) 
      {
         tmpAtomObj.setAtInf(strLine.substr(11,15));        
        
         strX = strLine.substr(30,8); strX = removeSpace(strX);
         strY = strLine.substr(38,8); strY = removeSpace(strY);
         strZ = strLine.substr(46,8); strZ = removeSpace(strZ);
         tmpAtomObj.setPose( atof(strX.c_str()),atof(strY.c_str()),atof(strZ.c_str()) );
         
         strSubLine = strLine.substr(54,7); strSubLine = removeSpace(strSubLine);
         tmpAtomObj.setCharge( (delphi_real)atof(strSubLine.c_str()) );

         strSubLine = strLine.substr(61,7); strSubLine = removeSpace(strSubLine);
         tmpAtomObj.setRadius( (delphi_real)atof(strSubLine.c_str()) );
         
         vctapAtomPdb.push_back(tmpAtomObj);

         vctiAtomMediaNum.push_back(iObjectMediaNum);
         
         iAtomIndex += 1;       
      }
      else
         continue;         

   } // ---------- end of while (!ifPdbFileStream.eof())
}

//-----------------------------------------------------------------------//
void CIO::readMod4File(ifstream& ifPdbFileStream)
{
   int  iObjectType = 0;   // objecttype, 1-sphere, 2-cylinder
   delphi_real fInDielec = 0.0;   // repsintmp
   string strLine,strSubLine,strHeader,strX,strY,strZ;
   delphi_integer iObjectIndex = 0,iAtomIndex = 0;
      
   ifPdbFileStream.clear(); ifPdbFileStream.seekg(0); // rewind   

   CAtomPdb tmpAtomObj;

   while (!ifPdbFileStream.eof())
   {   
      getline(ifPdbFileStream,strLine); 
      
      // to ignore possible empty lines
      strSubLine = removeSpace(strLine); 
      if (0 == strSubLine.compare(0,1,"")) continue;
      
      strHeader = strLine.substr(0,6); strHeader = toUpperCase(strHeader);
         
      if (0 == strHeader.compare("MEDIA ")) 
         // # of media has been read during first file reading           
         continue;    
      else if (0 == strHeader.compare("OBJECT"))
      {
         strSubLine = strLine.substr(7,3);  strSubLine = removeSpace(strSubLine);
         iObjectIndex = atoi( strSubLine.c_str() );  
         
         strSubLine = strLine.substr(11,3); strSubLine = removeSpace(strSubLine);
         iObjectType = atoi( strSubLine.c_str() ); 
         
         strSubLine = strLine.substr(15,3); strSubLine = removeSpace(strSubLine);         
         iObjectMediaNum = atoi( strSubLine.c_str() );
          
         strSubLine = strLine.substr(19,8); strSubLine = removeSpace(strSubLine);         
         fInDielec = atof( strSubLine.c_str() );
         
         // entries of vctfMediaEps have been initialized to be 1/epkt
         vctfMediaEps.push_back(fInDielec/fEPKT);
         
         if (0 != iObjectType)
         {
            bOnlyMolecule = false;
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(strLine);
            vctstrObject.push_back(string(" "));
         }
         else
         {
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(string("is a molecule    0"));
            vctstrObject.push_back(string(" "));
            prgiObjectMediaNum.push_back(iObjectMediaNum);
         }
      }
      else if (0 == strHeader.compare("CRGDST"))
         throw CUnsupportedCRGDST(strLine);
      else if (0 == strHeader.compare("ATOM  ") || 0 == strHeader.compare("HETATM")) 
      {
         tmpAtomObj.setAtInf(strLine.substr(11,15));
         
         strX = strLine.substr(30,8); strX = removeSpace(strX);
         strY = strLine.substr(38,8); strY = removeSpace(strY);
         strZ = strLine.substr(46,8); strZ = removeSpace(strZ);
         tmpAtomObj.setPose( atof(strX.c_str()),atof(strY.c_str()),atof(strZ.c_str()) );
         
         strSubLine = strLine.substr(54,7); strSubLine = removeSpace(strSubLine);
         tmpAtomObj.setRadius( (delphi_real)atof(strSubLine.c_str()) );

         strSubLine = strLine.substr(61,8); strSubLine = removeSpace(strSubLine);
         tmpAtomObj.setCharge( (delphi_real)atof(strSubLine.c_str()) );
        
         vctapAtomPdb.push_back(tmpAtomObj);
        
         vctiAtomMediaNum.push_back(iObjectMediaNum);
         
         iAtomIndex += 1;       
      }
      else
         continue;         

   } // ---------- end of while (!ifPdbFileStream.eof())
}

//-----------------------------------------------------------------------//
void CIO::readPqr4File(ifstream & ifPdbFileStream)
{
   int  iObjectType = 0;   // objecttype, 1-sphere, 2-cylinder
   delphi_real fInDielec = 0.0;   // repsintmp
   string strLine,strSubLine,strHeader,strX,strY,strZ;
   delphi_integer iObjectIndex = 0,iAtomIndex = 0;
      
   ifPdbFileStream.clear(); ifPdbFileStream.seekg(0); // rewind   

   CAtomPdb tmpAtomObj;

   while (!ifPdbFileStream.eof())
   {   
      getline(ifPdbFileStream,strLine); 
      
      // to ignore possible empty lines
      strSubLine = removeSpace(strLine); 
      if (0 == strSubLine.compare(0,1,"")) continue;
      
      strHeader = strLine.substr(0,6); strHeader = toUpperCase(strHeader);
         
      if (0 == strHeader.compare("MEDIA ")) 
         // # of media has been read during first file reading           
         continue;    
      else if (0 == strHeader.compare("OBJECT"))
      {
         strSubLine = strLine.substr(7,3);  strSubLine = removeSpace(strSubLine);
         iObjectIndex = atoi( strSubLine.c_str() );  
         
         strSubLine = strLine.substr(11,3); strSubLine = removeSpace(strSubLine);
         iObjectType = atoi( strSubLine.c_str() ); 
         
         strSubLine = strLine.substr(15,3); strSubLine = removeSpace(strSubLine);         
         iObjectMediaNum = atoi( strSubLine.c_str() );
          
         strSubLine = strLine.substr(19,8); strSubLine = removeSpace(strSubLine);         
         fInDielec = atof( strSubLine.c_str() );
         
         // entries of vctfMediaEps have been initialized to be 1/epkt
         vctfMediaEps.push_back(fInDielec/fEPKT);
         
         if (0 != iObjectType)
         {
            bOnlyMolecule = false;
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(strLine);
            vctstrObject.push_back(string(" "));
         }
         else
         {
            // dataobject(iObjectIndex,1) and dataobject(iObjectIndex,2)
            vctstrObject.push_back(string("is a molecule    0"));
            vctstrObject.push_back(string(" "));
            prgiObjectMediaNum.push_back(iObjectMediaNum);
         }
      }
      else if (0 == strHeader.compare("CRGDST"))
         throw CUnsupportedCRGDST(strLine);
      else if (0 == strHeader.compare("ATOM  ") || 0 == strHeader.compare("HETATM")) 
      {
         tmpAtomObj.setAtInf(strLine.substr(11,15));
         
         strX = strLine.substr(30,8); strX = removeSpace(strX);
         strY = strLine.substr(38,8); strY = removeSpace(strY);
         strZ = strLine.substr(46,8); strZ = removeSpace(strZ);
         tmpAtomObj.setPose( atof(strX.c_str()),atof(strY.c_str()),atof(strZ.c_str()) );
         
         strSubLine = strLine.substr(54,8); strSubLine = removeSpace(strSubLine);
         tmpAtomObj.setCharge( (delphi_real)atof(strSubLine.c_str()) );

         strSubLine = strLine.substr(62,7); strSubLine = removeSpace(strSubLine);
         tmpAtomObj.setRadius( (delphi_real)atof(strSubLine.c_str()) );
         
         vctapAtomPdb.push_back(tmpAtomObj);
         
         vctiAtomMediaNum.push_back(iObjectMediaNum);
         
         iAtomIndex += 1;       
      }
      else
         continue;         

   } // ---------- end of while (!ifPdbFileStream.eof())
}

//----------------------------------------------------------------------//
void CIO::commVecPDB(const vector<string>& strCommPDB)
{
    string strLine,strSubLine,strHeader,strX,strY,strZ;
    delphi_integer iAtomIndex = 0;
    
    
    
    CAtomPdb tmpAtomObj;
    
    for (int i = 0; i < strCommPDB.size(); ++i)
    {
        
        strLine = strCommPDB[i];
        // to ignore possible empty lines
        strSubLine = removeSpace(strLine);
        if (0 == strSubLine.compare(0,1,"")) continue;
        
        strHeader = strLine.substr(0,6); strHeader = toUpperCase(strHeader);
        
        
        if (0 == strHeader.compare("ATOM  ") || 0 == strHeader.compare("HETATM"))
        {
            tmpAtomObj.setAtInf(strLine.substr(11,15));
            
            strX = strLine.substr(30,8); strX = removeSpace(strX);
            strY = strLine.substr(38,8); strY = removeSpace(strY);
            strZ = strLine.substr(46,8); strZ = removeSpace(strZ);
            tmpAtomObj.setPose( atof(strX.c_str()),atof(strY.c_str()),atof(strZ.c_str()) );
            
            strSubLine = strLine.substr(54,8); strSubLine = removeSpace(strSubLine);
            tmpAtomObj.setCharge( (delphi_real)atof(strSubLine.c_str()) );
            
            strSubLine = strLine.substr(62,7); strSubLine = removeSpace(strSubLine);
            tmpAtomObj.setRadius( (delphi_real)atof(strSubLine.c_str()) );
            
            vctapAtomPdb.push_back(tmpAtomObj);
            
            vctiAtomMediaNum.push_back(iObjectMediaNum);
            
            iAtomIndex += 1;       
        }
        else
            continue;
    }
}


//-----------------------------------------------------------------------//
void CIO::readUnformattedPdb(const string& strPdbFile, ifstream& ifPdbFileStream, bool& bPostProcess)
{
   int iUnformatID = 0;  // idfrm  
        
   delphi_real fDummy[5];
   
   char cLine[80]; // line

   ifPdbFileStream.clear(); ifPdbFileStream.seekg(0); // rewind   

   ifPdbFileStream.read(cLine,80); // header
     
   ifPdbFileStream.read( reinterpret_cast<char*>(&iUnformatID), sizeof(int) ); // idfrm
      
   // First reading of file to determine number of atoms
   switch (iUnformatID)
   {
      case 0:  // one form of unformatted file with       
         while (!ifPdbFileStream.eof())
         {
            // one coord, 2 reals
            ifPdbFileStream.read( reinterpret_cast<char*>(fDummy), 5*sizeof(delphi_real) ); 
            iAtomNum += 1;
         }           
                       
         break;            
         
      case 1: // another form of unformatted file
         ifPdbFileStream.read( reinterpret_cast<char*>(&iAtomNum),sizeof(delphi_integer) ); 
   
         break;
      default:
         throw CUnknownUnformattedPdb(strPdbFile); 
   }   
      
   cout << "number of atoms read in = " << iAtomNum << " unformatted" << endl;

   ifPdbFileStream.clear(); ifPdbFileStream.seekg(0); // rewind

   ifPdbFileStream.read(cLine,80); // header
     
   ifPdbFileStream.read( reinterpret_cast<char*>(&iUnformatID), sizeof(int) ); // idfrm

   CAtomPdb tmpAtomObj;
   
   switch (iUnformatID)
   {
      case 0:  // one form of unformatted file with       
         for (delphi_integer i = 0; i < iAtomNum; i++)
         {
            ifPdbFileStream.read( reinterpret_cast<char*>(fDummy), 5*sizeof(delphi_real) ); 
            
            tmpAtomObj.setPose(fDummy[0],fDummy[1],fDummy[2]); 
            tmpAtomObj.setRadius(fDummy[3]);
            tmpAtomObj.setCharge(fDummy[4]);   
            tmpAtomObj.setAtInf(" ");
            vctapAtomPdb.push_back(tmpAtomObj);
         } 

         bPostProcess = true;
                       
         break;            
      case 1: // another form of unformatted file
         ifPdbFileStream.read( reinterpret_cast<char*>(&iAtomNum),sizeof(delphi_integer) ); 

         for (delphi_integer i = 0; i < iAtomNum; i++)
         {
            ifPdbFileStream.read(cLine,80); string strLine(cLine);
            ifPdbFileStream.read( reinterpret_cast<char*>(fDummy), 5*sizeof(delphi_real) ); 

            tmpAtomObj.setPose(fDummy[0],fDummy[1],fDummy[2]); 
            tmpAtomObj.setRadius(fDummy[3]);
            tmpAtomObj.setCharge(fDummy[4]);   
            tmpAtomObj.setAtInf(strLine);
            vctapAtomPdb.push_back(tmpAtomObj);
         }
         
         bPostProcess = false;
         
         break;
   }   
}


//-----------------------------------------------------------------------//

void CIO::readPdbFile(const string& strPdbFile, const int& iPdbFormat, const bool& bPdbUnformatIn,
                      const vector<string>& strCommPDB)
{
    if (CommPDB == iPdbFormat) {
        
        bExistRadiiInfo = true;
        bool bPostProcess = true;
        vctfMediaEps.push_back(0.0); // medeps(0)
        
        commVecPDB(strCommPDB);
        
        iAtomNum = vctapAtomPdb.size();
        
        if (bPostProcess)
        {
            if (0 == vctfMediaEps.size()) vctfMediaEps.push_back(0.0); // medeps(0)
            
            cout << "You are not reading from an object file! \n";
            cout << "Assuming having only molecules, and one medium \n";
            
            vctfMediaEps.push_back(fDielec/fEPKT); // medeps(1)=repsin/epkt
            // dataobject(1,1) and dataobject(1,2)
            vctstrObject.push_back(string("is a molecule    0"));
            vctstrObject.push_back(string(" "));
            prgiObjectMediaNum.push_back(iObjectMediaNum);
        }
        
        iObjectNum = vctstrObject.size()/2;
        iMediaNum  = vctfMediaEps.size()-1;
        
        // regardless of being a molecule or not, iatmmed has a field to say which is its medium
        if (0 != prgiObjectMediaNum.size())
        {
            for (unsigned i=0; i < prgiObjectMediaNum.size(); i++)
                vctiAtomMediaNum.push_back(prgiObjectMediaNum[i]);
            
            prgiObjectMediaNum.clear();
        }
    }
    
    else
    {
      bExistRadiiInfo = false; // iatrad, whether there is radius info
       
      ifstream ifPdbFileStream;
   
      ifPdbFileStream.open(strPdbFile.c_str());
           
      if (!ifPdbFileStream.is_open()) throw CUnknownIOFile(strPdbFile);
        
      bool bPostProcess = true;
      
      int iFormatID = 0; // idfrm 
   
      if(!bPdbUnformatIn) // ---------- formatted pdb file ---------- //
      {
         string strLine, strSubLine, strHeader;  
     
         cout << "\nopening formatted file: " << strPdbFile << endl;
              
         // Now only delphi mod type pdb files have keyword "DELPHI PDB", other pdb files
         // such as standard, mod4 and pqr pdb files do not check the keyword any more 
         if (MODPDB == iPdbFormat) // mod pdb file
         {
            getline(ifPdbFileStream,strLine); // read the 1st line
            string strHeadFirst = strLine.substr(0,6); 
            strHeadFirst = removeSpace(strHeadFirst);
            strHeadFirst = toUpperCase(strHeadFirst);
            if (0 == strHeadFirst.compare(0,5,"MEDIA")) bPostProcess = false;
            
            if (string::npos != strLine.find("DELPHI")) // find keyword "DELPHI"
            {      
               if (string::npos == strLine.find("PDB")) // no keyword "PDB" is found
                  throw CUnknownDelphiPdb(strPdbFile);        
         
               cout << "Reading a Delphi-style pdb file \n";      
               
               getline(ifPdbFileStream,strLine); // read the 2nd line            
   
               size_t found = strLine.find_first_of('='); 
               strSubLine = strLine.substr(found+1,strLine.size()-found-1);
               strSubLine = removeSpace(strSubLine);
               iFormatID = atoi(strSubLine.c_str());
               
               cout << "Read Delphi Format Number = " << iFormatID << endl;  
               
               if (1 != iFormatID && 4 != iFormatID) 
                  throw CUnknownModPdb(strPdbFile,iFormatID);
            }
         }
         else
         {
   
            ifPdbFileStream.clear(); ifPdbFileStream.seekg(0); // rewind   
            cout << "No DELPHI keyword, assuming Delphi Format number = " << iFormatID << endl;
         
         }
   
          
          
         if (STDPDB == iPdbFormat) bExistRadiiInfo = false; // standard pdb
         else                      bExistRadiiInfo = true;  // others                  
   
         // medeps(0:nmediamax)
         vctfMediaEps.push_back(0.0); // medeps(0)
   
         switch (iPdbFormat)
         {
            case STDPDB: // standard pdb
               readStdPdbFile(ifPdbFileStream);
               break;
            case MODPDB: // mod pdb    
               if (1 == iFormatID) readModFile1(ifPdbFileStream);
               else readModFile4(ifPdbFileStream); // 4 == iFormatID    
               break;
            case PQRPDB: // pqr pdb
               readPqrFile(ifPdbFileStream);   
               break;      
            case MOD4PDB: // mod4 pdb
               readMod4File(ifPdbFileStream);  
               break;       
            case PQR4PDB: // pqr4 pdb
               readPqr4File(ifPdbFileStream);
               break;
         }
   
         iAtomNum = vctapAtomPdb.size();
         cout << "number of atoms read in = " << iAtomNum << " formatted \n";
      }   
      else // ---------- unformatted pdb file ---------- //
      {
         int iUnformatID = 0;  // idfrm  
             
         char cLine[80]; // line
         
         bExistRadiiInfo = true;
   
         ifPdbFileStream.read(cLine,80); // header
      
         string strHeader(cLine); // string to be converted from cLine for simpler manipulation
      
         if (string::npos != strHeader.find("DELPHI")) // find keyword "DELPHI"
         {      
            if (string::npos == strHeader.find("PDB")) // no keyword "PDB" is found
               throw CUnknownDelphiPdb(strPdbFile);        
         
            cout << "Reading a Delphi-style pdb file \n";  
         }
      
         ifPdbFileStream.read( reinterpret_cast<char*>(&iUnformatID), sizeof(int) ); // idfrm
         
         if (0 != iUnformatID && 1 != iUnformatID) throw CUnknownUnformattedPdb(strPdbFile); 
         
         readUnformattedPdb(strPdbFile,ifPdbFileStream,bPostProcess); 
      }
           
   #ifdef DEBUG_IO_PDB
      printPDB();
   #endif
   
      // Some post-processing
      if (bPostProcess) 
      {
         if (0 == vctfMediaEps.size()) vctfMediaEps.push_back(0.0); // medeps(0)
   
         cout << "You are not reading from an object file! \n";
         cout << "Assuming having only molecules, and one medium \n";
      
         vctfMediaEps.push_back(fDielec/fEPKT); // medeps(1)=repsin/epkt
         // dataobject(1,1) and dataobject(1,2)
         vctstrObject.push_back(string("is a molecule    0"));
         vctstrObject.push_back(string(" "));
         prgiObjectMediaNum.push_back(iObjectMediaNum);
      }
   
      iObjectNum = vctstrObject.size()/2;
      iMediaNum  = vctfMediaEps.size()-1;
   
      // regardless of being a molecule or not, iatmmed has a field to say which is its medium
      if (0 != prgiObjectMediaNum.size())
      {
         for (unsigned i=0; i < prgiObjectMediaNum.size(); i++)
            vctiAtomMediaNum.push_back(prgiObjectMediaNum[i]);
   
         prgiObjectMediaNum.clear();
      }
      
      ifPdbFileStream.close();
    }
   
}

//-----------------------------------------------------------------------//
void CIO::writeUnformatPdb(const string& strPdbFile)
{
   SGrid<delphi_real> gXYZ;
   delphi_real        fRadius,fCharge;

   ofstream ofPdbFileStream(strPdbFile.c_str(),ios::binary);
   
   for (delphi_integer iThisAtom = 0; iThisAtom < iAtomNum; iThisAtom++)
   {   
      gXYZ    = vctapAtomPdb[iThisAtom].getPose();
      fRadius = vctapAtomPdb[iThisAtom].getRadius();
      fCharge = vctapAtomPdb[iThisAtom].getCharge();

      ofPdbFileStream.write(reinterpret_cast<char*>(&gXYZ),    sizeof(gXYZ));
      ofPdbFileStream.write(reinterpret_cast<char*>(&fRadius), sizeof(fRadius));
      ofPdbFileStream.write(reinterpret_cast<char*>(&fCharge), sizeof(fCharge));  
   }

   ofPdbFileStream.close();
}

//-----------------------------------------------------------------------//
void CIO::writeModifiedPdb(const string& strPdbFile, const int& iModPdbFormatOut)
{  
   cout << "atomic coordinates, charges and radii written to file " << strPdbFile << endl << endl;
        
   if (0 == iAtomNum && 0 == iObjectNum) { CNoAtomObjectData warning; return;}

   string strLine = "ATOM  ";

   ofstream ofPdbFileStream(strPdbFile.c_str());        
   
   switch (iModPdbFormatOut)
   {
      case MODPDB: // mod pdb w/ FORMAT = 1   
         ofPdbFileStream << "DELPHI PDB FILE" << endl;  
         ofPdbFileStream << "FORMAT = 1" << "   " << iModPdbFormatOut << endl;
         ofPdbFileStream << "HEADER output from qdiff" << endl;
         ofPdbFileStream << "HEADER atom radii in columns 55-60" << endl;
         ofPdbFileStream << "HEADER atom charges in columns 61-67" << endl;
         
         for (delphi_integer iThisAtom = 0; iThisAtom < iAtomNum; iThisAtom++)
         {         
            ofPdbFileStream << setw(6)  << left  << "ATOM  ";                            // line(1:6)
            ofPdbFileStream << setw(5)  << right << iThisAtom+1;                         // line(7:11)
            ofPdbFileStream << setw(15) << left << vctapAtomPdb[iThisAtom].getAtInf();   // line(12:26)
            ofPdbFileStream << "    ";                                                   // line(27:30)

            SGrid<delphi_real> gXYZ = vctapAtomPdb[iThisAtom].getPose();
            ofPdbFileStream << fixed << setprecision(3);
            ofPdbFileStream << setw(8)  << right << gXYZ.nX;                             // line(31:38) f8.3
            ofPdbFileStream << setw(8)  << right << gXYZ.nY;                             // line(39:46) f8.3
            ofPdbFileStream << setw(8)  << right << gXYZ.nZ;                             // line(47:54) f8.3
            
            ofPdbFileStream << fixed << setprecision(2);
            ofPdbFileStream << setw(6)  << right << vctapAtomPdb[iThisAtom].getRadius(); // line(55:60) f6.2

            ofPdbFileStream << fixed << setprecision(3);
            ofPdbFileStream << setw(7)  << right << vctapAtomPdb[iThisAtom].getCharge(); // line(61:67) f7.3
            ofPdbFileStream << "\n";                
         }
         
         break;

      case PQRPDB: // pqr pdb
         ofPdbFileStream << "DELPHI PDB FILE" << endl;  
         ofPdbFileStream << "FORMAT = PQR" << endl;
         ofPdbFileStream << "HEADER output from qdiff" << endl;
         ofPdbFileStream << "HEADER atom charges in columns 55-61" << endl;
         ofPdbFileStream << "HEADER atom radii   in columns 62-68" << endl;
         
         for (delphi_integer iThisAtom = 0; iThisAtom < iAtomNum; iThisAtom++)
         {         
            ofPdbFileStream << setw(6)  << left  << "ATOM  ";                            // line(1:6)
            
            ofPdbFileStream << setw(5)  << right << iThisAtom+1;                         // line(7:11)
            
            ofPdbFileStream << setw(19) << left  << vctapAtomPdb[iThisAtom].getAtInf();  // line(12:30)
                        
            SGrid<delphi_real> gXYZ = vctapAtomPdb[iThisAtom].getPose();
            ofPdbFileStream << fixed << setprecision(3);
            ofPdbFileStream << setw(8)  << right << gXYZ.nX;                             // line(31:38) f8.3
            ofPdbFileStream << setw(8)  << right << gXYZ.nY;                             // line(39:46) f8.3
            ofPdbFileStream << setw(8)  << right << gXYZ.nZ;                             // line(47:54) f8.3
            
            ofPdbFileStream << fixed << setprecision(3);
            ofPdbFileStream << setw(7)  << right << vctapAtomPdb[iThisAtom].getCharge(); // line(55:61) f7.3
            ofPdbFileStream << setw(7)  << right << vctapAtomPdb[iThisAtom].getRadius(); // line(62:68) f7.3
            ofPdbFileStream << "\n";                
         }

         break;      

      case MOD4PDB: // mod4 pdb
         ofPdbFileStream << "DELPHI PDB FILE" << endl;  
         ofPdbFileStream << "FORMAT = 1" << "   " << iModPdbFormatOut << endl;
         ofPdbFileStream << "4 digits precison" << endl;
         ofPdbFileStream << "HEADER output from qdiff" << endl;
         ofPdbFileStream << "HEADER atom radii in columns 55-61" << endl;
         ofPdbFileStream << "HEADER atom charges in columns 62-69" << endl;
         
         for (delphi_integer iThisAtom = 0; iThisAtom < iAtomNum; iThisAtom++)
         {         
            ofPdbFileStream << setw(6)  << left  << "ATOM  ";                            // line(1:6)
            
            ofPdbFileStream << setw(5)  << right << iThisAtom+1;                         // line(7:11)
            
            ofPdbFileStream << setw(19) << left  << vctapAtomPdb[iThisAtom].getAtInf();  // line(12:30)
                        
            SGrid<delphi_real> gXYZ = vctapAtomPdb[iThisAtom].getPose();
            ofPdbFileStream << fixed << setprecision(3);
            ofPdbFileStream << setw(8)  << right << gXYZ.nX;                             // line(31:38) f8.3
            ofPdbFileStream << setw(8)  << right << gXYZ.nY;                             // line(39:46) f8.3
            ofPdbFileStream << setw(8)  << right << gXYZ.nZ;                             // line(47:54) f8.3
            
            ofPdbFileStream << fixed << setprecision(4);
            ofPdbFileStream << setw(7)  << right << vctapAtomPdb[iThisAtom].getRadius(); // line(55:61) f7.4
            ofPdbFileStream << setw(8)  << right << vctapAtomPdb[iThisAtom].getCharge(); // line(62:69) f8.4
            ofPdbFileStream << "\n";                
         }

         break;       

      case PQR4PDB: // pqr4 pdb
         ofPdbFileStream << "DELPHI PDB FILE" << endl;  
         ofPdbFileStream << "FORMAT = PQR" << endl;
         ofPdbFileStream << "4 digits precision" << endl;
         ofPdbFileStream << "HEADER output from qdiff" << endl;
         ofPdbFileStream << "HEADER atom charges in columns 55-62" << endl;
         ofPdbFileStream << "HEADER atom radii   in columns 63-69" << endl;
         
         for (delphi_integer iThisAtom = 0; iThisAtom < iAtomNum; iThisAtom++)
         {         
            ofPdbFileStream << setw(6)  << left  << "ATOM  ";                               // line(1:6)
            
            ofPdbFileStream << setw(5)  << right << iThisAtom+1;                            // line(7:11)
            
            ofPdbFileStream << setw(19) << left  << vctapAtomPdb[iThisAtom].getAtInf();     // line(12:30)
                        
            ofPdbFileStream << fixed << setprecision(3);
            ofPdbFileStream << setw(8)  << right << (vctapAtomPdb[iThisAtom].getPose()).nX; // line(31:38) f8.3
            ofPdbFileStream << setw(8)  << right << (vctapAtomPdb[iThisAtom].getPose()).nY; // line(39:46) f8.3
            ofPdbFileStream << setw(8)  << right << (vctapAtomPdb[iThisAtom].getPose()).nZ; // line(47:54) f8.3
            
            ofPdbFileStream << fixed << setprecision(4);
            ofPdbFileStream << setw(8)  << right << vctapAtomPdb[iThisAtom].getCharge();    // line(55:62) f7.4
            ofPdbFileStream << setw(7)  << right << vctapAtomPdb[iThisAtom].getRadius();    // line(63:69) f8.4
            ofPdbFileStream << "\n";                
         }

         break;
   }   
   
   ofPdbFileStream.close(); 
}


#ifdef DEBUG_IO_PDB
//-----------------------------------------------------------------------//
void CIO::printPDB()
{
   cout << endl;
   cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";
   cout << "+                                                        + \n";
   cout << "+                        ATOM PDB                        + \n";
   cout << "+                                                        + \n";
   cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n";

   for (delphi_integer i = 0; i < iAtomNum; i++)
   {
      cout << setw(4) << i+1 << " : ";

      cout << setw(8)  << left << (vctapAtomPdb[i].getPose()).nX << "  "   << setw(8) << left << (vctapAtomPdb[i].getPose()).nY << "  "   << setw(8)  << left << (vctapAtomPdb[i].getPose()).nZ << " || ";
           << setw(8)  << left << vctapAtomPdb[i].getRadius()    << " || " << setw(8) << left << vctapAtomPdb[i].getCharge()    << " || " << setw(15) << left << vctapAtomPdb[i].getAtInf()     << endl;
   }

   cout << "========================================================== \n";
}
#endif




