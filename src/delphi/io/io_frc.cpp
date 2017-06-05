/*
 * io_frc.cpp
 *
 *  Created on: Feb 4, 2014
 *      Author: chuan
 */

#include "io.h"

//-----------------------------------------------------------------------//
SGrid<delphi_real> CIO::readFrcFile(const string& strFrcFile, const SGrid<delphi_real>& fgOffCenter, const delphi_real& fScale)
{  
   SGrid<delphi_real> fgBoxCenter = {0.0,0.0,0.0};
 
   SGrid<delphi_real> gBoxCenterSum = {0.0,0.0,0.0}, gThisCenter = {0.0,0.0,0.0}; 

   delphi_integer iCenterNum = 0;
 
   ifstream ifFileStream;
 
   // open the file with name strFrcFile
   ifFileStream.open(strFrcFile.c_str());
   
   // if the file doesnt exists, exit 
   if (!ifFileStream.is_open()) throw CUnknownIOFile(strFrcFile);   

   string strLine, strSubLine;

   while (!ifFileStream.eof()) 
   {
      getline(ifFileStream,strLine); strLine = toUpperCase(strLine);
           
      if (0 != strLine.compare(0,6,"ATOM  ") && 0 != strLine.compare(0,6,"HETATM")) continue;
      
      iCenterNum += 1; 

      strSubLine = strLine.substr(30,24);
      
      gThisCenter.nX = atof( strSubLine.substr(0,8).c_str() );
      gThisCenter.nY = atof( strSubLine.substr(8,8).c_str() );
      gThisCenter.nZ = atof( strSubLine.substr(16,8).c_str() );
      
      gBoxCenterSum = gBoxCenterSum + gThisCenter;
      
      if (1.0e-6 > abs(fgOffCenter.nY-999.0) || 1.0e-6 > abs(fgOffCenter.nY-777.0)) break;
   } 

   if (0 < iCenterNum)
      fgBoxCenter = gBoxCenterSum/(delphi_real)iCenterNum;
   else
   {
      if (1.0e-6 > abs(fgOffCenter.nX-999.0)) CEmptyAtomsInFile warning(strFrcFile);
         
      if (1.0e-6 > abs(fgOffCenter.nX-777.0)) CEmptyAtomsInFile warning(strFrcFile);
         
      fgBoxCenter = fgBoxCenter + fgOffCenter/fScale;            
   }      

   ifFileStream.close();

   return fgBoxCenter;
}                      
