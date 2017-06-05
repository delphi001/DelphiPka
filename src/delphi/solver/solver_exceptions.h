#ifndef SOLVER_EXCEPTIONS_H_
#define SOLVER_EXCEPTIONS_H_

#include "../interface/exceptions.h"

class CUnknownBioModelSolver : public CException 
{
   public:
      CUnknownBioModelSolver(const string& strBioModel, const string& strNumSolver)
      {
         cerr << "BIOMODEL = " << strBioModel << " AND SOLVER = " << strNumSolver << " CANNOT BE PARIED.\n";
         cerr << "(CHECK THE PARAMETER FILE AND TRY AGAIN... \n";
      }
};

class CDirectEpsilonMap : public CException
{
   public:
   CDirectEpsilonMap(const int& iDirectEpsMap)
      {
         cerr << "idirectalg = " << iDirectEpsMap << " IN SETCRG THIS PART NEEDS TO BE CHANGED.\n";
         cerr << "(CHECK THE PARAMETER FILE AND TRY AGAIN... \n";
      }
};

class CUnknownBndyCondition : public CException
{
   public:
      CUnknownBndyCondition(const int& iBndyType)
      {
         cerr << "UNRECOGNIZED BOUNDARY CONDITION BNDCON = " << iBndyType << " (CHECK THE PARAMETER FILE AND TRY AGAIN)\n";
      }
};

class CSettingBndyError : public CException
{
   public:
      CSettingBndyError(const int& iBndyType)
      {
         cerr << "PROBLEM OCCURS WHEN SETTING BOUNDARY CONDITION BNDCON = " << iBndyType << " (PROGRAM EXITS...)\n";
      }
};

class CUnknownPhiiFile : public CException
{
   public:
      CUnknownPhiiFile(const string& strPhiiFile)
      {
         cerr << "NP POTL MAP FOR FOCUSSING BOUNDARY CONDITIONS (FILE " << strPhiiFile << " DOES NOT EXISITS)\n";
      }
};

class COutsideFocus : public CException
{
   public:
      COutsideFocus(const delphi_real& fScale1,const SGrid<delphi_real>& fgBoxCenter1,const delphi_real& fScale,const SGrid<delphi_real>& fgBoxCenter)
      {
         cerr << "PART OF NEW GRID LIES OUTSIDE OLD GRID\n";
         cerr << "CHECK SCALING OF BOTH GRIDS\n";
         cerr << "OLD SCALE:\n";
         cerr << "scale (grids/A)      : " << fScale1 << endl;
         cerr << "object centre at (A) : " << fgBoxCenter1.nX << " " << fgBoxCenter1.nY << " " << fgBoxCenter1.nZ << endl;
         cerr << "NEW SCALE:\n";
         cerr << "scale (grids/A)      : " << fScale << endl;
         cerr << "object centre at (A) : " << fgBoxCenter.nX << " " << fgBoxCenter.nY << " " <<fgBoxCenter.nZ << endl;
      }
};

class CZeorLinIterateNum : public CException
{
   public:
   CZeorLinIterateNum(const bool& bAutoConverge, const int& iLinIterateNum)
      {
         cerr << "AUTOC = " << bAutoConverge << " AND LINIT = " << iLinIterateNum << " ";
         cerr << "(PROGRAM DOES NOT KNOW WHAT TO DO AND THEREFORE EXITS...) \n";
      }
};

class CCrgedPtsInSolution : public CWarning
{
   public:
   CCrgedPtsInSolution(const delphi_integer& ico)
      {
         cerr << "OUT OF THEM, " << ico << " CHARGED POINT(S) ARE IN SOLUTION ";
         cerr << "(INCREASED RESOLUTION IS NEEDED) \n";
      }
};

class CSmallRelaxParam : public CWarning
{
   public:
   CSmallRelaxParam()
      {
         cerr << "CONVERGENCE IS MORE RELIABLE IF RELAXATION PARAMETER > 0.05 ";
         cerr << "(IT IS ADVISABLE TO INCREASE IT) \n";
      }
};


#endif /* SOLVER_EXCEPTIONS_H_ */
