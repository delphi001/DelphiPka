/**
 * @file interface_datamarshal.h
 * @brief interface IDataMarshal
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * This file declares the interface IDataMarshal, the base class for CDelphiDataMarshal.
 * Direct realizing an object of this base class is forbidden, but a point to this class can be used to point to
 * an instance of its derived class via polymorphism to provide unified access to various classes.
 */

#ifndef IDATAMARSHAL_H_
#define IDATAMARSHAL_H_

#include <iostream>
#include <fstream>
#include <string> // std::string

#include "environment.h"
#include "interface_exceptions.h"

class IDataMarshal
{
   protected:     
      /**
       * Function to read biomodel and solver names
       * @param strLineNoSpace line of the parameter file with no spaces or comments
       * @return success/false
       */
      bool getBiomodel(const string &strLineNoSpace);
      
      /**
       * Function to determine the presence of a qinclude command and read it
       * @param strLineNoSpace line of the parameter file with no spaces or comments
       * @return success/false
       */
      bool getQinclude(const string &strLineNoSpace);

      /**
       * Function that reads a line, and determines whether the line is a statement or a function.
       * @param strLineNoSpace line of the parameter file without spaces and comments, and then determines if the function
       * is a statement by the presence of an equal sign which is NOT inside brackets.
       * @return success/false
       */
      bool getParameter(const string &strLineNoSpace); 
      
      /**
       * Function to obtain the value(s) of a statement
       * @param strLineNoSpace line of the parameter file with no spaces or comments
       * @return success/false
       */
      virtual bool getStatement(string strLineNoSpace) = 0;
      
      /**
       * Function to obtain the value(s) from a function
       * @param strLineNoSpace line of the parameter file with no spaces or comments
       * @return success/false
       */
      virtual bool getFunction(string strLineNoSpace) = 0;
            
   public:        
      string strParamFile; ///< read-in parameter file (default: fort.10)
      string strBioModel;  ///< Bio-model to solve (default: PBE)
      string strNumSolver; ///< numerical solver (default: DelPhi)

      /**
       * constructor I (for regular delphi runs)
       */
      IDataMarshal()
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                IDataMarshal is constructed                   *\n";
         cout << "****************************************************************\n";
#endif

         strParamFile = "fort.10";
         strBioModel  = "PBE";
         strNumSolver = "DELPHI"; // default
      };

      /**
       * constructor II (for mcce runs)
       */
      IDataMarshal(int argc, char *argv[])
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                IDataMarshal is constructed                   *\n";
         cout << "****************************************************************\n";
#endif

         strParamFile = "fort.10";
         strBioModel  = "PBE";
         strNumSolver = "DELPHI"; // default
         
         if (1 < argc) strParamFile = argv[1]; // input parameter file name 
      };

      /**
       * destructor
       */
      virtual ~IDataMarshal()
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                 IDataMarshal is destroyed                    *\n";
         cout << "****************************************************************\n";
#endif
      };
     
      /**
       * Function to read the parameter file and executes the parsing.
       * @param[in] strFileName Name of the parameter file
       */
      void read(string strFileName);

      /**
       * Function to perform after-reading updates. Must be implemented by the derived classes.
       */
      virtual void updateParameters() = 0;
};

#endif // IDATAMARSHAL_H_
