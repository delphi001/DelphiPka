/**
 * @file interface_abstractmodule.h
 * @brief interface IAbstractModule
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * This file declares the interface IAbstractModule, the base class for classes such as CDelphiSolver.
 * Direct realizing an object of this base class is forbidden, but a point to this class can be used to point to
 * an instance of its derived class via polymorphism to provide unified access to various classes.
 */

#ifndef IABSTRACTMODULE_H_
#define IABSTRACTMODULE_H_

#include <iostream> 
#include <memory> 

#include "environment.h"
#include "interface_datacontainer.h"

using namespace std;

class IAbstractModule
{
   protected:
      shared_ptr<IDataContainer> pdc;

   public:
      /**
       * constructor
       */
      IAbstractModule(shared_ptr<IDataContainer> pdcIn)
      {
         pdc = pdcIn;

#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*              IAbstractModule is constructed                  *\n";
         cout << "****************************************************************\n";
#endif
      };
      
      /**
       * destructor. must be virtual so that an instance of a derived class can be deleted properly
       * through a pointer to this base class
       */
      virtual ~IAbstractModule()
      {
#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*               IAbstractModule is destroyed                   *\n";
         cout << "****************************************************************\n";
#endif
      };
      
      /**
       * pure abstract function. must be implemented in a derived class in order to take specific actions
       */
      virtual void run() = 0;
      
      /**
       * another pure abstract function allowing to perform additional validation checks for particular variables
       */
      virtual void validateInput() = 0;
};

#endif // IABSTRACTMODULE_H_
