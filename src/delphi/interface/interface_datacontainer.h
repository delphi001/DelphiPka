/**
 * @file interface_datacontainer.h
 * @brief interface IDataContainer
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * This file declares the interface IDataContainer, the base class for classes such as CDelphiData.
 * Direct realizing an object of this base class is forbidden, but a point to this class can be used to point to
 * an instance of its derived class via polymorphism to provide unified access to various data containers.
 */

#ifndef IDATACONTAINER_H_
#define IDATACONTAINER_H_

#include <iostream> 
#include <string.h>      // STL::string 
#include <vector>        // STL::vector
#include <map>           // STL::map
#include <boost/any.hpp> // boost::any

#include "environment.h"
#include "../io/io_datatype.h"
#include "interface_exceptions.h"

using namespace std;
using boost::any_cast;

typedef map<string, boost::any> DataMap;

//-----------------------------------------------------------------------//
class IDataContainer
{
   protected:   

      /**
       * data map containing variables to be read/modified among many other classes. Each entry contains a
       * key (string type) and associated value (boost::any type)
       */
      DataMap myData;
	               	                  
      /**
       * member function to set above data map. must be virtual since only the derived class knows the
       * contents to be written in the data container.
       */
	   virtual void setMap() = 0;

	public:

	   /**
	    * constructor
	    */
	   IDataContainer()
	   {
#ifdef DEBUG_OBJECT
	      cout << endl;
	      cout << "****************************************************************\n";
	      cout << "*               IDataContainer is constructed                  *\n";
	      cout << "****************************************************************\n";
#endif
	   };
	   
      /**
       * virtual destructor allowing an instance of a derived class can be deleted properly
       * through a pointer to this base class
       */
	   virtual ~IDataContainer()
	   {
	      DataMap().swap(myData);

#ifdef DEBUG_OBJECT
         cout << endl;
         cout << "****************************************************************\n";
         cout << "*                IDataContainer is destroyed                   *\n";
         cout << "****************************************************************\n";
#endif
	   };

      /**
       * Function to check if the user-specified key exists or not in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    true if the key exists, false otherwise
       */
	   bool keyExists(const string &strKey);
         
	   //----------read-only key content in the map
	   /**
	    * Template function to get a constant reference to an entry of the data container
	    *
	    * @param[in] strKey The key to be searched in the data container
	    * @return    a constant reference to an entry if the key is found. Otherwise, an exception is thrown.
	    */
	   template <class T> const T& getKey_constRef(const string& strKey);
	   	   
	   /**
	    * Template function to get a constant pointer pointing to the data of a vector-type entry in the data container
	    *
	    * @param[in] strKey The key to be searched in the data container
	    * @return    a constant pointer if the key is found. Otherwise, an exception is thrown.
	    */
	   template <class T> const T* getKey_constPtr(const string& strKey);
	   
	   /**
	    * Template function to get a constant 2D pointer pointing to the data of a vector-type entry in the data container
	    *
	    * @param[in] strKey The key to be searched in the data container
	    * @param[in] iRows Number of rows
	    * @param[in] iColumns Number of columns
	    * @return    a constant 2D pointer if the key is found. Otherwise, an exception is thrown.
	    */
	   template <class T> const T** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns);

      /**
       * Template function to get a constant 3D pointer pointing to the data of a vector-type entry in the data container
       *
	    * @param[in] strKey The key to be searched in the data container
       * @param[in] iRows Number of rows
       * @param[in] iColumns Number of columns
       * @param[in] iPages Number of pages
       * @return    a constant 3D pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> const T*** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages);

      /**
       * Template function to get a constant 4D pointer pointing to the data of a vector-type entry in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @param[in] iRows Number of rows
       * @param[in] iColumns Number of columns
       * @param[in] iPages Number of pages
       * @param[in] iSects Number of sections
       * @return    a constant 4D pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> const T**** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages,const int& iSects);

	   //----------read and write key content in the map
      /**
       * Template function to get a reference to an entry of the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    a reference to an entry if the key is found. Otherwise, an exception is thrown.
       */
      template <class T> T& getKey_Ref(const string& strKey);

      /**
       * Template function to get the value to an entry of the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    Value of an entry if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> T getKey_Val(const string& strKey);
	   
      /**
       * Template function to get a pointer pointing to the data of a vector-type entry in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @return    a pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> T* getKey_Ptr(const string& strKey);

      /**
       * Template function to get a 3D pointer pointing to the data of a vector-type entry in the data container
       *
       * @param[in] strKey The key to be searched in the data container
       * @param[in] iRows Number of rows
       * @param[in] iColumns Number of columns
       * @param[in] iPages Number of pages
       * @return    a 3D pointer if the key is found. Otherwise, an exception is thrown.
       */
	   template <class T> T*** getKey_Ptr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages);

      /**
       * Function to show contents in the data container
       *
       * @param[in] strMapFile The file name that the data container to be written into
       */
	   virtual void showMap(const string& strMapFile) = 0;

      /**
       * Function to reset data container by the values obtained from FORTRAN program
       *
       * @param[in] strF95File The file name storing values obtained from FORTRAN program
       */
	   virtual void reset(const string& strF95File) = 0;
};

#endif // IDATACONTAINER_H_
