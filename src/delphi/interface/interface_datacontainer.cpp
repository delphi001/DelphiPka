/*
 * interface_datacontainer.cpp
 *
 *  Created on: Feb 1, 2014
 *      Author: chuan
 */

#include "interface_datacontainer.h"
#include "interface_datacontainer_impls.cpp"

//-----------------------------------------------------------------------//
bool IDataContainer::keyExists(const string& strKey)
{
   DataMap::iterator it = myData.find(strKey);
   
   if (myData.end() == it) return false;
   
   return true; 
}  

//-----------------------------------------------------------------------//
template <class T> const T& IDataContainer::getKey_constRef(const string& strKey)
{  
   if (!keyExists(strKey)) throw CInexistentKey(strKey);
   
   const T * nRetConstPtr = any_cast<const T>(&myData[strKey]);

   const T& nConstRetRef = *nRetConstPtr;
   
   return nConstRetRef;
}

//-----------------------------------------------------------------------//
template <class T> const T * IDataContainer::getKey_constPtr(const string& strKey)
{  
   if (!keyExists(strKey)) throw CInexistentKey(strKey);
   
   const T * nConstRetPtr = any_cast<const T>(&myData[strKey]);
   
   return nConstRetPtr;
}

//-----------------------------------------------------------------------//
template <class T> const T ** IDataContainer::getKey_constPtr(const string& strKey, const int& iRows, const int& iColumns)
{
   if (!keyExists(strKey)) throw CInexistentKey(strKey);

   vector<T> * nConstVectorPtr = any_cast< vector<T> >(&myData[strKey]);

   if (nConstVectorPtr->size() != (unsigned int)iRows*iColumns) throw CMismatchSize(strKey);

   const T * nConstDataPtr = nConstVectorPtr->data();

   const T** prg2D = new const T * [iRows];

   for (int i = 0; i < iRows; i++)
      prg2D[i] = &nConstDataPtr[i*iColumns];

   return prg2D;
}

//-----------------------------------------------------------------------//
template <class T> const T *** IDataContainer::getKey_constPtr(const string& strKey, const int& iRows, const int& iColumns, const int& iPages)
{
   if (!keyExists(strKey)) throw CInexistentKey(strKey);

   vector<T> * nConstVectorPtr = any_cast< vector<T> >(&myData[strKey]);

   if (nConstVectorPtr->size() != (unsigned int)iRows*iColumns*iPages) throw CMismatchSize(strKey);

   T * nConstDataPtr = nConstVectorPtr->data();

   const T *** prg3D = new const T ** [iRows];

   for (int i = 0; i < iRows; i++)
   {
      prg3D[i] = new const T * [iColumns];

      for (int j = 0; j < iColumns; j++)
         prg3D[i][j] = &nConstDataPtr[i*iColumns*iPages+j*iPages];
   }

   return prg3D;
}

//-----------------------------------------------------------------------//
template <class T> T *** IDataContainer::getKey_Ptr(const string& strKey, const int& iRows, const int& iColumns, const int& iPages)
{
   if (!keyExists(strKey)) throw CInexistentKey(strKey);

   vector<T> * nConstVectorPtr = any_cast< vector<T> >(&myData[strKey]);

   if (nConstVectorPtr->size() != (unsigned int)iRows*iColumns*iPages) throw CMismatchSize(strKey);

   T * nConstDataPtr = nConstVectorPtr->data();

   T *** prg3D = new T ** [iRows];

   for (int i = 0; i < iRows; i++)
   {
      prg3D[i] = new T * [iColumns];

      for (int j = 0; j < iColumns; j++)
         prg3D[i][j] = &nConstDataPtr[i*iColumns*iPages+j*iPages];
   }

   return prg3D;
}

//-----------------------------------------------------------------------//
template <class T> const T **** IDataContainer::getKey_constPtr(const string& strKey, const int& iRows, const int& iColumns, const int& iPages, const int& iSects)
{
   if (!keyExists(strKey)) throw CInexistentKey(strKey);

   vector<T> * nConstVectorPtr = any_cast< vector<T> >(&myData[strKey]);

   if (nConstVectorPtr->size() != (unsigned int)iRows*iColumns*iPages*iSects) throw CMismatchSize(strKey);

   const T * nConstDataPtr = nConstVectorPtr->data();

   const T **** prg4D = new const T *** [iRows];

   for (int i = 0; i < iRows; i++)
   {
      prg4D[i] = new const T ** [iColumns];

      for (int j = 0; j < iColumns; j++)
      {
         prg4D[i][j] = new const T * [iPages];

         for (int k = 0; k < iPages; k++)

            prg4D[i][j][k] = &nConstDataPtr[i*iColumns*iPages*iSects+j*iPages*iSects+k*iSects];
      }
   }

   return prg4D;
}

//-----------------------------------------------------------------------//
template <class T> T& IDataContainer::getKey_Ref(const string& strKey)
{  
   if (!keyExists(strKey)) throw CInexistentKey(strKey);
   
   T * nRetPtr = any_cast<T>(&myData[strKey]);

   T& nRetRef = *nRetPtr;
   
   return nRetRef;
}

//-----------------------------------------------------------------------//
template <class T> T IDataContainer::getKey_Val(const string& strKey)
{  
   if (!keyExists(strKey)) throw CInexistentKey(strKey);
   
   T nRetVal = any_cast<T>(myData[strKey]);

   return nRetVal;
}

//-----------------------------------------------------------------------//
template <class T> T * IDataContainer::getKey_Ptr(const string& strKey)
{  
   if (!keyExists(strKey)) throw CInexistentKey(strKey);
   
   T * nRetPtr = any_cast<T>(&myData[strKey]);
     
   return nRetPtr;
}
