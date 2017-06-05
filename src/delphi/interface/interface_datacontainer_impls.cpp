/**
 * @file interface_datacontainer_impls.h
 * @brief concrete implementations of template member functions in IDataContainer
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * DO NOT COMPILE THIS FILE SEPERATELY! ONLY NEED THIS FILE TO BE INCLUDED IN interface_datacontainer.cpp!!
 */

//----------template <class T> const T& getKey_constRef(const string& strKey)
template const int& IDataContainer::getKey_constRef(const string&);
template const long int& IDataContainer::getKey_constRef(const string&);
template const float& IDataContainer::getKey_constRef(const string&);
template const double& IDataContainer::getKey_constRef(const string&);
template const long double& IDataContainer::getKey_constRef(const string&);
template const bool& IDataContainer::getKey_constRef(const string&);
template const string& IDataContainer::getKey_constRef(const string&);
template const SGrid<int>& IDataContainer::getKey_constRef(const string&);
template const SGrid<long int>& IDataContainer::getKey_constRef(const string&);
template const SGrid<float>& IDataContainer::getKey_constRef(const string&);
template const SGrid<double>& IDataContainer::getKey_constRef(const string&);
template const SGrid<long double>& IDataContainer::getKey_constRef(const string&);
template const SGridValue<int>& IDataContainer::getKey_constRef(const string&);
template const SGridValue<long int>& IDataContainer::getKey_constRef(const string&);
template const SGridValue<float>& IDataContainer::getKey_constRef(const string&);
template const SGridValue<double>& IDataContainer::getKey_constRef(const string&);
template const SGridValue<long double>& IDataContainer::getKey_constRef(const string&);
template const SExtrema<delphi_integer>& IDataContainer::getKey_constRef(const string&);
template const SExtrema<delphi_real>& IDataContainer::getKey_constRef(const string&);
template const CForce& IDataContainer::getKey_constRef(const string&);
template const CAtomPdb& IDataContainer::getKey_constRef(const string&);
template const vector<delphi_real>& IDataContainer::getKey_constRef(const string&);
template const vector<int>& IDataContainer::getKey_constRef(const string&);
template const vector<bool>& IDataContainer::getKey_constRef(const string&);
template const vector< SExtrema<delphi_real> >& IDataContainer::getKey_constRef(const string&);
template const vector< SGrid<delphi_integer> >& IDataContainer::getKey_constRef(const string&);
template const vector< SGrid<delphi_real> >& IDataContainer::getKey_constRef(const string&);
template const vector< SGridValue<delphi_real> >& IDataContainer::getKey_constRef(const string&);
template const vector< SGridValue<delphi_integer> >& IDataContainer::getKey_constRef(const string&);
template const vector< CAtomPdb >& IDataContainer::getKey_constRef(const string&);
template const vector< string >& IDataContainer::getKey_constRef(const string&);
template const vector< SDoubleGridValue >& IDataContainer::getKey_constRef(const string&);

//----------template <class T> const T* getKey_constPtr(const string& strKey)
template const int * IDataContainer::getKey_constPtr(const string&);
template const long int * IDataContainer::getKey_constPtr(const string&);
template const float * IDataContainer::getKey_constPtr(const string&);
template const double * IDataContainer::getKey_constPtr(const string&);
template const long double * IDataContainer::getKey_constPtr(const string&);
template const bool * IDataContainer::getKey_constPtr(const string&);
template const string * IDataContainer::getKey_constPtr(const string&);
template const SGrid<int> * IDataContainer::getKey_constPtr(const string&);
template const SGrid<long int> * IDataContainer::getKey_constPtr(const string&);
template const SGrid<float> * IDataContainer::getKey_constPtr(const string&);
template const SGrid<double> * IDataContainer::getKey_constPtr(const string&);
template const SGrid<long double> * IDataContainer::getKey_constPtr(const string&);
template const SGridValue<int> * IDataContainer::getKey_constPtr(const string&);
template const SGridValue<long int> * IDataContainer::getKey_constPtr(const string&);
template const SGridValue<float> * IDataContainer::getKey_constPtr(const string&);
template const SGridValue<double> * IDataContainer::getKey_constPtr(const string&);
template const SGridValue<long double> * IDataContainer::getKey_constPtr(const string&);
template const SExtrema<delphi_integer> * IDataContainer::getKey_constPtr(const string&);
template const SExtrema<delphi_real> * IDataContainer::getKey_constPtr(const string&);
template const CForce * IDataContainer::getKey_constPtr(const string&);
template const CAtomPdb * IDataContainer::getKey_constPtr(const string&);
template const vector<delphi_real>* IDataContainer::getKey_constPtr(const string&);
template const vector<int>* IDataContainer::getKey_constPtr(const string&);
template const vector<bool>* IDataContainer::getKey_constPtr(const string&);
template const vector< SExtrema<delphi_real> >* IDataContainer::getKey_constPtr(const string&);
template const vector< SGrid<delphi_integer> >* IDataContainer::getKey_constPtr(const string&);
template const vector< SGrid<delphi_real> >* IDataContainer::getKey_constPtr(const string&);
template const vector< SGridValue<delphi_real> >* IDataContainer::getKey_constPtr(const string&);
template const vector< SGridValue<delphi_integer> >* IDataContainer::getKey_constPtr(const string&);
template const vector< CAtomPdb >* IDataContainer::getKey_constPtr(const string&);
template const vector< string >* IDataContainer::getKey_constPtr(const string&);
template const vector< SDoubleGridValue >* IDataContainer::getKey_constPtr(const string&);

//----------template <class T> const T** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns)
template const SGrid<delphi_integer> ** IDataContainer::getKey_constPtr(const string&,const int&,const int&);

//----------template <class T> const T*** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages)
template const SGrid<delphi_integer> *** IDataContainer::getKey_constPtr(const string&,const int&,const int&,const int&);

//----------template <class T> const T**** getKey_constPtr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages,const int& iSects)
template const SGrid<delphi_integer> **** IDataContainer::getKey_constPtr(const string&,const int&,const int&,const int&,const int&);
template const delphi_real *** IDataContainer::getKey_constPtr(const string&,const int&,const int&,const int&);

//----------template <class T> T& getKey_Ref(const string& strKey)
template int& IDataContainer::getKey_Ref(const string&);
template long int& IDataContainer::getKey_Ref(const string&);
template float& IDataContainer::getKey_Ref(const string&);
template double& IDataContainer::getKey_Ref(const string&);
template long double& IDataContainer::getKey_Ref(const string&);
template bool& IDataContainer::getKey_Ref(const string&);
template string& IDataContainer::getKey_Ref(const string&);
template SGrid<int>& IDataContainer::getKey_Ref(const string&);
template SGrid<long int>& IDataContainer::getKey_Ref(const string&);
template SGrid<float>& IDataContainer::getKey_Ref(const string&);
template SGrid<double>& IDataContainer::getKey_Ref(const string&);
template SGrid<long double>& IDataContainer::getKey_Ref(const string&);
template SGridValue<int>& IDataContainer::getKey_Ref(const string&);
template SGridValue<long int>& IDataContainer::getKey_Ref(const string&);
template SGridValue<float>& IDataContainer::getKey_Ref(const string&);
template SGridValue<double>& IDataContainer::getKey_Ref(const string&);
template SGridValue<long double>& IDataContainer::getKey_Ref(const string&);
template SExtrema<delphi_integer>& IDataContainer::getKey_Ref(const string&);
template SExtrema<delphi_real>& IDataContainer::getKey_Ref(const string&);
template CForce& IDataContainer::getKey_Ref(const string&);
template CAtomPdb& IDataContainer::getKey_Ref(const string&);
template vector<delphi_real>& IDataContainer::getKey_Ref(const string&);
template vector<int>& IDataContainer::getKey_Ref(const string&);
template vector<bool>& IDataContainer::getKey_Ref(const string&);
template vector< SExtrema<delphi_real> >& IDataContainer::getKey_Ref(const string&);
template vector< SGrid<delphi_integer> >& IDataContainer::getKey_Ref(const string&);
template vector< SGrid<delphi_real> >& IDataContainer::getKey_Ref(const string&);
template vector< SGridValue<delphi_real> >& IDataContainer::getKey_Ref(const string&);
template vector< SGridValue<delphi_integer> >& IDataContainer::getKey_Ref(const string&);
template vector< CAtomPdb >& IDataContainer::getKey_Ref(const string&);
template vector< string >& IDataContainer::getKey_Ref(const string&);
template vector< SDoubleGridValue >& IDataContainer::getKey_Ref(const string&);

//----------template <class T> T getKey_Val(const string& strKey)
template int IDataContainer::getKey_Val(const string&);
template long int IDataContainer::getKey_Val(const string&);
template float IDataContainer::getKey_Val(const string&);
template double IDataContainer::getKey_Val(const string&);
template long double IDataContainer::getKey_Val(const string&);
template bool IDataContainer::getKey_Val(const string&);
template string IDataContainer::getKey_Val(const string&);
template SGrid<int> IDataContainer::getKey_Val(const string&);
template SGrid<long int> IDataContainer::getKey_Val(const string&);
template SGrid<float> IDataContainer::getKey_Val(const string&);
template SGrid<double> IDataContainer::getKey_Val(const string&);
template SGrid<long double> IDataContainer::getKey_Val(const string&);
template SGridValue<int> IDataContainer::getKey_Val(const string&);
template SGridValue<long int> IDataContainer::getKey_Val(const string&);
template SGridValue<float> IDataContainer::getKey_Val(const string&);
template SGridValue<double> IDataContainer::getKey_Val(const string&);
template SGridValue<long double> IDataContainer::getKey_Val(const string&);
template SExtrema<delphi_integer> IDataContainer::getKey_Val(const string&);
template SExtrema<delphi_real> IDataContainer::getKey_Val(const string&);
template CForce IDataContainer::getKey_Val(const string&);
template CAtomPdb IDataContainer::getKey_Val(const string&);
template vector<delphi_real> IDataContainer::getKey_Val(const string&);
template vector<int> IDataContainer::getKey_Val(const string&);
template vector<bool> IDataContainer::getKey_Val(const string&);
template vector< SExtrema<delphi_real> > IDataContainer::getKey_Val(const string&);
template vector< SGrid<delphi_integer> > IDataContainer::getKey_Val(const string&);
template vector< SGrid<delphi_real> > IDataContainer::getKey_Val(const string&);
template vector< SGridValue<delphi_real> > IDataContainer::getKey_Val(const string&);
template vector< SGridValue<delphi_integer> > IDataContainer::getKey_Val(const string&);
template vector< CAtomPdb > IDataContainer::getKey_Val(const string&);
template vector< string > IDataContainer::getKey_Val(const string&);
template vector< SDoubleGridValue > IDataContainer::getKey_Val(const string&);

//----------template <class T> T * getKey_Ptr(const string& strKey)
template int * IDataContainer::getKey_Ptr(const string&);
template long int * IDataContainer::getKey_Ptr(const string&);
template float * IDataContainer::getKey_Ptr(const string&);
template double * IDataContainer::getKey_Ptr(const string&);
template long double * IDataContainer::getKey_Ptr(const string&);
template bool * IDataContainer::getKey_Ptr(const string&);
template string * IDataContainer::getKey_Ptr(const string&);
template SGrid<int> * IDataContainer::getKey_Ptr(const string&);
template SGrid<long int> * IDataContainer::getKey_Ptr(const string&);
template SGrid<float> * IDataContainer::getKey_Ptr(const string&);
template SGrid<double> * IDataContainer::getKey_Ptr(const string&);
template SGrid<long double> * IDataContainer::getKey_Ptr(const string&);
template SGridValue<int> * IDataContainer::getKey_Ptr(const string&);
template SGridValue<long int> * IDataContainer::getKey_Ptr(const string&);
template SGridValue<float> * IDataContainer::getKey_Ptr(const string&);
template SGridValue<double> * IDataContainer::getKey_Ptr(const string&);
template SGridValue<long double> * IDataContainer::getKey_Ptr(const string&);
template SExtrema<delphi_integer> * IDataContainer::getKey_Ptr(const string&);
template SExtrema<delphi_real> * IDataContainer::getKey_Ptr(const string&);
template CForce * IDataContainer::getKey_Ptr(const string&);
template CAtomPdb * IDataContainer::getKey_Ptr(const string&);
template vector<delphi_real>* IDataContainer::getKey_Ptr(const string&);
template vector<int>* IDataContainer::getKey_Ptr(const string&);
template vector<bool>* IDataContainer::getKey_Ptr(const string&);
template vector< SExtrema<delphi_real> >* IDataContainer::getKey_Ptr(const string&);
template vector< SGrid<delphi_integer> >* IDataContainer::getKey_Ptr(const string&);
template vector< SGrid<delphi_real> >* IDataContainer::getKey_Ptr(const string&);
template vector< SGridValue<delphi_real> >* IDataContainer::getKey_Ptr(const string&);
template vector< SGridValue<delphi_integer> >* IDataContainer::getKey_Ptr(const string&);
template vector< CAtomPdb >* IDataContainer::getKey_Ptr(const string&);
template vector< string >* IDataContainer::getKey_Ptr(const string&);
template vector< SDoubleGridValue >* IDataContainer::getKey_Ptr(const string&);

//----------template <class T> const T*** getKey_Ptr(const string& strKey,const int& iRows,const int& iColumns,const int& iPages)
template delphi_real *** IDataContainer::getKey_Ptr(const string&,const int&,const int&,const int&);
template SGrid <delphi_integer> *** IDataContainer::getKey_Ptr(const string&,const int&,const int&,const int&);
