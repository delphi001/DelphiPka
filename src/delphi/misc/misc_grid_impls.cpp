/**
 * @file misc_grid_impls.h
 * @brief concrete implementations of template member functions used for grids
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * DO NOT COMPILE THIS FILE SEPERATELY! ONLY NEED THIS FILE TO BE INCLUDED IN .cpp file!!
 */

//Instantiations of overload operator <<
template ostream& operator<<(ostream &, SGrid<int> &);
template ostream& operator<<(ostream &, SGrid<long int> &);
template ostream& operator<<(ostream &, SGrid<float> &);
template ostream& operator<<(ostream &, SGrid<double> &);
template ostream& operator<<(ostream &, SGrid<long double> &);
template ostream& operator<<(ostream &, SGridValue<int> &);
template ostream& operator<<(ostream &, SGridValue<long int> &);
template ostream& operator<<(ostream &, SGridValue<float> &);
template ostream& operator<<(ostream &, SGridValue<double> &);
template ostream& operator<<(ostream &, SGridValue<long double> &);
template ostream& operator<<(ostream &, SExtrema<int> &);
template ostream& operator<<(ostream &, SExtrema<long int> &);
template ostream& operator<<(ostream &, SExtrema<float> &);
template ostream& operator<<(ostream &, SExtrema<double> &);
template ostream& operator<<(ostream &, SExtrema<long double> &);

//Instantiations of overload operator +
template SGrid <int> operator+(const SGrid <int> &, const SGrid <int> &);
template SGrid <int> operator+(const SGrid <int> &, const int &); 
template SGrid <int> operator+(const delphi_integer &, const SGrid <int> &);
template SGrid <long int> operator+(const SGrid <long int> &, const SGrid <long int> &);
template SGrid <long int> operator+(const SGrid <long int> &, const long int &); 
template SGrid <long int> operator+(const long int &, const SGrid <long int> &);
template SGrid <float> operator+(const SGrid <float> &, const SGrid <float> &);
template SGrid <float> operator+(const SGrid <float> &, const float &);
template SGrid <float> operator+(const float &, const SGrid <float> &);
template SGrid <double> operator+(const SGrid <double> &, const SGrid <double> &);
template SGrid <double> operator+(const SGrid <double> &, const double &);
template SGrid <double> operator+(const double &, const SGrid <double> &);
template SGrid <long double> operator+(const SGrid <long double> &, const SGrid <long double> &);
template SGrid <long double> operator+(const SGrid <long double> &, const long double &);
template SGrid <long double> operator+(const long double &, const SGrid <long double> &);

//Instantiations of overload operator -
template SGrid <int> operator-(const SGrid <int> &, const SGrid <int> &);
template SGrid <int> operator-(const SGrid <int> &, const int &);
template SGrid <int> operator-(const int &, const SGrid <int> &);
template SGrid <int> operator-(const SGrid <int> &);
template SGrid <long int> operator-(const SGrid <long int> &, const SGrid <long int> &);
template SGrid <long int> operator-(const SGrid <long int> &, const long int &);
template SGrid <long int> operator-(const long int &, const SGrid <long int> &);
template SGrid <long int> operator-(const SGrid <long int> &);
template SGrid <float> operator-(const SGrid <float> &, const SGrid <float> &);
template SGrid <float> operator-(const SGrid <float> &, const float &);
template SGrid <float> operator-(const float &, const SGrid <float> &);
template SGrid <float> operator-(const SGrid <float> &);
template SGrid <double> operator-(const SGrid <double> &, const SGrid <double> &);
template SGrid <double> operator-(const SGrid <double> &, const double &);
template SGrid <double> operator-(const double &, const SGrid <double> &);
template SGrid <double> operator-(const SGrid <double> &);
template SGrid <long double> operator-(const SGrid <long double> &, const SGrid <long double> &);
template SGrid <long double> operator-(const SGrid <long double> &, const long double &);
template SGrid <long double> operator-(const long double &, const SGrid <long double> &);
template SGrid <long double> operator-(const SGrid <long double> &);

//Instantiations of overload operator *
template SGrid <int> operator*(const SGrid <int> &, const SGrid <int> &);
template SGrid <int> operator*(const SGrid <int> &, const int &);
template SGrid <int> operator*(const int &, const SGrid <int> &);
template SGrid <long int> operator*(const SGrid <long int> &, const SGrid <long int> &);
template SGrid <long int> operator*(const SGrid <long int> &, const long int &);
template SGrid <long int> operator*(const long int &, const SGrid <long int> &);
template SGrid <float> operator*(const SGrid <float> &, const SGrid <float> &);
template SGrid <float> operator*(const SGrid <float> &, const float &);
template SGrid <float> operator*(const float &, const SGrid <float> &);
template SGrid <double> operator*(const SGrid <double> &, const SGrid <double> &);
template SGrid <double> operator*(const SGrid <double> &, const double &);
template SGrid <double> operator*(const double &, const SGrid <double> &);
template SGrid <long double> operator*(const SGrid <long double> &, const SGrid <long double> &);
template SGrid <long double> operator*(const SGrid <long double> &, const long double &);
template SGrid <long double> operator*(const long double &, const SGrid <long double> &);

//Instantiations of overload operator /
template SGrid <float> operator/(const SGrid <float> &, const float &);  
template SGrid <double> operator/(const SGrid <double> &, const double &);  
template SGrid <long double> operator/(const SGrid <long double> &, const long double &);  
template SGrid <float> operator/(const SGrid <int> &, const float &);  
template SGrid <double> operator/(const SGrid <int> &, const double &);  
template SGrid <long double> operator/(const SGrid <int> &, const long double &);
template SGrid <float> operator/(const SGrid <long int> &, const float &);  
template SGrid <double> operator/(const SGrid <long int> &, const double &);  
template SGrid <long double> operator/(const SGrid <long int> &, const long double &);  
template SGrid <int> operator/(const SGrid <int> &, const int &);
template SGrid <long int> operator/(const SGrid <long int> &, const long int &);

//Instantiations of dot product
template int optDot(const SGrid <int> &, const SGrid <int> &);
template long int optDot(const SGrid <long int> &, const SGrid <long int> &);
template float optDot(const SGrid <float> &, const SGrid <float> &);
template double optDot(const SGrid <double> &, const SGrid <double> &);
template long double optDot(const SGrid <long double> &, const SGrid <long double> &);

//Instantiations of cross product
template SGrid <int> optCross(const SGrid <int> &, const SGrid <int> &);     
template SGrid <long int> optCross(const SGrid <long int> &, const SGrid <long int> &);
template SGrid <float> optCross(const SGrid <float> &, const SGrid <float> &);
template SGrid <double> optCross(const SGrid <double> &, const SGrid <double> &);
template SGrid <long double> optCross(const SGrid <long double> &, const SGrid <long double> &);

//Instantiations of optABS()
template SGrid <int> optABS(const SGrid <int> &);
template SGrid <long int> optABS(const SGrid <long int> &);
template SGrid <float> optABS(const SGrid <float> &);   
template SGrid <double> optABS(const SGrid <double> &);   
template SGrid <long double> optABS(const SGrid <long double> &);   

//Instantiations of optCast()
template SGrid <int> optCast(const SGrid <float> &);
template SGrid <int> optCast(const SGrid <double> &);
template SGrid <int> optCast(const SGrid <long double> &);
template SGrid <long int> optCast(const SGrid <float> &);
template SGrid <long int> optCast(const SGrid <double> &);
template SGrid <long int> optCast(const SGrid <long double> &);
template SGrid <float> optCast(const SGrid <int> &);
template SGrid <double> optCast(const SGrid <int> &);
template SGrid <long double> optCast(const SGrid <int> &);
template SGrid <float> optCast(const SGrid <long int> &);
template SGrid <double> optCast(const SGrid <long int> &);
template SGrid <long double> optCast(const SGrid <long int> &);

//Instantiations of optMin()  
template SGrid <int> optMin(const int &, const SGrid <int> &);   
template SGrid <long int> optMin(const long int &, const SGrid <long int> &);   
template SGrid <int> optMin(const SGrid <int> &, const int &);
template SGrid <long int> optMin(const SGrid <long int> &, const long int &);
template SGrid <int> optMin(const SGrid <int> &, const SGrid <int> &);
template SGrid <long int> optMin(const SGrid <long int> &, const SGrid <long int> &);
template int optMin(const SGrid <int> & nVector);
template long int optMin(const SGrid <long int> & nVector);
template SGrid <float> optMin(const float &, const SGrid <float> &);
template SGrid <double> optMin(const double &, const SGrid <double> &);
template SGrid <long double> optMin(const long double &, const SGrid <long double> &);
template SGrid <float> optMin(const SGrid <float> &, const float &);
template SGrid <double> optMin(const SGrid <double> &, const double &);
template SGrid <long double> optMin(const SGrid <long double> &, const long double &);
template SGrid <float> optMin(const SGrid <float> &, const SGrid <float> &);
template SGrid <double> optMin(const SGrid <double> &, const SGrid <double> &);
template SGrid <long double> optMin(const SGrid <long double> &, const SGrid <long double> &);
template float optMin(const SGrid <float> & nVector);
template double optMin(const SGrid <double> & nVector);
template long double optMin(const SGrid <long double> & nVector);

//Instantiations of optMax()
template SGrid <int> optMax(const int &, const SGrid <int> &);
template SGrid <long int> optMax(const long int &, const SGrid <long int> &);
template SGrid <int> optMax(const SGrid <int> &, const int &);
template SGrid <long int> optMax(const SGrid <long int> &, const long int &); 
template SGrid <int> optMax(const SGrid <int> &, const SGrid <int> &);
template SGrid <long int> optMax(const SGrid <long int> &, const SGrid <long int> &);
template int optMax(const SGrid <int> &);
template long int optMax(const SGrid <long int> &);   
template SGrid <float> optMax(const float &, const SGrid <float> &);
template SGrid <double> optMax(const double &, const SGrid <double> &);
template SGrid <long double> optMax(const long double &, const SGrid <long double> &);
template SGrid <float> optMax(const SGrid <float> &, const float &);
template SGrid <double> optMax(const SGrid <double> &, const double &);
template SGrid <long double> optMax(const SGrid <long double> &, const long double &);
template SGrid <float> optMax(const SGrid <float> &, const SGrid <float> &);
template SGrid <double> optMax(const SGrid <double> &, const SGrid <double> &);
template SGrid <long double> optMax(const SGrid <long double> &, const SGrid <long double> &);
template float optMax(const SGrid <float> &);
template double optMax(const SGrid <double> &);
template long double optMax(const SGrid <long double> &);

//Instantiations of optMinSign()
template SGrid <float> optMinSign(const float &, const SGrid <float> &);
template SGrid <double> optMinSign(const double &, const SGrid <double> &);
template SGrid <long double> optMinSign(const long double &, const SGrid <long double> &);
template SGrid <float> optMinSign(const SGrid <float> &, const float &);
template SGrid <double> optMinSign(const SGrid <double> &, const double &);
template SGrid <long double> optMinSign(const SGrid <long double> &, const long double &);

//Instantiations of optMaxSign()
template SGrid <float> optMaxSign(const float &, const SGrid <float> &);
template SGrid <double> optMaxSign(const double &, const SGrid <double> &);
template SGrid <long double> optMaxSign(const long double &, const SGrid <long double> &);
template SGrid <float> optMaxSign(const SGrid <float> &, const float &);
template SGrid <double> optMaxSign(const SGrid <double> &, const double &);
template SGrid <long double> optMaxSign(const SGrid <long double> &, const long double &);
   
//Instantiations of optSubMin() 
template SGrid <float> optSubMin(const SGrid <float> &, const SGrid <float> &, const SGrid <float> &); 
template SGrid <double> optSubMin(const SGrid <double> &, const SGrid <double> &, const SGrid <double> &); 
template SGrid <long double> optSubMin(const SGrid <long double> &, const SGrid <long double> &, const SGrid <long double> &);
                                   
//Instantiations of optSubMax()                                 
template SGrid <float> optSubMax(const SGrid <float> &, const SGrid <float> &, const SGrid <float> &);
template SGrid <double> optSubMax(const SGrid <double> &, const SGrid <double> &, const SGrid <double> &);
template SGrid <long double> optSubMax(const SGrid <long double> &, const SGrid <long double> &, const SGrid <long double> &);
                                   
//Instantiations of optSum()     
template int optSum(const SGrid <int> &);   
template long int optSum(const SGrid <long int> &);   
template float optSum(const SGrid <float> &);
template double optSum(const SGrid <double> &);
template long double optSum(const SGrid <long double> &);
     
//Instantiations of optComp()    
template int optComp(const SGrid <int> &, const int &);
template long int optComp(const SGrid <long int> &, const int &);
template float optComp(const SGrid <float> &, const int &);   
template double optComp(const SGrid <double> &, const int &);   
template long double optComp(const SGrid <long double> &, const int &);   
   
//Instantiations of optORLT()  
template bool optORLT(const SGrid <int> &, const SGrid <int> &);    
template bool optORLT(const SGrid <long int> &, const SGrid <long int> &);       
template bool optORLT(const SGrid <int> &, const int &);   
template bool optORLT(const SGrid <long int> &, const long int &);   
template bool optORLT(const SGrid <float> &, const SGrid <float> &); 
template bool optORLT(const SGrid <double> &, const SGrid <double> &); 
template bool optORLT(const SGrid <long double> &, const SGrid <long double> &);
template bool optORLT(const SGrid <float> &, const float &);
template bool optORLT(const SGrid <double> &, const double &);
template bool optORLT(const SGrid <long double> &, const long double &);
   
//Instantiations of optANDLT()     
template bool optANDLT(const SGrid <int> &, const SGrid <int> &);
template bool optANDLT(const SGrid <long int> &, const SGrid <long int> &);
template bool optANDLT(const SGrid <int> &, const int &);
template bool optANDLT(const SGrid <long int> &, const long int &);
template bool optANDLT(const SGrid <float> &, const SGrid <float> &);
template bool optANDLT(const SGrid <double> &, const SGrid <double> &);
template bool optANDLT(const SGrid <long double> &, const SGrid <long double> &);
template bool optANDLT(const SGrid <float> &, const float &);
template bool optANDLT(const SGrid <double> &, const double &);
template bool optANDLT(const SGrid <long double> &, const long double &);
   
//Instantiations of optORLE()     
template bool optORLE(const SGrid <int> &, const SGrid <int> &);
template bool optORLE(const SGrid <long int> &, const SGrid <long int> &);
template bool optORLE(const SGrid <int> &, const int &);  
template bool optORLE(const SGrid <long int> &, const long int &);   
template bool optORLE(const SGrid <float> &, const SGrid <float> &);
template bool optORLE(const SGrid <double> &, const SGrid <double> &);
template bool optORLE(const SGrid <long double> &, const SGrid <long double> &);
template bool optORLE(const SGrid <float> &, const float &);
template bool optORLE(const SGrid <double> &, const double &);
template bool optORLE(const SGrid <long double> &, const long double &);
   
//Instantiations of optANDLE()    
template bool optANDLE(const SGrid <int> &, const SGrid <int> &);
template bool optANDLE(const SGrid <long int> &, const SGrid <long int> &);
template bool optANDLE(const SGrid <int> &, const int &);   
template bool optANDLE(const SGrid <long int> &, const long int &);   
template bool optANDLE(const SGrid <float> &, const SGrid <float> &);
template bool optANDLE(const SGrid <double> &, const SGrid <double> &);
template bool optANDLE(const SGrid <long double> &, const SGrid <long double> &);
template bool optANDLE(const SGrid <float> &, const float &);
template bool optANDLE(const SGrid <double> &, const double &);
template bool optANDLE(const SGrid <long double> &, const long double &);
   
//Instantiations of optORGT()   
template bool optORGT(const SGrid <int> &, const SGrid <int> &);
template bool optORGT(const SGrid <long int> &, const SGrid <long int> &);
template bool optORGT(const SGrid <int> &, const int &);
template bool optORGT(const SGrid <long int> &, const long int &);
template bool optORGT(const SGrid <float> &, const SGrid <float> &);
template bool optORGT(const SGrid <double> &, const SGrid <double> &);
template bool optORGT(const SGrid <long double> &, const SGrid <long double> &);
template bool optORGT(const SGrid <float> &, const float &);
template bool optORGT(const SGrid <double> &, const double &);
template bool optORGT(const SGrid <long double> &, const long double &);
   
//Instantiations of optORGE()
template bool optORGE(const SGrid <int> &, const SGrid <int> &);
template bool optORGE(const SGrid <long int> &, const SGrid <long int> &);
template bool optORGE(const SGrid <int> &, const int &);  
template bool optORGE(const SGrid <long int> &, const long int &);  
template bool optORGE(const SGrid <float> &, const SGrid <float> &);
template bool optORGE(const SGrid <double> &, const SGrid <double> &);
template bool optORGE(const SGrid <long double> &, const SGrid <long double> &);
template bool optORGE(const SGrid <float> &, const float &);
template bool optORGE(const SGrid <double> &, const double &);
template bool optORGE(const SGrid <long double> &, const long double &);
   
//Instantiations of optANDGT()   
template bool optANDGT(const SGrid <int> &, const SGrid <int> &);
template bool optANDGT(const SGrid <long int> &, const SGrid <long int> &);
template bool optANDGT(const SGrid <int> &, const int &);   
template bool optANDGT(const SGrid <long int> &, const long int &);   
template bool optANDGT(const SGrid <float> &, const SGrid <float> &);
template bool optANDGT(const SGrid <double> &, const SGrid <double> &);
template bool optANDGT(const SGrid <long double> &, const SGrid <long double> &);
template bool optANDGT(const SGrid <float> &, const float &);
template bool optANDGT(const SGrid <double> &, const double &);
template bool optANDGT(const SGrid <long double> &, const long double &);
   
//Instantiations of optANDGE()      
template bool optANDGE(const SGrid <int> &, const SGrid <int> &);
template bool optANDGE(const SGrid <long int> &, const SGrid <long int> &);
template bool optANDGE(const SGrid <int> &, const int &);
template bool optANDGE(const SGrid <long int> &, const long int &);
template bool optANDGE(const SGrid <float> &, const SGrid <float> &);
template bool optANDGE(const SGrid <double> &, const SGrid <double> &);
template bool optANDGE(const SGrid <long double> &, const SGrid <long double> &);
template bool optANDGE(const SGrid <float> &, const float &);
template bool optANDGE(const SGrid <double> &, const double &);
template bool optANDGE(const SGrid <long double> &, const long double &);

//Instantiations of overlaod operator !=()
template bool operator!=(const SGrid <int> &, const int &);
template bool operator!=(const SGrid <long int> &, const long int &);
template bool operator!=(const SGrid <int> &, const SGrid <int> &);
template bool operator!=(const SGrid <long int> &, const SGrid <long int> &);
template bool operator!=(const SGrid <float> &, const float &);
template bool operator!=(const SGrid <double> &, const double &);
template bool operator!=(const SGrid <long double> &, const long double &);
template bool operator!=(const SGrid <float> &, const SGrid <float> &);
template bool operator!=(const SGrid <double> &, const SGrid <double> &);
template bool operator!=(const SGrid <long double> &, const SGrid <long double> &);


       
