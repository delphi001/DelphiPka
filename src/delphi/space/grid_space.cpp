#include "space.h"


//int_coord(): covert 3 delphi_integer numbers to be a SGrid <delphi_integer> structure
SGrid <int> CDelphiSpace::int_coord(const int& a, const int& b, const int& c)
{
    SGrid <int>  d;
    d.nX=a;
    d.nY=b;
    d.nZ=c;

    return(d);
}


//coord(): covert 3 float numbers to be a SGrid <float> structure
SGrid <float> CDelphiSpace::coord(const float& a,  const float& b,  const float& c)
{
    SGrid <float>  d;
    d.nX=a;
    d.nY=b;
    d.nZ=c;

    return(d);
}



//Float2Int(): 1. convert SGrid <float> structure to SGrid <int> structure; 2. convert float to int


SGrid <int> CDelphiSpace::Float2Int( const SGrid <float>& a )
{
    SGrid <int> b;

    b.nX=int(a.nX);
    b.nY=int(a.nY);
    b.nZ=int(a.nZ);

    return(b);
}

int CDelphiSpace::Float2Int(const float& a){

    int b;
    b=int(a);
    return(b);

}

//Int2Float(): 1. convert SGrid <int> structure to SGrid <float> structure; 2. convert int to float
SGrid <float> CDelphiSpace::Int2Float( const SGrid <int>& a )
{
    SGrid <float> b;

    b.nX=float(a.nX);
    b.nY=float(a.nY);
    b.nZ=float(a.nZ);

    return(b);
}

float CDelphiSpace::Int2Float(const int& a){

    float b;
    b=float(a);
    return(b);

}




