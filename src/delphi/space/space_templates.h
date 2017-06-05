#ifndef SPACE_TEMPLATES_H
#define SPACE_TEMPLATES_H

#include "../interface/environment.h"

//template <class T> T ** get_pt2d(const delphi_integer& length, const delphi_integer& width)
template <class T> void get_pt2d(T ** & pt2d, const delphi_integer& length, const delphi_integer& width)
{
    delphi_integer i;
    //T ** pt2d;

    pt2d= new T *[length];
    for (i=0; i<length; i++)
    {
        pt2d[i] = new T [width]();

    }

    //return (pt2d);
}

//template <class T> T *** get_pt3d(const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight)
template <class T> void get_pt3d(T *** & pt3d, const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight)
//template <class T> T *** get_pt3d( delphi_integer& length,  delphi_integer& width,  delphi_integer& hight)
{
    delphi_integer i,j;
    //T *** pt3d;

    pt3d= new T **[length];
    for (i=0; i<length; i++)
    {
        pt3d[i] = new T* [width];
        for (j=0; j<width; j++)
        {
            pt3d[i][j]= new T[hight]();
        }
    }

    //return (pt3d);
}

//template <class T> T **** get_pt4d(const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight, const delphi_integer& page)
template <class T> void get_pt4d(T **** & pt4d, const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight, const delphi_integer& page)
//template <class T> T *** get_pt3d( delphi_integer& length,  delphi_integer& width,  delphi_integer& hight)
{
    delphi_integer i,j,k;
    //T **** pt4d;

    pt4d= new T ***[length];
    for (i=0; i<length; i++)
    {
        pt4d[i] = new T** [width];
        for (j=0; j<width; j++)
        {
            pt4d[i][j]= new T* [hight];
            for(k=0; k<hight; k++){
                 pt4d[i][j][k]= new T [page]();
            }
        }
    }

    //return (pt4d);
}



template <class T> void free_pt2d(T ** & pt2d, const delphi_integer& length, const delphi_integer& width)
{
    delphi_integer i;


    pt2d= new T *[length];
    for (i=0; i<length; i++)
    {
        //pt2d[i] = new T [width]();
        delete [] pt2d[i];
    }
    delete [] pt2d;
    pt2d=NULL;

}
 template <class T> void free_pt3d(T *** & pt3d, const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight)
//template <class T> T *** get_pt3d( delphi_integer& length,  delphi_integer& width,  delphi_integer& hight)
{
    delphi_integer i,j;
    //T *** pt3d;

    //pt3d= new T **[length];
    for (i=0; i<length; i++)
    {
        //pt3d[i] = new T* [width];
        for (j=0; j<width; j++)
        {
            //pt3d[i][j]= new T[hight]();
            delete [] pt3d[i][j];
        }
        delete [] pt3d[i];

    }
    delete [] pt3d;
    pt3d=NULL;
}



template <class T> void free_pt4d(T **** & pt4d, const delphi_integer& length, const delphi_integer& width, const delphi_integer& hight, const delphi_integer& page)
//template <class T> T *** get_pt3d( delphi_integer& length,  delphi_integer& width,  delphi_integer& hight)
{
   delphi_integer i,j,k;


    //pt4d= new T ***[length];
    for (i=0; i<length; i++)
    {
        //pt4d[i] = new T** [width];
        for (j=0; j<width; j++)
        {
            //pt4d[i][j]= new T* [hight];
            for(k=0; k<hight; k++){
                 //pt4d[i][j][k]= new T [page]();
                 delete [] pt4d[i][j][k];
            }
            delete []  pt4d[i][j];
        }
        delete [] pt4d[i];
    }
    delete [] pt4d;
    pt4d=NULL;
}



/*
// Move index:
//template <class T> T ** Move_index_2d(T ** ptr_old,const int& length, const int& width)
template <class T> void Move_index_2d(T ** & ptr_new, T ** & ptr_old,const int& length, const int& width)
{
    int i;
    //T ** ptr_new;
    ptr_new=new T *[length+1];
    for(i=0;i<=length-1;i++){
        ptr_new[i+1] = new T [width+1];
        ptr_new[i+1]=ptr_old[i]-1;
    }
    //return(ptr_new);

}

//template <class T> T *** Move_index_3d(T *** ptr_old,const int& length, const int& width, const int& height)
template <class T> void Move_index_3d(T *** & ptr_new, T *** & ptr_old,const int& length, const int& width, const int& height)
{
    int i,j;
    //T *** ptr_new;

    ptr_new=new T ** [length+1];
    for(i=0;i<=length-1;i++){
        ptr_new[i+1] = new T * [width+1];
        for(j=0;j<=width-1;j++){
            //cout << "i,j: " << i << j << endl;
            ptr_new[i+1][j+1] = new T [height+1];
            ptr_new[i+1][j+1]=ptr_old[i][j]-1;
        }
    }
    //return(ptr_new);

}

template <class T> void free_index_3d(T *** & ptr_new,const int& length, const int& width, const int& height)
{
    int i,j;
    //T *** ptr_new;


    for(i=2;i<=length-1;i++){

        for(j=2;j<=width-1;j++){
            //cout << "i,j: " << i << j << endl;
            cout << "### i,j: " << i << j << endl;
            delete [] ptr_new[i][j];

            //ptr_new[i+1][j+1]=ptr_old[i][j]-1;
        }

        delete [] ptr_new[i];

    }

    delete [] ptr_new;
    //return(ptr_new);

}
*/

#endif // SPACE_TEMPLATES_H
