#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "../space/space.h"

using namespace std;

void CDelphiSpace::VdwToMs_piece(bool& cycle_flag, const delphi_integer& j, const delphi_integer& jx, const delphi_integer& jy, const delphi_integer& jz, const delphi_real& prbrd2, const SGrid <delphi_real>& s123, delphi_integer * eps, delphi_integer * eps2, delphi_integer & iacv ,delphi_real & dmn)
{
    delphi_integer ii,liml,limu,iarv;
    delphi_real dist=0;

    cycle_flag=false;
    liml=iab1[jx][jy][jz];
    limu=iab2[jx][jy][jz];

    for(ii=liml; ii<=limu; ii++)
    {
        iarv= icume[ii];
        //cout << "ii,icume[ii],iarv: " << ii << " " << icume[ii] << " " << iarv << " " << endl;

//2011-05-17 Using operations on coord and int_coord type
//variables defined in module operators_on_coordinates
        dist=optDot((s123-expos[iarv]),(s123-expos[iarv]));
        //cout << "dist: " << dist << endl;
        if (dist<prbrd2)
        {
            eps[j]=-1;
            eps2[j]=-1;
            //cycle D200
            cycle_flag=true;
        }
        else if (dist<dmn)
        {
            iacv=iarv;
            dmn=dist;
        }// if
    }// do
    //cout << "liml,limu,eps[j],eps2[j],iacv,dmn: " << liml << " " << limu << " " << eps[j] << " " << eps2[j] << " " << iacv << " " << dmn << endl;
}
