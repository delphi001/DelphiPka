/*
 *
 * Author: Arghya Chakravorty (Argo)
 * writing Epsilon map for
 * 1> 2-dielectric model
 * 2> Gaussian module
 * 3> COnvolution module
 *
 * Written on 21-Sep, 2016
 *
 */

#include "io.h"

using namespace std;

// CIO::writeEpsMap(const delphi_integer& iAtomNumIn,const delphi_integer& iObjectNumIn,const delphi_integer& iGrid,
//                  const delphi_real& fScale,const SGrid<delphi_real>& fgBoxCenter,
//                  const vector< SGrid<delphi_integer> >& vctigEpsMap,const vector<bool>& vctbDielecMap,
//                  const string& strEpsFile)

void CIO::writeHomoEpsMap(const delphi_integer& iGrid,
                          const delphi_real&    repsout,
                          const delphi_real&    repsin,
                          const delphi_real& fScale,
                          const SGrid<delphi_real>& fgBoxCenter,
                          const vector<bool>& vctbDielecMap,
                          const string& strEpsFile)
{

  int ix, iy, iz, iw, l;

  delphi_real coeff=0.5291772108;
  delphi_real stepsize=1.0/fScale;
  SGrid<delphi_real> origin = (fgBoxCenter-stepsize*(iGrid-1)/2.0)/coeff;

  ofstream eps_out(strEpsFile.c_str());

  eps_out     << setw(10) << fixed << setprecision(6) << fScale << setw(6) << iGrid
              << setw(10) << setprecision(6) << fgBoxCenter.nX
              << setw(10) << setprecision(6) << fgBoxCenter.nY
              << setw(10) << setprecision(6) << fgBoxCenter.nZ
              << endl;
  eps_out << "Gaussian cube format Epsilon map for homogeneous dielectric distro." << endl;

  eps_out << fixed << setprecision(6);
  eps_out << setw(5) << right << 1 << setw(14) << right << origin.nX << setw(14) << right << origin.nY << setw(14) << right << origin.nZ << endl;
  eps_out << setw(5) << right << iGrid << setw(14) << right << stepsize/coeff << setw(14) << right << 0.0 << setw(14) << right << 0.0 << endl;
  eps_out << setw(5) << right << iGrid << setw(14) << right << 0.0 << setw(14) << right << stepsize/coeff << setw(14) << right << 0.0 << endl;
  eps_out << setw(5) << right << iGrid << setw(14) << right << 0.0 << setw(14) << right << 0.0 << setw(14) << right << stepsize/coeff << endl;
  eps_out << setw(5) << right << 1 << setw(14) << right << 0.0 << setw(14) << right << 0.0 << setw(14) << right << 0.0 << setw(14) << right << 0.0 << endl;

  eps_out.unsetf(ios_base::floatfield); // return to eps_out default notation

  eps_out.precision(5);


  for (ix = 1; ix <= iGrid; ix++)
  {
    for (iy = 1; iy <= iGrid; iy++)
    {
      l = 0;
      for (iz = 1; iz <= iGrid; iz++)
      {
        iw = (iz-1)*iGrid*iGrid + (iy-1)*iGrid + (ix-1);

        eps_out << setw(13) << scientific << (vctbDielecMap[iw]?repsout:repsin) << " " ;
        l++;
        if(l == 6)
        {
          eps_out << endl;
          l=0;
        }

      }
      eps_out << endl;
    }
  }

  eps_out.close();

}






void CIO::writeGaussEpsMap(const delphi_integer& iGrid,
                           const delphi_real&    repsout,
                           const delphi_real& fScale,
                           const SGrid<delphi_real>& fgBoxCenter,
                           const vector< SGrid<delphi_real> >& vctigEpsMap,
                           const string& strEpsFile)
{
  // cout << "Function called - Gaussian" << endl;

  delphi_real coeff=0.5291772108;
  delphi_real stepsize=1.0/fScale;
  delphi_real sum_eps = 0;

  SGrid<delphi_real> origin = (fgBoxCenter-stepsize*(iGrid-1)/2.0)/coeff;
  vector<delphi_real> grid_eps(iGrid*iGrid*iGrid);

  int i,j,k,l,d;

  for (i=1; i <= iGrid; i++)
  {
    for (j = 1; j <= iGrid; j++)
    {
      for ( k = 1; k <= iGrid; k++)
      {

        if ( !( i == iGrid || j == iGrid || k == iGrid || i == 1 || j == 1 || k == 1) )
        {
          sum_eps = vctigEpsMap[(i-1)*iGrid*iGrid + (j-1)*iGrid + (k-1)].nX +
                    vctigEpsMap[(i-2)*iGrid*iGrid + (j-1)*iGrid + (k-1)].nX +
                    vctigEpsMap[(i-1)*iGrid*iGrid + (j-1)*iGrid + (k-1)].nY +
                    vctigEpsMap[(i-1)*iGrid*iGrid + (j-2)*iGrid + (k-1)].nY +
                    vctigEpsMap[(i-1)*iGrid*iGrid + (j-1)*iGrid + (k-1)].nZ +
                    vctigEpsMap[(i-1)*iGrid*iGrid + (j-1)*iGrid + (k-2)].nZ;

          sum_eps /= 6;

          grid_eps[(i-1)*iGrid*iGrid + (j-1)*iGrid + (k-1)] = sum_eps;

        }
        else
        {
          grid_eps[(i-1)*iGrid*iGrid + (j-1)*iGrid + (k-1)] = repsout;
        }

        // cout << grid_eps[(i-1)*iGrid*iGrid + (j-1)*iGrid + (k-1)] << endl;

      }
    }
  }

  ofstream eps_out(strEpsFile.c_str());

  eps_out     << setw(10) << fixed << setprecision(6) << fScale << setw(6) << iGrid
              << setw(10) << setprecision(6) << fgBoxCenter.nX
              << setw(10) << setprecision(6) << fgBoxCenter.nY
              << setw(10) << setprecision(6) << fgBoxCenter.nZ
              << endl;
  eps_out << "Gaussian cube format Epsilon map for gaussian-smoothing mode" << endl;

  eps_out << fixed << setprecision(6);
  eps_out << setw(5) << right << 1 << setw(14) << right << origin.nX << setw(14) << right << origin.nY << setw(14) << right << origin.nZ << endl;
  eps_out << setw(5) << right << iGrid << setw(14) << right << stepsize/coeff << setw(14) << right << 0.0 << setw(14) << right << 0.0 << endl;
  eps_out << setw(5) << right << iGrid << setw(14) << right << 0.0 << setw(14) << right << stepsize/coeff << setw(14) << right << 0.0 << endl;
  eps_out << setw(5) << right << iGrid << setw(14) << right << 0.0 << setw(14) << right << 0.0 << setw(14) << right << stepsize/coeff << endl;
  eps_out << setw(5) << right << 1 << setw(14) << right << 0.0 << setw(14) << right << 0.0 << setw(14) << right << 0.0 << setw(14) << right << 0.0 << endl;

  eps_out.unsetf(ios_base::floatfield); // return to eps_out default notation

  eps_out.precision(5);

  vector<delphi_real>::iterator it = grid_eps.begin();

  l = 0;
  d = 0;
  while ( it != grid_eps.end())
  {
    eps_out << setw(13) << scientific << *it << " " ;
    it++;
    l++;
    d++;

    if (( l == 6) || ( d == iGrid ))
    {
      eps_out << endl;
      l = 0;
      if ( d == iGrid)
      {
        d = 0;
      }
    }

  }

  eps_out.close();
  grid_eps.resize(0);

}





void CIO::writeConvEpsMap(const delphi_integer& iGrid,
                          const delphi_real& fScale,
                          const SGrid<delphi_real>& fgBoxCenter,
                          vector<delphi_real>& vct_cepsmap,
                          const string& strEpsFile)
{

  int ix, iy, iz, iw, l,d;

  delphi_real coeff=0.5291772108;
  delphi_real stepsize=1.0/fScale;
  SGrid<delphi_real> origin = (fgBoxCenter-stepsize*(iGrid-1)/2.0)/coeff;

  ofstream eps_out(strEpsFile.c_str());

  eps_out     << setw(10) << fixed << setprecision(6) << fScale << setw(6) << iGrid
              << setw(10) << setprecision(6) << fgBoxCenter.nX
              << setw(10) << setprecision(6) << fgBoxCenter.nY
              << setw(10) << setprecision(6) << fgBoxCenter.nZ
              << endl;
  eps_out << "Gaussian cube format Epsilon map for convolution mode" << endl;

  eps_out << fixed << setprecision(6);
  eps_out << setw(5) << right << 1 << setw(14) << right << origin.nX << setw(14) << right << origin.nY << setw(14) << right << origin.nZ << endl;
  eps_out << setw(5) << right << iGrid << setw(14) << right << stepsize/coeff << setw(14) << right << 0.0 << setw(14) << right << 0.0 << endl;
  eps_out << setw(5) << right << iGrid << setw(14) << right << 0.0 << setw(14) << right << stepsize/coeff << setw(14) << right << 0.0 << endl;
  eps_out << setw(5) << right << iGrid << setw(14) << right << 0.0 << setw(14) << right << 0.0 << setw(14) << right << stepsize/coeff << endl;
  eps_out << setw(5) << right << 1 << setw(14) << right << 0.0 << setw(14) << right << 0.0 << setw(14) << right << 0.0 << setw(14) << right << 0.0 << endl;

  eps_out.unsetf(ios_base::floatfield); // return to eps_out default notation

  eps_out.precision(5);

  vector<delphi_real>::iterator it = vct_cepsmap.begin();

  l = 0;
  d = 0;
  while ( it != vct_cepsmap.end())
  {
    eps_out << setw(13) << scientific << *it << " " ;
    it++;
    l++;
    d++;

    if (( l == 6) || ( d == iGrid ))
    {
      eps_out << endl;
      l = 0;
      if ( d == iGrid)
      {
        d = 0;
      }
    }

  }

  eps_out.close();
}
