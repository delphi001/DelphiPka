#include <complex>
#include <iostream>
#include <cmath>
#include <valarray>
#include <algorithm>

#include "math.h"
#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "space_templates.h"
#include "space.h"

using namespace std;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;


void CDelphiSpace::Convolute ()
{

	//Input signal
    	//----------------------Change to HRhomap


    	int	nnt = pow(2,ceil(log(iGrid)/log(2)));
    	cout << "  FFT> Grid size enlarged to " << nnt << " for FFT/iFFT operations .." << endl;


    	//Creating GAUSSIAN kernel
    	CArray gkernel_X({0.0,0.0},nnt);
    	gkernel_X = gaussianKernel( nnt );

    	CArray fftg;
    	fftg = fft(gkernel_X);

	//Will pick up one ix,iy value and Convolute the z-HS
    //	int ix = (iGrid+1)*0.5 - 3;
    	int ix,iy,iz;

		CArray ffth({0.0,0.0} , nnt);
		CArray fftp({0.0,0.0} , nnt);
		CArray ffti({0.0,0.0} , nnt);

    	for ( ix = 1; ix <= iGrid; ix++) {
			for ( iy = 1; iy <=iGrid; iy ++ )
			{
				//Padding the input map row with zeroes
				//Final length = nnt == (2^n)
				CArray hs({1.0,0.0},nnt);

				for ( iz = 1; iz <= iGrid; iz++ ) {
					hs[iz - 1 + (nnt - iGrid + 1)/2]= {HRhomap[ix][iy][iz], 0.0 };
				}
				ffth = fft(hs);

				//Multiplying the FFTs
				fftp = ffth * fftg;

				//Inverse FFTs ----------------- should go to HRhomap ---- Check length
				ffti = ifft(fftp);

				ffti /= ffti[0].real();

				//update HRhomap
				for (  iz = 1; iz <= iGrid; iz++ ) {

					HRhomap[ix][iy][iz] = ffti[iz - 1 + (nnt - iGrid + 1)/2].real();
					//cout << HRhomap[ix][iy][iz].nZ << " ";
				}
				//cout << "" << endl;

		}//foreach iY

   	}//foreach ix
   	cout << "  FFT> Completed Z ..... " << endl;

   	for ( ix = 1; ix <= iGrid; ix++) {
			for ( iz = 1; iz <=iGrid; iz ++ )
			{

				CArray hs({1.0,0.0},nnt);
				for ( iy = 1; iy <= iGrid; iy++ ) {
					hs[iy - 1 + (nnt - iGrid + 1)/2]= {HRhomap[ix][iy][iz], 0.0 };
				}
				ffth = fft(hs);
				//Multiplying the FFTs
				fftp = ffth * fftg;
				//Inverse FFTs ----------------- should go to HRhomap ---- Check length
				ffti = ifft(fftp);
				ffti /= ffti[0].real();

				//update HRhomap
				for (  iy = 1; iy <= iGrid; iy++ ) {
					HRhomap[ix][iy][iz] = ffti[iy - 1 + (nnt - iGrid + 1)/2].real();
					//cout << (int) HRhomap[ix][iy][iz].nZ ;
				}
				//cout << "" << endl;
		}//foreach iZ
   	}//foreach ix
   	cout << "  FFT> Completed Y ..... " << endl;


   	for ( iy = 1; iy <= iGrid; iy++) {
			for ( iz = 1; iz <=iGrid; iz ++ )
			{
				//Padding the input map row with zeroes
				//Final length = nnt == (2^n)
				CArray hs({1.0,0.0},nnt);

				for ( ix = 1; ix <= iGrid; ix++ ) {
					hs[ix - 1 + (nnt - iGrid + 1)/2]= {HRhomap[ix][iy][iz], 0.0 };
				}

				ffth = fft(hs);

				//Multiplying the FFTs
				fftp = ffth * fftg;

				//Inverse FFTs ----------------- should go to HRhomap ---- Check length
				ffti = ifft(fftp);
				ffti /= ffti[0].real();

				//update HRhomap
				for (  ix = 1; ix <= iGrid; ix++ ) {

					HRhomap[ix][iy][iz]= ffti[ix - 1 + (nnt - iGrid + 1)/2].real();
					//cout << (int) HRhomap[ix][iy][iz].nZ ;
				}
				//cout << "" << endl;

		}//foreach iZ

   	}//foreach ix
   	cout << "  FFT> Completed X ..... " << endl;

   	fftg.resize(0);
   	ffth.resize(0);
   	ffti.resize(0);
   	fftp.resize(0);


}







//Create Gaussian discrete series to convolute
//------------ gaussian array -----------------
CArray CDelphiSpace::gaussianKernel( int N )
{

	double invscale = (1.0/fScale);
	int midpoint = N/2;

	cout << "  FFT> Midpoint of extended array = " << midpoint << endl;
	cout << "  FFT> Creating GAUSSIAN kernel with SIGMA = " << fksigma << endl;

	CArray gauss({0.0,0.0},N);

    for ( int j = 0; j < N; j++ ){

        double pwr = -1*(j - midpoint)*(j - midpoint)*invscale*invscale/(2*fksigma*fksigma);
        gauss[j] = (1/(sqrt(2*3.14)*fksigma))*exp(pwr);

    }

    //wrapping gauss
    CArray gauss2({0.0,0.0},N);

    for ( int j  = midpoint; j < 2*midpoint; j++ ) {
        gauss2[j - midpoint] = gauss[j];
    }
    for ( int j = 0; j < midpoint; j++ ) {
        gauss2[j + midpoint] = gauss[j];
    }



    //cout << "------------GAUSS2 ------------ ^^^ " << endl;
    gauss.resize(0);

    return gauss2;
}








CArray CDelphiSpace::fft ( CArray& x ) {

    // DFT
    unsigned long N = x.size(), k = N, n;
    double thetaT = 3.14159265358979323846264338328L / N;
    Complex phiT = Complex(cos(thetaT), sin(thetaT)), T;
    while (k > 1)
    {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0L;
        for (unsigned int l = 0; l < k; l++)
        {
            for (unsigned long a = l; a < N; a += n)
            {
                unsigned long b = a + k;
                Complex t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    // Decimate
    unsigned int m = (unsigned int)log2(N);
    for (unsigned int a = 0; a < N; a++)
    {
        unsigned int b = a;
        // Reverse bits
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a)
        {
            Complex t = x[a];
            x[a] = x[b];
            x[b] = t;
        }
    }

    return x;
}



CArray CDelphiSpace::ifft( CArray& x )
{
    unsigned long sz = x.size();
    CArray ifftx(sz);

    // conjugate the complex numbers
    x = x.apply(conj);

    // forward fft
    ifftx = fft( x );

    // conjugate the complex numbers again
    ifftx = ifftx.apply(std::conj);

    // scale the numbers
    ifftx /= ifftx.size();

    return ifftx;
}



//CArray CDelphiSpace::getMax ( CArray & x ) {
//
//    double m = 0.0;
//    for (int i = 0; i < x.size(); i++ ) {
//        m = max(m,x[i].real());
//    }
//    CArray maxArr({1/m,0.0},x.size());
//    return maxArr;
//}
