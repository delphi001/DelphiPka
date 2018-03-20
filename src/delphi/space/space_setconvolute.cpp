//#define DEBUG_DELPHI_SPACE

//ARGO: 2016-FEB-09 --> Including headers for I/O for VDW_ssurface grids
#include <iostream>
#include <fstream>
//


#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include <vector>


#include "../space/space.h"
#ifdef PARALLEL_OMP
#include <omp.h>
#endif

using namespace std;
//#####################################################################
//ATTENTION! This file is part of epsmakmod module.
// Do not compile separately!
//
// 2011-05-12 Parameters transfered via modules
//#####################################################################

/*------------------------- CONVOLUTE ----------------------------
creating a HEAVYSIDE RHO DISTRO --> CONVOLUTING IT WITH GAUSSIAN KERNEL
template of 1st part - till DO ATOM - is taken from setOut;
it is then merged with setGaussian routine and the calculations carry on thereafter.
---------------------------------------------------------------------*/

void CDelphiSpace::setConvolute()
{

//2011-05-12 Some local variables

    SGrid <delphi_real>  sqtemp[31],rad2aavtemp[31];
    SGrid <delphi_real> * sq=sqtemp+15;
    SGrid <delphi_real> * rad2aav=rad2aavtemp+15;
    SGrid <delphi_integer>* ioff;   ioff = NULL;


    delphi_real coeff,stepsize;
    delphi_real eps_front, eps_back, eps_left, eps_right, eps_top, eps_bott,eps_min,eps_max;
    delphi_real eps_diff = 0.6;
    vector<delphi_real> vct_cepsmap;
	vector<delphi_real> neighbor_eps;

    bool itobig,itest2;
    string strtmp,strtmp1;
    delphi_integer epsdim;
    delphi_integer longint,n;

    SGrid <delphi_real> ddist,dxyz,ddxyz,xn;
    SGrid <delphi_real> rad2av,fxn,vtemp;
    SGrid <delphi_real> origin;
    SGrid <delphi_integer> ismin,ismax,idist,ixyz,itest,ixn,i123;

    debug_space = false;
    bool hs_rhomap=true;
    bool conv_rhomap=true;

/**
 * here fRadPrb is not zero only if one wants to map the ext}//ed
 *surf. that has to be done with fRadPrb(1) if solute-solvent
 *interface is concerned imedia = medium number for a object
 *a non-zero entry in iEpsMap indicates an atom # plus 1 (to
 *properly treat re-entrant mid-points later (15 Aug 93)
 */
//2011-05-27 Declarations added due to IMPLICIT NONE
    //delphi_integer iboxt;
    delphi_integer iac,ibox,igrdc,iv,ix,iy,iz;
    delphi_integer limmax,lim;
    delphi_real dis2min2,dis2min1,distsq,dist,rad,rad2;
    //delphi_real rad4,rad5;
    delphi_real rad2a,radmax2,fRMid,radp2,radtest,radprobe;
    delphi_real radp,temp;
    delphi_real sigma = 0.93;
    delphi_real dentemp = exp(-1/(sigma*sigma));
    delphi_real den;
    delphi_real sigmatime,pi,radsq;
    int i,j,k;
    delphi_real eps_medium,rho_hs;

	//Central Y-Z PLANE
	int fx = (iGrid - 1)/2;


	//HEAVYSIDE EPSILON based on GAUSSIAN RHO[in].
	rho_hs = (repsout - fhvsd_eps)/(repsout-repsin);


    //inhomo = 1 ==> WATER
    if (!(inhomo==0 && logs))
	{

		if(debug_space) cout << "############ beginning of setConv: RUN1 ##############" << endl;

#ifdef VERBOSE
		cout << " Conv> Starting creation of Convoluted Epsilon Map: " << endl;
#endif
    cout << " Conv> Heavyside step will be set up where Rho[in] = " << rho_hs << " corresponding to Heavyside Epsilon of " << fhvsd_eps << endl;

		radprobe=0; // this radprobe seems not useful.
		epsdim=iNatom+iNObject+2;
		//iboxt=0;
		radmax2=0.0;
		fRMid=float((iGrid+1)/2);
		itest2=false;

		pi=3.14159265359;
		sigmatime=3.0;



#ifdef DEBUG_DELPHI_SPACE
		cout << " Conv> pdb radius 0 & iNatom-1: " << sDelPhiPDB[0].radius << " " << sDelPhiPDB[iNatom-1].radius << endl;
#endif // DEBUG_DELPHI_SPACE

		if (iNatom>0)
		{
	//2011-05-12 Label 691 removed
			for(ix=1; ix<=iNatom; ix++)
			{
	//2011-05-13 Changed to derive-type array sDelPhiPDB (pointers module)
				radmax2=max(radmax2,sDelPhiPDB[ix].radius);
			}// do

#ifdef DEBUG_DELPHI_SPACE
			cout << "radmax2: " << radmax2 << endl;
#endif // DEBUG_DELPHI_SPACE

	/**
	 * this is probably the best way to do it,depending on which
	 * surf. is desired
	 */
			temp=max(radprobe,fExternRadius);



#ifdef DEBUG_DELPHI_SPACE
			cout << "radprobe: " << radprobe << " fExternRadius: " << fExternRadius << endl;
#endif // DEBUG_DELPHI_SPACE


			lim=1+radmax2;
			radmax2=sigmatime*fScale*(radmax2*sigma+temp); //Gaussian: 3 sigma plus temp. sigma=sigma*radius.Now radmax is 3 sigma.
			limmax = 12;
			itobig=false;

			if(lim>limmax) itobig=true;
			igrdc=pow((2*lim+1),3);
			ioff = new SGrid <int> [igrdc];



			if (!itobig)
			{
				//2011-05-12 removed goto 7878 statement
				radtest= pow( (radmax2 + 0.5*sqrt(3.0)),2 );
				ibox=-1;

				//2011-05-12 Strange statement. May allocate or may not
				//allocate array that used later in the program
				//irrespectively of itobig value, thus moved array
				//allocation before if condition


				if(debug_space) cout << "radmax2: " << radmax2 << " radtest: " << radtest << " lim: " << lim << endl;
#ifdef PARALLEL_OMP
			#pragma omp for schedule(auto)
#endif

				for(ix=-lim; ix<=lim; ix++)
				{
					for(iy=-lim; iy<=lim; iy++)
					{
						for(iz=-lim; iz<=lim; iz++)
						{
							//2011-05-13 Replaced by faster operation
							//2011-05-13 Using operations on coord and
							//int_coord type variables defined in module
							//operators_on_coordinates
							idist=int_coord(ix,iy,iz);
							dist=float ( optDot(idist,idist) );
							ddist = dist + float(0.25) + optCast <delphi_real,delphi_integer> (idist);

							if ((dist<radtest)|| optORLT( ddist,radtest ))
							{
								ibox++;
								ioff[ibox]=idist;
								//cout << ibox <<": " << ioff[ibox] << idist << endl;

							}// if
						}// do
					}// do
				}// do
			}// if


		}// if
		if(debug_space) cout << ibox <<": " << ioff[ibox] << idist << endl;


		//set interiors in MOLECULES

		//ARGO: 2016-FEB-09: modified lines here
		//ORIGINIAL: 	if(itest2||itobig) cout <<"setout method 1 " << itest2 << " " << itobig << endl;
		if(itest2||itobig) cout << " Conv> Setout method 1 : itest2 = " << itest2 << " and itobig = " << itobig << endl;

		//DoATOMS:
		for( iv=1; iv<=iNatom; iv++)
		{
			/**
			 * restore values
			 */
			rad= sDelPhiPDB[iv].radius;
			rad= sigmatime*sigma*rad; // 3 sigma for Gaussian

			xn=xn2[iv];
			//cout << "iv,xn1[iv],xn2[iv]: " <<iv<< xn1[iv+1] << xn2[iv+1] << endl;

			//2011-05-13 Removed GOTO statement
			if (rad<1.e-6)
			{
				continue;
			}// if

			//fScale radius to grid
			rad=rad*fScale;
			//rad5=pow( (rad+0.5),2);
			radp=rad+fExternRadius*fScale;
			rad=rad+radprobe*fScale;
			//rad4=pow( (rad+0.5),2); // not used
			rad2=rad*rad;
			radp2=radp*radp;

			radsq=rad2/(sigmatime*sigmatime); // from Gaussian

			/**
			 * set dielectric map
			 * check if sphere sits within limits of box
			*/
			itest2=false;

			ismin=optCast <delphi_integer,delphi_real> (xn-radmax2-1.0);
			ismax=optCast <delphi_integer,delphi_real> (xn+radmax2+1.0);
			itest=ismin;
			ismin=optMin(ismin,iGrid);
			ismin=optMax(ismin,1);
			if(itest!=ismin) itest2=true;
			itest=ismax;
			ismax=optMin(ismax,iGrid);
			ismax=optMax(ismax,1);
			if(itest != ismax) itest2=true;



			if (itest2||itobig)   //slow method;
			{
				//2011-05-13 Seems redundant statement
				if(debug_space) cout << " Conv> ### slow method:"  << endl;
				if(debug_space) cout << " Conv> itest, itest2, itobig: " << itest << ", " << itest2 << ", " << itobig << endl;
				rad2a = rad2 - 0.25;



				for(iz=ismin.nZ; iz<=ismax.nZ; iz++)
				{
					for(iy=ismin.nY; iy<=ismax.nY; iy++)
					{
						for(ix=ismin.nX; ix<=ismax.nX; ix++)
						{

							ixyz=int_coord(ix,iy,iz);
							dxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
							distsq=optDot(dxyz,dxyz);
							dxyz=dxyz+distsq;

//	                        if (dxyz.nX<rad2a)
//	                        {
//	                            iepsmp[iz][iy][ix].nX=iv+1+iAtomMed[iv]*epsdim;
//	                        }// if
//
//	                        if (dxyz.nY<rad2a)
//	                        {
//	                            iepsmp[iz][iy][ix].nY=iv+1+iAtomMed[iv]*epsdim;
//	                        }// if
//
//	                        if (dxyz.nZ<rad2a)
//	                        {
//	                            iepsmp[iz][iy][ix].nZ=iv+1+iAtomMed[iv]*epsdim;
//	                        }// if

							if(distsq<radp2) {
//								idebmap[ix][iy][iz] =false;
								den=exp(-(distsq/(sigma*sigma*radsq)));
								ginit_rhomap[ix][iy][iz]=1-(1-ginit_rhomap[ix][iy][iz])*(1-den);
							}


						}// do
					}// do
				}// do


//				for(iz=ismin.nZ; iz<=ismax.nZ; iz++)
//				{
//					for(iy=ismin.nY; iy<=ismax.nY; iy++)
//					{
//						for(ix=ismin.nX; ix<=ismax.nX; ix++)
//						{
//							ixyz=int_coord(ix,iy,iz);
//							dxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
//							distsq=optDot(dxyz,dxyz);
//							dxyz=dxyz+distsq;
////
////	                        if (dxyz.nX<rad2a)
////	                        {
////	                            iepsmp[iz][iy][ix].nX=iv+1+iAtomMed[iv]*epsdim;
////	                        }// if
////
////	                        if (dxyz.nY<rad2a)
////	                        {
////	                            iepsmp[iz][iy][ix].nY=iv+1+iAtomMed[iv]*epsdim;
////	                        }// if
////
////	                        if (dxyz.nZ<rad2a)
////	                        {
////	                            iepsmp[iz][iy][ix].nZ=iv+1+iAtomMed[iv]*epsdim;
////	                        }// if
//
//							if(distsq<radp2) {
////	                        	idebmap[iz][iy][ix] =false;
//								den=exp(-(distsq/(sigma*sigma*radsq)));
//								ginit_rhomap[iz][iy][ix]=1-(1-ginit_rhomap[iz][iy][ix])*(1-den);
//							}
//
//
//						}// do
//					}// do
//				}// do

			} //if


			else  /**faster method;*/
			{
				//cout << "####faster method:" << endl;
				rad2a=rad2-0.25;

				ixn=optRound(xn);

				fxn=optCast <delphi_real,delphi_integer> (ixn)-xn;
				rad2av=rad2a-fxn;


				for(ix=-lim; ix<=lim; ix++)
				{
					vtemp= double(ix)+fxn;
					//sq[ix]=vtemp*vtemp;
					//rad2aav[ix]=rad2a-vtemp;
					sqtemp[ix+15]=vtemp*vtemp;
					rad2aavtemp[ix+15]=rad2a-vtemp;

				}// do

				//adjust inter-atom, different epsilon bgps+++04/2004 Walter

				if (iNMedia>1&&bOnlyMol)
				{
					if(debug_space)cout << " Conv> #### multi media:" << endl;
					for(i=0; i<=ibox; i++)
					{

						i123=ioff[i];
						ixyz=ixn+i123;
						ix=ixyz.nX;
						iy=ixyz.nY;
						iz=ixyz.nZ;
						//distsq=sq[i123.nX].nX+sq[i123.nY].nY+sq[i123.nZ].nZ;
						distsq = sqtemp[i123.nX+15].nX +sqtemp[i123.nY+15].nY + sqtemp[i123.nZ+15].nZ;

//							if (distsq<rad2aav[i123.nX].nX)
//	                    	if (distsq<rad2aavtemp[i123.nX+15].nX)
//	                    	{
//									//iac=(iEpsMap[ix][iy][iz].nX % epsdim)-1;
//									iac=(iepsmp[ix][iy][iz].nX % epsdim)-1;
//
//									if (iac==-1||iac>iNatom)
//									{
//										//iEpsMap[ix][iy][iz].nX=iv+1+iAtomMed[iv]*epsdim;
//										iepsmp[ix][iy][iz].nX=iv+1+iAtomMed[iv]*epsdim;
//									}
//									else
//									{
//										//2011-05-14 Using operations on coord and int_coord type variables defined in module operators_on_coordinates
//										ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
//										ddxyz.nX= ddxyz.nX+0.5;
//										dis2min1=optDot(ddxyz,ddxyz)-rad2;
//
//										ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn2[iac];
//										ddxyz.nX= ddxyz.nX+0.5;
//										dis2min2=optDot(ddxyz,ddxyz)-pow((sDelPhiPDB[iac].radius*fScale),2);
//
//										if (dis2min2>dis2min1) iac=iv;
//										//iEpsMap[ix][iy][iz].nX=iac+1+iAtomMed[iac]*epsdim;
//										iepsmp[ix][iy][iz].nX=iac+1+iAtomMed[iac]*epsdim;
//									}// if
//								}// if
//
//								if (distsq<rad2aav[i123.nY].nY)
//	                    		if (distsq<rad2aavtemp[i123.nY+15].nY)
//								{
//									//iac=(iEpsMap[ix][iy][iz].nY % epsdim)-1;
//									iac=(iepsmp[ix][iy][iz].nY % epsdim)-1;
//									if (iac==-1||iac>iNatom)
//									{
//										//iEpsMap[ix][iy][iz].nY=iv+1+iAtomMed[iv]*epsdim;
//										iepsmp[ix][iy][iz].nY=iv+1+iAtomMed[iv]*epsdim;
//									}
//									else
//									{
//										//2011-05-14 Using operations on coord and int_coord type variables defined in module operators_on_coordinates
//										ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
//										ddxyz.nY= ddxyz.nY+0.5;
//										dis2min1=optDot(ddxyz,ddxyz)-rad2;
//
//										ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn2[iac];
//										ddxyz.nY= ddxyz.nY+0.5;
//										dis2min2=optDot(ddxyz,ddxyz)-pow((sDelPhiPDB[iac].radius*fScale),2);
//
//										if (dis2min2>dis2min1) iac=iv;
//										//iEpsMap[ix][iy][iz].nY=iac+1+iAtomMed[iac]*epsdim;
//										iepsmp[ix][iy][iz].nY=iac+1+iAtomMed[iac]*epsdim;
//									}// if
//								}// if
//
//								if (distsq<rad2aav[i123.nZ].nZ)
//	                    		if (distsq<rad2aavtemp[i123.nZ+15].nZ)
//								{
//									//iac=(iEpsMap[ix][iy][iz].nZ%epsdim)-1;
//									iac=(iepsmp[ix][iy][iz].nZ%epsdim)-1;
//									if (iac==-1||iac>iNatom)
//									{
//										//iEpsMap[ix][iy][iz].nZ=iv+1+iAtomMed[iv]*epsdim;
//										iepsmp[ix][iy][iz].nZ=iv+1+iAtomMed[iv]*epsdim;
//									}
//									else
//									{
//
//										ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn;
//										ddxyz.nZ= ddxyz.nZ+0.5;
//										dis2min1=optDot(ddxyz,ddxyz)-rad2;
//
//										ddxyz=optCast <delphi_real,delphi_integer> (ixyz)-xn2[iac];
//										ddxyz.nZ=ddxyz.nZ+0.5;
//										dis2min2=optDot(ddxyz,ddxyz)-pow((sDelPhiPDB[iac].radius*fScale),2);
//
//										if (dis2min2>dis2min1) iac=iv;
//										//iEpsMap[ix][iy][iz].nZ=iac+1+iAtomMed[iac]*epsdim;
//										iepsmp[ix][iy][iz].nZ=iac+1+iAtomMed[iac]*epsdim;
//									}// if
//								}// if

//								if(distsq<radp2) bDebMap[ix][iy][iz]=false;


								if(distsq<radp2) {
//									idebmap[ix][iy][iz]=false;
									den=exp(-(distsq/(sigma*sigma*radsq)));
									ginit_rhomap[ix][iy][iz]=1-(1-ginit_rhomap[ix][iy][iz])*(1-den);
								}

					}// do
				}
				else
				{
					if(debug_space) cout << " Conv> ### Fast method + single media:" << endl;

					//cout << "ix,rad2aav[ix]" << ix <<" " << rad2aav[ix] << endl;
					//cout << "#### ibox: " << ibox << endl;

					//cout << "### omp start: " << endl;
					//omp_set_num_threads(4);
					//#pragma omp parallel for schedule(static) private(xn,rad,radp,rad2,radp2,itest2,ismin,ismax,itest,rad2a,ixn,fxn,rad2av,vtemp,rad2aavtemp,i123,ixyz,ix,iy,iz,distsq,sqtemp)
					//#pragma omp parallel for schedule(static) private(xn,rad,radp,rad2,radp2,itest2,ismin,ismax,itest,rad2a,ixn,fxn,rad2av,vtemp,rad2aavtemp,i123,ixyz,ix,iy,iz,distsq,sqtemp) shared(idebmap)

#ifdef PARALLEL_OMP
			#pragma omp for schedule(auto)
			//#pragma omp parallel for schedule(static) private(xn,rad,radp,rad2,radp2,itest2,ismin,ismax,itest,rad2a,ixn,fxn,rad2av,vtemp,rad2aavtemp,i123,ixyz,ix,iy,iz,distsq,sqtemp)
#endif

					for(i=0; i<=ibox; i++)
					{

						i123=ioff[i];
						ixyz=ixn+i123;
						ix=ixyz.nX;
						iy=ixyz.nY;
						iz=ixyz.nZ;
						distsq = sq[i123.nX].nX +sq[i123.nY].nY + sq[i123.nZ].nZ;
//						cout << "i123,ixn,ix,iy,iz,distsq: " << i123 << " " <<ixn << " " <<ix << " " <<iy << " " <<iz << " " <<distsq << endl;
//	                    if (distsq<rad2aav[i123.nX].nX)
//	                    {
//	                        //Lin Li: all indexes -1, because of C++ arrays' indexes start from 0;
//	                        //Lin Li: because iv start from 0, so +2:
//	                        //iEpsMap[ix-1][iy-1][iz-1].nX=iv+2+iAtomMed[iv]*epsdim;
//	                        iepsmp[iz][iy][ix].nX=iv+1+iAtomMed[iv]*epsdim;
//
//	                        //cout << "ix,iy,iz,iEpsMap[ix-1][iy-1][iz-1].nX: " << ix<< " " << iy << " " << iz << " " <<iEpsMap[ix-1][iy-1][iz-1].nX << endl;
//	                        //cout << "ix,iy,iz,iv,iAtomMed[iv],epsdim " << ix<< " " << iy << " " << iz << " " << iv << " " << iAtomMed[iv]<< " " << epsdim << endl;
//	                    }// if
//
//
//	                    if (distsq<rad2aav[i123.nY].nY)
//	                    {
//	                        //iEpsMap[ix-1][iy-1][iz-1].nY=iv+2+iAtomMed[iv]*epsdim;
//	                        iepsmp[iz][iy][ix].nY=iv+1+iAtomMed[iv]*epsdim;
//	                    }// if
//
//	                    if (distsq<rad2aav[i123.nZ].nZ)
//	                    {
//	                        //iEpsMap[ix-1][iy-1][iz-1].nZ=iv+2+iAtomMed[iv]*epsdim;
//	                        iepsmp[iz][iy][ix].nZ=iv+1+iAtomMed[iv]*epsdim;
//	                    }// if



							if (distsq<radp2) {
//		                    	idebmap[ix][iy][iz]=false;
								den=exp(-(distsq/(sigma*sigma*radsq)));
								ginit_rhomap[ix][iy][iz]=1-(1-ginit_rhomap[ix][iy][iz])*(1-den);
							}

					}// do


				}// if
			}// if
		}// do DoATOMS //end do of atoms;

		if (ioff != NULL)    delete [] ioff;

		/* By this point GINIT_RHOMAP has been created.
		 * Now using the calculation Argo made, for the original f=0.5 to place HS step, the density should exp (-9/16) = 0.56978 for fRw from Rvdw
		 * Hence if ginit > 0.56978 => HRhomap = 1
		 * else HRhomap = 0
		 */

		for ( i = 1; i <= iGrid; i++)
		{
			for (j = 1; j <= iGrid; j++)
			{
				for (k = 1; k <= iGrid; k++)
				{
					if (ginit_rhomap[i][j][k] > rho_hs) {
						HRhomap[i][j][k] = 0.0;
					} else {
						HRhomap[i][j][k] = 1.0;
					}

				}
			}
		}


		if (hs_rhomap)
    //if (false)
		{
			//ofstream ginit_out("ginit_rhomap.dat");
			ofstream hs_out("hs_rhomap.dat");



			for (j = 1; j <= iGrid; j++) {
				for (k = 1; k <= iGrid; k++) {

					//ginit_out <<  j << "\t" << k << "\t"  << ginit_rhomap[fx][j][k] << endl;
					hs_out <<  j << "\t" << k << "\t"  << HRhomap[fx][j][k] << endl;

				}
				//ginit_out << " " << endl;
				hs_out << " " << endl;
			}

			//ginit_out.close();
			hs_out.close();
		}


		/*
		 * HRhomap has been created.
		 * Will proceed to FFT/iFFT
		 */

		Convolute();

    if (conv_rhomap)
    //if (false)
		{
			std::string convFileName = "conv_rhomap" + std::to_string(inhomo) + ".dat";
			ofstream conv_out(convFileName);

			for (j = 1; j <= iGrid; j++) {
				for (k = 1; k <= iGrid; k++) {

					conv_out << j << "\t" << k << "\t"  << HRhomap[fx][j][k] << endl;

				}
				conv_out << " " << endl;
			}

			conv_out.close();
		}


	}// if !(iConvolute!=0 && inhomo=0 && logs)

	/*
	 * FFT/iFFT is done on HRhomap
	 * HRhomap is convoluted.
	 * Now assign EPS values using the formula
	 * E = E[in](1-Rho[out]) + E[out]Rho[out].
	 * From all the above codes, we have the convoluted HRhomap.
	 * So Rho[out] is with us.
	 * Depending on the medium now, EPS values will be assigned.
	 */

	if ( inhomo == 1 )
	{
		eps_medium = 1;
    	// In gaussian, when INHOMO == 1, the E_out = INDI
    	//eps_medium = repsin;
	}
	else if ( inhomo == 0 )
	{
		eps_medium = repsout;
	}

	cout << " Conv> Outside dielectric is set to " << eps_medium << endl;
  	if ( debug_space ) cout << " INHOMO = " << inhomo << endl;

		sgrid_rho_real.nX = eps_medium;
		sgrid_rho_real.nY = eps_medium;
		sgrid_rho_real.nZ = eps_medium;



	//assigning dielectric constants based on Rho[out] values
		for ( i = 1; i <= iGrid; i++)
		{
			for (j = 1; j <= iGrid; j++)
			{
				for (k = 1; k <= iGrid; k++)
				{
				  cepsmap[i][j][k] = repsin*(1 - HRhomap[i][j][k]) + eps_medium*HRhomap[i][j][k];
				  vct_cepsmap.push_back(cepsmap[i][j][k]);
				}
			}
		}



    /*
     * Begin Linear Interpolation
     * For now the method is very basic. Looping is also simplistic
     * Try to get a memory efficient way to do the interpolation.
     */

	for ( i = 1; i < iGrid; i++)
	{
		for (j = 1; j < iGrid; j++)
		{
			for (k = 1; k < iGrid; k++)
			{

				gepsmp2[i][j][k].nX = 0.5*(cepsmap[i][j][k] + cepsmap[i+1][j][k]);
				gepsmp2[i][j][k].nY = 0.5*(cepsmap[i][j][k] + cepsmap[i][j+1][k]);
				gepsmp2[i][j][k].nZ = 0.5*(cepsmap[i][j][k] + cepsmap[i][j][k+1]);

			}
		}
	}

	// The following part wasn't necessary given gepsmp2 was initialised with repsout
	for ( j = 1; j <= iGrid; j++ )
	{
		for (k = 1; k <= iGrid; k++)
		{
			gepsmp2[iGrid][j][k] = sgrid_rho_real;
		}
	}

	for ( i = 1; i <= iGrid; i++ )
	{
		for (k = 1; k <= iGrid; k++)
		{
			gepsmp2[i][iGrid][k] = sgrid_rho_real;
		}
	}
	for ( i = 1; i <= iGrid; i++ )
	{
		for (j = 1; j <= iGrid; j++)
		{
			gepsmp2[i][j][iGrid] = sgrid_rho_real;
		}
	}


// ARGO
/////////////// FROM GAUSSIAN - to set the BNDYGRIDPOINTS ////////////////
  iBoundNum=0;
  longint=0;
  cout << "Setting BNDY POINTS for Convolute module ... " << endl;
  for(i=2; i<iGrid; i++)
  {
      //for(j=2; j<=iGrid; j++)
      for(j=2; j<iGrid; j++)
      {
          //for(k=2; k<=iGrid; k++)
          for(k=2; k<iGrid; k++)
          {
              if(HRhomap[i][j][k] < (1 - dentemp))
              {

                  idebmap[k][j][i]=false; //Like iGaussian's changed method of generating bDebMap;

                  //see if the grid-point is is close to the region with eps = repsin
//                  eps_front = gepsmp2[i-1][j][k].nX;
//                  eps_back  = gepsmp2[i][j][k].nX;
//
//                  eps_left  = gepsmp2[i][j-1][k].nY;
//                  eps_right = gepsmp2[i][j][k].nY;
//
//                  eps_bott  = gepsmp2[i][j][k-1].nZ;
//                  eps_top   = gepsmp2[i][j][k].nZ;
//				
//				  neighbor_eps.push_back(eps_front);	
//				  neighbor_eps.push_back(eps_back);	
//				  neighbor_eps.push_back(eps_left);	
//				  neighbor_eps.push_back(eps_right);	
//				  neighbor_eps.push_back(eps_top);	
//				  neighbor_eps.push_back(eps_bott);	
//                  
//				  sort(neighbor_eps.begin(),neighbor_eps.end());
//				  if (neighbor_eps.back()/neighbor_eps.front() > eps_diff )
//				  {
//                  	longint += 1;
//				  }
//				  
//				  neighbor_eps.resize(0);
				   //see if the epslon distribution is flat or not
                   eps_min=HRhomap[i][j][k];
                   eps_max=HRhomap[i][j][k];
                  
                   eps_min=min(eps_min,HRhomap[i-1][j][k]);
                   eps_min=min(eps_min,HRhomap[i][j-1][k]);
                   eps_min=min(eps_min,HRhomap[i][j][k-1]);
                  
                   eps_max=max(eps_max,HRhomap[i-1][j][k]);
                   eps_max=max(eps_max,HRhomap[i][j-1][k]);
                   eps_max=max(eps_max,HRhomap[i][j][k-1]);
                  
                   if(eps_max/eps_min > eps_diff) longint=longint+1;

              }//if
          }//do
      }//do
  }//do

  iBoundNum=longint;
  //allocate(ibgrd[iBoundNum]);
  if(debug_space) cout << "###### iboundnum: " <<iBoundNum << endl;
  if(debug_space)cout << "inhomo: " << inhomo << endl;
  ibgrd_v.assign(iBoundNum, sgrid_temp_int);
  ibgrd=&ibgrd_v[-1];
  
  cout << " IBOUNDNUM FOR CONVOLUTE = " << iBoundNum << endl;
  n=0;
  //for(i=2; i<=iGrid; i++)
  for(i=2; i<iGrid; i++)
  {
      //for(j=2; j<=iGrid; j++)
      for(j=2; j<iGrid; j++)
      {
          //for(k=2; k<=iGrid; k++)
          for(k=2; k<iGrid; k++)
          {

              if(HRhomap[i][j][k] < (1 - dentemp))
              {
//                  //see if the grid-point is is close to the region with eps = repsin
//                  eps_front = gepsmp2[i-1][j][k].nX;
//                  eps_back  = gepsmp2[i][j][k].nX;
//
//                  eps_left  = gepsmp2[i][j-1][k].nY;
//                  eps_right = gepsmp2[i][j][k].nY;
//
//                  eps_bott  = gepsmp2[i][j][k-1].nZ;
//                  eps_top   = gepsmp2[i][j][k].nZ;
//				  
//				  neighbor_eps.push_back(eps_front);	
//				  neighbor_eps.push_back(eps_back);	
//				  neighbor_eps.push_back(eps_left);	
//				  neighbor_eps.push_back(eps_right);	
//				  neighbor_eps.push_back(eps_top);	
//				  neighbor_eps.push_back(eps_bott);	
//				 
//				  sort(neighbor_eps.begin(),neighbor_eps.end()); 
//				  if (neighbor_eps.back()/neighbor_eps.front() > eps_diff )
//                  {
//                    n=n+1;
//                    ibgrd[n].nX=i;
//                    ibgrd[n].nY=j;
//                    ibgrd[n].nZ=k;
//                  }
//
//				  neighbor_eps.resize(0);

                   //see if the epslon distribution is flat or not
                   eps_min=HRhomap[i][j][k];
                   eps_max=HRhomap[i][j][k];
                  
                   eps_min=min(eps_min,HRhomap[i-1][j][k]);
                   eps_min=min(eps_min,HRhomap[i][j-1][k]);
                   eps_min=min(eps_min,HRhomap[i][j][k-1]);
                  
                   eps_max=max(eps_max,HRhomap[i-1][j][k]);
                   eps_max=max(eps_max,HRhomap[i][j-1][k]);
                   eps_max=max(eps_max,HRhomap[i][j][k-1]);
                  
                   if(eps_max/eps_min > eps_diff)
                   {
                  
                       n=n+1;
                       ibgrd[n].nX=i;
                       ibgrd[n].nY=j;
                       ibgrd[n].nZ=k;
                   }

              }//if
          }//do
      }//do
  }//do


  if (inhomo == 0 && bEpsOut)
	{

    unique_ptr<CIO> pio(new CIO()); // smart unique_ptr
    // pio->writeEpsMap(iAtomNum,iObjectNum,iGrid,fScale,fgBoxCenter,prgigEpsMap,prgbDielecMap,strEpsFile);
    pio->writeConvEpsMap(iGrid,fScale,fgBoxCenter,vct_cepsmap,strEpsFile);
    pio.reset();
    vct_cepsmap.resize(0);
	}



#ifdef VERBOSE
    cout <<" Conv> Finished creating Convoluted Epsilon Map " << endl;
#endif
    if ( debug_space ) cout <<"############# Done with SetConvolute part ############\n" << endl;



    debug_space=false;

    return;

}// void setout;
