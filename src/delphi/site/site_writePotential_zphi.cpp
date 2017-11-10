/*
 * site_writePotential_zphi.cpp
 *
 *  Created on: Feb 17, 2016
 *      Author: Argo
 */

#include "site.h"
#include "math.h"
#include <vector>
#include <complex>

typedef vector<delphi_real> potentials;
typedef vector<potentials> cluster;

void CSite::writeZetaPhiFile()
{

   size_t MAXWIDTH = 45;
   cout << left << setw(MAXWIDTH) << " Surface> Writing Surface Potential values now to " << strZetaPhiFile << endl;


    string tempFileName = "temp." + strZetaPhiFile;
    ofstream zetaPhi_file;
    zetaPhi_file.open(strZetaPhiFile);

    int i,j,k;			//ARGO the grid indices from the zetaSite-file written in the space module
    delphi_real px,py,pz;	//AGO the real-time coords for the zetaSite-file coordinates from space module


    if ( bZetaPhiOut )
    {
	    potentials surfPots;

	    std::ifstream infile(("temp." + strZetaPhiFile).c_str());
	    while ( infile >> px >> py >> pz >> i >> j >> k )
	    {
				surfPots.push_back(phimap[k][j][i]);
	    }

			infile.close();

	    //cout << " Surf> Read " << surfPots.size() << " lines" << endl;


	    /* KSplit method starts here
	     *
	     *
	     *
	    long int N = surfPots.size();
	    int k = kclusters;
	    cout << " Surf> Distributing the surface potentials in " << k << " clusters" << endl;

	    cluster kSplit;
	    potentials B;
	    potentials head_potentials;
	    delphi_real artificial_rnd = 23456;//FOR ERASING ELEMENTS LATER

	    //creating the first CLUSTER B1
	    for ( potentials::iterator itp = surfPots.begin(); itp != surfPots.end(); itp++)
	    {
	        B.push_back(*itp);
	    }
	    kSplit.push_back(B);

	    head_potentials.push_back(B[0]);

	    for ( int l = 1; l < k ; l++)
	    {
	        delphi_real h = 0,h1;
	        delphi_real next_head = 0.0;
	        int ll = 0;

	        while (ll < l) {
	            potentials present_B = kSplit[ll];
	            cout << " Size of present B = " << present_B.size() << endl;
	            delphi_real present_head = head_potentials[ll];
	            potentials::iterator itB = present_B.begin();
	             while ( itB != present_B.end())
	             {
	                 h1 = abs( (*itB) - present_head);
	                 if ( h1 > h)
	                 {
	                     h = h1;
	                     next_head = *itB;
	                 }

	                 itB++;
	            }
	            ll++;
	        }
	        cout << " Max  h = " << h <<  endl;
	        cout << " Number of clusters at 1 = " << kSplit.size() << endl;

	        potentials next_B;
	        next_B.push_back(next_head);
	        cout << " Next head = " << next_head << endl;
	        head_potentials.push_back(next_B[0]);

	        //time to move elements to the newly created cluster
	        ll = 0;
	        while (ll < l)
	        {
	            //cout << "Rellocating elements to new cluster" << endl;
	            potentials present_B = kSplit[ll];
	            potentials::iterator itB = present_B.begin();
	            delphi_real present_head = head_potentials[ll];

	            while ( itB != present_B.end() )
	            {
	                if ( (abs(*itB - present_head) > abs(*itB - next_head)) )
	                {
	                    next_B.push_back(*itB); //rellocation to the new cluster!!! OOOHHH!!!
	                    *itB = artificial_rnd;
	                }
	                itB++;
	            }

	            //now remove elements that are ARTIFICIAL_RND
	            potentials temp_B;
	            itB = present_B.begin();
	            while ( itB != present_B.end())
	            {
	                if ( *itB != artificial_rnd )
	                {
	                    temp_B.push_back(*itB);
	                }
	                itB++;
	            }

	            present_B.resize(0);
	            kSplit[ll] = temp_B;
	            temp_B.resize(0);

	            ll++;
	        }

	        //erase the first element of next_B as it is redundant
	        next_B.erase(next_B.begin());

	        kSplit.push_back(next_B);
	        cout << " Number of clusters at 2 = " << kSplit.size() << endl;
	        cout << "--------------------------------------------" << endl;

	    }

	    // ofstream cluster_OUT;
	    // cluster_OUT.open("/Users/arghyachakravorty/Desktop/zeta_potential_project/spheres/cluster.dat");

	    delphi_real avg_SurfPot = 0;
	    delphi_real simple_avg_SurfPot = 0;
	    simple_avg_SurfPot = meanPotential(surfPots);

	    for ( int l = 0; l < k; l++ )
	    {
	        potentials present_B = kSplit[l];
	        long int B_size = present_B.size();

	        //cout << (double)B_size << endl;
	        avg_SurfPot += (B_size * meanPotential(present_B))/N;

	        // potentials::iterator itB2 = present_B.begin();
					//
	        // while (itB2 != present_B.end())
	        // {
	        //     cluster_OUT << head_potentials[l] << "\t" << *itB2 << "\t" << l+1 << endl;
	        //     itB2++;
	        // }
	        // cluster_OUT << " " << endl;
	    }
	    //cluster_OUT.close();

	    cout << " Surf> Weighted Average Potential = " << avg_SurfPot << " kT/e" << endl;
	    cout << " Surf> Simple Average Potential   = " << simple_avg_SurfPot << " kT/e" << endl;
	    cout << "--------------------------------------------" << endl;
	    *
	    *
	    *
	    */

	    delphi_real simple_avg_SurfPot = meanPotential(surfPots);
	    cout << left << setw(MAXWIDTH) << " Surface> Average of the Potentials on the surface" << " : " << simple_avg_SurfPot << " kT/e" << endl;

	    infile.open(("temp." + strZetaPhiFile).c_str());
	    //zetaPhi_file << "#REMARK WEIGHTED AVERAGE SURFACE POTENTIAL = " << avg_SurfPot << " kT/e" << endl;
	    zetaPhi_file << "#REMARK SIMPLE   AVERAGE SURFACE POTENTIAL = " << simple_avg_SurfPot << " kT/e" << endl;
	    //zetaPhi_file << "#REMARK NUMBER OF K-CLUSTERS      = " << k << endl;

	    while ( infile >> px >> py >> pz >> i >> j >> k )
	    {
	   	zetaPhi_file << setw(10) << setprecision(6) << px
	   		     << setw(10) << setprecision(6) << py
	   		     << setw(10) << setprecision(6) << pz
	   		     << setw(15) << setprecision(6) << phimap[k][j][i]
	   		     << endl;
	    }

	    zetaPhi_file.close();
   }
    remove(("temp." + strZetaPhiFile).c_str());

    //Time to remove zeta_sites.dat file. Not needed anyway hereafter

}

delphi_real CSite::meanPotential ( potentials& pots )
{
	delphi_real mn = 0;
	potentials::iterator itpm = pots.begin();

	while ( itpm != pots.end() )
	{
		mn += *itpm;
		itpm++;
	}

	mn /= pots.size();
	return mn;
}
