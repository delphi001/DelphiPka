//
//  topology_import.cpp
//  DelPhiPKA
//
//  Created by Lin Wang on 10/22/14.
//  Copyright (c) 2014 Lin. All rights reserved.
//

#include "topology_import.h"
#include <fstream>
#include <string>


void CTologyImport::string_to_topology(const string& line, TOPOLOGY& topology_temp)
{
    topology_temp.resname = line.substr(3,3);
    topology_temp.atom = line.substr(7,4);
    topology_temp.orbital = line.substr(12,5);
    topology_temp.conf = line.substr(18,4);
    
    remove_all_whitespace_string(topology_temp.resname);
    remove_all_whitespace_string(topology_temp.atom);
    remove_all_whitespace_string(topology_temp.orbital);
    remove_all_whitespace_string(topology_temp.conf);
    
    if(line.length()>24 && !line.substr(23,4).empty()) {
        topology_temp.bondatom1 = line.substr(23,4);
        remove_all_whitespace_string(topology_temp.bondatom1);
    }
    else topology_temp.bondatom1.clear();
    
    if(line.length()>29 && !line.substr(28,4).empty()) {
        topology_temp.bondatom2 = line.substr(28,4);
        remove_all_whitespace_string(topology_temp.bondatom2);
    }
    else topology_temp.bondatom2.clear();
    
    if(line.length()>34 && !line.substr(33,4).empty()) {
        topology_temp.bondatom3 = line.substr(33,4);
        remove_all_whitespace_string(topology_temp.bondatom3);
    }
    else topology_temp.bondatom3.clear();
    
    if(line.length()>39 && !line.substr(38,4).empty()) {
        topology_temp.bondatom4 = line.substr(38,4);
        remove_all_whitespace_string(topology_temp.bondatom4);
    }
    else topology_temp.bondatom4.clear();
}

void CTologyImport::topology_import()
{
    
    if(HofGLU!= "OE1" && HofGLU!= "OE2")    HofGLU = "OE2";
    if(HofASP!= "OD1" && HofASP!= "OD2")    HofASP = "OD2";
    
    ifstream param(parameter_input);
    
    
    
    if(!param)
    {
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
            
            cout << "[ERROR] The topology file can not be found. " << endl;
            
#ifdef MPI_PARALLEL
        }
        MPI::Finalize();
#endif
        
        exit(0);
    }
    
#ifdef MPI_PARALLEL
    if(id==0) {
#endif
        cout << "Starting import the topology parameters from " << parameter_input << "  ..." << endl;

#ifdef MPI_PARALLEL
    }
#endif
    
    while(getline(param, line))
    {
        if(line.substr(0,1) == "$") {
            
            
        //------------------- For choose the Hydrogen acceptor ATOM in GLU and ASP residue ---------------//
        //------------------- OD1/OD2 of ASP, OE1/OE2 of GLU. Choose one of them -------------------------//
        //------------------- The corresponding toplogy information needs modification -------------------//
            
            if(HofASP == "OD1") {
                if(line == "$  ASP  OD2  sp2   SD   CG   HD2") {
                    line = "$  ASP  OD2  sp2   SD   CG";
                }
                if(line == "$  ASP  HD2  s     SD   OD2")   continue;
            }
            else if (HofASP == "OD2") {
                if(line == "$  ASP  OD1  sp2   SD   CG   HD1") {
                    line = "$  ASP  OD1  sp2   SD   CG";
                }
                if(line == "$  ASP  HD1  s     SD   OD1")   continue;
                
            }
            
            if(HofGLU == "OE1") {
                if(line == "$  GLU  OE2  sp2   SD   CD   HE2") {
                    line = "$  GLU  OE2  sp2   SD   CD";
                }
                if(line == "$  GLU  HE2  s     SD   OE2")   continue;
            }
            else if (HofGLU == "OE2") {
                if(line == "$  GLU  OE1  sp2   SD   CD   HE1") {
                    line = "$  GLU  OE1  sp2   SD   CD";
                }
                if(line == "$  GLU  HE1  s     SD   OE1")   continue;
            }
        //-------------------------------------------------------------------------------------------------//
            
            string_to_topology(line, topology_temp);
            key = topology_temp.resname + " " + topology_temp.atom;
            topology_map.insert (pair<string, TOPOLOGY>(key, topology_temp));
            //	cout << topology_temp << endl;  // Print a complete topology_map line by line on screen
        }
        
        else if(line.substr(0,1) == "*") {
            key = line.substr(7,5);
            remove_leadtail_whitespace_string(key);
            
            linetmp = line.substr(13,7);
            remove_leadtail_whitespace_string(linetmp);
            val = stof(linetmp);
            
            pkamap.insert({key,val});
        }
        
        else    continue;
    }

#ifdef MPI_PARALLEL
    if(id==0) {
#endif
        cout << " topology parameter imported successfully..." << endl;
        cout << endl;
        
#ifdef MPI_PARALLEL
    }
#endif

}


void CTologyImport::chargeinfo_import()
{
    
    ifstream crg_info(crg_file);
    
    if(!crg_info)
    {
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
            
            cout << "[ERROR] The charge file can not be found. " << endl;
            
#ifdef MPI_PARALLEL
        }
        MPI::Finalize();
#endif
        
        exit(0);
    }
    
#ifdef MPI_PARALLEL
    if(id==0) {
#endif
    
        cout << "Starting import the charge parameters from " << crg_file << "  ..." << endl;

        
#ifdef MPI_PARALLEL
    }
#endif
    
    bool bCorrectCrgFile = false;
    
    while(getline(crg_info, line))
    {

        remove_leadtail_whitespace_string(line);

        if(line.substr(0,1)!="!" && !line.empty()) {
            if(line == "atom__resnumbc_charge_") {
                
                bCorrectCrgFile = true;
                continue;
            }
            else if (line.substr(4,2)=="  " && line.substr(14,1)==" ") {
                string resname = line.substr(6,3);
                remove_all_whitespace_string(resname);

                string resatm = line.substr(0,4);
                remove_all_whitespace_string(resatm);
                
                key = resname + " " + resatm;
                val = stof(line.substr(15,7));

                crgmap.insert({key,val});
            }
            
        }
        
    }

    if (bCorrectCrgFile) {
        
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
        
            cout << " charge parameter imported successfully..." << endl;
            cout << endl;
            
#ifdef MPI_PARALLEL
        }
#endif
        
    }
    else {
        
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
            
            cout << " Incorrect charge file. Please check the charge file again...";
            cout << endl;
            
#ifdef MPI_PARALLEL
        }
        MPI::Finalize();
#endif
        exit(0);
        
    }
}

void CTologyImport::sizeinfo_import()
{
    ifstream siz_info(siz_file);
    
    if(!siz_info)
    {
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
            
            cout << "[ERROR] The size file can not be found. " << endl;
            
#ifdef MPI_PARALLEL
        }
        MPI::Finalize();
#endif
        
        exit(0);
    }

#ifdef MPI_PARALLEL
    if(id==0) {
#endif
        
        cout << "Starting import the size parameters from " << siz_file << "  ..." << endl;

#ifdef MPI_PARALLEL
    }
#endif
    
    bool bCorrectSizFile = 0;
    
    while(getline(siz_info, line))
    {
        remove_leadtail_whitespace_string(line);
        
        if(line.substr(0,1)!="!" && !line.empty()) {
            if(line == "atom__res_radius_") {
                bCorrectSizFile = 1;
                continue;
            }
            else if (line.substr(4,2)=="  " && line.substr(9,1)==" ") {
                string resname = line.substr(6,3);
                remove_all_whitespace_string(resname);
                
                string resatm = line.substr(0,4);
                remove_all_whitespace_string(resatm);
                
                key = resname + " " + resatm;
                val = stof(line.substr(12,5));
                
                sizmap.insert({key,val});
            }
            
        }
    }
    if (bCorrectSizFile) {
        
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
            
            cout << " size parameter imported successfully..." << endl;
            cout << endl;
            
#ifdef MPI_PARALLEL
        }
#endif
        
    }
    else {
        
#ifdef MPI_PARALLEL
        if(id==0) {
#endif
            
            cout << " Incorrect size file. Please check the size file again...";
            cout << endl;
            
#ifdef MPI_PARALLEL
        }
        MPI::Finalize();
#endif
        exit(0);
        
    }
}

void CTologyImport::run() {
    
#ifdef MPI_PARALLEL
    id = MPI::COMM_WORLD.Get_rank();
#endif
    
    ///////////++++++++++++  Import the topology file  ++++++++++++/////////////

    topology_import();
 
    ///////////++++++++++++  Import the charge file  ++++++++++++/////////////
   
    chargeinfo_import();
    
    ///////////++++++++++++  Import the size file  ++++++++++++/////////////

    sizeinfo_import();

#ifdef MPI_PARALLEL
    if(id==0) {
#endif
//        cout << "  HofASP = " << HofASP << "   HofGLU = " << HofGLU << endl;
        cout << "Done." << endl;
        cout << endl;

#ifdef MPI_PARALLEL
    }
#endif
    
}