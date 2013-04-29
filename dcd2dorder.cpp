//    MOLDY Tools 
//
//    Copyright (C) 2007,2008 Nava Whiteford, Justine Taylor
//
//
//    This file is part of MOLDY Tools.
//
//    MOLDY Tools is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    MOLDY Tools is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with MOLDY Tools.  If not, see <http://www.gnu.org/licenses/>.

#include <stdlib.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "vec.h"
#include "pdb.h"
#include "dcd.h"
#include "stringify.h"
#include "box.h"
#include "dorder.h"

using namespace std;

int main(int argc, char** argv)
{

  if(argc<4) {
     cout << " " << endl; 
     cout << "MOLDY Tools  Copyright (C) 2007,2008  Nava Whiteford, Justine Taylor" << endl;
     cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
     cout << "This is free software, and you are welcome to redistribute it" << endl;
     cout << "under certain conditions. See GPL.txt" << endl;
     cout << " " << endl;
     cout << "Usage: dcd2dorder <dcdfile> <pdbfile> <dictionary file>" << endl;
     cout << " " << endl;
     cout << "This program calculates the Deuterium order parameters down the lipid tails." << endl;
     cout << "The dictionary file defines which atoms are used for the calculation and is defined" << endl;
     cout << "fully in the README file." << endl;
     cout << " " << endl;
     exit(0);
  }

  DCDFile mydcd(argv[1]);
  PDBFile mypdb(argv[2]);
  string dictionary(argv[3]);

  map<string,double> sum_total_sda_sn1;
  map<string,double> sumsq_total_sda_sn1;
  map<string,int>    count_total_sda_sn1;

  map<string,double> sum_total_sda_sn2;
  map<string,double> sumsq_total_sda_sn2;
  map<string,int>    count_total_sda_sn2;
  
  map<string,double> sum_total_sdb_sn1;
  map<string,double> sumsq_total_sdb_sn1;
  map<string,int>    count_total_sdb_sn1;

  map<string,double> sum_total_sdb_sn2;
  map<string,double> sumsq_total_sdb_sn2;
  map<string,int>    count_total_sdb_sn2;
  
  map<string,double> sum_total_sdc_sn1;
  map<string,double> sumsq_total_sdc_sn1;
  map<string,int>    count_total_sdc_sn1;

  map<string,double> sum_total_sdc_sn2;
  map<string,double> sumsq_total_sdc_sn2;
  map<string,int>    count_total_sdc_sn2;

  double slice_length=0;
  bool first_frame = true;

  vector<int>    resnum_first;
  vector<string> sn1_c_names_first;
  vector<string> sn2_c_names_first;

  for(;mydcd.next_frame();) {
    DCDFrame f = mydcd.get_frame();

    DOrder dorder(f,mypdb,dictionary);
    dorder.dorder_calc();
    
    
    const vector<int> &resnums = dorder.get_residue_numbers();
    if(first_frame == true) {
      // I need some information from dorder, so I grab it from the first frame here
      resnum_first =  dorder.get_residue_numbers();
      sn1_c_names_first = dorder.get_sn1_c_names();
      sn2_c_names_first = dorder.get_sn2_c_names();

      first_frame=false;
    }

    vector<int>::const_iterator i_r=resnums.begin();

    cout << "FRAME ";
    cout << resnums.size() << endl;
    for(;i_r != resnums.end();i_r++) {
      vector<string>::const_iterator names  = dorder.get_sn1_c_names().begin();
      map<string, double > &angles_sn1_sda  = dorder.get_angles_sn1_sda((*i_r));
      map<string, double > &angles_sn1_sdb  = dorder.get_angles_sn1_sdb((*i_r));
      map<string, double > &angles_sn1_sdc  = dorder.get_angles_sn1_sdc((*i_r));
      
      for(;names != dorder.get_sn1_c_names().end();names++) {
        cout << (*i_r) << " SN1 " << (*names) << " ";
        cout << angles_sn1_sda[(*names)] << " ";
        cout << angles_sn1_sdb[(*names)] << " ";
        cout << angles_sn1_sdc[(*names)] << endl;
      
        if(angles_sn1_sda[(*names)] != NAN) { sumsq_total_sda_sn1[(*names)] += angles_sn1_sda[(*names)]*angles_sn1_sda[(*names)];  sum_total_sda_sn1[(*names)] += angles_sn1_sda[(*names)]; count_total_sda_sn1[(*names)]++;}
        if(angles_sn1_sdb[(*names)] != NAN) { sumsq_total_sdb_sn1[(*names)] += angles_sn1_sdb[(*names)]*angles_sn1_sdb[(*names)];  sum_total_sdb_sn1[(*names)] += angles_sn1_sdb[(*names)]; count_total_sdb_sn1[(*names)]++;}
        if(angles_sn1_sdc[(*names)] != NAN) { sumsq_total_sdc_sn1[(*names)] += angles_sn1_sdc[(*names)]*angles_sn1_sdc[(*names)];  sum_total_sdc_sn1[(*names)] += angles_sn1_sdc[(*names)]; count_total_sdc_sn1[(*names)]++;}
	
      }

      names = dorder.get_sn2_c_names().begin();
      map<string, double > &angles_sn2_sda  = dorder.get_angles_sn2_sda((*i_r));
      map<string, double > &angles_sn2_sdb  = dorder.get_angles_sn2_sdb((*i_r));
      map<string, double > &angles_sn2_sdc  = dorder.get_angles_sn2_sdc((*i_r));
      
      for(;names != dorder.get_sn2_c_names().end();names++) {
        cout << (*i_r) << " SN2 " << (*names) << " ";
	cout << angles_sn2_sda[(*names)] << " ";
        cout << angles_sn2_sdb[(*names)] << " ";
        cout << angles_sn2_sdc[(*names)] << endl;
       

        if(angles_sn2_sda[(*names)] != NAN) { sumsq_total_sda_sn2[(*names)] += angles_sn2_sda[(*names)]*angles_sn2_sda[(*names)];  sum_total_sda_sn2[(*names)] += angles_sn2_sda[(*names)]; count_total_sda_sn2[(*names)]++;}
        if(angles_sn2_sdb[(*names)] != NAN) { sumsq_total_sdb_sn2[(*names)] += angles_sn2_sdb[(*names)]*angles_sn2_sdb[(*names)];  sum_total_sdb_sn2[(*names)] += angles_sn2_sdb[(*names)]; count_total_sdb_sn2[(*names)]++;}
        if(angles_sn2_sdc[(*names)] != NAN) { sumsq_total_sdc_sn2[(*names)] += angles_sn2_sdc[(*names)]*angles_sn2_sdc[(*names)];  sum_total_sdc_sn2[(*names)] += angles_sn2_sdc[(*names)]; count_total_sdc_sn2[(*names)]++;}

      }
    }

  }

  // System averages and standard deviation
  
  cout << "FRAME AVERAGE" << endl;
  vector<string>::const_iterator names  = sn1_c_names_first.begin();
  for(;names !=sn1_c_names_first.end();names++) {
    cout << "SN1 " << (*names) << " ";
    
    if(count_total_sda_sn1[(*names)] != 0) cout << sum_total_sda_sn1[(*names)]/count_total_sda_sn1[(*names)] << " "; else cout << "NAN ";
    if(count_total_sdb_sn1[(*names)] != 0) cout << sum_total_sdb_sn1[(*names)]/count_total_sdb_sn1[(*names)] << " "; else cout << "NAN ";
    if(count_total_sdc_sn1[(*names)] != 0) cout << sum_total_sdc_sn1[(*names)]/count_total_sdc_sn1[(*names)] << endl; else cout << "NAN" << endl;

  }

  names  = sn2_c_names_first.begin();
  for(;names != sn2_c_names_first.end();names++) {
    cout << "SN2 " << (*names) << " ";

    if(count_total_sda_sn2[(*names)] != 0) cout << sum_total_sda_sn2[(*names)]/count_total_sda_sn2[(*names)] << " "; else cout << "NAN ";
    if(count_total_sdb_sn2[(*names)] != 0) cout << sum_total_sdb_sn2[(*names)]/count_total_sdb_sn2[(*names)] << " "; else cout << "NAN ";
    if(count_total_sdc_sn2[(*names)] != 0) cout << sum_total_sdc_sn2[(*names)]/count_total_sdc_sn2[(*names)] << endl; else cout << "NAN" << endl;
  }

  cout << "FRAME SD" << endl;
  names  = sn1_c_names_first.begin();
  for(;names != sn1_c_names_first.end();names++) {
    cout << "SN1 " << (*names) << " ";


    if(count_total_sda_sn1[(*names)] != 0) cout << sqrt((sumsq_total_sda_sn1[(*names)]-(sum_total_sda_sn1[(*names)]*sum_total_sda_sn1[(*names)]/count_total_sda_sn1[(*names)]))/(count_total_sda_sn1[(*names)]-1)) << " "; else cout << "NAN "; 
    if(count_total_sdb_sn1[(*names)] != 0) cout << sqrt((sumsq_total_sdb_sn1[(*names)]-(sum_total_sdb_sn1[(*names)]*sum_total_sdb_sn1[(*names)]/count_total_sdb_sn1[(*names)]))/(count_total_sdb_sn1[(*names)]-1)) << " "; else cout << "NAN "; 
    if(count_total_sdc_sn1[(*names)] != 0) cout << sqrt((sumsq_total_sdc_sn1[(*names)]-(sum_total_sdc_sn1[(*names)]*sum_total_sdc_sn1[(*names)]/count_total_sdc_sn1[(*names)]))/(count_total_sdc_sn1[(*names)]-1)) << endl; else cout << "NAN" << endl; 
  }

  names  = sn2_c_names_first.begin();
  for(;names != sn2_c_names_first.end();names++) {
    cout << "SN2 " << (*names) << " ";

    if(count_total_sda_sn2[(*names)] != 0) cout << sqrt((sumsq_total_sda_sn2[(*names)]-(sum_total_sda_sn2[(*names)]*sum_total_sda_sn2[(*names)]/count_total_sda_sn2[(*names)]))/(count_total_sda_sn2[(*names)]-1)) << " "; else cout << "NAN "; 
    if(count_total_sdb_sn2[(*names)] != 0) cout << sqrt((sumsq_total_sdb_sn2[(*names)]-(sum_total_sdb_sn2[(*names)]*sum_total_sdb_sn2[(*names)]/count_total_sdb_sn2[(*names)]))/(count_total_sdb_sn2[(*names)]-1)) << " "; else cout << "NAN "; 
    if(count_total_sdc_sn2[(*names)] != 0) cout << sqrt((sumsq_total_sdc_sn2[(*names)]-(sum_total_sdc_sn2[(*names)]*sum_total_sdc_sn2[(*names)]/count_total_sdc_sn2[(*names)]))/(count_total_sdc_sn2[(*names)]-1)) << endl;else cout << "NAN" << endl; 
  }


  return 0;
}
 
