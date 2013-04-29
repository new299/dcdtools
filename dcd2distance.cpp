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
#include "distance.h"

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
     cout << "Usage: dcd2Distance <dcdfile> <pdbfile> <atom1> <atom2>" << endl;
     cout << " " << endl;
     cout << "This program calculates the distance between two atoms, from atom1 to atom2." << endl;
     cout << " " << endl;
     exit(0);
  }
  
  DCDFile mydcd(argv[1]);
  PDBFile mypdb(argv[2]);
   
  double slice_length=0;
  for(;mydcd.next_frame();) {
    DCDFrame f = mydcd.get_frame();

    GenericVec genericvec(f,mypdb,string(argv[3]),string(argv[4]));
    genericvec.genericvec_calc();
    const vector<int>    &resnums = genericvec.get_residue_numbers();
    const vector<double> &meas_dist  = genericvec.get_meas_dist();

    vector<int>::const_iterator i_r=resnums.begin();
    vector<double>::const_iterator i_a=meas_dist.begin();

    cout << "# FRAME ";
    cout << resnums.size() << endl;
    cout << "# From: " << string(argv[3]) << endl; 
    cout << "# To: " << string(argv[4]) << endl;
    for(;i_r != resnums.end();i_r++,i_a++) {
      cout << (*i_r) << " " << (*i_a) << endl;
    }

  }

  return 0;
}
