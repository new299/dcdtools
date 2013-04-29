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
#include "stringify.h"
#include "box.h"
#include "slicecount.h"

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
     cout << "Usage: dcd2slicecount <dcdfile filtered> <dcdfile box> <slice direction> <num slices> [totalonly]" << endl;
     cout << " " << endl;
     cout << "The filtered DCD contains only the atoms to be counted, while the box DCD is the entire" << endl;
     cout << "simulation cell.  Both should be centered around the membrane." << endl;
     cout << "Slice direction should be one of x,y or z." << endl;
     cout << "The [totalonly] keyword will reported the sum total number of atoms in each slice for all frames." << endl;
     cout << " " << endl;
     exit(0);
  }
 
  bool totalonly=false;
  if(argc>5) totalonly=true;
  SliceCount::Direction slice_direction;
  if(string(argv[3]) == string("x")) slice_direction = SliceCount::x;
  if(string(argv[3]) == string("y")) slice_direction = SliceCount::y;
  if(string(argv[3]) == string("z")) slice_direction = SliceCount::z;
  
  int slices = convertTo<int>(argv[4]);
 
  DCDFile dcd_box(argv[2]);
  DCDFile dcd_use(argv[1]);

  int n=0;
 
  // Calculate aggregate box
  
  dcd_box.next_frame();
  Box b = dcd_box.get_frame().get_box();
  
  for(;dcd_box.next_frame();) {
    Box b_new = dcd_box.get_frame().get_box();
    
    b.extend(b_new);
  }


  vector<int> total_count(slices,0);

  double slice_length=0;
  for(;dcd_use.next_frame();) {
    DCDFrame f = dcd_use.get_frame();

    SliceCount s(b,f,slice_direction,slices);
    vector<int> current_slices = s.get_slices_vector();
    slice_length = s.get_slice_length();

    for(int i=0;i<total_count.size();i++) {
      total_count[i] += current_slices[i];
    }

    if(!totalonly) cout << s << endl;
  }

  cout << "# Slice width: " << slice_length << endl;
  cout << "# Slice direction: " << string(argv[3]) << endl;
  cout << "# Total slice count: " << endl;
  for(vector<int>::const_iterator i = total_count.begin();i != total_count.end();i++) {
    cout << (*i) << endl;
  }

  return 0;
}
