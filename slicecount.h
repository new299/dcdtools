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


#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "vec.h"
#include "dcd.h"
#include "box.h"

using namespace std;

class SliceCount {
 
public:

  enum Direction {
    x, y, z
  };
 
private:

  Box bounds;
  DCDFrame frame;
  Direction slice_direction;
  int slices;
  double slice_length;
  vector<int> slices_vector;


public:
  
  
  SliceCount (Box b,DCDFrame f,Direction d,int s) : bounds(b), frame(f), slice_direction(d), slices(s) {
    slicecount_calc();
  }

  double get_slice_length() const {
    return slice_length;
  }

  bool slicecount_calc() {
    
    double length;

    if(slice_direction == SliceCount::x) length = (bounds.x_max-bounds.x_min);
    if(slice_direction == SliceCount::y) length = (bounds.y_max-bounds.y_min);
    if(slice_direction == SliceCount::z) length = (bounds.z_max-bounds.z_min);

    slice_length = length/slices;
    

    vector<int> slicecounts(slices,0);

    for(vector<DCDAtom>::const_iterator i = frame.get_atoms().begin();i != frame.get_atoms().end();i++) {
      if(slice_direction == SliceCount::x) {
        int c = int( floor(((*i).get_position().get_x()-bounds.x_min)/slice_length) );
        if(c == slices) c--;
        slicecounts[c]++;
      }

      if(slice_direction == SliceCount::y) {
        int c = int( floor(((*i).get_position().get_y()-bounds.y_min)/slice_length) );
        if(c == slices) c--;
        slicecounts[c]++;
      }

      if(slice_direction == SliceCount::z) {
        int c = int( floor(((*i).get_position().get_z()-bounds.z_min)/slice_length) );
        if(c == slices) c--;
        slicecounts[c]++;
      }
    }

    slices_vector = slicecounts;

    return true;
  }

  vector<int> get_slices_vector() const {
    return slices_vector;
  }
};

inline std::ostream& operator<<(std::ostream& out,const SliceCount &rhs) {
  vector<int> v = rhs.get_slices_vector();

  for(vector<int>::iterator i = v.begin();i != v.end();i++) {
    out << (*i) << endl;
  } 

  return out;
}
