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

#ifndef PSIM_BOX_H
#define PSIM_BOX_H

class Box {
  typedef double _prec;

public:
  _prec x_min,x_max;
  _prec y_min,y_max;
  _prec z_min,z_max;


  Box(_prec x_min_in,
      _prec x_max_in,
      _prec y_min_in,
      _prec y_max_in,
      _prec z_min_in,
      _prec z_max_in) :
      x_min(x_min_in),
      x_max(x_max_in),
      y_min(y_min_in),
      y_max(y_max_in),
      z_min(z_min_in),
      z_max(z_max_in)
      {
  }

  void extend(const Box &b_new) {
    if(b_new.x_min < x_min) x_min = b_new.x_min;
    if(b_new.x_max > x_max) x_max = b_new.x_max;
    if(b_new.y_min < y_min) y_min = b_new.y_min;
    if(b_new.y_max > y_max) y_max = b_new.y_max;
    if(b_new.z_min < z_min) z_min = b_new.z_min;
    if(b_new.z_max > z_max) z_max = b_new.z_max;
  }
};

inline std::ostream& operator<<(std::ostream& out, const Box &rhs) {

  out << "x_min: " << rhs.x_min;
  out << " x_max: " << rhs.x_max;
  out << " y_min: " << rhs.y_min;
  out << " y_max: " << rhs.y_max;
  out << " z_min: " << rhs.z_min;
  out << " z_max: " << rhs.z_max;
  
  return out;
}

#endif
