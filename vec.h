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


#ifndef PSIM_VEC_H
#define PSIM_VEC_H

#include <iostream>
#include <math.h>

using namespace std;

/// \brief Class to represent a 3D vector
template <class _prec> class Vec {
public:
  _prec x;
  _prec y;
  _prec z;
  bool valid;

  Vec() {
    x=0;
    y=0;
    z=0;
    valid=false;
  }

  Vec(_prec x_in,_prec y_in,_prec z_in) {
    x = x_in;
    y = y_in;
    z = z_in;
  }

  Vec(const Vec<_prec> &v) {
    x=v.x;
    y=v.y;
    z=v.z;
    valid=v.valid;
  }

  void clear() {
    x=0;
    y=0;
    z=0;
    valid=false;
  }

  void randomize(_prec rand_max) {
    x = (((double) rand())/INT_MAX)*rand_max-(rand_max/2);
    y = (((double) rand())/INT_MAX)*rand_max-(rand_max/2);
    z = (((double) rand())/INT_MAX)*rand_max-(rand_max/2);
  }

  Vec<_prec> &operator= (const Vec<_prec> &rhs) {
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    valid = rhs.valid;

    return *this;
  }

  Vec<_prec> &operator= (const _prec rhs) {
    x = rhs;
    y = rhs;
    z = rhs;

    return *this;
  }

  _prec dot(const Vec<_prec> &rhs) const {
    return (x*rhs.x) + (y*rhs.y) + (z*rhs.z);
  }

  _prec len() const {
    return (exp(0.5*log((x*x) + (y*y) + (z*z))));
  }

  _prec sum() {
    return x + y + z;
  }

  Vec<_prec> operator+ (const Vec<_prec> &rhs) const {
    Vec<_prec> out(*this);

    out.x += rhs.x;
    out.y += rhs.y;
    out.z += rhs.z;
    return out;
  }

  Vec<_prec> operator- (const Vec<_prec> &rhs) const {
    Vec<_prec> out(*this);

    out.x -= rhs.x;
    out.y -= rhs.y;
    out.z -= rhs.z;
    return out;
  }

  Vec<_prec> operator* (const _prec rhs) const {
    Vec<_prec> out(*this);

    out.x = x*rhs;
    out.y = y*rhs;
    out.z = z*rhs;
    return out;
  }

  Vec<_prec> operator* (const Vec<_prec> &rhs) const {
    Vec<_prec> out(*this);
    
    out.x = x*rhs.x;
    out.y = y*rhs.y;
    out.z = z*rhs.z;
    return out;
  }

  Vec<_prec> operator/ (const _prec rhs) const {
    Vec<_prec> out(*this);

    out.x = x/rhs;
    out.y = y/rhs;
    out.z = z/rhs;
    return out;
  }

  Vec<_prec> operator/ (const Vec<_prec> &rhs) const {
    Vec<_prec> out(*this);

    out.x = x/rhs.x;
    out.y = y/rhs.y;
    out.z = z/rhs.z;
    return out;
  }

  Vec<_prec> operator-= (const _prec rhs) {
    x -= rhs;
    y -= rhs;
    z -= rhs;
  }

  Vec<_prec> cross(const Vec<_prec> &c) const {
    Vec<_prec> p;
    p.x = (y * c.z) - (z * c.y);
    p.y = (z * c.x) - (x * c.z);
    p.z = (x * c.y) - (y * c.x);
    p.set_valid();
    return p;
  }

  _prec angle(const Vec<_prec> &c) {
    _prec cos_val, sin_val, sinsq;

    cos_val = dot(c)/(len()*c.len());
    sinsq = 1.0-cos_val*cos_val;
    if(sinsq < 0.0) sinsq = 0.0;
    sin_val = sqrt(sinsq);

    return(atan2(sin_val,cos_val) * 57.29578); // pi radians
  }

  _prec distance(const Vec<_prec> &c) {
    Vec<_prec> dist;
    dist = (*this) - c;
    dist.valid = true;

    _prec radial = sqrt(dist.dot(dist));

    return radial; 
  }

  _prec mag() {
    return sqrt(x*x + y*y + z*z);  
  }

  Vec<double> make_double() const {
    Vec<double> v;
    v.set_x(x);
    v.set_y(y);
    v.set_z(z);
    v.valid = valid;
  
    return v;
  }

  void set_vec(_prec x_in,_prec y_in,_prec z_in) {
    x = x_in;
    y = y_in;
    z = z_in;
  }

  void set_x(_prec x_in) {
    x = x_in;
  }

  void set_y(_prec y_in) {
    y = y_in;
  }

  void set_z(_prec z_in) {
    z = z_in;
  }

  void set_valid() {
    valid = true;
  }

  void set_invalid() {
    valid = false;
  }

  bool get_valid() {
    return valid;
  }

  _prec get_x() const { return x; }
  _prec get_y() const { return y; }
  _prec get_z() const { return z; }

};

template <class _prec>
inline Vec<_prec> operator+ (_prec lhs, const Vec<_prec> &rhs) {
  Vec<_prec> out(lhs);

  out.x += rhs.x;
  out.y += rhs.y;
  out.z += rhs.z;
  return out;
}

template <class _prec>
inline Vec<_prec> operator- (_prec lhs, const Vec<_prec> &rhs) {
  Vec<_prec> out(lhs);

  out.x -= rhs.x;
  out.y -= rhs.y;
  out.z -= rhs.z;
  return out;
}

template <class _prec>
inline Vec<_prec> operator* (_prec lhs, const Vec<_prec> &rhs) {
  Vec<_prec> out;

  out.x = rhs.x * lhs;
  out.y = rhs.y * lhs;
  out.z = rhs.z * lhs;
  return out;
}

template <class _prec>
inline std::ostream& operator<<(std::ostream& out, const Vec<_prec> &rhs) {
 
  out << "x: " << rhs.x << " y: " << rhs.y << " z: " << rhs.z;
  if(rhs.valid==false) out << " (not tagged valid)";
  return out;  
}

template <class _prec>
_prec dihedral_angle(Vec<_prec> p1,Vec<_prec> p2,Vec<_prec> p3,Vec<_prec> p4) {
   
  Vec<_prec> v21 = p2 - p1;
  Vec<_prec> v23 = p2 - p3;
  Vec<_prec> v43 = p4 - p3;

  Vec<_prec> p   = v21.cross(v23);
  Vec<_prec> q   = v23.cross(v43);

  _prec theta = p.angle(q);
  Vec<_prec> r= p.cross(q);
  if(v23.dot(r) < 0) theta = -theta;

  return theta;
}

#endif
