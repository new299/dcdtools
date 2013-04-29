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


#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>

#ifndef STRINGIFY_H
#define STRINGIFY_H

class BadConversion : public std::runtime_error {
public:
 BadConversion(const std::string& s)
      : std::runtime_error(s) {}
};

template<class _type>
inline std::string stringify(_type x)
{
  std::ostringstream o;
  if (!(o << std::fixed << x))
    throw BadConversion("stringify()");
  return o.str();
}

template<class _type>
inline _type convertTo(const std::string& s)
{
  std::istringstream i(s);
  _type x;
  if (!(i >> x))
    throw BadConversion("convertTo(\"" + s + "\")");
  return x;
}

inline std::string string_find_and_replace(const std::string s1_in,const std::string s2,std::string s3) {
  std:string s1(s1_in);
  if(s2 == "") throw std::out_of_range("string_find_and_replace: find empty string");

  for(string::size_type i = s1.find(s2,0);i != string::npos;i = s1.find(s2,i)) {
    s1.erase(i,s2.length());
    s1.insert(i,s3);
  }

  return s1;
}

#endif
