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


#ifndef PROTTOOLS_DCD_H
#define PROTTOOLS_DCD_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "vec.h"
#include "box.h"

using namespace std;

typedef int int_32b;
typedef unsigned int unsigned_int_32b;

typedef float real_32b;

///\brief reverses the endianness of a string
template <typename T> inline
void reverse_endian(T& t){
  char* res = reinterpret_cast<char*>(&t);
  char *temp = new char[sizeof(T)];
  for(int n=0;n<sizeof(T);n++) temp[sizeof(T)-n] = res[n];
  for(int n=0;n<sizeof(T);n++) res[n] = temp[n];
  delete[] temp; //free(temp);
}


///\brief corrects the endian (calls reverse_endian if reverse is true)
template <typename T> inline
void correct_endian(T& t,bool reverse) {
  if(!reverse) return; else reverse_endian(t);
}

class DCDAtom {

public:

  DCDAtom() {
    position = Vec<real_32b>(0,0,0);
  }

  DCDAtom(Vec<real_32b> v) {
    position = v;
  }

  const Vec<real_32b> &get_position() const {
    return position;
  }

private:
  Vec<real_32b> position;
};

///\brief A Class to represent a DCD Frame
class DCDFrame {

protected:
  vector<DCDAtom> atoms; 

public:
  typedef vector<DCDAtom> _atoms_type;
  const vector<DCDAtom> &get_atoms() const {
    return atoms;
  }

  inline void push_back(DCDAtom a) {
    atoms.push_back(a);
  }

  inline void clear() {
    atoms.clear();
  }

  Box get_box() {
    
    double min_x = (*atoms.begin()).get_position().get_x();
    double max_x = (*atoms.begin()).get_position().get_x();
    double min_y = (*atoms.begin()).get_position().get_y();
    double max_y = (*atoms.begin()).get_position().get_y();
    double min_z = (*atoms.begin()).get_position().get_z();
    double max_z = (*atoms.begin()).get_position().get_z();


    for(vector<DCDAtom>::iterator i = atoms.begin();i != atoms.end();i++) {
      if(min_x > (*i).get_position().get_x()) min_x = (*i).get_position().get_x();
      if(max_x < (*i).get_position().get_x()) max_x = (*i).get_position().get_x();
      if(min_y > (*i).get_position().get_y()) min_y = (*i).get_position().get_y();
      if(max_y < (*i).get_position().get_y()) max_y = (*i).get_position().get_y();
      if(min_z > (*i).get_position().get_z()) min_z = (*i).get_position().get_z();
      if(max_z < (*i).get_position().get_z()) max_z = (*i).get_position().get_z();
    }

    return Box(min_x,max_x,min_y,max_y,min_z,max_z);

  }
};


///\brief A Class to represent a DCD File
class DCDFile {
private:

  DCDFrame current_frame;
  int_32b free_indexes;
  int_32b unit_cell_flag;
  int_32b atom_count;
  int_32b id;
  int file_ok;
  bool native_endian;
  ifstream *dcd_ptr;
  ostream &err;

public:

  DCDFile(string filename,ostream &err_in=std::cerr) : err(err_in) {
    load(filename);
  }

  DCDFile(const char* filename,ostream &err_in=std::cerr) : err(err_in) {
    load(filename);
  }

  bool detect_dcd(const char* filename,bool &native_endian, ostream &err) {
    ifstream dcd(filename);
    if(!dcd.is_open()) {
      err << "Could not open input file (detect_dcd)" << endl;
      return false;
    }

    char cord_str[5];
    dcd.read(reinterpret_cast<char *>(&id),sizeof(int_32b));
    dcd.read(cord_str,4);
    cord_str[4] = 0;
    dcd.close();

    int_32b id_reverse_endian = id;
    reverse_endian(id_reverse_endian);
   
    file_ok = true;
    if(string(cord_str) != "CORD") file_ok = false;
    if(id == 84)                   native_endian=true; else
    if(id_reverse_endian == 84)    native_endian=false; else file_ok=false;

    return file_ok;
  }
    
  bool read_coordinates_array(vector<real_32b> &coordinates, ostream &err) {
    ifstream &dcd = *dcd_ptr;

    int_32b c_atom_count;
      
    // Coordinate
    dcd.read(reinterpret_cast<char*>(&c_atom_count),sizeof(int_32b));
    correct_endian(c_atom_count,!native_endian);
    if((c_atom_count/(sizeof(int_32b))) != atom_count) { err << "ERROR: First atom counts do not match " << c_atom_count/(sizeof(int_32b)) << "!=" << atom_count << " (read_coordinates_array)" << endl; dcd.close(); return false;}

    for(int c=0;c<atom_count;c++) {
      real_32b coord;
      dcd.read(reinterpret_cast<char*>(&coord),sizeof(real_32b));
      correct_endian(coord,!native_endian);
      coordinates.push_back(coord);
    }

    dcd.read(reinterpret_cast<char*>(&c_atom_count),sizeof(int_32b));
    correct_endian(c_atom_count,!native_endian);
    if((c_atom_count/(sizeof(int_32b))) != atom_count) { err << "ERROR: Second atom counts do not match " << c_atom_count/(sizeof(int_32b)) << "!=" << atom_count << " (load_dcd)" << endl; dcd.close(); return false;}

    return true;
  }
  
  bool load(string filename) {
    load(filename.c_str());
  }

  bool load(const char* filename) {

    detect_dcd(filename,native_endian,err);
    
    dcd_ptr = new ifstream(filename);
    ifstream &dcd = *dcd_ptr;
    if(!dcd.is_open()) {
      err << "ERROR: Could not open input file (load_dcd)" << endl;
      return false;
    }

    // free_indexes
    dcd.seekg(40,ios::beg);
    dcd.read(reinterpret_cast<char*>(&free_indexes),sizeof(int_32b));
    correct_endian(free_indexes,!native_endian);
    
    // unit cell flag
    dcd.seekg(48,ios::beg);
    dcd.read(reinterpret_cast<char*>(&unit_cell_flag),sizeof(int_32b));
    correct_endian(unit_cell_flag,!native_endian);

    // skip header
    int_32b id;
    dcd.seekg(88,ios::beg);
    dcd.read(reinterpret_cast<char*>(&id),sizeof(int_32b));
    correct_endian(id,!native_endian);
    if(id != 84) {err << "ERROR: File id error in namd dcd file (load_dcd)" << endl; dcd.close(); return false;}

    // skip titles
    unsigned_int_32b l,m,n;

    dcd.read(reinterpret_cast<char*>(&l),sizeof(int_32b));
    correct_endian(l,!native_endian);

    dcd.read(reinterpret_cast<char*>(&m),sizeof(int_32b));
    correct_endian(m,!native_endian);
   
    string comment;
    for(unsigned int i=0,j=0;i<m;i++,j += 81) {
      char tempstr[81];
      dcd.read(tempstr,80);
      comment += tempstr + '\n';
    }

    dcd.read(reinterpret_cast<char*>(&n),sizeof(int_32b));
    correct_endian(n,!native_endian);

    if(((l-4)%80) != 0) { err << "ERROR: File Comment field not correct size (load_dcd)" << endl; dcd.close(); return false;}
    if(n != l) { err << "ERROR: Flag0 incorrect (load_dcd)" << endl; dcd.close(); return false;} 

    dcd.read(reinterpret_cast<char*>(&n),sizeof(int_32b));
    correct_endian(n,!native_endian);
    if(n != 4) { err << "ERROR: Flag1 incorrect (load_dcd)" << endl; dcd.close(); return false;}

    dcd.read(reinterpret_cast<char*>(&atom_count),sizeof(int_32b));
    correct_endian(atom_count,!native_endian);

    dcd.read(reinterpret_cast<char*>(&n),sizeof(int_32b));
    correct_endian(n,!native_endian);
    if(n != 4) { err << "ERROR: Flag2 incorrect (load_dcd)" << endl; dcd.close(); return false;}

    if(free_indexes > 0) dcd.seekg(4*(atom_count - free_indexes+2),ios::cur);

      
      
  }

  bool next_frame() {
    ifstream &dcd = *dcd_ptr;
    current_frame.clear();

    if(!dcd.eof()) {
      vector<real_32b> x_coordinates;
      vector<real_32b> y_coordinates;
      vector<real_32b> z_coordinates;

      bool res=true;

      // right now I'm throwing away unit cell information
      if(unit_cell_flag == 1) {
        int_32b unitcellsize;

        dcd.read(reinterpret_cast<char*>(&unitcellsize),sizeof(int_32b));
        correct_endian(atom_count,!native_endian);
        if(unitcellsize != 48) {if(dcd.eof()) return false; else err << "ERROR: Unit cell size incorrect" << endl;}

	char unitcell[48];
        dcd.read(reinterpret_cast<char*>(&unitcell),sizeof(int_32b)*12);
        
	dcd.read(reinterpret_cast<char*>(&unitcellsize),sizeof(int_32b));
        correct_endian(atom_count,!native_endian);
        if(unitcellsize != 48) {if(dcd.eof()) return false; else err << "ERROR: Unit cell size incorrect" << endl;}
        if(dcd.eof()) return false;
      }

      int_32b temp;
      res = read_coordinates_array(x_coordinates,err);
      if(res) res = read_coordinates_array(y_coordinates,err);
      if(res) res = read_coordinates_array(z_coordinates,err);

      vector<real_32b>::iterator i_x = x_coordinates.begin();
      vector<real_32b>::iterator i_y = y_coordinates.begin();
      vector<real_32b>::iterator i_z = z_coordinates.begin();

      for(int n=0;(i_x != x_coordinates.end()) || (i_y != y_coordinates.end()) || (i_z != z_coordinates.end());n++) {
	
	DCDAtom atom(Vec<real_32b>((*i_x),(*i_y),(*i_z)));
	current_frame.push_back(atom);
	i_x++;
        i_y++;
        i_z++;
      }
      return true;
    } else return false;
  }

  DCDFrame &get_frame() {
    return current_frame;
  }

};

#endif
