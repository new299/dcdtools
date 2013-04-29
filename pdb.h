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


#ifndef PDBPARSE_H
#define PDBPARSE_H

#include "vec.h"
#include "stringify.h"
#include <math.h>

using namespace std;

class PDBAtom {
private:
  int    atom_number;
  string atom_tag;
  string atom_name;
  string alt_location;
  string residue_type;
  string chain_id;
  int    residue_sequence_number;
  string insertion_code;
  double occupancy;
  double temperature_factor;
  string element;
  string charge;
  string remainder;
  string blankarea;

  Vec<double> position;

public:
 
  inline const string &get_atom_tag() const {
    return atom_tag;
  }
 
  inline const string &get_atom_name() const {
    return atom_name;
  }

  inline const string &get_alt_location() const {
    return alt_location;
  }

  inline const string &get_residue_type() const {
    return residue_type;
  }

  inline const string &get_chain_id() const {
    return chain_id;
  }

  inline int get_residue_sequence_number() const {
    return residue_sequence_number;
  }

  inline const string &get_insertion_code() const {
    return insertion_code;
  }

  inline double get_occupancy() const {
    return occupancy;
  }

  inline double get_temperature_factor() const {
    return temperature_factor;
  }

  inline const string &get_element() const {
    return element;
  }

  inline const string &get_charge() const {
    return charge;
  }

  inline const Vec<double> &get_position() const {
    return position;
  }

  inline const string &get_remainder() const {
    return remainder;
  }

  inline const string &get_blankarea() const {
    return blankarea;
  }

  void load_line(string tag, std::istream &in) {
    string current;
    char   temp[500];

    atom_tag = tag;

    // 7-11
    in.get(temp,6);
    try { atom_number = convertTo<int>(string(temp)); } 
    catch(...) { atom_number = 0; }

    // 12
    // BLANK CHARACTER
    in.get(temp,2);

    // 13-16
    // FIXED FOUR CHAR STRING    
    in.get(temp,5);
    atom_name = temp;
    
    // 17
    // FIXED TWO CHAR STRING
    in.get(temp,2);
    alt_location = temp;

    // 18-20
    // FIXED THREE CHAR STRING
    in.get(temp,4);
    residue_type = temp;

    // 21
    // BLANK CHAR
    in.get(temp,2);

    // 22
    // FIXED ONE CHAR STRING
    in.get(temp,2);
    chain_id = temp;

    // 23-26
    in.get(temp,5);
    try { residue_sequence_number = convertTo<int>(string(temp)); }
    catch(...) { residue_sequence_number = 0; }

    // 27
    in.get(temp,2);
    insertion_code = temp;
    
    // 28-30
    // 3 BLANK CHAR
    in.get(temp,4);

    // 31-38
    in.get(temp,9);
    try { this->position.set_x(convertTo<double>(string(temp))); }
    catch(...) { this->position.set_x(0); }

    // 39-46
    in.get(temp,9);
    try { this->position.set_y(convertTo<double>(string(temp))); }
    catch(...) { this->position.set_y(0); }

    // 47-54
    in.get(temp,9);
    try { this->position.set_z(convertTo<double>(string(temp))); }
    catch(...) { this->position.set_z(0); }

    // 55-60
    in.get(temp,7);
    try { occupancy = convertTo<double>(string(temp)); }
    catch(...) { occupancy = 0; }

    // 61-66
    in.get(temp,7);
    try { temperature_factor = convertTo<double>(string(temp)); }
    catch(...) { temperature_factor = 0; }


    // 67-66
    // BLANK SPACES
    in.get(temp,11);
    blankarea = temp;
    for(int n=blankarea.size();n<11;n++) {
      blankarea += " ";
    }

    // 77-78
    in.get(temp,3);
    if(!in.eof()) element = temp;
	     else element = "  ";

    // 79-80
    in.get(temp,3);
    if(!in.eof()) charge = temp;
	     else charge = "  ";

    for(;!in.eof();) {
      in.get(temp,1);
      remainder += temp[0];
    }
  }


};

//// \brief Parses and stores a pdb file
class PDBFile {
private:
  vector<PDBAtom> atoms;

public:
  typedef vector<PDBAtom> _atoms_type;

  PDBFile(string filename) {
    load(filename);
  }

  PDBFile(const char *filename) {
    load(string(filename));
  }


  void load(string filename) {
    load(filename.c_str());
  }

  void load(const char *filename) {
    ifstream infile(filename);

    // strip starting line
    string str1;
    getline(infile, str1);
    char temp[50];

    for(;!infile.eof();) {
      string str;
      getline(infile, str);
      istringstream iss(str);

      if(!iss.eof()) {
	string tag;
	iss.get(temp,7);
	tag = temp;

	if((tag.compare("ATOM  ") == 0) || (tag.compare("HETATM") == 0)) {
	  PDBAtom p;
	  p.load_line(tag,iss);
	  atoms.push_back(p);
	}
      }
    }
  }

  const vector<PDBAtom> &get_atoms() {
    return atoms;
  }
};

#endif
