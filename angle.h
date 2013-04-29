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
#include "pdb.h"

using namespace std;

inline double rad2deg(double rad) {
  return rad*(180/3.14159265);
}


///\brief A Class to represent a DCD2GenericVec
class GenericVec {

public:
  string atom1;
  string atom2;
  string atom3;


public:
  DCDFrame mydcdframe;
  PDBFile mypdbfile;
  vector<int> residue_numbers;
  vector<double> angles;
  ostream &err;

  
  
  GenericVec(DCDFrame f,PDBFile p,string atom_name1,string atom_name2,string atom_name3,ostream &err_in=std::cerr) : mydcdframe(f), mypdbfile(p), atom1(atom_name1), atom2(atom_name2), atom3(atom_name3), err(err_in) {
  
  }

  void genericvec_calc() {
    int cur_resid=-1;
    bool cur_output = false;

    bool toatom_clear = true;
    bool fromatom_clear = true;
    bool toatom2_clear = true;
    

    DCDAtom to_atom;
    DCDAtom from_atom;
    DCDAtom to_atom2;

    DCDFrame::_atoms_type::const_iterator i_d = mydcdframe.get_atoms().begin();
    PDBFile::_atoms_type::const_iterator i_p = mypdbfile.get_atoms().begin();
    for(;(i_p != mypdbfile.get_atoms().end()) && (i_d != mydcdframe.get_atoms().end());i_p++,i_d++) {
       if((*i_p).get_residue_sequence_number() != cur_resid) {cur_output=false; toatom_clear=true; fromatom_clear=true; toatom2_clear=true; cur_resid=(*i_p).get_residue_sequence_number();}

         
       // trim all whitespaces
       string name =  (*i_p).get_atom_name();
       for(int i = 0;i < name.length();) {
         if(name[i] == ' ') {name.erase(i,1);} else i++;
       }

       if(name == atom1) {toatom_clear = false; to_atom = (*i_d);}
       if(name == atom2)  {fromatom_clear  = false; from_atom = (*i_d);}
       if(name == atom3)  {toatom2_clear  = false; to_atom2 = (*i_d);}

       if(toatom2_clear==false && toatom_clear==false && fromatom_clear==false && cur_output==false) {
         cur_output = true;
       
	 Vec<double> our_vector = to_atom.get_position().make_double() - from_atom.get_position().make_double();
         Vec<double> our_vector2 = to_atom2.get_position().make_double() - from_atom.get_position().make_double();

         double angle = rad2deg(acos((our_vector.dot(our_vector2))/(our_vector.mag()*our_vector2.mag())));
	 
	 residue_numbers.push_back(cur_resid);
	 angles.push_back(angle);

       }

    }
  }

  const vector<int> &get_residue_numbers() {
    return residue_numbers;
  }

  const vector<double> &get_angles() {
    return angles;
  }

};
