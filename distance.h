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


public:
  DCDFrame mydcdframe;
  PDBFile mypdbfile;
  vector<int> residue_numbers;
  vector<double> meas_dist;
  ostream &err;

  
  
  GenericVec(DCDFrame f,PDBFile p,string atom_name1,string atom_name2,ostream &err_in=std::cerr) : mydcdframe(f), mypdbfile(p), atom1(atom_name1), atom2(atom_name2), err(err_in) {
  
  }

  void genericvec_calc() {
    int cur_resid=-1;
    bool cur_output = false;

    bool toatom_clear = true;
    bool fromatom_clear = true;


    DCDAtom to_atom;
    DCDAtom from_atom;

    DCDFrame::_atoms_type::const_iterator i_d = mydcdframe.get_atoms().begin();
    PDBFile::_atoms_type::const_iterator i_p = mypdbfile.get_atoms().begin();
    for(;(i_p != mypdbfile.get_atoms().end()) && (i_d != mydcdframe.get_atoms().end());i_p++,i_d++) {
       if((*i_p).get_residue_sequence_number() != cur_resid) {cur_output=false; toatom_clear=true; fromatom_clear=true; cur_resid=(*i_p).get_residue_sequence_number();}

         
       // trim all whitespaces
       string name =  (*i_p).get_atom_name();
       for(int i = 0;i < name.length();) {
         if(name[i] == ' ') {name.erase(i,1);} else i++;
       }

       if(name == atom2) {toatom_clear = false; to_atom = (*i_d);}
       if(name == atom1)  {fromatom_clear  = false; from_atom = (*i_d);}

       if(toatom_clear==false && fromatom_clear==false && cur_output==false) {
         cur_output = true;
       
	 Vec<double> our_vector = to_atom.get_position().make_double() - from_atom.get_position().make_double();

         double distance = our_vector.mag();
	 
	 residue_numbers.push_back(cur_resid);
	 meas_dist.push_back(distance);

       }

    }
  }

  const vector<int> &get_residue_numbers() {
    return residue_numbers;
  }

  const vector<double> &get_meas_dist() {
    return meas_dist;
  }

};
