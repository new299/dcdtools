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
#include <map>
#include "dcd.h"
#include "pdb.h"
#include <iostream>

using namespace std;

template <class _first,class _second,class _third>
class triple {

public:

_first first;
_second second;
_third third;

  triple() {
  }

  triple(_first first_in,_second second_in,_third third_in) : first(first_in), second(second_in), third(third_in) {
  }

};

inline double rad2deg(double rad) {
  return rad*(180/3.14159265);
}

///\brief A Class to represent a DOrder
class DOrder {

public:
  DCDFrame mydcdframe;
  PDBFile mypdbfile;
  
  vector<int>                   residue_numbers;
  map<int, map<string,double> > angles_sn1_sda;
  map<int, map<string,double> > angles_sn1_sdb;
  map<int, map<string,double> > angles_sn1_sdc;
  
  map<int, map<string,Vec<float> > > sn1_c_positions;
  map<int, map<string,Vec<float> > > sn1_h1_positions;
  map<int, map<string,Vec<float> > > sn1_h2_positions;
  map<int, map<string,Vec<float> > > sn1_h3_positions;

  map<int, map<string,double> > angles_sn2_sda;
  map<int, map<string,double> > angles_sn2_sdb;
  map<int, map<string,double> > angles_sn2_sdc;
  
  map<int, map<string,Vec<float> > > sn2_c_positions;
  map<int, map<string,Vec<float> > > sn2_h1_positions;
  map<int, map<string,Vec<float> > > sn2_h2_positions;
  map<int, map<string,Vec<float> > > sn2_h3_positions;

  ostream &err;
  string dictionary_file;

  vector<string> sn1_c_names;
  vector<string> sn1_h1_names;
  vector<string> sn1_h2_names;
  vector<string> sn1_h3_names;

  vector<string> sn2_c_names;
  vector<string> sn2_h1_names;
  vector<string> sn2_h2_names;
  vector<string> sn2_h3_names;
    
  map<string, triple<string, string, string> > sn1_c_to_h;
  map<string, triple<string, string, string> > sn2_c_to_h;

  DOrder(DCDFrame f,PDBFile p,string dictionary_file_in, ostream &err_in=std::cerr) : mydcdframe(f), mypdbfile(p), dictionary_file(dictionary_file_in), err(err_in) {
    
    load_dictionary(dictionary_file);
  
  }
  
  void process_line(ifstream &infile,vector<string> &names) {
    string str;
    
    getline(infile, str);
    istringstream iss(str);

    // get sn1_c_names
    for(int n=0;!iss.eof();n++) {
      string s;
      iss >> s;
      names.push_back(s);
    }
  }

  void load_dictionary(string filename) {
    ifstream infile(filename.c_str());

    process_line(infile,sn1_c_names);
    process_line(infile,sn1_h1_names);
    process_line(infile,sn1_h2_names);
    process_line(infile,sn1_h3_names);

    process_line(infile,sn2_c_names);
    process_line(infile,sn2_h1_names);
    process_line(infile,sn2_h2_names);
    process_line(infile,sn2_h3_names);

    // sn1: create c to h mapping
    for(int i = 0;i < sn1_c_names.size();i++) {
      sn1_c_to_h[sn1_c_names[i]] = triple<string,string,string>(sn1_h1_names[i],sn1_h2_names[i],sn1_h3_names[i]);
    }

    // sn2: create c to h mapping
    for(int i = 0;i < sn2_c_names.size();i++) {
      sn2_c_to_h[sn2_c_names[i]] = triple<string,string,string>(sn2_h1_names[i],sn2_h2_names[i],sn2_h3_names[i]);
    }

  }

  double calc_angle(Vec<float> c_position,Vec<float> h_position) {
    // calculate SDA angle here
    Vec<double> refvec;
    refvec = Vec<double>(0,0,1);

    Vec<double> vecA = h_position.make_double() - c_position.make_double();
    // double angle = rad2deg(acos(refvec.dot(vecA) / (refvec.mag()*vecA.mag())));
    double angle = refvec.dot(vecA) / (refvec.mag()*vecA.mag());

    double out = (3*angle*angle)-1;

    return out;
  }

  void dorder_calc() {
    int cur_resid=-1;

    DCDFrame::_atoms_type::const_iterator i_d = mydcdframe.get_atoms().begin();
    PDBFile::_atoms_type::const_iterator i_p = mypdbfile.get_atoms().begin();
    cur_resid=(*i_p).get_residue_sequence_number();
    bool resok=false;
    for(;(i_p != mypdbfile.get_atoms().end()) && (i_d != mydcdframe.get_atoms().end());i_p++,i_d++) {
    
      if((*i_p).get_residue_sequence_number() != cur_resid) {
        if(resok == true) residue_numbers.push_back(cur_resid);
        cur_resid=(*i_p).get_residue_sequence_number();
        resok=false;
      }

       // trim leading spaces
       string name =  (*i_p).get_atom_name();
       for(bool del=true;del == true;) {
         del = false;
         if(name[0] == ' ') {del=true; name.erase(0,1);}
       }

       // lookup name
       if(find(sn1_c_names.begin(),sn1_c_names.end(),name)   != sn1_c_names.end())  {resok=true; sn1_c_positions [cur_resid][name] = (*i_d).get_position(); } else
       if(find(sn1_h1_names.begin(),sn1_h1_names.end(),name) != sn1_h1_names.end()) {            sn1_h1_positions[cur_resid][name] = (*i_d).get_position(); } else
       if(find(sn1_h2_names.begin(),sn1_h2_names.end(),name) != sn1_h2_names.end()) {            sn1_h2_positions[cur_resid][name] = (*i_d).get_position(); } else
       if(find(sn1_h3_names.begin(),sn1_h3_names.end(),name) != sn1_h3_names.end()) {            sn1_h3_positions[cur_resid][name] = (*i_d).get_position(); } else
       if(find(sn2_c_names.begin(),sn2_c_names.end(),name)   != sn2_c_names.end())  {resok=true; sn2_c_positions [cur_resid][name] = (*i_d).get_position(); } else
       if(find(sn2_h1_names.begin(),sn2_h1_names.end(),name) != sn2_h1_names.end()) {            sn2_h1_positions[cur_resid][name] = (*i_d).get_position(); } else
       if(find(sn2_h2_names.begin(),sn2_h2_names.end(),name) != sn2_h2_names.end()) {            sn2_h2_positions[cur_resid][name] = (*i_d).get_position(); } else
       if(find(sn2_h3_names.begin(),sn2_h3_names.end(),name) != sn2_h3_names.end()) {            sn2_h3_positions[cur_resid][name] = (*i_d).get_position(); }
    }

    if(resok == true) residue_numbers.push_back(cur_resid);

    // frame complete, calculate angles

    for(vector<int>::iterator i = residue_numbers.begin();i != residue_numbers.end();i++) {
      
      // SN1
      for(vector<string>::iterator sn1_c_names_i= sn1_c_names.begin();sn1_c_names_i != sn1_c_names.end();sn1_c_names_i++) {
	Vec<float> c_position = sn1_c_positions[(*i)][(*sn1_c_names_i)];

        if(sn1_c_to_h[(*sn1_c_names_i)].first  != string("NA")) {
	  Vec<float> h1_position = sn1_h1_positions[(*i)][sn1_c_to_h[(*sn1_c_names_i)].first];
	  angles_sn1_sda[(*i)][(*sn1_c_names_i)] = calc_angle(c_position,h1_position); 
        } else angles_sn1_sda[(*i)][(*sn1_c_names_i)] = std::numeric_limits<double>::quiet_NaN();

	if(sn1_c_to_h[(*sn1_c_names_i)].second != string("NA")) {
	  Vec<float> h2_position = sn1_h2_positions[(*i)][sn1_c_to_h[(*sn1_c_names_i)].second];
          angles_sn1_sdb[(*i)][(*sn1_c_names_i)] = calc_angle(c_position,h2_position);
	} else angles_sn1_sdb[(*i)][(*sn1_c_names_i)] = std::numeric_limits<double>::quiet_NaN();

	if(sn1_c_to_h[(*sn1_c_names_i)].third != string("NA")) {
	  Vec<float> h3_position = sn1_h3_positions[(*i)][sn1_c_to_h[(*sn1_c_names_i)].third];
          angles_sn1_sdc[(*i)][(*sn1_c_names_i)] = calc_angle(c_position,h3_position);
        } else angles_sn1_sdc[(*i)][(*sn1_c_names_i)] = std::numeric_limits<double>::quiet_NaN();
      }
   
      // SN2
      for(vector<string>::iterator sn2_c_names_i= sn2_c_names.begin();sn2_c_names_i != sn2_c_names.end();sn2_c_names_i++) {
	Vec<float> c_position = sn2_c_positions[(*i)][(*sn2_c_names_i)];

        if(sn2_c_to_h[(*sn2_c_names_i)].first  != string("NA")) {
	  Vec<float> h1_position = sn2_h1_positions[(*i)][sn2_c_to_h[(*sn2_c_names_i)].first];
          angles_sn2_sda[(*i)][(*sn2_c_names_i)] = calc_angle(c_position,h1_position); 
        } else angles_sn2_sda[(*i)][(*sn2_c_names_i)] = std::numeric_limits<double>::quiet_NaN();

	if(sn2_c_to_h[(*sn2_c_names_i)].second != string("NA")) {
	  Vec<float> h2_position = sn2_h2_positions[(*i)][sn2_c_to_h[(*sn2_c_names_i)].second];
          angles_sn2_sdb[(*i)][(*sn2_c_names_i)] = calc_angle(c_position,h2_position);
        } else angles_sn2_sdb[(*i)][(*sn2_c_names_i)] = std::numeric_limits<double>::quiet_NaN();

	if(sn2_c_to_h[(*sn2_c_names_i)].third != string("NA")) {
	  Vec<float> h3_position = sn2_h3_positions[(*i)][sn2_c_to_h[(*sn2_c_names_i)].third];
          angles_sn2_sdc[(*i)][(*sn2_c_names_i)] = calc_angle(c_position,h3_position);
        } else angles_sn2_sdc[(*i)][(*sn2_c_names_i)] = std::numeric_limits<double>::quiet_NaN();
      }
 
    }
  }

  const vector<int> &get_residue_numbers() {
    return residue_numbers;
  }

  const vector<string> &get_sn1_c_names() {
    return sn1_c_names;
  }
  
  const vector<string> &get_sn2_c_names() {
    return sn2_c_names;
  }

  map<string, double > &get_angles_sn1_sda(int residue) {
    return angles_sn1_sda[residue];
  }

  map<string, double > &get_angles_sn1_sdb(int residue) {
    return angles_sn1_sdb[residue];
  }
  
  map<string, double > &get_angles_sn1_sdc(int residue) {
    return angles_sn1_sdc[residue];
  }
  
  map<string, double > &get_angles_sn2_sda(int residue) {
    return angles_sn2_sda[residue];
  }

  map<string, double > &get_angles_sn2_sdb(int residue) {
    return angles_sn2_sdb[residue];
  }

  map<string, double > &get_angles_sn2_sdc(int residue) {
    return angles_sn2_sdc[residue];
  }

};
