#    MOLDY Tools 
#
#    Copyright (C) 2007,2008 Nava Whiteford, Justine Taylor
#
#
#    This file is part of MOLDY Tools.
#
#    MOLDY Tools is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MOLDY Tools is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MOLDY Tools.  If not, see <http://www.gnu.org/licenses/>.
#
#
RM = /bin/rm -f

C++ = g++
CCFLAGS = -g -I. 

LDFLAGS = -lpthread 

.PHONY: clean

all: dcd2slicecount dcd2Distance dcd2Angle dcd2anyvec dcd2anyvec_positive dcd2dorder

dcd2slicecount: dcd2slicecount.cpp box.h slicecount.h dcd.h pdb.h stringify.h vec.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) dcd2slicecount.cpp -o dcd2slicecount $(LDFLAGS)

dcd2Distance: dcd2distance.cpp distance.h dcd.h pdb.h stringify.h vec.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) dcd2distance.cpp -o dcd2Distance $(LDFLAGS)
 
dcd2Angle: dcd2angle.cpp angle.h dcd.h pdb.h stringify.h vec.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) dcd2angle.cpp -o dcd2Angle $(LDFLAGS)

dcd2anyvec: dcd2anyvec.cpp genericvec.h dcd.h pdb.h stringify.h vec.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) dcd2anyvec.cpp -o dcd2anyvec $(LDFLAGS)

dcd2anyvec_positive: dcd2anyvec_positive.cpp genericvec_positive.h dcd.h pdb.h stringify.h vec.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) dcd2anyvec_positive.cpp -o dcd2anyvec_positive $(LDFLAGS)

dcd2dorder: dcd2dorder.cpp dorder.h dcd.h pdb.h stringify.h vec.h
	@echo --- $@ ---
	$(C++) $(CCFLAGS) dcd2dorder.cpp -o dcd2dorder $(LDFLAGS)

clean:
	@echo --- $@ ---
	-$(RM) *.o core


