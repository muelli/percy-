##  Percy++
##  Copyright 2007 Ian Goldberg <iang@cs.uwaterloo.ca>
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of version 2 of the GNU General Public License as
##  published by the Free Software Foundation.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  There is a copy of the GNU General Public License in the COPYING file
##  packaged with this plugin; if you cannot find it, write to the Free
##  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
##  02111-1307  USA

CXXFLAGS=-Wall -g -O2 -I/usr/local/include/NTL
LDLIBS=-lntl -lgmp

TARGETS=pirserver pirclient splitdatabase

RUNFILES=database database.* out.client out.real

CLIENT_O=percyclient.o percyparams.o recover.o percyio.o rr_roots.o ZZ_pXY.o
SERVER_O=percyserver.o percyparams.o datastore.o percyio.o

all: $(TARGETS)

pirserver: pirserver.o $(SERVER_O)
	g++ -o $@ $^ $(LDLIBS)

pirclient: pirclient.o $(CLIENT_O)
	g++ -o $@ $^ $(LDLIBS)

splitdatabase: splitdatabase.o percyio.o
	g++ -o $@ $^ $(LDLIBS)

clean:
	-rm -f *.o

veryclean: clean
	-rm -f $(TARGETS)

distclean: veryclean
	-rm -f $(RUNFILES)
