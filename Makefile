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

CLIENT_O=percyclient.o percyparams.o recover.o percyio.o rr_roots.o ZZ_pXY.o \
		gf28.o
SERVER_O=percyserver.o percyparams.o datastore.o percyio.o gf28.o
SRCS=$(subst .o,.cc,$(CLIENT_O) $(SERVER_O) pirclient.o pirserver.o splitdatabase.o percyio.o)

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

depend:
	makedepend -Y -- $(CXXFLAGS) -- $(SRCS) 2>/dev/null

# DO NOT DELETE

percyclient.o: /usr/local/include/NTL/vec_vec_ZZ_p.h
percyclient.o: /usr/local/include/NTL/ZZ_pX.h recover.h
percyclient.o: /usr/local/include/NTL/vec_ZZ_p.h percyresult.h gf28.h
percyclient.o: percyclient.h percyparams.h
percyparams.o: percyparams.h percyio.h /usr/local/include/NTL/ZZ.h
recover.o: /usr/local/include/NTL/mat_ZZ_p.h /usr/local/include/NTL/ZZ_pX.h
recover.o: rr_roots.h ZZ_pXY.h /usr/local/include/NTL/ZZ_pXFactoring.h
recover.o: recover.h /usr/local/include/NTL/vec_ZZ_p.h percyresult.h gf28.h
percyio.o: percyio.h /usr/local/include/NTL/ZZ.h
rr_roots.o: ZZ_pXY.h /usr/local/include/NTL/ZZ_pX.h
rr_roots.o: /usr/local/include/NTL/ZZ_pXFactoring.h
ZZ_pXY.o: ZZ_pXY.h /usr/local/include/NTL/ZZ_pX.h
ZZ_pXY.o: /usr/local/include/NTL/ZZ_pXFactoring.h
gf28.o: gf28.h
percyserver.o: /usr/local/include/NTL/vec_ZZ_p.h percyserver.h datastore.h
percyserver.o: /usr/local/include/NTL/ZZ.h percyparams.h
percyparams.o: percyparams.h percyio.h /usr/local/include/NTL/ZZ.h
datastore.o: datastore.h /usr/local/include/NTL/ZZ.h percyparams.h percyio.h
percyio.o: percyio.h /usr/local/include/NTL/ZZ.h
gf28.o: gf28.h
pirclient.o: pstream.h /usr/local/include/NTL/ZZ_p.h percyclient.h
pirclient.o: /usr/local/include/NTL/vec_vec_ZZ_p.h percyresult.h
pirclient.o: percyparams.h gf28.h config.h version.h
pirserver.o: datastore.h /usr/local/include/NTL/ZZ.h percyparams.h
pirserver.o: percyserver.h /usr/local/include/NTL/vec_ZZ_p.h config.h
pirserver.o: version.h
splitdatabase.o: percyio.h /usr/local/include/NTL/ZZ.h config.h version.h
percyio.o: percyio.h /usr/local/include/NTL/ZZ.h
