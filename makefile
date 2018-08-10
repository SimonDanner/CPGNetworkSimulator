cc=clang++
LINK=ln
BUILDDIR = ./build
SRCDIR = ./src
INCDIR = ./include
INSTALLDIR = /bin
PROFILE =  #-fprofile-use #-fprofile-generate
BOOSTINCLUDEDIR = /home/simon/boost_1_67_0/
CPPFLAGS = -I$(BOOSTINCLUDEDIR) -I$(INCDIR) -std=c++14 -I$(SRCDIR)  -Ofast -march=native  $(PROFILE) $(ADDCPPFLAGS) -fPIC
LINKFLAGS = -lyaml-cpp
OUTFILE = a.out


all: $(BUILDDIR)/Network.o $(BUILDDIR)/CPGNetworkSimulator.o $(BUILDDIR)/typedefs.o executable python

python:  $(BUILDDIR)/Network.o $(BUILDDIR)/CPGNetworkSimulator.o $(BUILDDIR)/typedefs.o $(SRCDIR)/CPGNetworkSimulatorPY.cpp
	$(cc) $(CPPFLAGS) -undefined dynamic_lookup `python3 -m pybind11 --includes`  -shared $(SRCDIR)/CPGNetworkSimulatorPY.cpp $(BUILDDIR)/CPGNetworkSimulator.o  $(BUILDDIR)/typedefs.o $(BUILDDIR)/Network.o -o $(BUILDDIR)/CPGNetworkSimulator`python3-config --extension-suffix`

executable:  $(BUILDDIR)/Network.o $(BUILDDIR)/CPGNetworkSimulator.o $(BUILDDIR)/Main.o $(BUILDDIR)/typedefs.o
	$(cc)  $(BUILDDIR)/CPGNetworkSimulator.o $(BUILDDIR)/Main.o $(BUILDDIR)/typedefs.o $(BUILDDIR)/Network.o $(CPPFLAGS) $(LINKFLAGS)  -o $(BUILDDIR)/$(OUTFILE)  

$(BUILDDIR)/CPGNetworkSimulator.o: $(SRCDIR)/CPGNetworkSimulator.cpp $(INCDIR)/CPGNetworkSimulator.hpp $(INCDIR)/typedefs.h $(INCDIR)/Network.hpp
	$(cc) -c $(SRCDIR)/CPGNetworkSimulator.cpp -o $(BUILDDIR)/CPGNetworkSimulator.o $(CPPFLAGS)

$(BUILDDIR)/Network.o: $(SRCDIR)/Network.cpp $(INCDIR)/Network.hpp $(INCDIR)/typedefs.h
	$(cc) -c $(SRCDIR)/Network.cpp -o $(BUILDDIR)/Network.o $(CPPFLAGS)

$(BUILDDIR)/typedefs.o: $(SRCDIR)/typedefs.cpp $(INCDIR)/typedefs.h
	$(cc) -c  $(SRCDIR)/typedefs.cpp  -o $(BUILDDIR)/typedefs.o $(CPPFLAGS)

$(BUILDDIR)/Main.o: $(SRCDIR)/Main.cpp $(INCDIR)/Network.hpp $(INCDIR)/CPGNetworkSimulator.hpp
	$(cc) -c $(SRCDIR)/Main.cpp  -o $(BUILDDIR)/Main.o $(CPPFLAGS)


clean:
	rm $(BUILDDIR)/*.o


