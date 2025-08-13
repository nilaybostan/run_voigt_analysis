# Compiler and ROOT flags
CXX = g++
CXXFLAGS = -O2 -Wall -fPIC $(shell root-config --cflags)

# Default target
all: libRoccoR.so

libRoccoR_Dict.cxx : RoccoR.h Linkdef.h
	rootcling -f -rmf libRoccoR.rootmap -rml libRoccoR.so $@ $^
libRoccoR.so : libRoccoR_Dict.cxx RoccoR.cc
	g++ -shared -o $@ $(CXXFLAGS) $^ `root-config --libs`

# Clean everything
clean:
	rm -f *.o *.so *Dict.* *.pcm *.rootmap
