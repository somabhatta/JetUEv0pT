# Build of the per-event analyzer jetv0ana (ROOT, FastJet + fjcontrib
# ConstituentSubtractor, PYTHIA 8).
#
# The production binary was built on stbc (el9) with the LCG_110 gcc16 view:
#   source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_110 x86_64-el9-gcc16-opt
#   make
# Use the same view and flags to reproduce the production trees bit by bit.
# (JEWEL itself must NOT be built with gcc16, see jewel/build_jewel.sh.)

LCG_VIEW ?= /cvmfs/sft.cern.ch/lcg/views/LCG_110/x86_64-el9-gcc16-opt

CXX      := g++
CXXFLAGS := -O2 $(shell root-config --cflags) $(shell fastjet-config --cxxflags) $(shell pythia8-config --cxxflags) -DWITH_PYTHIA8 -std=c++23
LIBS     := $(shell root-config --libs) $(shell fastjet-config --libs --plugins) -L$(LCG_VIEW)/lib -lConstituentSubtractor $(shell pythia8-config --libs)

all: jetv0ana

jetv0ana: jetv0ana.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f jetv0ana

.PHONY: all clean
