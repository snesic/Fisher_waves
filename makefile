CXX = g++
CXXFLAGS = -Wno-deprecated -w -O4 -Wall

DIFF = ./sdiff
PRE = ./
MAJOR = 1
MINOR = 0

%.o:           	%.cpp
		$(CXX) $(CXXFLAGS) -c $*.cpp

everything:    	fisher_waves

fisher_waves_obj = fisher_waves.o simulation.o splitting_routines.o random_no_generators.o initial_conditions.o read_write_msg.o
fisher_waves:    $(fisher_waves_obj)
		$(CXX) -o $@ $(fisher_waves_obj) -L. -lm -lfftw3 $(CXXFLAGS)



fisher_waves.o: fisher_waves.cpp simulation.h read_write_msg.h initial_conditions.h

simulation.o: simulation.cpp splitting_routines.h

splitting_routines.o: splitting_routines.cpp random_no_generators.h

random_no_generators.o: random_no_generators.cpp

read_write_msg.o: read_write_msg.cpp

initial_conditions.o: initial_conditions.cpp

instabilities.txx:   	fisher_waves
		$(PRE) fisher_waves > fisher_waves.txx
		$(DIFF) fisher_waves.txt fisher_waves.txx


