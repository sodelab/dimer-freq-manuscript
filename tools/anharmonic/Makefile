# ----------------------------------------------- #
# ** Makefile for mbCO2-vib routine **
# ----------------------------------------------- #


CPP = g++
#CPPFLAGS = -O0 -fopenmp -llapack -lm -lmbCO2 -Wl,-rpath,/Users/oosode/Documents/Research/co2-dimer-freq/tools/library
CPPFLAGS = -O0 -fopenmp -llapack -lm -lmbCO2
USERLIB = -L../library

SRC = $(wildcard *.cpp)
#INC = $(wildcard *.h)
OBJ = $(SRC:.cpp=.o)

install: $(OBJ)
	$(CPP) $(USERLIB) $(CPPFLAGS) $(OBJ) -o anharmonic.x

clean: 
	rm -f *.o anharmonic.x

%.o:%.cpp
	$(CPP) $(CPPFLAGS) $(USERLIB) -c $<
