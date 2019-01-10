# ----------------------------------------------- #
# ** Makefile for mbCO2-vib routine **
# ----------------------------------------------- #


CPP = g++
CPPFLAGS = -O0 -fopenmp -llapack -lm -lmbCO2
USERLIB = -L../library

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

install: $(OBJ)
	$(CPP) $(OBJ) $(USERLIB) $(CPPFLAGS) -o anharmonic.x

clean: 
	rm -f *.o anharmonic.x

%.o:%.cpp
	$(CPP) -c $<
