# ----------------------------------------------- #
# ** Makefile for mbCO2-vib routine **
# ----------------------------------------------- #


CPP = g++
CPPFLAGS = -O0 -lm -lmbCO2
USERLIB = -L../library

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

install: $(OBJ)
	$(CPP) $(USERLIB) $(CPPFLAGS) $(OBJ) -o eval-pot.x

clean: 
	rm -f *.o eval-pot.x

%.o:%.cpp
	$(CPP) $(CPPFLAGS) $(USERLIB) -c $<
