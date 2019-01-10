# ----------------------------------------------- #
# ** Makefile for mbCO2-evaluate routine **
# ----------------------------------------------- #


CPP = g++
CPPFLAGS = -O0 -lm -lmbCO2
USERLIB = -L../library

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

install: $(OBJ)
	$(CPP) $(OBJ) $(USERLIB) $(CPPFLAGS) -o eval-pot.x

clean: 
	rm -f *.o eval-pot.x

%.o:%.cpp
	$(CPP) -c $<
