# ----------------------------------------------- #
# ** Makefile for mbCO2 routine **
# ----------------------------------------------- #


OPT = optimize

CPP = gcc
OP = -O0
CPPFLAGS =
USERLIB = -L../library
LIB = -lstdc++ -lgsl -lgslcblas -fopenmp -lmbCO2

SRC = $(wildcard *.cpp)
INC = $(wildcard *.h)
OBJ = $(SRC:.cpp=.o)

all: $(OPT).x 

$(OPT).x: $(OPT).cpp
	$(CPP) -c io-xyz.cpp -o io-xyz.o
	$(CPP) $(LIB) -c $(OPT).cpp -o $(OPT).o
	$(CPP) $(CPPFLAGS) $(OP) io-xyz.o $(OPT).o $(LIB) $(USERLIB) -o $(OPT).x

clean: 
	rm -f *.o *.x

%.o:%.cpp
	$(CPP) $(CPPFLAGS) $(OP) -c $<
