# ----------------------------------------------- #
# ** Makefile for mbCO2 routine **
# ----------------------------------------------- #


OPT = optimize

CPP = g++
OP = -O0
CPPFLAGS =
USERLIB = -L../library
INC = -I/opt/local/include
#LIB = -lstdc++ -lgsl -lgslcblas -fopenmp -lmbCO2
LIB = -L/opt/local/lib -lstdc++ -lgsl -lgslcblas -lmbCO2

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

all: $(OPT).x 

$(OPT).x: $(OPT).cpp
	$(CPP) -c io-xyz.cpp -o io-xyz.o
	$(CPP) $(INC) $(LIB) -c $(OPT).cpp -o $(OPT).o
	$(CPP) $(INC) $(CPPFLAGS) $(OP) io-xyz.o $(OPT).o $(LIB) $(USERLIB) -o $(OPT).x

clean: 
	rm -f *.o *.x

%.o:%.cpp
	$(CPP) $(CPPFLAGS) $(OP) -c $<
