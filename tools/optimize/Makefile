# ----------------------------------------------- #
# ** Makefile for mbCO2-optimize routine **
# ----------------------------------------------- #


OPT = optimize

CPP = g++
OP = -O0
USERLIB = -L../library
LIB = -lstdc++ -lgsl -lgslcblas -lmbCO2

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

install: $(OBJ)
	$(CPP) $(OBJ) $(INC) $(LIB) $(USERLIB) -o $(OPT).x

clean: 
	rm -f *.o *.x

%.o:%.cpp
	$(CPP) -c $<
