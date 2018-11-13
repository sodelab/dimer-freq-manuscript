# ----------------------------------------------- #
# ** Makefile for mbCO2 routine **
# ----------------------------------------------- #

MBCO2 = libmbCO2

CPP = gcc
OP = -O0
CPPFLAGS = -fPIC 
SONAME = -install_name
INC = -I/usr/local/include
LIB = -lstdc++

all: $(MBCO2).so $(MBCO2).a 

$(MBCO2).so: mbCO2.cpp
	$(CPP) $(CPPFLAGS) -c x1b.cpp -o x1b.o
	$(CPP) $(CPPFLAGS) -c x2b.cpp -o x2b.o
	$(CPP) $(CPPFLAGS) -c poly-1b.cpp -o poly-1b.o
	$(CPP) $(CPPFLAGS) -c poly-2b.cpp -o poly-2b.o
	$(CPP) $(CPPFLAGS) -c sapt-s.cpp -o sapt-s.o
	$(CPP) $(CPPFLAGS) -c mbCO2.cpp -o mbCO2.o
	$(CPP) $(OP) $(CPPFLAGS) -shared -o $(MBCO2)-NITROGEN.so mbCO2.o x1b.o x2b.o poly-1b.o poly-2b.o sapt-s.o $(LIB) 

$(MBCO2).a: mbCO2.cpp
	ar -rv $(MBCO2).a mbCO2.o x1b.o x2b.o poly-1b.o poly-2b.o sapt-s.o

clean: 
	rm -f *.o *.so *.a

%.o:%.cpp
	$(CPP) $(CPPFLAGS) $(OP) -c $<
