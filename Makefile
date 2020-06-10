CPP= g++
CFLAGS= -O3
SRC_CPP = aux.cpp instance.cpp lagrangean.cpp main.cpp
OBJECTS_CPP = $(subst .cpp,.o, $(SRC_CPP))

all: set_covering

set_covering: $(OBJECTS_CPP) expknap.o
	$(CPP) $(CFLAGS) -o set_covering $(OBJECTS_CPP) expknap.o

%.o: %.cpp
	$(CPP) $(CFLAGS) -o $@ -c $<

expknap.o: expknap.c
	gcc -c expknap.c -lm

clean:
	rm -f set_covering $(OBJECTS_CPP)