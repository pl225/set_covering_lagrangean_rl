CPP= g++
CFLAGS= -O3
SRC_CPP = aux.cpp instance.cpp
OBJECTS_CPP = $(subst .cpp,.o, $(SRC_CPP))

all: set_covering

set_covering: $(OBJECTS_CPP)
	$(CPP) $(CFLAGS) -o set_covering $(OBJECTS_CPP)

%.o: %.cpp
	$(CPP) $(CFLAGS) -o $@ -c $<

clean:
	rm -f set_covering $(OBJECTS_CPP)