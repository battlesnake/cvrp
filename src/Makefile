CC=g++
CFLAGS=-c -Wall -std=c++11 -g -Wextra
LDFLAGS=-g
LIBS=

ifeq ($(OMP),y)
CFLAGS+=-fopenmp
LDFLAGS+=-fopenmp
LIBS+=-lpthread
else
CFLAGS+=-Wno-unknown-pragmas
LDFLAGS+=-Wno-unknown-pragmas
endif

ifeq ($(O),y)
CFLAGS+=-O2 -flto
LDFLAGS+=-O2 -flto
else
CFLAGS+=-O0
LDFLAGS+=
endif

SOURCES=main.cpp \
	cvrp_idataModel.cpp \
	cvrp_dataModel.cpp \
	cvrp_vehicleTrip.cpp \
	cvrp_solutionModel.cpp \
	cvrp_solutionFinder.cpp \
	cvrp_util.cpp \
	jsoncpp.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=cvrp

all: $(SOURCES) $(EXECUTABLE)

test:
	+$(MAKE) clean
	+$(MAKE) $(EXECUTABLE)
	/bin/time ./$(EXECUTABLE) ../data/data.json

clean :
	rm -f $(EXECUTABLE) *.o
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

jsoncpp.o: CFLAGS+=-w

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
