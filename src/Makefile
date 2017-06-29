CXXFLAGS:=-c -Wall -std=c++14 -Wextra -ffunction-sections -fdata-sections -march=native
LDFLAGS:=-Wl,--gc-sections
LIBS=

OMP?=y
O?=y

ifeq ($(OMP),y)
CXXFLAGS+=-fopenmp
LDFLAGS+=-fopenmp
LIBS+=-lpthread
else
CXXFLAGS+=-Wno-unknown-pragmas
LDFLAGS+=-Wno-unknown-pragmas
endif

ifeq ($(O),y)
CXXFLAGS+=-O2 -flto -s
LDFLAGS+=-O2 -flto -s
else
CXXFLAGS+=-O0 -g
LDFLAGS+=-g
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

all: $(EXECUTABLE)

test:
	+$(MAKE) clean
	+$(MAKE) $(EXECUTABLE)
	/bin/time ./$(EXECUTABLE) ../data/data.json

bench:
	@echo 'Testing with OpenMP enabled'
	+@$(MAKE) -s --no-print-directory HIDE_PROGRESS=y BENCH=y O=y OMP=y test
	@echo 'Testing with OpenMP disabled'
	+@$(MAKE) -s --no-print-directory HIDE_PROGRESS=y BENCH=y O=y OMP=n test

clean :
	rm -f $(EXECUTABLE) *.o

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

jsoncpp.o: CXXFLAGS+=-w
