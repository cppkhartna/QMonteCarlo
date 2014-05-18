MPICXX = mpiCC
CXX = g++
CXXFLAGS = -pg -lm -lstdc++ -std=c++11 -std=gnu++11 -fPIC -O3
CFLAGS = -pg -lm -fPIC

SRCS = mtwist/mtwist.c mtwist/randistrs.c math.cpp QMolecule.cpp QDMC.cpp

OBJS = $(SRCS:.cpp=.o)
OBJS = $(SRCS:.c=.o)

default: QDMC

%.o: %.cpp %.hpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

QDMC: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ 

PDMC: $(OBJS)
	$(MPICXX) $(CXXFLAGS) $^ -o $@ 

qlib: $(OBJS)
	$(CXX) $(CXXFLAGS) -shared -Wl,-soname,qlib.so $^ -o qlib.so

ifneq (clean, $(MAKECMDGOALS))
-include deps.mk
endif

deps.mk: $(SRCS) 
	$(CXX) -MM $^ > $@

run: QDMC
	./QDMC
clean:
	rm -rf *.o 
	rm -rf mtwist/*.o 
