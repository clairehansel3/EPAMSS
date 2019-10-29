CXX   = mpicxx
FLAGS = -Wall -Wextra -std=c++17 -fsanitize=undefined
OBJS  = ../../source/LambertW/LambertW.o ../../source/input.o \
				../../source/main.o ../../source/moliere.o ../../source/scattering.o \
				../../source/solver.o ../../source/statistics.o ../../source/random.o
LIBS  = /usr/local/lib/libboost_mpi.dylib
EXEC  = epamss

$(EXEC): $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o $(EXEC)

%.o: %.cxx
	$(CXX) $(FLAGS) -c $^ -o $@

%.o: %.cc
	$(CXX) $(FLAGS) -c $^ -o $@

.PHONY: clean

clean:
	rm $(OBJS) $(EXEC)
