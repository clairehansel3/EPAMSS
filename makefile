CXX   = mpicxx
FLAGS = -Wall -Wextra -std=c++17 #-Ofast -flto -ffast-math -DNDEBUG #\
				-I /u/home/c/claireha/prefix/include
OBJS  = source/LambertW/LambertW.o source/input.o source/main.o \
				source/moliere.o source/scattering.o source/solver.o \
				source/statistics.o source/random.o
LIBS  = -lboost_mpi #/u/home/c/claireha/prefix/lib/libboost_mpi.a
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
