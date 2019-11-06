CXX   = mpicxx
FLAGS = -Wall -Wextra -std=c++17 -O3 -flto -ffast-math -DNDEBUG \
				-I /u/home/c/claireha/prefix/include \
				-L /u/home/c/claireha/prefix/lib/
OBJS  = source/LambertW/LambertW.o source/input.o source/main.o \
				source/moliere.o source/scattering.o source/solver.o \
				source/statistics.o source/random.o
LIBS  = -lboost_mpi -lboost_serialization
EXEC  = epamss

$(EXEC): $(OBJS)
	$(CXX) $(FLAGS) $(OBJS) $(LIBS) -o $(EXEC)

%.o: %.cxx
	$(CXX) $(FLAGS) -c $^ -o $@

%.o: %.cc
	$(CXX) $(FLAGS) -c $^ -o $@

.PHONY: clean

clean:
	-rm $(OBJS) $(EXEC)
