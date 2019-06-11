all:
	g++-9 -Wall -Wextra -std=c++17 -O3 -flto -ffast-math -DNDEBUG -I \
	../../source/ -I /usr/local/include ../../source/moliere.cxx \
	../../source/random.cxx moliere.cxx /usr/local/lib/libboost_mpi.dylib \
	-o moliere
