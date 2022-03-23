CXX=clang++
#CXX=mcxx

MPICXX=I_MPI_CXX=$(CXX) OMPI_CXX=$(CXX) mpicxx
MPILINK=$(MPICXX)

#CXXFLAGS+=-mavx
#CXXFLAGS+=--ompss-2
#LDFLAGS+=--ompss-2

#CXXFLAGS+=-fompss-2
#LDFLAGS+=-fompss-2

# Use tampi?
#CXXFLAGS+=-DUSE_TAMPI
#LDFLAGS+=-ltampi-c

CXXFLAGS+=-O0
#CXXFLAGS+=-O3

CXXFLAGS+=-g
LDFLAGS+=-g
#CXXFLAGS+=-k
#CXXFLAGS+=--line-markers

CXXFLAGS+=-fsanitize=address
LDFLAGS+=-fsanitize=address

#CXXFLAGS+=-Wno-attributes

SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
EXE=miniMD

all: $(EXE)

# Compile
%.o: %.cpp
	$(MPICXX) $(CXXFLAGS) -c -o $@ $^

# Link
#$(EXE): $(OBJ)
#	$(MPILINK) $(LDFLAGS) -o $@ $^
$(EXE): ljs.o box.o input.o force.o safe.o packbuf.o comm.o neigh.o thermo.o integrate.o
	$(MPILINK) $(LDFLAGS) -o $@ $^

clean:
	rm -f *.o mcxx_* $(EXE)
