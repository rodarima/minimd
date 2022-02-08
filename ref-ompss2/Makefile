CXX=I_MPI_CXX=mcxx OMPI_CXX=mcxx mpicxx
LINK=$(CXX)

CXXFLAGS+=--ompss-2 -mavx
LDFLAGS+=--ompss-2

# Use tampi?
CXXFLAGS+=-DUSE_TAMPI
LDFLAGS+=-ltampi-c

#CXXFLAGS+=-O0
CXXFLAGS+=-O2

CXXFLAGS+=-g
LDFLAGS+=-g
#CXXFLAGS+=-k
CXXFLAGS+=--line-markers

CXXFLAGS+=-fsanitize=address
LDFLAGS+=-fsanitize=address

CXXFLAGS+=-Wno-attributes

SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
EXE=miniMD

all: $(EXE)

# Compile
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

# Link
$(EXE): $(OBJ)
	$(LINK) $(LDFLAGS) -o $@ $^

clean:
	rm -f *.o mcxx_* $(EXE)
