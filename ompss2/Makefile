CXX=I_MPI_CXX=mcxx mpicxx
LINK=$(CXX)

CXXFLAGS=--ompss-2 -O3 -mavx -DUSE_TAMPI

LDFLAGS=--ompss-2 -O3 -ltampi

CXXFLAGS+=-g
LDFLAGS+=-g
#CXXFLAGS+=-k
CXXFLAGS+=--line-markers

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
