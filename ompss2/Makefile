CXX=I_MPI_CXX=mcxx mpicxx
LINK=$(CXX)

CXXFLAGS=--ompss-2 -O3 -DMPICH_IGNORE_CXX_SEEK -DNOCHUNK -mavx \
      -UUSE_SIMD -DUSE_TAMPI -DUSE_TASKS -DPRECISION=2

LDFLAGS=--ompss-2 -O3 -ltampi

CXXFLAGS+=-g
LDFLAGS+=-g
#CXXFLAGS+=-k
CXXFLAGS+=--line-markers

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
	rm -f *.o $(EXE)
