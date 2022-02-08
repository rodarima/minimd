CXX=mpicxx
LINK=$(CXX)

#CXXFLAGS+=-O0
CXXFLAGS+=-O2 -g

CXXFLAGS+=-fsanitize=address
LDFLAGS+=-fsanitize=address

CXXFLAGS+=-Wno-attributes

CXXFLAGS+=-DDONT_CHECK_DELTAX

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
