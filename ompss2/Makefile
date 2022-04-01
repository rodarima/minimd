CC=clang
#CC=mcxx

MPICC=I_MPI_CC=$(CC) OMPI_CC=$(CC) mpicc
MPILINK=$(MPICC)


#CFLAGS+=-mavx
#CFLAGS+=--ompss-2
#LDFLAGS+=--ompss-2

CFLAGS+=-fompss-2
LDFLAGS+=-fompss-2

LIBS+=-lm

# Use tampi?
CFLAGS+=-DUSE_TAMPI
LIBS+=-lstdc++ -l:libtampi-c.a

#CFLAGS+=-O0
CFLAGS+=-O3

CFLAGS+=-g
LDFLAGS+=-g
#CFLAGS+=-k
#CFLAGS+=--line-markers

CFLAGS+=-std=c11

#CFLAGS+=-fsanitize=address
#LDFLAGS+=-fsanitize=address

#CFLAGS+=-Wno-attributes

SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
EXE=miniMD

all: $(EXE)

# Compile
%.o: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $^

# Link
#$(EXE): $(OBJ)
#	$(MPILINK) $(LDFLAGS) -o $@ $^
$(EXE): ljs.o box.o input.o force.o safe.o packbuf.o comm.o neigh.o thermo.o integrate.o
	$(MPILINK) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o mcxx_* $(EXE)
