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
LIBS+=-lstdc++ -l:libtampi-c.a
#LIBS+=-lGPI2 -ltagaspi

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

SRC=box.c \
  comm.c \
  comm_send.c \
  comm_recv.c \
  comm_pack.c \
  comm_unpack.c \
  comm_wait.c \
  force.c \
  input.c \
  integrate.c \
  ljs.c \
  neigh.c \
  packbuf.c \
  packbuf_mpi.c \
  ref.c \
  safe.c \
  thermo.c

#  packbuf_gaspi.c
OBJ=$(SRC:.c=.o)
EXE=miniMD

all: $(EXE)

#.depend: $(SRC)
#	rm -f "$@"
#	$(MPICC) $(CFLAGS) -MM $^ -MF "$@"
#
#include .depend

# Compile
%.o: %.c
	$(MPICC) $(CFLAGS) -c -o $@ $^

$(EXE): $(OBJ)
	$(MPILINK) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(OBJ) $(EXE)
