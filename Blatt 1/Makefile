CC = gcc
CFLAGS = -std=c99 -Wall -D_GNU_SOURCE
.SUFFIXES: .o .c
.c.o: ; $(CC) -c $(CFLAGS) $*.c

OBJ = fields.o blas.o blatt1.o

blatt0: $(OBJ)
	$(CC) $(CFLAGS) -o blatt1 $(OBJ) -lm

fields.o:	fields.h
blas.o:		blas.h

clean:
	rm -rf *.o blatt1
 
