ALL = cismer 
OBJ = seq.o 
INCLUDE = seq.h 
CFLAGS = -g -Wall 

all: ${ALL}
clean:
	-rm ${ALL}
	-rm ${OBJ}

cismer: seq.o
	gcc -O2 ${CFLAGS} -lm seq.o -o cismer cismer.c
cismer.o: ${INCLUDE} cismer.c
	gcc -c ${CFLAGS} cismer.c
seq.o: ${INCLUDE} seq.c
	gcc -c ${CFLAGS} seq.c

