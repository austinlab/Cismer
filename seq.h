#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <ctype.h>

#define CHUNK_SIZE 2048

typedef struct Seq 	
{	
	char *annotation;	
	char *sequence; 
} Seq;

struct Seq *getSeq(FILE *fasta_stream);
unsigned int count_seqs(char *fasta_file);
char *rc(char *seq, int length);
